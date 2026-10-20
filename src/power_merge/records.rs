//! Ordered FASTA reporting with persistent worker state. Only encoded result
//! bytes cross back to the writer, never a dense document vector per record.
use anyhow::{ensure, Context, Result};
use deko::read::AnyDecoder;
use rayon::prelude::*;
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Read, Write},
    path::Path,
    sync::mpsc::{self, Receiver, SyncSender},
    thread,
};

pub(super) type Record = (Vec<u8>, Vec<u8>);
type Batch = Vec<Record>;

#[derive(Default)]
struct Output {
    chunks: Vec<Vec<u8>>,
    active: usize,
}

pub(super) fn query_file<S: Send>(
    path: &Path,
    batch_size: usize,
    init: impl Fn() -> S,
    query: impl Fn(&mut S, &[Record], &mut Vec<u8>) -> Result<()> + Sync,
) -> Result<()> {
    ensure!(batch_size > 0, "record batch size must be positive");
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let input: Box<dyn Read + Send> = if crate::is_compressed_path(path) {
        Box::new(AnyDecoder::new(file))
    } else {
        Box::new(file)
    };
    let stdout = std::io::stdout();
    // StdoutLock cannot cross threads. Stdout itself locks only the large
    // buffered writes, while the dedicated writer overlaps lookup and parsing.
    let mut out = BufWriter::with_capacity(1024 * 1024, stdout);
    query_reader_chunks(BufReader::new(input), &mut out, batch_size, init, query)?;
    out.flush()?;
    Ok(())
}

fn query_reader_chunks<S: Send>(
    reader: impl BufRead + Send,
    out: &mut (impl Write + Send),
    batch_size: usize,
    init: impl Fn() -> S,
    query: impl Fn(&mut S, &[Record], &mut Vec<u8>) -> Result<()> + Sync,
) -> Result<()> {
    ensure!(batch_size > 0, "record batch size must be positive");
    writeln!(out, "query_id\tvalid_kmers\thits")?;
    let workers = rayon::current_num_threads();
    // Two queued batches bound read-ahead/output memory. Recycle channels only
    // carry these bounded in-flight buffers; nobody waits for a recycled buffer.
    let (input_tx, input_rx) = mpsc::sync_channel(2);
    let (input_recycle_tx, input_recycle_rx) = mpsc::channel();
    let (output_tx, output_rx) = mpsc::sync_channel::<Output>(2);
    let (output_recycle_tx, output_recycle_rx) = mpsc::channel();
    let records = thread::scope(|scope| -> Result<usize> {
        let parser =
            scope.spawn(move || read_batches(reader, batch_size, input_tx, input_recycle_rx));
        let writer = scope.spawn(move || -> Result<()> {
            for output in output_rx {
                for chunk in &output.chunks[..output.active] {
                    out.write_all(chunk)?;
                }
                let _ = output_recycle_tx.send(output);
            }
            Ok(())
        });
        // Stay on the caller's Rayon pool (including custom/single-thread pools).
        // Dedicated parser/writer threads never consume lookup worker slots.
        let compute = compute_batches(
            input_rx,
            input_recycle_tx,
            output_tx,
            output_recycle_rx,
            workers,
            init,
            query,
        );
        let writing = writer
            .join()
            .map_err(|_| anyhow::anyhow!("record writer panicked"))?;
        let parsing = parser
            .join()
            .map_err(|_| anyhow::anyhow!("record parser panicked"))?;
        writing?;
        compute?;
        parsing
    })?;
    eprintln!(
        "benchmark_records={} batch_size={} workers={}",
        records, batch_size, workers
    );
    Ok(())
}

fn read_batches(
    mut reader: impl BufRead,
    batch_size: usize,
    send: SyncSender<Batch>,
    recycle: Receiver<Batch>,
) -> Result<usize> {
    let mut batch = Vec::with_capacity(batch_size);
    let mut spare_records = Vec::new();
    let mut current: Option<(Vec<u8>, Vec<u8>)> = None;
    let mut line = Vec::new();
    let mut records = 0usize;
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        while line.last().is_some_and(|b| b.is_ascii_whitespace()) {
            line.pop();
        }
        if line.is_empty() {
            continue;
        }
        if line[0] == b'>' {
            if let Some(record) = current.take() {
                batch.push(record);
                records += 1;
                if batch.len() == batch_size {
                    if send.send(batch).is_err() {
                        return Ok(records);
                    }
                    spare_records = recycle.try_recv().unwrap_or_default();
                    batch = Vec::with_capacity(batch_size);
                }
            }
            let id = line[1..]
                .split(|b| b.is_ascii_whitespace())
                .next()
                .unwrap_or_default();
            ensure!(!id.is_empty(), "empty FASTA identifier");
            std::str::from_utf8(id).context("invalid FASTA identifier")?;
            let (mut name, mut seq) = spare_records.pop().unwrap_or_default();
            name.clear();
            name.extend_from_slice(id);
            seq.clear();
            current = Some((name, seq));
        } else {
            let record = current.as_mut().context("--records requires FASTA")?;
            ensure!(
                !line.iter().any(|b| b.is_ascii_whitespace()),
                "whitespace inside FASTA sequence"
            );
            record.1.extend_from_slice(&line);
        }
    }
    if let Some(record) = current {
        batch.push(record);
        records += 1;
    }
    ensure!(records > 0, "no FASTA records");
    if !batch.is_empty() {
        let _ = send.send(batch);
    }
    Ok(records)
}

fn compute_batches<S: Send>(
    input: Receiver<Batch>,
    input_recycle: mpsc::Sender<Batch>,
    output: SyncSender<Output>,
    output_recycle: Receiver<Output>,
    threads: usize,
    init: impl Fn() -> S,
    query: impl Fn(&mut S, &[Record], &mut Vec<u8>) -> Result<()> + Sync,
) -> Result<()> {
    let mut states = Vec::new();
    while let Ok(mut batch) = input.recv() {
        // If parsing is ahead, amortise one Rayon dispatch over two already
        // available batches. Never wait for more work or exceed this bound.
        let first_len = batch.len();
        let mut second = input.try_recv().ok();
        if let Some(next) = &mut second {
            batch.append(next);
        }
        let chunk_size = batch.len().div_ceil(threads.min(batch.len()));
        let active = batch.len().div_ceil(chunk_size);
        while states.len() < active {
            states.push(init());
        }
        let mut encoded = output_recycle.try_recv().unwrap_or_default();
        encoded.active = active;
        encoded
            .chunks
            .resize_with(encoded.chunks.len().max(active), Vec::new);
        states[..active]
            .par_iter_mut()
            .zip(encoded.chunks[..active].par_iter_mut())
            .zip(batch.par_chunks(chunk_size))
            .try_for_each(|((state, bytes), records)| -> Result<()> {
                bytes.clear();
                query(state, records, bytes)
            })?;
        // Ordered chunks and FIFO batches preserve exact input record order.
        if let Some(mut spare) = second {
            spare.extend(batch.drain(first_len..));
            let _ = input_recycle.send(spare);
        }
        let _ = input_recycle.send(batch);
        if output.send(encoded).is_err() {
            break;
        }
    }
    Ok(())
}

#[cfg(test)]
fn query_reader<S: Send>(
    reader: impl BufRead + Send,
    out: &mut (impl Write + Send),
    batch_size: usize,
    init: impl Fn() -> S,
    query: impl Fn(&mut S, &[u8], &mut Vec<u8>) -> Result<()> + Sync,
) -> Result<()> {
    query_reader_chunks(reader, out, batch_size, init, |state, records, bytes| {
        for (id, seq) in records {
            bytes.extend_from_slice(id);
            bytes.push(b'\t');
            query(state, seq, bytes)?;
        }
        Ok(())
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ordered_records_reuse_state_across_batches_and_short_final_chunks() {
        let input = b">a description\r\nAC\r\nGT\r\n>b\nT\n>c\n\n>d\nAAA\n>e\nC\n";
        let expected = b"query_id\tvalid_kmers\thits\na\t4\t\nb\t1\t\nc\t0\t\nd\t3\t\ne\t1\t\n";
        for threads in [1, 2, 4] {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| {
                    for batch in [1, 2, 3, 4, 8] {
                        let mut output = Vec::new();
                        query_reader(
                            &input[..],
                            &mut output,
                            batch,
                            Vec::<u8>::new,
                            |state, seq, out| {
                                state.clear();
                                state.extend_from_slice(seq);
                                writeln!(out, "{}\t", state.len())?;
                                Ok(())
                            },
                        )
                        .unwrap();
                        assert_eq!(output, expected);
                    }
                });
        }
    }

    #[test]
    fn invalid_record_inputs_are_rejected() {
        for (input, batch) in [
            (&b""[..], 2),
            (&b"ACGT\n"[..], 2),
            (&b">\nACGT\n"[..], 2),
            (&b">a\nAC GT\n"[..], 2),
            (&b">a\nACGT\n"[..], 0),
        ] {
            assert!(query_reader(input, &mut Vec::new(), batch, || (), |_, _, _| Ok(())).is_err());
        }
    }

    #[test]
    fn query_and_writer_errors_cancel_bounded_pipeline_without_deadlock() {
        struct BrokenWriter;
        impl Write for BrokenWriter {
            fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
                if bytes.starts_with(b"query_id") {
                    Ok(bytes.len())
                } else {
                    Err(std::io::Error::new(
                        std::io::ErrorKind::BrokenPipe,
                        "test closed output",
                    ))
                }
            }
            fn flush(&mut self) -> std::io::Result<()> {
                Ok(())
            }
        }
        let input = ">r\nACGT\n".repeat(257);
        assert!(query_reader(
            input.as_bytes(),
            &mut BrokenWriter,
            1,
            || (),
            |_, _, out| {
                writeln!(out, "1\t0:1")?;
                Ok(())
            }
        )
        .is_err());
        let mut output = Vec::new();
        assert!(query_reader(
            input.as_bytes(),
            &mut output,
            1,
            || (),
            |_, _, _| { anyhow::bail!("test lookup error") }
        )
        .unwrap_err()
        .to_string()
        .contains("test lookup error"));
    }

    #[test]
    fn long_multibatch_output_is_ordered_with_bounded_buffer_reuse() {
        let mut input = String::new();
        let mut expected = "query_id\tvalid_kmers\thits\n".to_owned();
        for i in 0..1003 {
            input.push_str(&format!(">r{i}\n{}\n", "A".repeat(i % 97)));
            expected.push_str(&format!("r{i}\t{}\t\n", i % 97));
        }
        rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap()
            .install(|| {
                let mut output = Vec::new();
                query_reader(
                    input.as_bytes(),
                    &mut output,
                    13,
                    || (),
                    |_, seq, out| {
                        writeln!(out, "{}\t", seq.len())?;
                        Ok(())
                    },
                )
                .unwrap();
                assert_eq!(output, expected.as_bytes());
            });
    }
}
