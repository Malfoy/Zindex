//! Canonical SimdMinimizers FASTAs with bounded, ordered parallel processing.
//! Each selected position is emitted once, including repeated seed values.
use anyhow::{ensure, Result};
use flate2::{write::GzEncoder, Compression};
use packed_seq::{PackedSeqVec, SeqVec};
use rayon::prelude::*;
use serde::Serialize;
use std::{
    io::{Read, Write},
    sync::mpsc::{self, SyncSender},
    thread,
};

#[derive(Clone, Copy, Debug)]
pub struct Config {
    /// Seed length, NOT the span of the complete minimizer window.
    pub k: usize,
    /// Number of candidate k-mers per minimizer window.
    pub w: usize,
    pub chunk_bases: usize,
    pub batch_chunks: usize,
    /// None: plain FASTA. Some(0..=9): concatenated gzip members.
    pub gzip_level: Option<u32>,
}

impl Config {
    pub fn validate(self) -> Result<()> {
        ensure!((1..=64).contains(&self.k), "k must be in 1..=64");
        ensure!((1..=65535).contains(&self.w), "window must be in 1..=65535");
        ensure!(
            (self.k + self.w - 1) % 2 == 1,
            "canonical selection requires k + window - 1 to be odd"
        );
        ensure!(
            (1..=268_435_456).contains(&self.chunk_bases),
            "chunk-bases must be in 1..=268435456"
        );
        ensure!(self.batch_chunks > 0, "batch-chunks must be positive");
        ensure!(
            self.gzip_level.is_none_or(|level| level <= 9),
            "gzip level must be in 0..=9"
        );
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Default, Serialize, PartialEq, Eq)]
pub struct Stats {
    pub records: u64,
    pub eligible_kmer_positions: u64,
    pub selected_positions: u64,
}

impl std::ops::AddAssign for Stats {
    fn add_assign(&mut self, other: Self) {
        self.records += other.records;
        self.eligible_kmer_positions += other.eligible_kmer_positions;
        self.selected_positions += other.selected_positions;
    }
}

struct Job {
    seq: Vec<u8>,
    record: u64,
    base: usize,
    owned_start: usize,
    owned_end: usize,
}

#[derive(Default)]
struct Encoded {
    bytes: Vec<u8>,
    selected: u64,
}

#[derive(Default)]
struct Worker {
    packed: PackedSeqVec,
    positions: Vec<u32>,
    marks: Vec<u64>,
    plain: Vec<u8>,
}

#[inline]
fn is_dna(base: u8) -> bool {
    matches!(base, b'A' | b'C' | b'G' | b'T')
}

#[inline]
fn complement(base: u8) -> u8 {
    // Only called for uppercase ACGT. The ASCII complement mapping is an XOR.
    base ^ if base & 2 == 0 { 21 } else { 4 }
}

fn encode_job(job: &Job, worker: &mut Worker, out: &mut Encoded, config: Config) -> Result<()> {
    worker.plain.clear();
    out.selected = 0;
    let mut offset = 0;
    let mut canonical = [0u8; 64];
    let mut integer = itoa::Buffer::new();
    for run in job.seq.split(|&b| !is_dna(b)) {
        if run.len() >= config.k + config.w - 1 {
            worker.packed.clear();
            worker.packed.push_ascii(run);
            worker.positions.clear();
            let _ = simd_minimizers::canonical_minimizers(config.k, config.w)
                .run(worker.packed.as_slice(), &mut worker.positions);
            // Canonical minimizer positions need not be monotone and may recur
            // non-adjacently. A reusable bitset provides exact sorted uniqueness
            // in linear time, without comparison sorting or changing tie rules.
            let positions = run.len() - config.k + 1;
            worker.marks.resize(positions.div_ceil(64), 0);
            worker.marks.fill(0);
            for &position in &worker.positions {
                let position = position as usize;
                let absolute = job.base + offset + position;
                if absolute >= job.owned_start && absolute < job.owned_end {
                    worker.marks[position / 64] |= 1 << (position % 64);
                }
            }
            for (block, &marked) in worker.marks.iter().enumerate() {
                let mut bits = marked;
                while bits != 0 {
                    let position = block * 64 + bits.trailing_zeros() as usize;
                    bits &= bits - 1;
                    let word = &run[position..position + config.k];
                    let reverse = word
                        .iter()
                        .copied()
                        .cmp(word.iter().rev().map(|&b| complement(b)))
                        .is_gt();
                    let word = if reverse {
                        for (out, &base) in canonical[..config.k].iter_mut().zip(word.iter().rev())
                        {
                            *out = complement(base);
                        }
                        &canonical[..config.k]
                    } else {
                        word
                    };
                    worker.plain.extend_from_slice(b">r");
                    worker
                        .plain
                        .extend_from_slice(integer.format(job.record).as_bytes());
                    worker.plain.extend_from_slice(b"_p");
                    worker
                        .plain
                        .extend_from_slice(integer.format(job.base + offset + position).as_bytes());
                    worker.plain.push(b'\n');
                    worker.plain.extend_from_slice(word);
                    worker.plain.push(b'\n');
                    out.selected += 1;
                }
            }
        }
        offset += run.len() + 1;
    }
    out.bytes.clear();
    if let Some(level) = config.gzip_level {
        if !worker.plain.is_empty() {
            // Each independent chunk is compressed by its Rayon worker. Ordered
            // concatenated members are a standard gzip stream; compression is
            // parallel even for one very long input chromosome.
            let mut gzip = GzEncoder::new(&mut out.bytes, Compression::new(level));
            gzip.write_all(&worker.plain)?;
            gzip.finish()?;
        }
    } else {
        std::mem::swap(&mut worker.plain, &mut out.bytes);
    }
    Ok(())
}

struct Producer {
    config: Config,
    send: SyncSender<Vec<Job>>,
    batch: Vec<Job>,
    seq: Vec<u8>,
    base: usize,
    next: usize,
    length: usize,
    run_length: usize,
    stats: Stats,
}

impl Producer {
    fn push_job(&mut self, end: usize, context_end: usize) -> Result<()> {
        if end == self.next {
            return Ok(());
        }
        self.batch.push(Job {
            seq: self.seq[..context_end - self.base].to_vec(),
            record: self.stats.records,
            base: self.base,
            owned_start: self.next,
            owned_end: end,
        });
        self.next = end;
        if self.batch.len() == self.config.batch_chunks {
            let batch = std::mem::replace(
                &mut self.batch,
                Vec::with_capacity(self.config.batch_chunks),
            );
            self.send
                .send(batch)
                .map_err(|_| anyhow::anyhow!("minimizer consumer closed"))?;
        }
        Ok(())
    }

    fn finish_record(&mut self) -> Result<()> {
        self.push_job(self.length, self.length)?;
        self.seq.clear();
        self.base = 0;
        self.next = 0;
        self.length = 0;
        self.run_length = 0;
        Ok(())
    }

    fn base(&mut self, base: u8) -> Result<()> {
        let base = base.to_ascii_uppercase();
        self.seq.push(base);
        self.length += 1;
        self.run_length = if is_dna(base) { self.run_length + 1 } else { 0 };
        self.stats.eligible_kmer_positions += u64::from(self.run_length >= self.config.k);
        let right = self.config.k + self.config.w - 2;
        if self.length - self.next >= self.config.chunk_bases + right {
            let end = self.next + self.config.chunk_bases;
            self.push_job(end, end + right)?;
            // Every window that can select an owned start is present: w-1
            // bases before it, k+w-2 after the owned interval. Ownership, not
            // chunk-local position, prevents duplicates at chunk boundaries.
            let keep = end.saturating_sub(self.config.w - 1);
            self.seq.copy_within(keep - self.base.., 0);
            self.seq.truncate(self.length - keep);
            self.base = keep;
        }
        Ok(())
    }
}

fn produce(mut reader: impl Read, config: Config, send: SyncSender<Vec<Job>>) -> Result<Stats> {
    let mut producer = Producer {
        config,
        send,
        batch: Vec::with_capacity(config.batch_chunks),
        seq: Vec::with_capacity(config.chunk_bases + config.k + 2 * config.w),
        base: 0,
        next: 0,
        length: 0,
        run_length: 0,
        stats: Stats::default(),
    };
    let mut block = [0u8; 65536];
    let (mut line_start, mut header) = (true, false);
    loop {
        let n = reader.read(&mut block)?;
        if n == 0 {
            break;
        }
        for &base in &block[..n] {
            if header {
                if base == b'\n' {
                    header = false;
                    line_start = true;
                }
                continue;
            }
            if line_start && base == b'>' {
                producer.finish_record()?;
                producer.stats.records += 1;
                header = true;
                line_start = false;
                continue;
            }
            line_start = base == b'\n';
            if base.is_ascii_whitespace() {
                continue;
            }
            ensure!(
                producer.stats.records > 0,
                "sequence before FASTA header (FASTA input required)"
            );
            producer.base(base)?;
        }
    }
    producer.finish_record()?;
    if !producer.batch.is_empty() {
        producer
            .send
            .send(producer.batch)
            .map_err(|_| anyhow::anyhow!("minimizer consumer closed"))?;
    }
    Ok(producer.stats)
}

/// Uses the caller's Rayon pool. Parser and writer have dedicated OS threads;
/// bounded channels apply backpressure. Neither full chromosomes nor complete
/// selected-position lists/output files are retained. Input/output errors close
/// channels before joining threads, so cancellation cannot deadlock the parser.
pub fn generate(
    reader: impl Read + Send,
    mut writer: impl Write + Send,
    config: Config,
) -> Result<Stats> {
    config.validate()?;
    let (input_send, input_receive) = mpsc::sync_channel(1);
    let (output_send, output_receive) = mpsc::sync_channel::<Vec<Encoded>>(1);
    let (recycle_send, recycle_receive) = mpsc::channel();
    thread::scope(|scope| -> Result<Stats> {
        let parser = scope.spawn(move || produce(reader, config, input_send));
        let writing = scope.spawn(move || -> Result<u64> {
            let mut selected = 0;
            let mut nonempty = false;
            for output in output_receive {
                for chunk in &output {
                    writer.write_all(&chunk.bytes)?;
                    selected += chunk.selected;
                    nonempty |= !chunk.bytes.is_empty();
                }
                let _ = recycle_send.send(output);
            }
            if !nonempty {
                if let Some(level) = config.gzip_level {
                    GzEncoder::new(&mut writer, Compression::new(level)).finish()?;
                }
            }
            writer.flush()?;
            Ok(selected)
        });
        let compute = (|| -> Result<()> {
            let mut states: Vec<_> = (0..config.batch_chunks)
                .map(|_| Worker::default())
                .collect();
            while let Ok(jobs) = input_receive.recv() {
                let mut output = recycle_receive.try_recv().unwrap_or_default();
                output.resize_with(jobs.len(), Encoded::default);
                jobs.par_iter()
                    .zip(states.par_iter_mut())
                    .zip(output.par_iter_mut())
                    .try_for_each(|((job, state), out)| encode_job(job, state, out, config))?;
                output_send
                    .send(output)
                    .map_err(|_| anyhow::anyhow!("FASTA writer closed"))?;
            }
            Ok(())
        })();
        drop(input_receive);
        drop(output_send);
        let written = writing
            .join()
            .map_err(|_| anyhow::anyhow!("FASTA writer panicked"))?;
        let parsed = parser
            .join()
            .map_err(|_| anyhow::anyhow!("FASTA parser panicked"))?;
        let selected = written?;
        compute?;
        let mut stats = parsed?;
        stats.selected_positions = selected;
        Ok(stats)
    })
}
