//! Bounded sorting/merging of construction-only abandoned-key records.
//! No feature keys are externally sorted; records remain 32-bit hashes with
//! construction-group IDs and origin masks, never expanded per dataset on disk.
use anyhow::{ensure, Result};
use std::{
    cmp::Reverse,
    collections::BinaryHeap,
    fs::{self, File, OpenOptions},
    io::{BufReader, BufWriter, Read, Write},
    path::{Path, PathBuf},
};

pub(super) type Record = (u32, u32, u64);
const SORT_RECORDS: usize = 65_536; // 1 MiB of records per active construction job.
const GROUP_FAN_IN: usize = 8;
const GLOBAL_FAN_IN: usize = 64;
const BUFFER_BYTES: usize = 16 * 1024;

fn read_record(reader: &mut impl Read) -> Result<Option<Record>> {
    let mut bytes = [0u8; 16];
    if reader.read(&mut bytes[..1])? == 0 {
        return Ok(None);
    }
    reader.read_exact(&mut bytes[1..])?;
    Ok(Some((
        u32::from_le_bytes(bytes[..4].try_into().unwrap()),
        u32::from_le_bytes(bytes[4..8].try_into().unwrap()),
        u64::from_le_bytes(bytes[8..].try_into().unwrap()),
    )))
}

fn write_record(writer: &mut impl Write, (hash, group, mask): Record) -> Result<()> {
    let mut bytes = [0u8; 16];
    bytes[..4].copy_from_slice(&hash.to_le_bytes());
    bytes[4..8].copy_from_slice(&group.to_le_bytes());
    bytes[8..].copy_from_slice(&mask.to_le_bytes());
    writer.write_all(&bytes)?;
    Ok(())
}

fn merge(paths: &[PathBuf], mut emit: impl FnMut(Record) -> Result<()>) -> Result<()> {
    let mut readers = paths
        .iter()
        .map(|path| {
            ensure!(
                fs::metadata(path)?.len() % 16 == 0,
                "truncated abandoned-key run: {}",
                path.display()
            );
            Ok(BufReader::with_capacity(BUFFER_BYTES, File::open(path)?))
        })
        .collect::<Result<Vec<_>>>()?;
    let mut heap = BinaryHeap::with_capacity(paths.len());
    for (id, reader) in readers.iter_mut().enumerate() {
        if let Some(record) = read_record(reader)? {
            heap.push(Reverse((record, id)));
        }
    }
    let mut previous = None;
    while let Some(Reverse((record, id))) = heap.pop() {
        ensure!(
            previous.is_none_or(|old| old <= record),
            "unsorted abandoned-key run"
        );
        emit(record)?;
        previous = Some(record);
        if let Some(next) = read_record(&mut readers[id])? {
            ensure!(next >= record, "unsorted abandoned-key run");
            heap.push(Reverse((next, id)));
        }
    }
    Ok(())
}

fn writer(path: &Path) -> Result<BufWriter<File>> {
    Ok(BufWriter::with_capacity(
        BUFFER_BYTES,
        OpenOptions::new().write(true).create_new(true).open(path)?,
    ))
}

fn reduce(mut paths: Vec<PathBuf>, prefix: &Path, fan_in: usize) -> Result<Vec<PathBuf>> {
    ensure!(fan_in >= 2, "invalid abandoned-key merge fan-in");
    let mut pass = 0;
    while paths.len() > fan_in {
        let mut next = Vec::with_capacity(paths.len().div_ceil(fan_in));
        for (id, chunk) in paths.chunks(fan_in).enumerate() {
            if chunk.len() == 1 {
                next.push(chunk[0].clone());
                continue;
            }
            let path = prefix.with_extension(format!("merge-{pass}-{id}"));
            let mut out = writer(&path)?;
            merge(chunk, |record| write_record(&mut out, record))?;
            out.flush()?;
            drop(out);
            for input in chunk {
                fs::remove_file(input)?;
            }
            next.push(path);
        }
        paths = next;
        pass += 1;
    }
    Ok(paths)
}

pub(super) fn write_sorted(path: &Path, records: impl Iterator<Item = Record>) -> Result<()> {
    write_sorted_limited(path, records, SORT_RECORDS, GROUP_FAN_IN)
}

fn write_sorted_limited(
    path: &Path,
    mut records: impl Iterator<Item = Record>,
    limit: usize,
    fan_in: usize,
) -> Result<()> {
    ensure!(limit > 0, "invalid abandoned-key sort buffer");
    let mut buffer = Vec::with_capacity(limit);
    let mut runs = Vec::new();
    loop {
        buffer.clear();
        buffer.extend(records.by_ref().take(limit));
        if buffer.is_empty() {
            break;
        }
        buffer.sort_unstable();
        let run = path.with_extension(format!("sort-{}", runs.len()));
        let mut out = writer(&run)?;
        for &record in &buffer {
            write_record(&mut out, record)?;
        }
        out.flush()?;
        runs.push(run);
    }
    drop(buffer);
    let runs = reduce(runs, path, fan_in)?;
    if runs.len() == 1 {
        fs::rename(&runs[0], path)?;
        return Ok(());
    }
    let mut out = writer(path)?;
    merge(&runs, |record| write_record(&mut out, record))?;
    out.flush()?;
    drop(out);
    for run in runs {
        fs::remove_file(run)?;
    }
    Ok(())
}

/// Consumes only the supplied temporary files. Memory/file-descriptor use of
/// record buffers and the merge heap is bounded independently of record count.
pub(super) fn visit_sorted(
    paths: Vec<PathBuf>,
    directory: &Path,
    emit: impl FnMut(Record) -> Result<()>,
) -> Result<()> {
    visit_sorted_limited(paths, directory, GLOBAL_FAN_IN, emit)
}

fn visit_sorted_limited(
    paths: Vec<PathBuf>,
    directory: &Path,
    fan_in: usize,
    emit: impl FnMut(Record) -> Result<()>,
) -> Result<()> {
    let paths = reduce(paths, &directory.join("global"), fan_in)?;
    merge(&paths, emit)?;
    for path in paths {
        fs::remove_file(path)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};
    static NEXT: AtomicU64 = AtomicU64::new(0);

    #[test]
    fn record_io_handles_short_reads_and_reports_each_truncation_boundary() {
        struct ShortReader<'a>(&'a [u8]);
        impl Read for ShortReader<'_> {
            fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
                let len = out.len().min(3).min(self.0.len());
                out[..len].copy_from_slice(&self.0[..len]);
                self.0 = &self.0[len..];
                Ok(len)
            }
        }
        for record in [
            (0, 0, 0),
            (u32::MAX, u32::MAX, u64::MAX),
            (17, 31, 1u64 << 63),
        ] {
            let mut bytes = Vec::new();
            write_record(&mut bytes, record).unwrap();
            assert_eq!(read_record(&mut ShortReader(&bytes)).unwrap(), Some(record));
            assert_eq!(read_record(&mut ShortReader(&[])).unwrap(), None);
            for end in 1..16 {
                assert!(read_record(&mut ShortReader(&bytes[..end])).is_err());
            }
        }
    }

    #[test]
    fn record_io_propagates_injected_read_and_storage_full_errors() {
        struct Broken;
        impl Read for Broken {
            fn read(&mut self, _: &mut [u8]) -> std::io::Result<usize> {
                Err(std::io::Error::other("injected read error"))
            }
        }
        struct FullDisk(usize);
        impl Write for FullDisk {
            fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
                if self.0 == 0 {
                    return Err(std::io::Error::other("injected storage full"));
                }
                let n = self.0.min(bytes.len()).min(3);
                self.0 -= n;
                Ok(n)
            }
            fn flush(&mut self) -> std::io::Result<()> {
                Ok(())
            }
        }
        assert!(read_record(&mut Broken)
            .unwrap_err()
            .to_string()
            .contains("injected read error"));
        for available in 0..16 {
            assert!(write_record(&mut FullDisk(available), (1, 2, 3))
                .unwrap_err()
                .to_string()
                .contains("injected storage full"));
        }
        write_record(&mut FullDisk(16), (1, 2, 3)).unwrap();
    }

    #[test]
    fn failed_emission_or_unsorted_run_does_not_delete_input_runs() {
        let directory = std::env::temp_dir().join(format!(
            "zor-abandoned-failure-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&directory).unwrap();
        let a = directory.join("a");
        let b = directory.join("b");
        write_sorted(&a, [(1, 2, 3), (3, 4, 5)].into_iter()).unwrap();
        write_sorted(&b, [(2, 3, 4)].into_iter()).unwrap();
        let before = [fs::read(&a).unwrap(), fs::read(&b).unwrap()];
        assert!(
            visit_sorted(vec![a.clone(), b.clone()], &directory, |_| anyhow::bail!(
                "injected emitter error"
            ))
            .unwrap_err()
            .to_string()
            .contains("injected emitter error")
        );
        assert_eq!(fs::read(&a).unwrap(), before[0]);
        assert_eq!(fs::read(&b).unwrap(), before[1]);
        let mut bytes = Vec::new();
        write_record(&mut bytes, (9, 0, 1)).unwrap();
        write_record(&mut bytes, (1, 0, 1)).unwrap();
        fs::write(&a, &bytes).unwrap();
        assert!(visit_sorted(vec![a.clone()], &directory, |_| Ok(())).is_err());
        assert_eq!(fs::read(&a).unwrap(), bytes);
        fs::remove_file(a).unwrap();
        fs::remove_file(b).unwrap();
        fs::remove_dir(directory).unwrap();
    }

    #[test]
    fn tiny_buffers_force_multiple_sort_and_merge_passes_without_losing_records() {
        let directory = std::env::temp_dir().join(format!(
            "zor-abandoned-runs-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&directory).unwrap();
        let mut expected = Vec::new();
        let mut paths = Vec::new();
        for group in 0..17 {
            let records: Vec<_> = (0..53)
                .map(|i| ((i * 17 + group) % 23, group, 1u64 << (i % 7)))
                .rev()
                .collect();
            expected.extend_from_slice(&records);
            let path = directory.join(format!("group-{group}.records"));
            write_sorted_limited(&path, records.into_iter(), 3, 2).unwrap();
            paths.push(path);
        }
        expected.sort_unstable();
        let mut actual = Vec::new();
        visit_sorted_limited(paths, &directory, 2, |record| {
            actual.push(record);
            Ok(())
        })
        .unwrap();
        assert_eq!(actual, expected);
        assert_eq!(fs::read_dir(&directory).unwrap().count(), 0);
        fs::remove_dir(&directory).unwrap();
    }

    #[test]
    fn empty_and_truncated_runs_are_handled_without_silent_data_loss() {
        let directory = std::env::temp_dir().join(format!(
            "zor-abandoned-runs-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&directory).unwrap();
        let path = directory.join("empty.records");
        write_sorted(&path, std::iter::empty()).unwrap();
        visit_sorted(vec![path.clone()], &directory, |_| panic!("empty run")).unwrap();
        fs::write(&path, [0; 15]).unwrap();
        assert!(visit_sorted(vec![path.clone()], &directory, |_| Ok(())).is_err());
        fs::remove_file(path).unwrap();
        fs::remove_dir(&directory).unwrap();
    }
}
