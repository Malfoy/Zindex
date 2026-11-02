//! Independently decodable Zstd frames with an optional skippable seek table.
//! The decoded bank stream is unchanged; legacy streaming readers skip the table.
use anyhow::{ensure, Context, Result};
use rayon::prelude::*;
use std::{
    fs::File,
    io::{self, BufWriter, Read, Seek, SeekFrom, Write},
    path::Path,
};

const CHUNK: usize = 4 * 1024 * 1024;
const HEADER: usize = 40;
const SKIP: u32 = 0x184D2A50;
const INDEX: &[u8; 8] = b"ZORFRM1\0";
const FOOTER: &[u8; 8] = b"ZORSEEK1";
const MAX_TABLE: usize = 64 * 1024 * 1024;

#[derive(Clone, Debug)]
struct Frame {
    compressed: u64,
    compressed_len: usize,
    raw: usize,
    raw_len: usize,
}
fn checksum(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf29ce484222325, |h, &b| {
        (h ^ b as u64).wrapping_mul(0x100000001b3)
    })
}
fn integer(bytes: &[u8]) -> u64 {
    u64::from_le_bytes(bytes.try_into().unwrap())
}

pub(super) struct Writer {
    output: BufWriter<File>,
    pending: Vec<u8>,
    frames: Vec<Frame>,
    compressed: u64,
    raw: usize,
}
impl Writer {
    pub(super) fn new(path: &Path) -> Result<Self> {
        Ok(Self {
            output: BufWriter::new(File::create(path)?),
            pending: Vec::with_capacity(CHUNK),
            frames: Vec::new(),
            compressed: 0,
            raw: 0,
        })
    }
    fn frame(&mut self) -> io::Result<()> {
        if self.pending.is_empty() {
            return Ok(());
        }
        let mut encoder = zstd::bulk::Compressor::new(crate::DEFAULT_INDEX_ZSTD_LEVEL)?;
        encoder.set_parameter(zstd::zstd_safe::CParameter::ChecksumFlag(true))?;
        let encoded = encoder.compress(&self.pending)?;
        self.output.write_all(&encoded)?;
        self.frames.push(Frame {
            compressed: self.compressed,
            compressed_len: encoded.len(),
            raw: self.raw,
            raw_len: self.pending.len(),
        });
        self.compressed += encoded.len() as u64;
        self.raw += self.pending.len();
        self.pending.clear();
        Ok(())
    }
    pub(super) fn finish(mut self) -> Result<()> {
        self.frame()?;
        ensure!(!self.frames.is_empty(), "empty bank stream");
        let len = self
            .frames
            .len()
            .checked_mul(32)
            .and_then(|n| n.checked_add(40))
            .context("seek table overflow")?;
        ensure!(len <= MAX_TABLE, "bank seek table too large");
        let mut table = Vec::with_capacity(len);
        table.extend_from_slice(INDEX);
        table.extend_from_slice(&(self.frames.len() as u64).to_le_bytes());
        for f in self.frames {
            for n in [
                f.compressed,
                f.compressed_len as u64,
                f.raw as u64,
                f.raw_len as u64,
            ] {
                table.extend_from_slice(&n.to_le_bytes());
            }
        }
        table.extend_from_slice(&checksum(&table).to_le_bytes());
        table.extend_from_slice(FOOTER);
        table.extend_from_slice(&(len as u64).to_le_bytes());
        self.output.write_all(&SKIP.to_le_bytes())?;
        self.output.write_all(&(len as u32).to_le_bytes())?;
        self.output.write_all(&table)?;
        self.output.flush()?;
        Ok(())
    }
}
impl Write for Writer {
    fn write(&mut self, mut bytes: &[u8]) -> io::Result<usize> {
        let total = bytes.len();
        while !bytes.is_empty() {
            let n = bytes.len().min(CHUNK - self.pending.len());
            self.pending.extend_from_slice(&bytes[..n]);
            bytes = &bytes[n..];
            if self.pending.len() == CHUNK {
                self.frame()?;
            }
        }
        Ok(total)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.frame()?;
        self.output.flush()
    }
}

fn table(path: &Path, raw_total: usize) -> Result<Option<Vec<Frame>>> {
    let mut file = File::open(path)?;
    let length = file.metadata()?.len();
    if length < 16 {
        return Ok(None);
    }
    file.seek(SeekFrom::End(-16))?;
    let mut footer = [0; 16];
    file.read_exact(&mut footer)?;
    if &footer[..8] != FOOTER {
        return Ok(None);
    }
    let len = usize::try_from(integer(&footer[8..])).context("seek table length overflow")?;
    ensure!(
        (40..=MAX_TABLE).contains(&len) && (len - 40) % 32 == 0 && len as u64 + 8 <= length,
        "invalid bank seek table length"
    );
    let start = length - len as u64 - 8;
    file.seek(SeekFrom::Start(start))?;
    let mut skip = [0; 8];
    file.read_exact(&mut skip)?;
    ensure!(
        u32::from_le_bytes(skip[..4].try_into().unwrap()) == SKIP
            && u32::from_le_bytes(skip[4..].try_into().unwrap()) as usize == len,
        "invalid skippable bank index"
    );
    let mut bytes = vec![0; len];
    file.read_exact(&mut bytes)?;
    ensure!(
        &bytes[..8] == INDEX
            && &bytes[len - 16..len - 8] == FOOTER
            && integer(&bytes[len - 8..]) == len as u64,
        "invalid bank seek table tag"
    );
    ensure!(
        integer(&bytes[len - 24..len - 16]) == checksum(&bytes[..len - 24]),
        "bank seek table checksum mismatch"
    );
    let count = usize::try_from(integer(&bytes[8..16])).context("frame count overflow")?;
    ensure!(
        count == (len - 40) / 32 && count != 0 && count <= raw_total,
        "invalid bank frame count"
    );
    let mut frames = Vec::with_capacity(count);
    let (mut compressed, mut raw) = (0u64, 0usize);
    for values in bytes[16..len - 24].chunks_exact(32) {
        let f = Frame {
            compressed: integer(&values[..8]),
            compressed_len: usize::try_from(integer(&values[8..16]))?,
            raw: usize::try_from(integer(&values[16..24]))?,
            raw_len: usize::try_from(integer(&values[24..32]))?,
        };
        ensure!(
            f.compressed == compressed
                && f.raw == raw
                && (1..=CHUNK).contains(&f.raw_len)
                && f.compressed_len > 0
                && f.compressed_len <= f.raw_len + f.raw_len / 128 + 1024,
            "invalid bank frame geometry"
        );
        compressed = compressed
            .checked_add(f.compressed_len as u64)
            .context("frame offset overflow")?;
        raw = raw.checked_add(f.raw_len).context("raw offset overflow")?;
        ensure!(
            compressed <= start && raw <= raw_total,
            "bank frame exceeds payload"
        );
        frames.push(f);
    }
    ensure!(
        compressed == start && raw == raw_total && frames[0].raw_len >= HEADER,
        "incomplete bank frame coverage"
    );
    Ok(Some(frames))
}

pub(super) fn zeroed(bytes: usize) -> Vec<u8> {
    if bytes < CHUNK {
        vec![0; bytes]
    } else {
        (0..bytes)
            .into_par_iter()
            .with_min_len(1024 * 1024)
            .map(|_| 0u8)
            .collect()
    }
}
fn encoded(file: &mut File, f: &Frame, bytes: &mut Vec<u8>) -> Result<()> {
    file.seek(SeekFrom::Start(f.compressed))?;
    bytes.resize(f.compressed_len, 0);
    file.read_exact(bytes)?;
    ensure!(
        zstd::zstd_safe::find_frame_compressed_size(&bytes)
            .map_err(|e| anyhow::anyhow!("zstd frame: {e}"))?
            == bytes.len(),
        "frame contains trailing data"
    );
    ensure!(
        zstd::zstd_safe::get_frame_content_size(&bytes)
            .map_err(|e| anyhow::anyhow!("zstd size: {e:?}"))?
            == Some(f.raw_len as u64),
        "frame content size differs from seek table"
    );
    Ok(())
}

pub(super) fn read_parallel(
    path: &Path,
    payload: usize,
) -> Result<Option<([u8; HEADER], Vec<u8>)>> {
    let total = payload.checked_add(HEADER).context("bank size overflow")?;
    let Some(frames) = table(path, total)? else {
        return Ok(None);
    };
    if frames.len() == 1 {
        return Ok(None);
    }
    let mut data = zeroed(payload);
    let mut header = [0; HEADER];
    let first_len = frames[0].raw_len - HEADER;
    let (first, mut tail) = data.split_at_mut(first_len);
    let mut jobs = Vec::with_capacity(frames.len() - 1);
    for frame in &frames[1..] {
        let (part, rest) = tail.split_at_mut(frame.raw_len);
        jobs.push((frame, part));
        tail = rest;
    }
    let (a, b) = rayon::join(
        || -> Result<()> {
            let mut file = File::open(path)?;
            let mut input = Vec::new();
            encoded(&mut file, &frames[0], &mut input)?;
            let bytes = zstd::bulk::decompress(&input, frames[0].raw_len)?;
            ensure!(bytes.len() == frames[0].raw_len, "short first frame");
            header.copy_from_slice(&bytes[..HEADER]);
            first.copy_from_slice(&bytes[HEADER..]);
            Ok(())
        },
        || {
            let chunk = jobs.len().div_ceil(rayon::current_num_threads()).max(1);
            jobs.par_chunks_mut(chunk)
                .try_for_each(|jobs| -> Result<()> {
                    // Keep one descriptor, compressed buffer and decoder per
                    // parallel job. Reallocating them for every frame causes
                    // avoidable page faults and allocator contention.
                    let mut file = File::open(path)?;
                    let mut input = Vec::new();
                    let mut decoder = zstd::bulk::Decompressor::new()?;
                    for (frame, out) in jobs {
                        encoded(&mut file, frame, &mut input)?;
                        ensure!(
                            decoder.decompress_to_buffer(&input, &mut **out)? == out.len(),
                            "short bank frame"
                        );
                    }
                    Ok(())
                })
        },
    );
    a?;
    b?;
    Ok(Some((header, data)))
}

#[cfg(test)]
mod tests {
    use super::*;
    struct Temp(std::path::PathBuf);
    impl Drop for Temp {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }
    fn temp() -> Temp {
        let p = std::env::temp_dir().join(format!(
            "zor-frame-test-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir(&p).unwrap();
        Temp(p)
    }
    #[test]
    fn parallel_frames_match_streaming_and_validate_corruption() {
        let dir = temp();
        let path = dir.0.join("bank.zst");
        let raw: Vec<u8> = (0..CHUNK * 2 + 917)
            .map(|i| crate::splitmix64(i as u64) as u8)
            .collect();
        let mut writer = Writer::new(&path).unwrap();
        for part in raw.chunks(70001) {
            writer.write_all(part).unwrap();
        }
        writer.finish().unwrap();
        let mut decoded = Vec::new();
        crate::open_index_reader(&path)
            .unwrap()
            .read_to_end(&mut decoded)
            .unwrap();
        assert_eq!(decoded, raw);
        let (header, data) = read_parallel(&path, raw.len() - HEADER).unwrap().unwrap();
        assert_eq!(header.as_slice(), &raw[..HEADER]);
        assert_eq!(data, &raw[HEADER..]);
        assert!(read_parallel(&path, raw.len() - HEADER - 1).is_err());
        let original = std::fs::read(&path).unwrap();
        let mut bad = original.clone();
        let n = bad.len();
        bad[n - 25] ^= 1;
        std::fs::write(&path, &bad).unwrap();
        assert!(read_parallel(&path, raw.len() - HEADER).is_err());
        let mut bad = original.clone();
        bad[100] ^= 0x80;
        std::fs::write(&path, &bad).unwrap();
        assert!(read_parallel(&path, raw.len() - HEADER).is_err());
        std::fs::write(&path, &original[..original.len() - 19]).unwrap();
        let mut reader = crate::open_index_reader(&path).unwrap();
        assert!(reader.read_to_end(&mut Vec::new()).is_err());
    }
    #[test]
    fn short_frames_flush_and_legacy_streams_remain_readable() {
        let dir = temp();
        let path = dir.0.join("bank.zst");
        let raw: Vec<u8> = (0..10000).map(|i| (i % 251) as u8).collect();
        let mut writer = Writer::new(&path).unwrap();
        writer.write_all(&raw[..317]).unwrap();
        writer.flush().unwrap();
        writer.write_all(&raw[317..]).unwrap();
        writer.finish().unwrap();
        let (h, d) = read_parallel(&path, raw.len() - HEADER).unwrap().unwrap();
        assert_eq!([h.as_slice(), &d].concat(), raw);
        let mut writer = crate::create_index_writer(&path).unwrap();
        writer.write_all(&raw).unwrap();
        crate::finish_index_writer(writer).unwrap();
        assert!(read_parallel(&path, raw.len() - HEADER).unwrap().is_none());
    }
}
