use minimizer_fasta::{generate, Config, Stats};
use packed_seq::{PackedSeqVec, SeqVec};
use std::{
    io::{self, Read, Write},
    sync::atomic::{AtomicUsize, Ordering},
};

fn config(k: usize, w: usize, chunk_bases: usize, threads: usize) -> Config {
    Config {
        k,
        w,
        chunk_bases,
        batch_chunks: threads,
        gzip_level: None,
    }
}

fn dna(n: usize, mut seed: u64) -> Vec<u8> {
    (0..n)
        .map(|_| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            b"ACGT"[(seed & 3) as usize]
        })
        .collect()
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => unreachable!(),
    }
}

// Direct reference copied from the original benchmark generator's algorithm:
// full-run SIMD selection, comparison sort+dedup, ASCII lexical canonicalization.
fn reference(records: &[Vec<u8>], k: usize, w: usize) -> (Vec<u8>, Stats) {
    let mut out = Vec::new();
    let mut stats = Stats {
        records: records.len() as u64,
        ..Stats::default()
    };
    for (id, record) in records.iter().enumerate() {
        let mut offset = 0;
        for run in record.split(|b| !b"ACGT".contains(b)) {
            stats.eligible_kmer_positions += run.len().saturating_sub(k - 1) as u64;
            if run.len() >= k + w - 1 {
                let packed = PackedSeqVec::from_ascii(run);
                let mut positions =
                    simd_minimizers::canonical_minimizer_positions(packed.as_slice(), k, w);
                positions.sort_unstable();
                positions.dedup();
                for p in positions {
                    let p = p as usize;
                    let word = &run[p..p + k];
                    let reverse: Vec<_> = word.iter().rev().map(|&b| complement(b)).collect();
                    writeln!(out, ">r{}_p{}", id + 1, offset + p).unwrap();
                    out.extend_from_slice(word.min(reverse.as_slice()));
                    out.push(b'\n');
                    stats.selected_positions += 1;
                }
            }
            offset += run.len() + 1;
        }
    }
    (out, stats)
}

fn fasta(records: &[Vec<u8>]) -> Vec<u8> {
    let mut out = Vec::new();
    for (id, record) in records.iter().enumerate() {
        writeln!(out, ">record{id}").unwrap();
        out.extend_from_slice(record);
        out.push(b'\n');
    }
    out
}

fn run(input: &[u8], conf: Config, threads: usize) -> (Vec<u8>, Stats) {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    let mut output = Vec::new();
    let stats = pool.install(|| generate(input, &mut output, conf)).unwrap();
    (output, stats)
}

#[test]
fn every_chunk_boundary_matches_full_run_reference_including_repeated_ties() {
    for seed in 1..=3 {
        let records = vec![
            dna(1007, seed),
            b"ACGT".repeat(131),
            vec![b'A'; 431],
            [
                dna(101, seed + 1),
                b"NNNR".to_vec(),
                dna(279, seed + 2),
                b"NN".to_vec(),
            ]
            .concat(),
            vec![],
            dna(48, seed),
            dna(49, seed),
        ];
        for (k, w) in [(31, 19), (21, 9), (32, 18), (7, 1), (1, 3)] {
            let expected = reference(&records, k, w);
            for chunk in [1, 7, 31, 127, 1024, 262144] {
                assert_eq!(
                    run(&fasta(&records), config(k, w, chunk, 4), 4),
                    expected,
                    "seed={seed}, k={k}, w={w}, chunk={chunk}"
                );
            }
        }
    }
}

#[test]
fn gzip_is_standard_multimember_and_deterministic_across_thread_counts() {
    let records = vec![dna(20000, 63), vec![], dna(300, 42)];
    let expected = reference(&records, 31, 19);
    let input = fasta(&records);
    let mut conf = config(31, 19, 1000, 1);
    conf.gzip_level = Some(1);
    let one = run(&input, conf, 1);
    conf.batch_chunks = 4;
    let four = run(&input, conf, 4);
    assert_eq!(one, four);
    let mut plain = Vec::new();
    flate2::read::MultiGzDecoder::new(one.0.as_slice())
        .read_to_end(&mut plain)
        .unwrap();
    assert_eq!((plain, one.1), expected);
}

#[test]
fn whitespace_case_crlf_and_short_runs_match_original_semantics() {
    let input = format!(
        ">first\r\n{}\r\n>empty\n>third\n{}\n",
        "acgt ".repeat(40),
        "a".repeat(48)
    );
    let records = vec![b"ACGT".repeat(40), vec![], vec![b'A'; 48]];
    assert_eq!(
        run(input.as_bytes(), config(31, 19, 59, 4), 4),
        reference(&records, 31, 19)
    );
}

#[test]
fn single_long_unwrapped_record_does_not_require_record_sized_buffers() {
    let records = vec![dna(300_000, 42)];
    assert_eq!(
        run(&fasta(&records), config(31, 19, 4096, 4), 4),
        reference(&records, 31, 19)
    );
}

#[test]
fn empty_plain_and_gzip_outputs_are_valid() {
    assert_eq!(
        run(b"", config(31, 19, 100, 1), 1),
        (vec![], Stats::default())
    );
    let mut conf = config(31, 19, 100, 1);
    conf.gzip_level = Some(1);
    let (gzip, stats) = run(b">empty\n", conf, 1);
    let mut decoded = Vec::new();
    flate2::read::MultiGzDecoder::new(gzip.as_slice())
        .read_to_end(&mut decoded)
        .unwrap();
    assert!(decoded.is_empty());
    assert_eq!(stats.records, 1);
}

#[test]
fn invalid_parameters_and_non_fasta_input_fail() {
    let conf = config(31, 19, 100, 1);
    for invalid in [
        Config { k: 0, ..conf },
        Config { w: 0, ..conf },
        Config { w: 18, ..conf },
        Config {
            chunk_bases: 0,
            ..conf
        },
        Config {
            batch_chunks: 0,
            ..conf
        },
        Config {
            gzip_level: Some(10),
            ..conf
        },
    ] {
        assert!(invalid.validate().is_err());
    }
    assert!(generate(&b"ACGT\n"[..], Vec::new(), conf).is_err());
    assert!(generate(&b"@read\nACGT\n+\n!!!!\n"[..], Vec::new(), conf).is_err());
}

struct BadReader;
impl Read for BadReader {
    fn read(&mut self, _: &mut [u8]) -> io::Result<usize> {
        Err(io::Error::other("broken input"))
    }
}
struct BadWriter;
impl Write for BadWriter {
    fn write(&mut self, _: &[u8]) -> io::Result<usize> {
        Err(io::Error::other("broken output"))
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

#[test]
fn reader_and_writer_errors_cancel_bounded_pipeline() {
    let conf = config(31, 19, 100, 4);
    assert!(generate(BadReader, Vec::new(), conf).is_err());
    let input = fasta(&[dna(100_000, 42)]);
    assert!(generate(input.as_slice(), BadWriter, conf).is_err());
}

static NEXT: AtomicUsize = AtomicUsize::new(0);
struct Fixture(std::path::PathBuf);
impl Fixture {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "minimizer-fasta-test-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn write(&self, name: &str, bytes: &[u8]) -> std::path::PathBuf {
        let path = self.0.join(name);
        std::fs::write(&path, bytes).unwrap();
        path
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

#[test]
fn cli_gzip_input_output_collection_and_no_overwrite() {
    let fixture = Fixture::new();
    let records = vec![dna(10000, 42), dna(77, 63)];
    let data = fasta(&records);
    let expected = reference(&records, 31, 19);
    let mut gzip = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
    gzip.write_all(&data).unwrap();
    let gzip = gzip.finish().unwrap();
    // Two gzip members, magic-detected even without a .gz extension.
    let input = fixture.write(
        "input.bin",
        &[
            gzip,
            vec![
                31, 139, 8, 0, 0, 0, 0, 0, 4, 255, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            ],
        ]
        .concat(),
    );
    let output = fixture.0.join("out.fa.gz");
    let result = std::process::Command::new(env!("CARGO_BIN_EXE_minimizer-fasta"))
        .args(["--threads", "4", "--chunk-bases", "121"])
        .arg(&input)
        .arg(&output)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let stats: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(stats["selected_positions"], expected.1.selected_positions);
    let mut decoded = Vec::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(&output).unwrap())
        .read_to_end(&mut decoded)
        .unwrap();
    assert_eq!(decoded, expected.0);
    let before = std::fs::read(&output).unwrap();
    assert!(
        !std::process::Command::new(env!("CARGO_BIN_EXE_minimizer-fasta"))
            .arg(&input)
            .arg(&output)
            .output()
            .unwrap()
            .status
            .success()
    );
    assert_eq!(std::fs::read(&output).unwrap(), before);
    let fof = fixture.write(
        "input.fof",
        format!("a\t{}\nb\t{}\n", input.display(), input.display()).as_bytes(),
    );
    let directory = fixture.0.join("collection");
    let result = std::process::Command::new(env!("CARGO_BIN_EXE_minimizer-fasta"))
        .args(["--threads", "4", "--gzip", "--fof"])
        .arg(&fof)
        .arg("--output-dir")
        .arg(&directory)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    for name in ["00000000.minimizers.fa.gz", "00000001.minimizers.fa.gz"] {
        let mut decoded = Vec::new();
        flate2::read::MultiGzDecoder::new(std::fs::File::open(directory.join(name)).unwrap())
            .read_to_end(&mut decoded)
            .unwrap();
        assert_eq!(decoded, expected.0);
    }
    assert_eq!(
        std::fs::read_to_string(directory.join("minimizers.fof"))
            .unwrap()
            .lines()
            .count(),
        2
    );
}
