use minimizer_fasta::{generate, Config, MinimizerScheme, Stats};
use packed_seq::{PackedSeqVec, SeqVec};
use std::{
    io::{self, Read, Write},
    sync::atomic::{AtomicUsize, Ordering},
};

fn config(k: usize, w: usize, chunk_bases: usize, threads: usize) -> Config {
    Config {
        k,
        w,
        minimizer_scheme: minimizer_fasta::MinimizerScheme::Simd,
        seed: 69,
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

fn reference_ocm(records: &[Vec<u8>], k: usize, w: usize, seed: u64) -> (Vec<u8>, Stats) {
    let mut output = Vec::new();
    let mut stats = Stats {
        records: records.len() as u64,
        ..Stats::default()
    };
    for (id, record) in records.iter().enumerate() {
        let mut offset = 0;
        for run in record.split(|b| !b"ACGT".contains(b)) {
            stats.eligible_kmer_positions += run.len().saturating_sub(k - 1) as u64;
            let mut positions = std::collections::BTreeSet::new();
            ocm_minimizers::scan_positions(run, k + w - 1, k, seed, &mut |p, _| {
                positions.insert(p);
            });
            for p in positions {
                let word = &run[p..p + k];
                let reverse: Vec<_> = word.iter().rev().map(|&b| complement(b)).collect();
                writeln!(output, ">r{}_p{}", id + 1, offset + p).unwrap();
                output.extend_from_slice(word.min(reverse.as_slice()));
                output.push(b'\n');
                stats.selected_positions += 1;
            }
            offset += run.len() + 1;
        }
    }
    (output, stats)
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
    let actual = run(input.as_bytes(), config(31, 19, 59, 4), 4);
    let expected = reference(&records, 31, 19);
    assert_eq!(actual.1, expected.1);
    assert_eq!(
        String::from_utf8(actual.0).unwrap(),
        String::from_utf8(expected.0).unwrap()
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

#[test]
fn ocm_chunks_match_whole_record_positions_across_parameters_ties_and_strands() {
    let random = dna(2307, 71);
    let reverse: Vec<_> = random.iter().rev().map(|&b| complement(b)).collect();
    let records = vec![
        random,
        reverse,
        b"ACGT".repeat(221),
        vec![b'A'; 501],
        [dna(101, 19), b"NNRYN".to_vec(), dna(291, 21)].concat(),
        vec![],
        dna(18, 20),
    ];
    for (k, w) in [(31, 19), (21, 9), (12, 8), (7, 1), (1, 3)] {
        for seed in [0, 69, u64::MAX] {
            let expected = reference_ocm(&records, k, w, seed);
            for chunk in [1, 7, 31, 127, 262144] {
                let conf = Config {
                    minimizer_scheme: MinimizerScheme::OpenClosedMod,
                    seed,
                    ..config(k, w, chunk, 2)
                };
                assert_eq!(
                    run(&fasta(&records), conf, 2),
                    expected,
                    "OCM k={k}, w={w}, seed={seed}, chunk={chunk}"
                );
            }
        }
    }
}

#[test]
fn ocm_output_and_gzip_bytes_are_thread_deterministic() {
    let records = vec![dna(40_000, 731), vec![], b"ACGT".repeat(700)];
    let expected = reference_ocm(&records, 31, 19, 69);
    let input = fasta(&records);
    let mut conf = Config {
        minimizer_scheme: MinimizerScheme::OpenClosedMod,
        gzip_level: Some(1),
        ..config(31, 19, 701, 1)
    };
    let one = run(&input, conf, 1);
    conf.batch_chunks = 4;
    let many = run(&input, conf, 4);
    assert_eq!(one, many);
    let mut plain = Vec::new();
    flate2::read::MultiGzDecoder::new(one.0.as_slice())
        .read_to_end(&mut plain)
        .unwrap();
    assert_eq!((plain, one.1), expected);
}

#[test]
fn both_schemes_read_all_compressions_by_magic_in_library_and_cli() {
    use deko::{
        write::{AnyEncoder, Compression},
        Format,
    };
    let fixture = Fixture::new();
    let records = vec![dna(1003, 42), vec![], b"ACGT".repeat(109)];
    let input = fasta(&records);
    for scheme in [MinimizerScheme::Simd, MinimizerScheme::OpenClosedMod] {
        let scheme_name = if scheme == MinimizerScheme::Simd {
            "simd"
        } else {
            "open-closed-mod"
        };
        let expected = if scheme == MinimizerScheme::Simd {
            reference(&records, 21, 9)
        } else {
            reference_ocm(&records, 21, 9, 17)
        };
        for (id, format) in [
            Format::Verbatim,
            Format::Gz,
            Format::Xz,
            Format::Zstd,
            Format::Bz,
        ]
        .into_iter()
        .enumerate()
        {
            let mut encoder = AnyEncoder::new(Vec::new(), format, Compression::Default).unwrap();
            encoder.write_all(&input).unwrap();
            let encoded = encoder.finish().unwrap();
            let conf = Config {
                minimizer_scheme: scheme,
                seed: 17,
                ..config(21, 9, 73, 2)
            };
            assert_eq!(
                run(&encoded, conf, 2),
                expected,
                "scheme={scheme_name}, format={format:?}"
            );
            // No extension: decompression must use magic, never the filename.
            let path = fixture.write(&format!("input-{scheme_name}-{id}.bin"), &encoded);
            let output = fixture.0.join(format!("output-{scheme_name}-{id}.fa"));
            let result = std::process::Command::new(env!("CARGO_BIN_EXE_minimizer-fasta"))
                .args([
                    "--minimizer-scheme",
                    scheme_name,
                    "--seed",
                    "17",
                    "--k",
                    "21",
                    "--window",
                    "9",
                    "--threads",
                    "2",
                    "--chunk-bases",
                    "73",
                ])
                .arg(&path)
                .arg(&output)
                .output()
                .unwrap();
            assert!(
                result.status.success(),
                "{}",
                String::from_utf8_lossy(&result.stderr)
            );
            assert_eq!(std::fs::read(output).unwrap(), expected.0);
            let summary: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
            assert_eq!(summary["minimizer_scheme"], scheme_name);
            assert_eq!(summary["seed"], 17);
            assert_eq!(summary["selected_positions"], expected.1.selected_positions);
        }
    }
}

#[test]
fn ocm_collection_and_compressed_stdin_match_single_file_selection() {
    use std::process::{Command, Stdio};
    let fixture = Fixture::new();
    let records = vec![dna(1003, 77), vec![], dna(70, 81)];
    let input = fasta(&records);
    let expected = reference_ocm(&records, 31, 19, 69);
    let mut encoder = deko::write::AnyEncoder::new(
        Vec::new(),
        deko::Format::Zstd,
        deko::write::Compression::Default,
    )
    .unwrap();
    encoder.write_all(&input).unwrap();
    let compressed = encoder.finish().unwrap();
    let a = fixture.write("input.data", &compressed);
    let b = fixture.write("input.fa", &input);
    let fof = fixture.write(
        "inputs.fof",
        format!("a\t{}\nb\t{}\n", a.display(), b.display()).as_bytes(),
    );
    let directory = fixture.0.join("ocm-collection");
    let result = Command::new(env!("CARGO_BIN_EXE_minimizer-fasta"))
        .args([
            "--minimizer-scheme",
            "ocm",
            "--threads",
            "2",
            "--chunk-bases",
            "71",
            "--fof",
        ])
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
    for name in ["00000000.minimizers.fa", "00000001.minimizers.fa"] {
        assert_eq!(std::fs::read(directory.join(name)).unwrap(), expected.0);
    }
    let mut child = Command::new(env!("CARGO_BIN_EXE_minimizer-fasta"))
        .args([
            "--minimizer-scheme",
            "open-closed-mod",
            "--threads",
            "2",
            "-",
            "-",
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(&compressed).unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, expected.0);
    let summary: serde_json::Value = serde_json::from_slice(&output.stderr).unwrap();
    assert_eq!(summary["minimizer_scheme"], "open-closed-mod");
}

#[test]
fn compressed_read_failures_and_ocm_invalid_parameters_are_errors() {
    use deko::{
        write::{AnyEncoder, Compression},
        Format,
    };
    let data = fasta(&[dna(1000, 54)]);
    for format in [Format::Gz, Format::Xz, Format::Zstd] {
        let mut encoder = AnyEncoder::new(Vec::new(), format, Compression::Default).unwrap();
        encoder.write_all(&data).unwrap();
        let bytes = encoder.finish().unwrap();
        for end in [bytes.len() / 2, bytes.len() - 1] {
            for scheme in [MinimizerScheme::Simd, MinimizerScheme::OpenClosedMod] {
                let conf = Config {
                    minimizer_scheme: scheme,
                    ..config(31, 19, 73, 2)
                };
                assert!(
                    generate(&bytes[..end], Vec::new(), conf).is_err(),
                    "accepted truncated {format:?} {end}"
                );
            }
        }
    }
    let invalid = Config {
        minimizer_scheme: MinimizerScheme::OpenClosedMod,
        ..config(32, 18, 73, 2)
    };
    assert!(invalid
        .validate()
        .unwrap_err()
        .to_string()
        .contains("1..=31"));
    let valid = Config {
        minimizer_scheme: MinimizerScheme::Simd,
        ..invalid
    };
    assert!(valid.validate().is_ok());
    let records = vec![b"ACGT".repeat(500), vec![], dna(31, 54)];
    let conf = Config {
        minimizer_scheme: MinimizerScheme::OpenClosedMod,
        ..config(31, 19, 73, 2)
    };
    assert!(generate(fasta(&records).as_slice(), BadWriter, conf).is_err());
}

#[test]
fn helicase_headers_empty_records_whitespace_and_reader_boundaries_keep_coordinates() {
    let sequence = dna(3001, 419);
    let mut input = Vec::new();
    input.extend_from_slice(b">");
    input.extend_from_slice(&vec![b'h'; 150_000]);
    input.extend_from_slice(b"\r\n");
    for chunk in sequence.chunks(37) {
        for &base in chunk {
            input.push(base.to_ascii_lowercase());
            input.extend_from_slice(b" \t");
        }
        input.extend_from_slice(b"\r\n");
    }
    for _ in 0..130 {
        input.extend_from_slice(b">empty >description\r\n");
    }
    input.extend_from_slice(b">last\nacgtNNRY");
    input.extend_from_slice(&sequence);
    // Deliberately no final newline.
    let mut records = vec![sequence.clone()];
    records.extend(std::iter::repeat_n(Vec::new(), 130));
    records.push([b"ACGTNNRY".to_vec(), sequence].concat());
    for scheme in [MinimizerScheme::Simd, MinimizerScheme::OpenClosedMod] {
        let expected = if scheme == MinimizerScheme::Simd {
            reference(&records, 31, 19)
        } else {
            reference_ocm(&records, 31, 19, 69)
        };
        let conf = Config {
            minimizer_scheme: scheme,
            ..config(31, 19, 127, 2)
        };
        assert_eq!(run(&input, conf, 2), expected);
    }
}

#[test]
fn ocm_cli_rejects_invalid_scheme_parameters_before_creating_output() {
    let fixture = Fixture::new();
    let input = fixture.write("input.fa", &fasta(&[dna(80, 12)]));
    for (id, args) in [
        vec!["--minimizer-scheme", "unknown"],
        vec![
            "--minimizer-scheme",
            "open-closed-mod",
            "--k",
            "32",
            "--window",
            "18",
        ],
        vec!["--minimizer-scheme", "ocm", "--window", "18"],
        vec!["--minimizer-scheme", "ocm", "--seed", "-1"],
    ]
    .into_iter()
    .enumerate()
    {
        let output = fixture.0.join(format!("invalid-{id}.fa"));
        let result = std::process::Command::new(env!("CARGO_BIN_EXE_minimizer-fasta"))
            .arg(&input)
            .arg(&output)
            .args(args)
            .output()
            .unwrap();
        assert!(!result.status.success());
        assert!(!output.exists());
    }
}
