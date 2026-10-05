use std::{
    fs,
    path::{Path, PathBuf},
    process::{Command, Output},
    sync::atomic::{AtomicU64, Ordering},
};

static NEXT: AtomicU64 = AtomicU64::new(0);
struct Fixture(PathBuf);
impl Fixture {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "zor-findere-test-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn path(&self, name: &str) -> String {
        self.0.join(name).to_str().unwrap().to_owned()
    }
    fn write(&self, name: &str, value: impl AsRef<[u8]>) -> String {
        let path = self.path(name);
        fs::write(&path, value).unwrap();
        path
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn binaries() -> [&'static str; 2] {
    [
        env!("CARGO_BIN_EXE_bloomindex"),
        env!("CARGO_BIN_EXE_zorindex"),
    ]
}
fn call(bin: &str, args: &[&str]) -> Output {
    Command::new(bin)
        .args(["--threads", "1"])
        .args(args)
        .output()
        .unwrap()
}
fn ok(bin: &str, args: &[&str]) -> String {
    let output = call(bin, args);
    assert!(
        output.status.success(),
        "{bin} {args:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}
fn is_zor(bin: &str) -> bool {
    Path::new(bin).file_name().unwrap() == "zorindex"
}
fn build(bin: &str, fof: &str, output: &str, k: &str, z: Option<&str>, extra: &[&str]) {
    let mut args = vec![
        "build",
        "--fof",
        fof,
        "--output",
        output,
        "--index-mode",
        "findere",
        "--k",
        k,
    ];
    if let Some(z) = z {
        args.extend(["--findere-z", z]);
    }
    if is_zor(bin) {
        args.extend(["--slot-scale", "8", "--fingerprint-bits", "16"]);
    } else {
        args.extend(["--bit-count", "65536"]);
    }
    args.extend(extra);
    ok(bin, &args);
}
fn matches(output: &str, color: &str) -> u64 {
    output
        .lines()
        .find_map(|line| {
            let mut fields = line.split('\t');
            if fields.next() == Some(color) {
                Some(fields.next().unwrap().parse().unwrap())
            } else {
                None
            }
        })
        .unwrap()
}
fn dna(n: usize, seed: u64) -> String {
    (0..n)
        .map(|i| b"ACGT"[(zorindex::splitmix64(i as u64 ^ seed) & 3) as usize] as char)
        .collect()
}
fn rc(seq: &str) -> String {
    seq.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            _ => 'A',
        })
        .collect()
}

#[test]
fn default_custom_zero_and_long_k_persist_and_query() {
    let f = Fixture::new();
    let seq = dna(150, 17);
    let reference = f.write("reference.fa", format!(">r\n{seq}\n"));
    let fof = f.write("inputs.fof", format!("sample\t{reference}\n"));
    let query = f.write(
        "query.fa",
        format!(
            ">forward\n{seq}\n>reverse\n{}\n>short\nACGT\n>ambiguous\nACGTNACGT\n",
            rc(&seq)
        ),
    );
    for bin in binaries() {
        for (k, z, stored_z) in [
            ("31", None, 10),
            ("9", Some("3"), 3),
            ("9", Some("0"), 0),
            ("41", None, 10),
        ] {
            let index = f.path("test.idx");
            build(bin, &fof, &index, k, z, &[]);
            let info = ok(bin, &["info", "--index", &index]);
            assert!(info.contains("index_mode\tfindere\n"));
            assert!(info.contains(&format!("findere_z\t{stored_z}\n")));
            assert!(info.contains(&format!(
                "indexed_k\t{}\n",
                k.parse::<usize>().unwrap() - stored_z
            )));
            let result = ok(bin, &["query", "--index", &index, "--query", &query]);
            let expected = 2 * (seq.len() - k.parse::<usize>().unwrap() + 1) as u64;
            assert_eq!(matches(&result, "sample"), expected, "{result}");
            assert!(
                result.contains(&format!("total_features={expected}\t")),
                "{result}"
            );
        }
    }
}

#[test]
fn append_repack_union_and_compressed_fastq_keep_findere_semantics() {
    let f = Fixture::new();
    let a = dna(300, 71);
    let b = dna(80, 93);
    let pa = f.write("a.fa", format!(">a\n{a}\n"));
    let pb = f.write("b.fa", format!(">b\n{b}\n"));
    let fof = f.write("a.fof", format!("a\t{pa}\n"));
    let append = f.write("b.fof", format!("b\t{pb}\n"));
    let fastq = format!(
        "@b\n{b}\n+\n{}\n@short\nACGTNACGT\n+\nIIIIIIIII\n",
        "I".repeat(b.len())
    );
    let query = f.write(
        "query.fq.zst",
        zstd::bulk::compress(fastq.as_bytes(), 1).unwrap(),
    );
    for bin in binaries() {
        let index = f.path("initial.idx");
        let appended = f.path("appended.idx");
        let extra = if is_zor(bin) {
            vec!["--union-graph", "--stack-compression", "lz4"]
        } else {
            vec![]
        };
        build(bin, &fof, &index, "31", Some("7"), &extra);
        let qa = f.write("qa.fa", format!(">a\n{a}\n"));
        let output = ok(bin, &["query", "--index", &index, "--query", &qa]);
        assert_eq!(matches(&output, "a"), (a.len() - 31 + 1) as u64);
        if is_zor(bin) {
            let ignored = ok(
                bin,
                &["query", "--index", &index, "--query", &qa, "--ignore-union"],
            );
            assert_eq!(matches(&ignored, "a"), matches(&output, "a"));
        }
        ok(
            bin,
            &[
                "append", "--index", &index, "--fof", &append, "--output", &appended,
            ],
        );
        let info = ok(bin, &["info", "--index", &appended]);
        assert!(info.contains("findere_z\t7\n"));
        let output = ok(bin, &["query", "--index", &appended, "--query", &query]);
        assert_eq!(matches(&output, "b"), (b.len() - 31 + 1) as u64);
        if is_zor(bin) {
            let repacked = f.path("repacked.idx");
            ok(
                bin,
                &[
                    "repack",
                    "--index",
                    &appended,
                    "--output",
                    &repacked,
                    "--stack-compression",
                    "lz4",
                ],
            );
            let output = ok(
                bin,
                &[
                    "query",
                    "--index",
                    &repacked,
                    "--query",
                    &query,
                    "--decode-stack",
                ],
            );
            assert_eq!(matches(&output, "b"), (b.len() - 31 + 1) as u64);
            assert!(ok(bin, &["info", "--index", &repacked]).contains("findere_z\t7\n"));
        }
    }
}

#[test]
fn rejects_invalid_parameters_and_supports_z_alias() {
    let f = Fixture::new();
    let seq = dna(90, 1);
    let reference = f.write("r.fa", format!(">r\n{seq}\n"));
    let fof = f.write("r.fof", format!("r\t{reference}\n"));
    let index = f.path("r.idx");
    for bin in binaries() {
        for (k, z) in [
            ("0", "0"),
            ("10", "10"),
            ("9", "10"),
            ("42", "10"),
            ("256", "225"),
        ] {
            let output = call(
                bin,
                &[
                    "build",
                    "--fof",
                    &fof,
                    "--output",
                    &index,
                    "--index-mode",
                    "findere",
                    "--k",
                    k,
                    "--findere-z",
                    z,
                ],
            );
            assert!(!output.status.success());
            assert!(
                !Path::new(&index).exists(),
                "invalid build should not write an index"
            );
        }
        ok(
            bin,
            &[
                "build",
                "--fof",
                &fof,
                "--output",
                &index,
                "--index-mode",
                "findere",
                "--k",
                "31",
                "--z",
                "4",
            ],
        );
        assert!(ok(bin, &["info", "--index", &index]).contains("findere_z\t4\n"));
        fs::remove_file(&index).unwrap();
    }
}

#[test]
fn shorter_features_in_different_datasets_do_not_make_a_kmer_match() {
    let f = Fixture::new();
    let mut rows = String::new();
    for (name, seq) in [("a", "AAATT"), ("b", "AACCC"), ("c", "ACGGG")] {
        let path = f.write(&format!("{name}.fa"), format!(">{name}\n{seq}\n"));
        rows.push_str(&format!("{name}\t{path}\n"));
    }
    let fof = f.write("input.fof", rows);
    let query = f.write("query.fa", ">q\nAAACG\n");
    for bin in binaries() {
        let index = f.path("test.idx");
        build(bin, &fof, &index, "5", Some("2"), &[]);
        let result = ok(bin, &["query", "--index", &index, "--query", &query]);
        assert!(result.contains("total_features=1\t"));
        for color in ["a", "b", "c"] {
            assert_eq!(matches(&result, color), 0);
        }
    }
}

#[test]
fn query_results_and_probe_counts_match_across_thread_counts() {
    let f = Fixture::new();
    let a = dna(500, 17);
    let b = dna(450, 93);
    let pa = f.write("a.fa", format!(">a\n{a}\n"));
    let pb = f.write("b.fa", format!(">b\n{b}\n"));
    let fof = f.write("input.fof", format!("a\t{pa}\nb\t{pb}\n"));
    let long = format!(
        ">long\n{}\n>broken\n{}N{}\n",
        a.repeat(300),
        &b[..17],
        b.repeat(300)
    );
    let short = format!(">short_a\n{a}\n>short_b\n{b}\n>tiny\nACGT\n").repeat(400);
    let plain = f.write("long.fa", long);
    let compressed = f.write(
        "short.fa.zst",
        zstd::bulk::compress(short.as_bytes(), 1).unwrap(),
    );
    for bin in binaries() {
        for z in ["0", "3", "10"] {
            let index = f.path("test.idx");
            let extra = if is_zor(bin) {
                vec!["--union-graph", "--stack-compression", "lz4"]
            } else {
                vec![]
            };
            build(bin, &fof, &index, "31", Some(z), &extra);
            for query in [&plain, &compressed] {
                let mut expected = None;
                for threads in ["1", "2", "4"] {
                    let output = Command::new(bin)
                        .args([
                            "--threads",
                            threads,
                            "query",
                            "--index",
                            &index,
                            "--query",
                            query,
                        ])
                        .output()
                        .unwrap();
                    assert!(
                        output.status.success(),
                        "{}",
                        String::from_utf8_lossy(&output.stderr)
                    );
                    let text = String::from_utf8(output.stdout).unwrap();
                    let normalized: Vec<_> = text
                        .lines()
                        .map(|line| line.split("\telapsed_s=").next().unwrap().to_owned())
                        .collect();
                    assert_eq!(
                        expected.get_or_insert_with(|| normalized.clone()),
                        &normalized,
                        "thread-dependent output: {bin}, z={z}, threads={threads}"
                    );
                }
            }
        }
    }
}
