//! Small correctness fixtures, not timing benchmarks. Expected k-mer membership
//! is computed independently using canonical DNA strings, not the index parser.
use std::{
    collections::HashSet,
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
            "zor-options-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn path(&self, name: &str) -> String {
        self.0.join(name).to_str().unwrap().into()
    }
    fn write(&self, name: &str, bytes: impl AsRef<[u8]>) -> String {
        let path = self.path(name);
        fs::write(&path, bytes).unwrap();
        path
    }
    fn no_staging(&self) {
        assert!(fs::read_dir(&self.0).unwrap().all(|e| !e
            .unwrap()
            .file_name()
            .to_string_lossy()
            .starts_with(".zor-pow2-building-")));
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn call(args: &[&str], threads: &str) -> Output {
    Command::new(env!("CARGO_BIN_EXE_zorindex"))
        .args(["--threads", threads])
        .args(args)
        .output()
        .unwrap()
}
fn ok(args: &[&str], threads: &str) -> String {
    let output = call(args, threads);
    assert!(
        output.status.success(),
        "{args:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}
fn build(fof: &str, index: &str, extra: &[&str]) {
    let mut args = vec!["build", "--fof", fof, "--output", index, "--k", "31"];
    args.extend_from_slice(extra);
    ok(&args, "2");
}
fn dna(n: usize, seed: u64) -> String {
    (0..n)
        .map(|i| b"ACGT"[(zorindex::splitmix64(seed ^ i as u64) & 3) as usize] as char)
        .collect()
}
fn rc(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => 'N',
        })
        .collect()
}
fn canonical(sequence: &str) -> String {
    let forward = sequence.to_ascii_uppercase();
    forward.clone().min(rc(&forward))
}
fn exact_kmers(sequence: &str) -> Vec<String> {
    sequence
        .as_bytes()
        .windows(31)
        .filter(|w| w.iter().all(|b| b"ACGTacgt".contains(b)))
        .map(|w| canonical(std::str::from_utf8(w).unwrap()))
        .collect()
}
fn manifest(index: &str) -> serde_json::Value {
    serde_json::from_slice(&fs::read(Path::new(index).join("manifest.json")).unwrap()).unwrap()
}
fn scores(output: &str) -> Vec<u64> {
    output
        .lines()
        .filter(|l| !l.starts_with('#'))
        .skip(1)
        .map(|l| l.split('\t').nth(1).unwrap().parse().unwrap())
        .collect()
}
fn record_rows(output: &str, colors: usize) -> Vec<(u64, Vec<u64>)> {
    assert!(
        output.starts_with("query_id\tvalid_kmers\thits\n"),
        "{output}"
    );
    output
        .lines()
        .skip(1)
        .map(|line| {
            let fields: Vec<_> = line.split('\t').collect();
            let total = fields[1].parse().unwrap();
            let mut counts = vec![0; colors];
            if !fields[2].is_empty() {
                for pair in fields[2].split(',') {
                    let (id, score) = pair.split_once(':').unwrap();
                    counts[id.parse::<usize>().unwrap()] = score.parse().unwrap();
                }
            }
            (total, counts)
        })
        .collect()
}

#[test]
fn canonical_string_oracle_checks_membership_boundaries_and_empty_documents() {
    let f = Fixture::new();
    let shared = dna(210, 37);
    let documents = [
        format!("{shared}NNN{}", dna(90, 11)),
        format!("{}NNN{}", rc(&shared), dna(90, 12)),
        dna(370, 13),
        "ACGTNacgt".into(),
        "NNNNNN".into(),
    ];
    let sets: Vec<HashSet<_>> = documents
        .iter()
        .map(|s| exact_kmers(s).into_iter().collect())
        .collect();
    let mut fof = String::new();
    for (id, document) in documents.iter().enumerate() {
        let path = f.write(&format!("d{id}.fa"), format!(">a\n{document}\n"));
        fof.push_str(&format!("d{id}\t{path}\n"));
    }
    let fof = f.write("inputs.fof", fof);
    // Cross-record windows must never be created; lowercase, ambiguity, repeated
    // keys, reverse complements, and records shorter than k are all deliberate.
    let queries = [
        shared.to_ascii_lowercase(),
        rc(&shared),
        format!("{}NN{}", &shared[..61], &documents[2][..75]),
        dna(150, 999),
        "ACGT".into(),
        "NNNN".into(),
        String::new(),
        format!("{}N{}N{}", &shared[..31], &shared[..31], &shared[..31]),
    ];
    let query = f.write(
        "query.fa",
        queries
            .iter()
            .enumerate()
            .map(|(i, s)| format!(">q{i}\r\n{s}\r\n"))
            .collect::<String>(),
    );
    let expected: Vec<_> = queries
        .iter()
        .map(|s| {
            let keys = exact_kmers(s);
            let counts: Vec<_> = sets
                .iter()
                .map(|set| keys.iter().filter(|k| set.contains(*k)).count() as u64)
                .collect();
            (keys.len() as u64, counts)
        })
        .collect();
    for (case, bits) in ["1", "2", "4", "8", "16"].into_iter().enumerate() {
        let index = f.path(&format!("oracle-{bits}"));
        let planner = ["streaming", "all-pairs", "streaming-unpruned"][case % 3];
        build(
            &fof,
            &index,
            &[
                "--fingerprint-bits",
                bits,
                "--merge-planner",
                planner,
                "--min-subindex-documents",
                "1",
            ],
        );
        let mut previous = None;
        for (threads, batch) in [("1", "1"), ("2", "3")] {
            let output = ok(
                &[
                    "query",
                    "--index",
                    &index,
                    "--query",
                    &query,
                    "--records",
                    "--record-batch-size",
                    batch,
                ],
                threads,
            );
            if let Some(old) = &previous {
                assert_eq!(old, &output);
            }
            let actual = record_rows(&output, documents.len());
            assert_eq!(actual.len(), expected.len());
            for ((total, counts), (want_total, exact)) in actual.iter().zip(&expected) {
                assert_eq!(total, want_total);
                for (&got, &present) in counts.iter().zip(exact) {
                    // An approximate index may add false positives, never lose
                    // actual membership when abandoned-key storage is enabled.
                    assert!(
                        got >= present && got <= *total,
                        "fp{bits}: {got} vs exact {present}/{total}"
                    );
                }
            }
            previous = Some(output);
        }
        let aggregate = ok(&["query", "--index", &index, "--query", &query], "2");
        let actual = scores(&aggregate);
        for (id, &got) in actual.iter().enumerate() {
            let exact: u64 = expected.iter().map(|(_, c)| c[id]).sum();
            assert!(got >= exact);
        }
    }
    f.no_staging();
}

#[test]
fn constructed_singleton_false_positive_rates_follow_fingerprint_width() {
    let f = Fixture::new();
    let sequence = dna(1054, 123);
    let reference = f.write("ref.fa", format!(">ref\n{sequence}\n"));
    let fof = f.write("inputs.fof", format!("reference\t{reference}\n"));
    let present: HashSet<_> = exact_kmers(&sequence).into_iter().collect();
    let negative: Vec<_> = (0..8192)
        .map(|i| dna(31, 10_000 + i))
        .filter(|s| !present.contains(&canonical(s)))
        .collect();
    assert_eq!(negative.len(), 8192);
    let query = f.write(
        "absent.fa",
        negative
            .iter()
            .enumerate()
            .map(|(i, s)| format!(">absent{i}\n{s}\n"))
            .collect::<String>(),
    );
    for bits in [1u32, 2, 4, 8, 16] {
        let index = f.path(&format!("fpr-{bits}"));
        build(
            &fof,
            &index,
            &["--fingerprint-bits", &bits.to_string(), "--merge-disabled"],
        );
        let output = ok(
            &[
                "query",
                "--index",
                &index,
                "--query",
                &query,
                "--ignore-abandoned",
            ],
            "2",
        );
        let got = scores(&output)[0] as f64;
        let p = 2f64.powi(-(bits as i32));
        let mean = negative.len() as f64 * p;
        // Fixed inputs, conservative tolerance: not a performance experiment
        // or a fragile assertion of exactly n*p random false positives.
        let tolerance = 8.0 * (mean * (1.0 - p)).sqrt() + 8.0;
        assert!(
            (got - mean).abs() <= tolerance,
            "fp{bits}: observed={got}, expected={mean}, tolerance={tolerance}"
        );
    }
}

#[test]
fn merged_colored_filters_survive_promotion_across_all_widths_and_storage_modes() {
    let f = Fixture::new();
    let mut fof = String::new();
    let mut references = Vec::new();
    for (family, length) in [480usize, 1900, 4200].into_iter().enumerate() {
        let sequence = dna(length, family as u64 + 97);
        let path = f.write(&format!("family{family}.fa"), format!(">ref\n{sequence}\n"));
        for copy in 0..4 {
            fof.push_str(&format!("f{family}-{copy}\t{path}\n"));
        }
        references.push(path);
    }
    let fof = f.write("inputs.fof", fof);
    for bits in ["1", "2", "4", "8", "16"] {
        for (mode, storage) in [
            ("ram", None),
            ("packed-disk", Some("--low-disk=false")),
            ("disk", Some("--disk-backed")),
        ] {
            for (cap, minimum) in [("1", "1"), ("2", "4"), ("3", "128")] {
                let index = f.path(&format!("promotion-{bits}-{mode}-{cap}-{minimum}"));
                let mut args = vec![
                    "--fingerprint-bits",
                    bits,
                    "--max-subindexes",
                    cap,
                    "--min-subindex-documents",
                    minimum,
                ];
                args.extend(storage);
                build(&fof, &index, &args);
                let m = manifest(&index);
                let banks = m["banks"].as_array().unwrap();
                assert!(banks.len() <= cap.parse::<usize>().unwrap());
                let groups: usize = banks
                    .iter()
                    .map(|b| b["groups"].as_array().unwrap().len())
                    .sum();
                if bits.parse::<usize>().unwrap() > 2 {
                    assert!(groups < 12, "expected actual color merges for fp{bits}");
                } else {
                    // With one bit per origin, two singleton 2-bit payloads
                    // already cost the same as their merged 4-bit payload.
                    // Strict bit savings therefore forbid starting a merge.
                    assert_eq!(groups, 12);
                }
                let mut ids = Vec::new();
                for bank in banks {
                    assert!(bank["slots"].as_u64().unwrap().is_power_of_two());
                    let members: Vec<usize> = bank["groups"]
                        .as_array()
                        .unwrap()
                        .iter()
                        .flat_map(|g| g["members"].as_array().unwrap())
                        .map(|id| id.as_u64().unwrap() as usize)
                        .collect();
                    assert!(members.len() >= minimum.parse::<usize>().unwrap().min(12));
                    ids.extend(members);
                }
                ids.sort_unstable();
                assert_eq!(ids, (0..12).collect::<Vec<_>>());
                for (family, reference) in references.iter().enumerate() {
                    let output = ok(&["query", "--index", &index, "--query", reference], "2");
                    let counts = scores(&output);
                    let length = [480usize, 1900, 4200][family];
                    for &got in &counts[family * 4..family * 4 + 4] {
                        assert_eq!(got, (length - 30) as u64);
                    }
                }
                f.no_staging();
            }
        }
    }
}

#[test]
fn failed_builds_clean_staging_and_never_overwrite_existing_output() {
    let f = Fixture::new();
    let good = f.write("ref.fa", format!(">ref\n{}\n", dna(100, 17)));
    let missing = f.path("missing.fa");
    let bad = f.write("corrupt.fa.zst", [0x28, 0xb5, 0x2f, 0xfd, 1, 2, 3]);
    for (id, input) in [&missing, &bad].into_iter().enumerate() {
        let fof = f.write("failed.fof", format!("good\t{good}\nbad\t{input}\n"));
        for mode in [None, Some("--low-disk=false"), Some("--disk-backed")] {
            let index = f.path(&format!("failed-{id}-{mode:?}"));
            let mut args = vec!["build", "--fof", &fof, "--output", &index, "--k", "31"];
            args.extend(mode);
            assert!(!call(&args, "2").status.success());
            assert!(!Path::new(&index).exists());
            f.no_staging();
        }
    }
    let fof = f.write("good.fof", format!("good\t{good}\n"));
    let existing = f.path("existing");
    fs::create_dir(&existing).unwrap();
    fs::write(Path::new(&existing).join("sentinel"), b"must survive").unwrap();
    assert!(!call(
        &["build", "--fof", &fof, "--output", &existing, "--k", "31"],
        "1"
    )
    .status
    .success());
    assert_eq!(
        fs::read(Path::new(&existing).join("sentinel")).unwrap(),
        b"must survive"
    );
    let parent_file = f.write("not-a-directory", "keep");
    let output = format!("{parent_file}/index");
    assert!(!call(
        &["build", "--fof", &fof, "--output", &output, "--k", "31"],
        "1"
    )
    .status
    .success());
    assert_eq!(fs::read(parent_file).unwrap(), b"keep");
    f.no_staging();
}

#[test]
fn truncated_bank_and_abandoned_streams_are_rejected_without_changing_files() {
    let f = Fixture::new();
    let reference = f.write("ref.fa", format!(">ref\n{}\n", dna(1054, 34)));
    let fof = f.write("inputs.fof", format!("a\t{reference}\n"));
    let index = f.path("index");
    build(&fof, &index, &["--fingerprint-bits", "2"]);
    let m = manifest(&index);
    for file in [
        m["banks"][0]["file"].as_str().unwrap(),
        m["abandoned_file"].as_str().unwrap(),
    ] {
        let path = Path::new(&index).join(file);
        let original = fs::read(&path).unwrap();
        for keep in [0, 4, original.len() / 2] {
            let truncated = &original[..keep];
            fs::write(&path, truncated).unwrap();
            assert!(
                !call(&["info", "--index", &index], "1").status.success(),
                "accepted {file} truncated to {keep}"
            );
            assert_eq!(fs::read(&path).unwrap(), truncated);
        }
        fs::write(&path, original).unwrap();
    }
    ok(&["info", "--index", &index], "1");
    let output = ok(&["query", "--index", &index, "--query", &reference], "2");
    assert_eq!(scores(&output), [1024]);
}

#[test]
fn invalid_option_boundaries_fail_before_publishing_an_index() {
    let f = Fixture::new();
    let reference = f.write("ref.fa", format!(">ref\n{}\n", dna(100, 17)));
    let fof = f.write("inputs.fof", format!("a\t{reference}\n"));
    for (id, extra) in [
        vec!["--fingerprint-bits", "0"],
        vec!["--fingerprint-bits", "3"],
        vec!["--fingerprint-bits", "32"],
        vec!["--hashes", "1"],
        vec!["--hashes", "33"],
        vec!["--segment-length", "0"],
        vec!["--segment-length", "3"],
        vec!["--segment-length", "8192"],
        vec!["--merge-sketch-size", "3"],
        vec!["--merge-max-colors", "0"],
        vec!["--fingerprint-bits", "4", "--merge-max-colors", "61"],
        vec!["--merge-min-saving", "1"],
        vec!["--merge-plan-margin", "0.5"],
        vec!["--max-subindexes", "0"],
        vec!["--min-subindex-documents", "0"],
        vec!["--build-memory-gib", "0"],
        vec!["--slot-scale", "NaN"],
        vec!["--slot-scale", "inf"],
        vec!["--tie-scan", "0"],
    ]
    .into_iter()
    .enumerate()
    {
        let index = f.path(&format!("invalid-{id}"));
        let mut args = vec!["build", "--fof", &fof, "--output", &index, "--k", "31"];
        args.extend(extra);
        assert!(
            !call(&args, "1").status.success(),
            "accepted invalid options: {args:?}"
        );
        assert!(!Path::new(&index).exists());
        f.no_staging();
    }
}

#[test]
fn maximal_color_payloads_keep_both_origins_distinct_at_all_mergeable_widths() {
    let f = Fixture::new();
    let a = f.write("a.fa", format!(">a\n{}\n", dna(700, 811)));
    let b = f.write("b.fa", format!(">b\n{}\n", dna(700, 923)));
    for bits in [4usize, 8, 16] {
        let cap = 64 - bits;
        let mut inputs: String = (0..cap).map(|i| format!("a{i:02}\t{a}\n")).collect();
        inputs.push_str(&format!("b0\t{b}\nb1\t{b}\n"));
        let fof = f.write("inputs.fof", inputs);
        for planner in ["streaming", "all-pairs"] {
            let index = f.path(&format!("full-payload-{bits}-{planner}"));
            build(
                &fof,
                &index,
                &[
                    "--fingerprint-bits",
                    &bits.to_string(),
                    "--merge-max-colors",
                    &cap.to_string(),
                    "--merge-planner",
                    planner,
                ],
            );
            let m = manifest(&index);
            assert!(m["banks"]
                .as_array()
                .unwrap()
                .iter()
                .flat_map(|b| b["groups"].as_array().unwrap())
                .any(|g| g["width"] == 64 && g["members"].as_array().unwrap().len() == cap));
            let counts_a = scores(&ok(&["query", "--index", &index, "--query", &a], "2"));
            let counts_b = scores(&ok(&["query", "--index", &index, "--query", &b], "1"));
            assert!(counts_a[..cap].iter().all(|&n| n == 670));
            assert!(counts_b[cap..].iter().all(|&n| n == 670));
            // Allow legitimate random false positives, but catch origin-mask
            // aliasing/broadcast that attributes the whole query to every origin.
            assert!(counts_a[cap..].iter().all(|&n| n < 150));
            assert!(counts_b[..cap].iter().all(|&n| n < 150));
        }
    }
}

#[test]
fn malformed_manifest_fields_fail_cleanly_before_loading_payloads() {
    let f = Fixture::new();
    let a = f.write("a.fa", format!(">a\n{}\n", dna(400, 1)));
    let b = f.write("b.fa", format!(">b\n{}\n", dna(400, 2)));
    let fof = f.write("inputs.fof", format!("a\t{a}\nb\t{b}\n"));
    let index = f.path("index");
    build(&fof, &index, &["--fingerprint-bits", "2"]);
    let original = manifest(&index);
    let path = Path::new(&index).join("manifest.json");
    for field in [
        "version",
        "hashes",
        "fingerprint_bits",
        "k",
        "bank-file",
        "slots",
        "row-width",
        "group-width",
        "origin-id",
        "duplicate-origin",
        "empty-colors",
    ] {
        let mut m = original.clone();
        match field {
            "version" => m["version"] = 999.into(),
            "hashes" => m["hashes"] = 0.into(),
            "fingerprint_bits" => m["fingerprint_bits"] = 3.into(),
            "k" => m["k"] = 0.into(),
            "bank-file" => m["banks"][0]["file"] = "../outside.zst".into(),
            "slots" => m["banks"][0]["slots"] = 3.into(),
            "row-width" => m["banks"][0]["row_bits"] = u64::MAX.into(),
            "group-width" => m["banks"][0]["groups"][0]["width"] = 0.into(),
            "origin-id" => m["banks"][0]["groups"][0]["members"][0] = 99.into(),
            "duplicate-origin" => m["banks"][0]["groups"][1]["members"][0] = 0.into(),
            "empty-colors" => m["colors"] = serde_json::json!([]),
            _ => unreachable!(),
        }
        let bytes = serde_json::to_vec(&m).unwrap();
        fs::write(&path, &bytes).unwrap();
        let output = call(&["info", "--index", &index], "1");
        assert!(!output.status.success(), "accepted invalid {field}");
        assert_ne!(
            output.status.code(),
            Some(101),
            "panicked for {field}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(fs::read(&path).unwrap(), bytes);
    }
    fs::write(path, serde_json::to_vec(&original).unwrap()).unwrap();
    ok(&["info", "--index", &index], "1");
}

#[test]
fn invalid_query_options_leave_index_bytes_unchanged() {
    let f = Fixture::new();
    let reference = f.write("ref.fa", format!(">ref\n{}\n", dna(100, 7)));
    let fof = f.write("inputs.fof", format!("a\t{reference}\n"));
    let index = f.path("index");
    build(&fof, &index, &[]);
    let before: Vec<_> = fs::read_dir(&index)
        .unwrap()
        .map(|e| {
            let path = e.unwrap().path();
            let bytes = fs::read(&path).unwrap();
            (path, bytes)
        })
        .collect();
    for extra in [
        vec!["--min-ratio=NaN"],
        vec!["--min-ratio=inf"],
        vec!["--min-ratio=-0.1"],
        vec!["--min-ratio=1.1"],
        vec!["--records", "--record-batch-size", "0"],
    ] {
        let mut args = vec!["query", "--index", &index, "--query", &reference];
        args.extend(extra);
        let output = call(&args, "1");
        assert!(!output.status.success(), "accepted {args:?}");
        assert_ne!(output.status.code(), Some(101));
    }
    assert!(
        !call(&["query", "--index", &index, "--query", &reference], "0")
            .status
            .success()
    );
    for (path, bytes) in before {
        assert_eq!(fs::read(path).unwrap(), bytes);
    }
}

#[test]
fn duplicate_name_append_fails_without_touching_plain_indexes() {
    let f = Fixture::new();
    let reference = f.write("ref.fa", format!(">ref\n{}\n", dna(100, 7)));
    let fof = f.write("inputs.fof", format!("a\t{reference}\n"));
    for bits in ["1", "2", "4", "8", "16"] {
        let index = f.path(&format!("plain-{bits}"));
        build(&fof, &index, &["--plain", "--fingerprint-bits", bits]);
        let before = fs::read(&index).unwrap();
        let output = call(&["append", "--index", &index, "--fof", &fof], "2");
        assert!(!output.status.success());
        assert!(String::from_utf8_lossy(&output.stderr).contains("duplicate color name"));
        assert_eq!(fs::read(&index).unwrap(), before);
    }
}

#[test]
fn all_empty_feature_collections_and_empty_queries_are_consistent() {
    let f = Fixture::new();
    let a = f.write("empty.fa", ">empty\n\n");
    let b = f.write("short.fa", ">short\nacgtNNNN\n");
    let fof = f.write("inputs.fof", format!("a\t{a}\nb\t{b}\n"));
    for bits in ["1", "2", "4", "8", "16"] {
        for (mode, flag) in [
            ("ram", None),
            ("packed-disk", Some("--low-disk=false")),
            ("disk", Some("--disk-backed")),
        ] {
            let index = f.path(&format!("empty-{bits}-{mode}"));
            let mut args = vec!["--fingerprint-bits", bits];
            args.extend(flag);
            build(&fof, &index, &args);
            let m = manifest(&index);
            assert!(m["colors"]
                .as_array()
                .unwrap()
                .iter()
                .all(|c| c["unique_kmers"] == 0));
            for input in [&a, &b] {
                assert_eq!(
                    scores(&ok(&["query", "--index", &index, "--query", input], "2")),
                    [0, 0]
                );
                let rows = record_rows(
                    &ok(
                        &["query", "--index", &index, "--query", input, "--records"],
                        "2",
                    ),
                    2,
                );
                assert_eq!(rows, [(0, vec![0, 0])]);
            }
            f.no_staging();
        }
    }
    let empty_fof = f.write("empty.fof", "");
    let index = f.path("no-inputs");
    assert!(!call(
        &["build", "--fof", &empty_fof, "--output", &index, "--k", "31"],
        "1"
    )
    .status
    .success());
    assert!(!Path::new(&index).exists());
    f.no_staging();
}
