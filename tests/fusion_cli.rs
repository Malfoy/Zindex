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
            "zor-index-merge-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn path(&self, name: &str) -> String {
        self.0.join(name).to_str().unwrap().to_owned()
    }
    fn write(&self, name: &str, bytes: impl AsRef<[u8]>) -> String {
        let path = self.path(name);
        fs::write(&path, bytes).unwrap();
        path
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn call(args: &[&str]) -> Output {
    let format =
        if args.first() == Some(&"query") && args.iter().any(|a| ["--records", "-r"].contains(a)) {
            vec!["--output-format", "compact"]
        } else {
            vec![]
        };
    Command::new(env!("CARGO_BIN_EXE_zorindex"))
        .args(["-t", "2"])
        .args(args)
        .args(format)
        .output()
        .unwrap()
}
fn ok(args: &[&str]) -> String {
    let result = call(args);
    assert!(
        result.status.success(),
        "{args:?}: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    String::from_utf8(result.stdout).unwrap()
}
fn dna(n: usize, seed: u64) -> String {
    (0..n)
        .map(|i| b"ACGT"[(zorindex::splitmix64(i as u64 ^ seed) & 3) as usize] as char)
        .collect()
}
fn build(fof: &str, output: &str, extra: &[&str]) {
    let mut args = vec![
        "build",
        "-f",
        fof,
        "-o",
        output,
        "-k",
        if extra.windows(2).any(|pair| pair == ["-i", "minimizers"]) {
            "49"
        } else {
            "31"
        },
        "-D",
        "1",
        "-S",
        "0.3",
    ];
    args.extend_from_slice(extra);
    ok(&args);
}
fn manifest(index: &str) -> serde_json::Value {
    serde_json::from_slice(&fs::read(Path::new(index).join("manifest.json")).unwrap()).unwrap()
}
fn write_manifest(index: &str, value: &serde_json::Value) {
    fs::write(
        Path::new(index).join("manifest.json"),
        serde_json::to_vec_pretty(value).unwrap(),
    )
    .unwrap();
}
fn rows(output: &str) -> Vec<String> {
    output
        .lines()
        .filter(|s| {
            !s.starts_with('#') && !s.starts_with("color\t") && !s.starts_with("query_id\t")
        })
        .map(str::to_owned)
        .collect()
}
fn records(output: &str) -> Vec<(String, u64, Vec<(usize, u64)>)> {
    rows(output)
        .into_iter()
        .map(|s| {
            let fields: Vec<_> = s.split('\t').collect();
            (
                fields[0].to_owned(),
                fields[1].parse().unwrap(),
                fields[2]
                    .split(',')
                    .filter(|s| !s.is_empty())
                    .map(|s| {
                        let (id, count) = s.split_once(':').unwrap();
                        (id.parse().unwrap(), count.parse().unwrap())
                    })
                    .collect(),
            )
        })
        .collect()
}

#[test]
fn fusion_preserves_every_query_mode_with_different_planners_and_partition_layouts() {
    let modes: &[&[&str]] = &[
        &[],
        &["-b", "1"],
        &["-b", "16"],
        &["-i", "minimizers", "-m", "31"],
        &["-i", "minimizers", "-m", "31", "-M", "simd"],
        &["-i", "modimizers", "-d", "4"],
        &["-i", "findere"],
        &["-p", "##_#_####_#_#_###_#_#_####_#_##"],
        &["-r", "-l", "31", "-u", "93"],
        &["-i", "minimizers", "-m", "31", "-r"],
        &["-c"],
    ];
    for mode in modes {
        let f = Fixture::new();
        let common = dna(4500, 77);
        let a = f.write("a.fa", format!(">a\n{common}\n"));
        let b = f.write("b.fa", format!(">b\n{}\n", dna(900, 78)));
        let first_fof = f.write("first.fof", format!("shared\t{a}\nleft\t{b}\n"));
        let second_fof = f.write("second.fof", format!("shared\t{a}\nright\t{b}\n"));
        let first = f.path("first");
        let second = f.path("second");
        let merged = f.path("merged");
        let twice = f.path("twice");
        let mut args = mode.to_vec();
        args.extend(["-P", "full"]);
        build(&first_fof, &first, &args);
        args.truncate(mode.len());
        args.extend(["-P", "streaming", "-X", "1KiB"]);
        build(&second_fof, &second, &args);
        let query = f.write(
            "queries.fa",
            format!(
                ">a\n{}\n>b\n{}\n>negative\n{}\n>reverse\n{}\n>short\nACGT\n",
                &common[..600],
                dna(600, 78),
                dna(600, 79),
                common[..600]
                    .bytes()
                    .rev()
                    .map(|b| match b {
                        b'A' => 'T',
                        b'C' => 'G',
                        b'G' => 'C',
                        b'T' => 'A',
                        _ => unreachable!(),
                    })
                    .collect::<String>()
            ),
        );
        let left = manifest(&first);
        let right = manifest(&second);
        let query_args = |index: &str, records: bool, memory: &str| {
            let mut args = vec!["query", "-i", index, "-q", &query, "-M", memory];
            if records {
                args.push("-r");
            }
            ok(&args)
        };
        let first_output = query_args(&first, false, "all");
        let mut table = first_output.lines().filter(|line| !line.starts_with('#'));
        assert_eq!(table.next(), Some("color\tmatches\tratio\tunique_kmers"));
        assert!(table.all(|line| line.split('\t').count() == 4));
        let mut expected = rows(&first_output);
        expected.extend(rows(&query_args(&second, false, "all")));
        fs::remove_file(&a).unwrap();
        fs::remove_file(&b).unwrap();
        ok(&["fusion", "-i", &first, &second, "-o", &merged]);
        let m = manifest(&merged);
        assert_eq!(m["version"], 10);
        assert_eq!(m["colors"].as_array().unwrap().len(), 4);
        assert_eq!(m["source_document_ids"], serde_json::json!([0, 1, 2, 3]));
        // Copying the bank streams must not alter even absent-key false positives.
        for (input, source) in [(&first, &left), (&second, &right)] {
            for bank in source["banks"].as_array().unwrap() {
                let name = bank["file"].as_str().unwrap();
                let component = if input == &first { 1 } else { 2 };
                assert_eq!(
                    fs::read(Path::new(input).join(name)).unwrap(),
                    fs::read(Path::new(&merged).join(format!("component-{component}-{name}")))
                        .unwrap()
                );
            }
        }
        for memory in ["all", "subindex"] {
            assert_eq!(
                rows(&query_args(&merged, false, memory)),
                expected,
                "{mode:?}, {memory}"
            );
            if !mode.contains(&"findere") {
                let mut expected_records = records(&query_args(&first, true, memory));
                for (left, right) in expected_records
                    .iter_mut()
                    .zip(records(&query_args(&second, true, memory)))
                {
                    assert_eq!((left.0.as_str(), left.1), (right.0.as_str(), right.1));
                    left.2
                        .extend(right.2.into_iter().map(|(id, count)| (id + 2, count)));
                }
                assert_eq!(
                    records(&query_args(&merged, true, memory)),
                    expected_records,
                    "{mode:?}, {memory}"
                );
            }
        }
        if mode.is_empty() {
            let serial = f.path("serial");
            let result = Command::new(env!("CARGO_BIN_EXE_zorindex"))
                .args(["-t", "1", "fusion", "-i", &first, &second, "-o", &serial])
                .output()
                .unwrap();
            assert!(
                result.status.success(),
                "{}",
                String::from_utf8_lossy(&result.stderr)
            );
            assert_eq!(manifest(&serial)["banks"], m["banks"]);
            assert_eq!(
                fs::read(Path::new(&serial).join("abandoned.zormap.zst")).unwrap(),
                fs::read(Path::new(&merged).join("abandoned.zormap.zst")).unwrap()
            );
            assert_eq!(rows(&query_args(&serial, false, "subindex")), expected);
            assert!(left["abandoned_hashes"].as_u64().unwrap() > 0);
            assert!(right["banks"]
                .as_array()
                .unwrap()
                .iter()
                .any(|b| b.get("routing").is_some()));
            assert!(
                m["abandoned_hashes"].as_u64().unwrap()
                    < left["abandoned_hashes"].as_u64().unwrap()
                        + right["abandoned_hashes"].as_u64().unwrap()
            );
        }
        ok(&["fusion", "--index", &merged, &first, "--output", &twice]);
        fs::remove_dir_all(&first).unwrap();
        fs::remove_dir_all(&second).unwrap();
        fs::remove_dir_all(&merged).unwrap();
        expected.extend(expected[..2].to_vec());
        for memory in ["all", "subindex"] {
            assert_eq!(rows(&query_args(&twice, false, memory)), expected);
        }
        let info = ok(&["info", "-i", &twice]);
        assert!(info.contains("components\t3\n"));
        assert!(info.contains("merge_planner\tpreserved\n"));
    }
}

#[test]
fn incompatible_or_corrupt_inputs_leave_no_output_and_do_not_modify_sources() {
    let f = Fixture::new();
    let a = f.write("a.fa", format!(">a\n{}\n", dna(1200, 12)));
    let fof = f.write("inputs.fof", format!("a\t{a}\n"));
    let first = f.path("first");
    let second = f.path("second");
    let output = f.path("output");
    build(&fof, &first, &[]);
    build(&fof, &second, &[]);
    let original = manifest(&second);
    for (field, value) in [
        ("k", serde_json::json!(29)),
        ("seed", serde_json::json!(70)),
        ("minimizer_scheme", serde_json::json!("simd")),
        ("modimizer_sampling", serde_json::json!(8)),
        ("findere_z", serde_json::json!(9)),
        ("homopolymer_compression", serde_json::json!(true)),
    ] {
        let mut changed = original.clone();
        changed[field] = value;
        write_manifest(&second, &changed);
        let result = call(&["fusion", "-i", &first, &second, "-o", &output]);
        assert!(!result.status.success(), "{field}");
        assert!(
            String::from_utf8_lossy(&result.stderr)
                .contains(&format!("incompatible indexes: {field} differs")),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert!(!Path::new(&output).exists());
    }
    write_manifest(&second, &original);
    for (flag, value, field) in [
        ("-H", "2", "hashes"),
        ("-b", "4", "fingerprint_bits"),
        ("-L", "32", "segment_length"),
    ] {
        let different = f.path(field);
        build(&fof, &different, &[flag, value]);
        let result = call(&["fusion", "-i", &first, &different, "-o", &output]);
        assert!(!result.status.success());
        assert!(String::from_utf8_lossy(&result.stderr)
            .contains(&format!("incompatible indexes: {field} differs")));
        assert!(!Path::new(&output).exists());
    }
    let before = fs::read(Path::new(&first).join("manifest.json")).unwrap();
    let result = call(&["fusion", "-i", &first, &second, "-o", &first]);
    assert!(!result.status.success());
    assert_eq!(
        before,
        fs::read(Path::new(&first).join("manifest.json")).unwrap()
    );
    let nested = Path::new(&first).join("nested");
    assert!(!call(&[
        "fusion",
        "-i",
        &first,
        &second,
        "-o",
        nested.to_str().unwrap()
    ])
    .status
    .success());
    let bank_path = Path::new(&second).join(original["banks"][0]["file"].as_str().unwrap());
    let original_bank = fs::read(&bank_path).unwrap();
    fs::write(&bank_path, &original_bank[..original_bank.len() / 2]).unwrap();
    assert!(!call(&["fusion", "-i", &first, &second, "-o", &output])
        .status
        .success());
    assert!(!Path::new(&output).exists());
    fs::write(bank_path, original_bank).unwrap();
    fs::write(Path::new(&second).join("abandoned.zormap.zst"), b"broken").unwrap();
    assert!(!call(&["fusion", "-i", &first, &second, "-o", &output])
        .status
        .success());
    assert!(!Path::new(&output).exists());
    assert!(!fs::read_dir(&f.0).unwrap().any(|e| e
        .unwrap()
        .file_name()
        .to_string_lossy()
        .starts_with(".zor-fusing-")));
}

#[test]
fn fusion_preserves_skipped_input_provenance_and_disabled_recovery() {
    let f = Fixture::new();
    let a = f.write("a.fa", format!(">a\n{}\n", dna(900, 55)));
    let fof = f.write("inputs.fof", format!("a\t{a}\n"));
    let first = f.path("first");
    let second = f.path("second");
    let merged = f.path("merged");
    build(&fof, &first, &["-R"]);
    build(&fof, &second, &["-R"]);
    for index in [&first, &second] {
        let mut m = manifest(index);
        m["input_document_count"] = serde_json::json!(2);
        m["source_document_ids"] = serde_json::json!([1]);
        m["skipped_documents"] = serde_json::json!([{
            "source_document_id": 0, "name": "skipped", "path": "large.fa", "stage": "extraction",
            "estimated_workspace_bytes": 2000, "memory_budget_bytes": 1000, "reason": "too large"
        }]);
        write_manifest(index, &m);
    }
    ok(&["fusion", "-i", &first, &second, "-o", &merged]);
    let m = manifest(&merged);
    assert_eq!(m["input_document_count"], 4);
    assert_eq!(m["source_document_ids"], serde_json::json!([1, 3]));
    assert_eq!(m["skipped_documents"][0]["source_document_id"], 0);
    assert_eq!(m["skipped_documents"][1]["source_document_id"], 2);
    assert_eq!(m["abandoned_recovery"], false);
    assert!(!Path::new(&merged).join("abandoned.zormap.zst").exists());
    assert!(Path::new(&merged).join("skipped-documents.json").exists());
    assert!(ok(&["info", "-i", &merged]).contains("collection_status\tincomplete\n"));
    let mut expected = rows(&ok(&["query", "-i", &first, "-q", &a]));
    expected.extend(rows(&ok(&["query", "-i", &second, "-q", &a])));
    assert_eq!(rows(&ok(&["query", "-i", &merged, "-q", &a])), expected);
    let recovery = f.path("recovery");
    build(&fof, &recovery, &[]);
    let failed = f.path("failed");
    let result = call(&["fusion", "-i", &first, &recovery, "-o", &failed]);
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("abandoned_recovery differs"));
    assert!(!Path::new(&failed).exists());
}

#[test]
fn nway_fusion_matches_pairwise_tree_and_preserves_order_and_threads() {
    let f = Fixture::new();
    let mut inputs = Vec::new();
    let mut original_manifests = Vec::new();
    let mut original_banks = Vec::new();
    for i in 0..4 {
        let sequence = f.write(
            &format!("{i}.fa"),
            format!(">s{i}\n{}\n", dna(2000, 300 + i % 2)),
        );
        let fof = f.write(&format!("{i}.fof"), format!("s{i}\t{sequence}\n"));
        let index = f.path(&format!("index{i}"));
        build(&fof, &index, &["-i", "minimizers", "-m", "31", "-H", "8"]);
        let source = manifest(&index);
        assert_eq!(
            source["column_order_method"],
            "union-sketch-fingerprints-lexicographic-v1"
        );
        for bank in source["banks"].as_array().unwrap() {
            original_banks
                .push(fs::read(Path::new(&index).join(bank["file"].as_str().unwrap())).unwrap());
        }
        original_manifests.push(source);
        inputs.push(index);
    }
    let query = f.write(
        "q.fa",
        format!(">q0\n{}\n>q1\n{}\n", dna(1000, 300), dna(1000, 999)),
    );
    let fused = f.path("fused");
    let output = Command::new(env!("CARGO_BIN_EXE_zorindex"))
        .args([
            "-t", "4", "fusion", "-i", &inputs[0], &inputs[1], &inputs[2], &inputs[3], "-o", &fused,
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(String::from_utf8_lossy(&output.stderr).contains("worker_threads=4"));
    let combined = manifest(&fused);
    // Fusion preserves each input layout and must not claim a new global ordering.
    assert!(combined["column_order_method"].is_null());
    assert_eq!(
        combined["banks"].as_array().unwrap().len(),
        original_banks.len()
    );
    for (bank, original) in combined["banks"]
        .as_array()
        .unwrap()
        .iter()
        .zip(original_banks)
    {
        assert_eq!(
            fs::read(Path::new(&fused).join(bank["file"].as_str().unwrap())).unwrap(),
            original
        );
    }
    for (input, original) in inputs.iter().zip(original_manifests) {
        assert_eq!(manifest(input), original);
    }
    let ab = f.path("ab");
    let cd = f.path("cd");
    let tree = f.path("tree");
    ok(&["fusion", "-i", &inputs[0], &inputs[1], "-o", &ab]);
    ok(&["fusion", "-i", &inputs[2], &inputs[3], "-o", &cd]);
    ok(&["fusion", "-i", &ab, &cd, "-o", &tree]);
    assert_eq!(
        manifest(&fused)["colors"]
            .as_array()
            .unwrap()
            .iter()
            .map(|c| c["name"].as_str().unwrap())
            .collect::<Vec<_>>(),
        vec!["s0", "s1", "s2", "s3"]
    );
    assert_eq!(
        manifest(&fused)["source_document_ids"],
        serde_json::json!([0, 1, 2, 3])
    );
    assert_eq!(manifest(&fused)["banks"], manifest(&tree)["banks"]);
    assert_eq!(
        fs::read(Path::new(&fused).join("abandoned.zormap.zst")).unwrap(),
        fs::read(Path::new(&tree).join("abandoned.zormap.zst")).unwrap()
    );
    let expected = ok(&["query", "-i", &tree, "-q", &query, "-r"]);
    for input in &inputs {
        fs::remove_dir_all(input).unwrap();
    }
    for memory in ["all", "subindex"] {
        assert_eq!(
            records(&ok(&[
                "query", "-i", &fused, "-q", &query, "-r", "-M", memory
            ])),
            records(&expected)
        );
    }
}

#[test]
fn nway_fusion_checks_every_input_and_cleans_failed_staging() {
    let f = Fixture::new();
    let sequence = f.write("a.fa", format!(">a\n{}\n", dna(1000, 1)));
    let fof = f.write("a.fof", format!("a\t{sequence}\n"));
    let a = f.path("a");
    let b = f.path("b");
    build(&fof, &a, &[]);
    build(&fof, &b, &["-H", "8"]);
    let output = f.path("output");
    assert!(!call(&["fusion", "-i", &a, &a, &b, "-o", &output])
        .status
        .success());
    assert!(!Path::new(&output).exists());
    assert!(!call(&["fusion", "-i", &a, "-o", &output]).status.success());
    let m = manifest(&b);
    let bank = m["banks"][0]["file"].as_str().unwrap();
    fs::write(Path::new(&b).join(bank), b"broken").unwrap();
    assert!(!call(&["fusion", "-i", &b, &b, &b, "-o", &output])
        .status
        .success());
    assert!(!Path::new(&output).exists());
    assert!(!fs::read_dir(&f.0).unwrap().any(|p| p
        .unwrap()
        .file_name()
        .to_string_lossy()
        .starts_with(".zor-fusing-")));
}
