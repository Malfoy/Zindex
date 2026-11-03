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
            "zor-pow2-test-{}-{}",
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
fn call(args: &[&str], threads: &str) -> Output {
    Command::new(env!("CARGO_BIN_EXE_zorindex"))
        .args(["--threads", threads])
        .args(args)
        .output()
        .unwrap()
}
fn ok(args: &[&str], threads: &str) -> String {
    let result = call(args, threads);
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
fn stat(output: &str, field: &str) -> u64 {
    output
        .lines()
        .find_map(|line| {
            line.strip_prefix(&format!("{field}\t"))
                .and_then(|s| s.parse().ok())
        })
        .unwrap()
}
fn build(fof: &str, output: &str, extra: &[&str]) {
    let mut args = vec![
        "build",
        "--fof",
        fof,
        "--output",
        output,
        "--power-of-two-merge",
    ];
    for (flag, value) in [
        ("--k", "31"),
        ("--fingerprint-bits", "16"),
        ("--slot-scale", "1.0"),
        ("--min-subindex-documents", "1"),
    ] {
        if !extra.contains(&flag) {
            args.extend([flag, value]);
        }
    }
    args.extend(extra);
    ok(&args, "4");
}

#[test]
fn naive_ocm_keeps_singletons_and_first_capacity_with_identical_query_engines() {
    let f = Fixture::new();
    let large = dna(30_000, 901);
    let small = dna(600, 902);
    let a = f.write("large.fa", format!(">large\n{large}\n"));
    let b = f.write("small.fa", format!(">small\n{small}\n"));
    let fof = f.write("naive.fof", format!("large\t{a}\nsmall\t{b}\n"));
    let index = f.path("naive-index");
    ok(
        &[
            "build",
            "--fof",
            &fof,
            "--output",
            &index,
            "--k",
            "31",
            "--index-mode",
            "minimizers",
            "--window",
            "19",
            "--minimizer-scheme",
            "open-closed-mod",
            "--hashes",
            "8",
            "--fingerprint-bits",
            "8",
            "--merge-planner",
            "naive",
            "--merge-plan-margin",
            "0",
            "--max-subindexes",
            "1",
        ],
        "4",
    );
    let m: serde_json::Value =
        serde_json::from_slice(&fs::read(Path::new(&index).join("manifest.json")).unwrap())
            .unwrap();
    assert_eq!(m["merge_planner"], "naive");
    assert_eq!(m["planning_pair_evaluations"], 0);
    assert_eq!(m["column_order"], "width");
    let subindexes = m["banks"].as_array().unwrap();
    assert_eq!(subindexes.len(), 1);
    let groups = subindexes[0]["groups"].as_array().unwrap();
    assert_eq!(groups.len(), 2);
    assert_eq!(groups[0]["members"], serde_json::json!([0]));
    assert_eq!(groups[1]["members"], serde_json::json!([1]));
    assert!(groups.iter().all(|g| g["width"] == 8));
    let query = f.write(
        "queries.fa",
        format!(
            ">positive-large\n{}\n>positive-small\n{small}\n>negative\n{}\n",
            &large[..1000],
            dna(1000, 903)
        ),
    );
    let mut baseline = None;
    for threads in ["1", "4"] {
        for memory in ["subindex", "all"] {
            let result = ok(
                &[
                    "query",
                    "--index",
                    &index,
                    "--query",
                    &query,
                    "--records",
                    "--min-ratio",
                    "0.8",
                    "--query-memory",
                    memory,
                ],
                threads,
            );
            if let Some(expected) = &baseline {
                assert_eq!(&result, expected);
            } else {
                baseline = Some(result);
            }
        }
    }
    // The reverse order creates the small capacity before the larger one;
    // --max-subindexes must not trigger a global reassignment in naive mode.
    let reverse = f.write("reverse.fof", format!("small\t{b}\nlarge\t{a}\n"));
    let reverse_index = f.path("reverse-index");
    ok(
        &[
            "build",
            "--fof",
            &reverse,
            "--output",
            &reverse_index,
            "--k",
            "21",
            "--window",
            "29",
            "--index-mode",
            "minimizers",
            "--hashes",
            "8",
            "--fingerprint-bits",
            "8",
            "--merge-planner",
            "naive",
            "--max-subindexes",
            "1",
        ],
        "4",
    );
    let m: serde_json::Value =
        serde_json::from_slice(&fs::read(Path::new(&reverse_index).join("manifest.json")).unwrap())
            .unwrap();
    assert_eq!(m["banks"].as_array().unwrap().len(), 2);
    // The public name full is accepted and retains all-pairs on-disk metadata.
    let full = f.path("full-index");
    ok(
        &[
            "build",
            "--fof",
            &fof,
            "--output",
            &full,
            "--k",
            "21",
            "--window",
            "29",
            "--index-mode",
            "minimizers",
            "--merge-planner",
            "full",
        ],
        "4",
    );
    let m: serde_json::Value =
        serde_json::from_slice(&fs::read(Path::new(&full).join("manifest.json")).unwrap()).unwrap();
    assert_eq!(m["merge_planner"], "all-pairs");
}

#[test]
fn short_default_command_matches_explicit_optimized_settings() {
    let f = Fixture::new();
    let shared = dna(1100, 42);
    let mut inputs = String::new();
    for id in 0..4 {
        let private = dna(100, id + 69);
        let path = f.write(
            &format!("default-{id}.fa"),
            format!(">shared\n{shared}\n>private\n{private}\n"),
        );
        inputs.push_str(&format!("{id}\t{path}\n"));
    }
    let fof = f.write("defaults.fof", inputs);
    let query = f.write("query.fa", format!(">query\n{shared}\n"));
    let explicit = [
        "--power-of-two-merge",
        "--merge-planner",
        "streaming",
        "--merge-column-order",
        "compression",
        "--hashes",
        "4",
        "--slot-scale",
        "1.0",
        "--cycle-break",
        "no-heuristic",
        "--tie-scan",
        "1",
    ];
    let mut indexes = Vec::new();
    let mut outputs = Vec::new();
    for (label, settings) in [("default", &[][..]), ("explicit", &explicit[..])] {
        let index = f.path(label);
        let mut args = vec![
            "build",
            "--fof",
            &fof,
            "--output",
            &index,
            "--k",
            "31",
            "--index-mode",
            "minimizers",
            "--window",
            "19",
            "--fingerprint-bits",
            "16",
        ];
        args.extend(settings);
        ok(&args, "4");
        let manifest: serde_json::Value =
            serde_json::from_slice(&fs::read(Path::new(&index).join("manifest.json")).unwrap())
                .unwrap();
        assert_eq!(manifest["merge_planner"], "streaming");
        assert_eq!(manifest["version"], 8);
        assert_eq!(manifest["segment_length"], 64);
        assert_eq!(manifest["min_subindex_documents"], 128);
        assert_eq!(manifest["column_order"], "compression");
        assert_eq!(
            manifest["construction_cell_storage"],
            "two-pass-temporary-disk-packed-tiled-packing"
        );
        assert_eq!(manifest["abandoned_recovery"], true);
        assert_eq!(manifest["rejected_groups"], 0);
        let result = ok(&["query", "--index", &index, "--query", &query], "4");
        outputs.push(
            result
                .lines()
                .filter(|line| !line.starts_with('#'))
                .collect::<Vec<_>>()
                .join("\n"),
        );
        indexes.push(index);
    }
    assert_eq!(outputs[0], outputs[1]);
    for entry in fs::read_dir(&indexes[0]).unwrap() {
        let name = entry.unwrap().file_name();
        if name == "manifest.json" {
            continue;
        }
        assert_eq!(
            fs::read(Path::new(&indexes[0]).join(&name)).unwrap(),
            fs::read(Path::new(&indexes[1]).join(&name)).unwrap(),
            "{name:?}"
        );
    }

    let dropped = f.path("drop-default");
    ok(
        &[
            "build",
            "--fof",
            &fof,
            "--output",
            &dropped,
            "--k",
            "31",
            "--drop-abandoned",
        ],
        "4",
    );
    let m: serde_json::Value =
        serde_json::from_slice(&fs::read(Path::new(&dropped).join("manifest.json")).unwrap())
            .unwrap();
    assert_eq!(m["abandoned_recovery"], false);
    assert!(!Path::new(&dropped).join("abandoned.zormap.zst").exists());
}

#[test]
fn merged_origins_and_serialized_size_classes_query_correctly() {
    let f = Fixture::new();
    let common = dna(700, 17);
    let mut fof = String::new();
    let mut privates = Vec::new();
    for id in 0..4 {
        let private = dna(32, 100 + id);
        let path = f.write(
            &format!("{id}.fa"),
            format!(">common\n{common}\n>private\n{private}\n"),
        );
        fof.push_str(&format!("c{id}\t{path}\n"));
        privates.push(private);
    }
    let big = f.write("big.fa", format!(">big\n{}\n", dna(4000, 91)));
    fof.push_str(&format!("big\t{big}\n"));
    let inputs = f.write("inputs.fof", fof);
    let index = f.path("merged");
    build(&inputs, &index, &[]);
    let info = ok(&["info", "--index", &index], "1");
    assert!(
        info.contains("construction_cell_storage\ttwo-pass-temporary-disk-packed-tiled-packing")
    );
    assert!(info.contains("merge_planner\tstreaming"));
    assert!(info.contains("merge_column_order\tcompression"));
    assert!(!Path::new(&index).join("cells").exists());
    assert!(!Path::new(&index).join("features").exists());
    assert!(fs::read_dir(&f.0).unwrap().all(|entry| !entry
        .unwrap()
        .file_name()
        .to_string_lossy()
        .starts_with(".zor-pow2-building-")));
    assert_eq!(stat(&info, "colors"), 5);
    assert!(stat(&info, "groups") < 5, "{info}");
    assert!(stat(&info, "banks") >= 2, "{info}");
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(Path::new(&index).join("manifest.json")).unwrap())
            .unwrap();
    let mut last = 0;
    for bank in manifest["banks"].as_array().unwrap() {
        let slots = bank["slots"].as_u64().unwrap();
        assert!(slots.is_power_of_two() && slots > last);
        last = slots;
        assert!(Path::new(&index)
            .join(bank["file"].as_str().unwrap())
            .exists());
    }
    let query = f.write("common.fa", format!(">q\n{common}\n"));
    for threads in ["1", "4"] {
        let result = ok(&["query", "--index", &index, "--query", &query], threads);
        for id in 0..4 {
            assert_eq!(matches(&result, &format!("c{id}")), 670, "{result}");
        }
    }
    for (owner, seq) in privates.iter().enumerate() {
        let query = f.write("private.fa", format!(">q\n{seq}\n"));
        let result = ok(&["query", "--index", &index, "--query", &query], "4");
        for id in 0..4 {
            assert_eq!(
                matches(&result, &format!("c{id}")),
                if id == owner { 2 } else { 0 },
                "{result}"
            );
        }
    }
    let singleton = f.path("singletons");
    build(&inputs, &singleton, &["--merge-disabled"]);
    let single_info = ok(&["info", "--index", &singleton], "1");
    assert_eq!(stat(&single_info, "groups"), 5);
    assert!(stat(&info, "structure_bytes") < stat(&single_info, "structure_bytes"));
    // Existing output must never be overwritten.
    let result = call(
        &[
            "build",
            "--fof",
            &inputs,
            "--output",
            &index,
            "--k",
            "31",
            "--power-of-two-merge",
        ],
        "1",
    );
    assert!(!result.status.success());
}

#[test]
fn findere_parallel_queries_keep_boundaries_and_compressed_input() {
    let f = Fixture::new();
    let a = dna(100_000, 41);
    let reference = f.write("a.fa", format!(">a\n{a}\n"));
    let inputs = f.write("inputs.fof", format!("a\t{reference}\nb\t{reference}\n"));
    let query_text = format!(
        ">forward\n{a}\n>broken\n{}N{}\n>short\nACGT\n",
        &a[..80],
        &a[80..160]
    );
    let query = f.write(
        "q.fa.zst",
        zstd::bulk::compress(query_text.as_bytes(), 1).unwrap(),
    );
    for (k, z) in [(31, 10), (41, 10), (31, 0)] {
        let index = f.path(&format!("merged-{k}-{z}"));
        build(
            &inputs,
            &index,
            &[
                "--index-mode",
                "findere",
                "--k",
                &k.to_string(),
                "--findere-z",
                &z.to_string(),
            ],
        );
        let info = ok(&["info", "--index", &index], "1");
        assert_eq!(stat(&info, "findere_z"), z as u64);
        let one = ok(&["query", "--index", &index, "--query", &query], "1");
        let four = ok(&["query", "--index", &index, "--query", &query], "4");
        assert_eq!(
            one.lines().skip(1).collect::<Vec<_>>(),
            four.lines().skip(1).collect::<Vec<_>>()
        );
        let expected = (a.len() - k + 1 + 2 * (80 - k + 1)) as u64;
        assert_eq!(matches(&one, "a"), expected);
        assert_eq!(matches(&one, "b"), expected);
    }
}

#[test]
fn expensive_exact_fallback_does_not_revert_the_planned_merge() {
    let f = Fixture::new();
    let seq = dna(700, 13);
    let reference = f.write("a.fa", format!(">a\n{seq}\n"));
    let inputs = f.write("inputs.fof", format!("a\t{reference}\nb\t{reference}\n"));
    for planner in ["all-pairs", "streaming"] {
        let index = f.path(planner);
        build(
            &inputs,
            &index,
            &["--slot-scale", "0.01", "--merge-planner", planner],
        );
        let info = ok(&["info", "--index", &index], "1");
        assert!(stat(&info, "fallback_bytes") > 0, "{info}");
        assert_eq!(stat(&info, "rejected_groups"), 0, "{info}");
        assert_eq!(stat(&info, "groups"), 1, "{info}");
        assert_eq!(stat(&info, "planned_groups"), 1, "{info}");
        assert!(info.contains("merge_reverts_enabled\tfalse"));
        let result = ok(&["query", "--index", &index, "--query", &reference], "4");
        assert_eq!(matches(&result, "a"), 670);
        assert_eq!(matches(&result, "b"), 670);
    }
}

#[test]
fn streaming_planner_preserves_input_order_and_all_origins() {
    let f = Fixture::new();
    let seq = dna(700, 888);
    let reference = f.write("a.fa", format!(">a\n{seq}\n"));
    let inputs = f.write(
        "inputs.fof",
        format!("z\t{reference}\na\t{reference}\nb\t{reference}\n"),
    );
    for fp in ["8", "16"] {
        let index = f.path(&format!("streaming-{fp}"));
        build(
            &inputs,
            &index,
            &["--merge-planner", "streaming", "--fingerprint-bits", fp],
        );
        let info = ok(&["info", "--index", &index], "1");
        assert!(info.contains("merge_planner\tstreaming"));
        assert_eq!(stat(&info, "groups"), 1);
        assert_eq!(stat(&info, "planning_pair_evaluations"), 2);
        assert_eq!(stat(&info, "rejected_groups"), 0);
        let query = ok(&["query", "--index", &index, "--query", &reference], "4");
        for name in ["z", "a", "b"] {
            assert_eq!(matches(&query, name), 670);
        }
        let names: Vec<_> = query
            .lines()
            .skip(2)
            .map(|line| line.split('\t').next().unwrap())
            .collect();
        assert_eq!(names, ["z", "a", "b"]);
    }
}

#[test]
fn global_abandoned_recovery_can_be_omitted_or_skipped_without_changing_filters() {
    let f = Fixture::new();
    // 1024 distinct 31-mers, at scale 1.0: exercise cycle breaking rather than
    // artificially changing filter sizing for this recovery test.
    let reference = f.write("a.fa", format!(">a\n{}\n", dna(1054, 13)));
    let inputs = f.write("inputs.fof", format!("a\t{reference}\nb\t{reference}\n"));
    let recovery = f.path("recovery");
    let dropped = f.path("dropped");
    build(&inputs, &recovery, &[]);
    build(&inputs, &dropped, &["--drop-abandoned"]);
    let info = ok(&["info", "--index", &recovery], "1");
    let drop_info = ok(&["info", "--index", &dropped], "1");
    assert!(info.contains("abandoned_recovery\ttrue"));
    assert!(info.contains("abandoned_hash_bits\t32"));
    assert!(stat(&info, "abandoned_keys") > 0, "{info}");
    assert!(stat(&info, "abandoned_unique_hashes") > 0);
    assert!(drop_info.contains("abandoned_recovery\tfalse"));
    assert_eq!(stat(&drop_info, "fallback_bytes"), 0);
    assert_eq!(stat(&drop_info, "abandoned_map_ram_bytes"), 0);
    assert_eq!(
        stat(&info, "abandoned_keys"),
        stat(&drop_info, "abandoned_keys")
    );
    assert!(!Path::new(&dropped).join("abandoned.zormap.zst").exists());
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(Path::new(&recovery).join("manifest.json")).unwrap())
            .unwrap();
    assert_eq!(manifest["version"], 8);
    for bank in manifest["banks"].as_array().unwrap() {
        let name = bank["file"].as_str().unwrap();
        assert_eq!(
            fs::read(Path::new(&recovery).join(name)).unwrap(),
            fs::read(Path::new(&dropped).join(name)).unwrap()
        );
    }
    let recovered = ok(&["query", "--index", &recovery, "--query", &reference], "4");
    let ignored = ok(
        &[
            "query",
            "--index",
            &recovery,
            "--query",
            &reference,
            "--ignore-abandoned",
        ],
        "4",
    );
    let raw = ok(&["query", "--index", &dropped, "--query", &reference], "4");
    for color in ["a", "b"] {
        assert_eq!(matches(&recovered, color), 1024, "{recovered}");
        assert!(matches(&raw, color) < 1024, "{raw}");
        assert_eq!(matches(&raw, color), matches(&ignored, color));
    }
    // Skipping recovery does not even open the global map file.
    fs::remove_file(Path::new(&recovery).join("abandoned.zormap.zst")).unwrap();
    let skipped = ok(
        &[
            "query",
            "--index",
            &recovery,
            "--query",
            &reference,
            "--ignore-abandoned",
        ],
        "1",
    );
    for color in ["a", "b"] {
        assert_eq!(matches(&skipped, color), matches(&raw, color));
    }
    assert!(
        !call(&["query", "--index", &recovery, "--query", &reference], "1")
            .status
            .success()
    );
}

#[test]
fn findere_recovers_abandoned_features_with_one_global_map_in_parallel() {
    let f = Fixture::new();
    // 1024 indexed 21-mers at the unchanged 1.0 sizing factor.
    let reference = f.write("a.fa", format!(">a\n{}\n", dna(1044, 13)));
    let inputs = f.write("inputs.fof", format!("a\t{reference}\nb\t{reference}\n"));
    let index = f.path("findere-recovery");
    build(
        &inputs,
        &index,
        &["--index-mode", "findere", "--findere-z", "10"],
    );
    let info = ok(&["info", "--index", &index], "1");
    assert!(stat(&info, "abandoned_keys") > 0, "{info}");
    for threads in ["1", "4"] {
        let result = ok(
            &["query", "--index", &index, "--query", &reference],
            threads,
        );
        for color in ["a", "b"] {
            assert_eq!(matches(&result, color), 1014, "{result}");
        }
    }
}

#[test]
fn pruned_streaming_builds_identical_banks_and_queries_to_original() {
    let f = Fixture::new();
    let mut fof = String::new();
    for id in 0..24 {
        let size = [64, 256, 1024, 4096][id % 4];
        let seq = dna(size, (id % 3) as u64);
        let path = f.write(&format!("g{id}.fa"), format!(">r\n{seq}\n"));
        fof.push_str(&format!("g{id}\t{path}\n"));
    }
    let inputs = f.write("inputs.fof", fof);
    let query = f.write("query.fa", format!(">r\n{}\n", dna(4096, 0)));
    for bits in ["8", "16"] {
        for order in ["width", "compression"] {
            let old = f.path(&format!("old-{bits}-{order}"));
            let new = f.path(&format!("new-{bits}-{order}"));
            for (path, planner) in [(&old, "streaming-unpruned"), (&new, "streaming")] {
                build(
                    &inputs,
                    path,
                    &[
                        "--merge-planner",
                        planner,
                        "--fingerprint-bits",
                        bits,
                        "--merge-column-order",
                        order,
                        "--merge-max-colors",
                        "3",
                    ],
                );
            }
            let load = |path: &str| -> serde_json::Value {
                serde_json::from_slice(&fs::read(Path::new(path).join("manifest.json")).unwrap())
                    .unwrap()
            };
            let a = load(&old);
            let b = load(&new);
            assert_eq!(a["slot_scale"], 1.0);
            assert_eq!(a["banks"], b["banks"]);
            assert_eq!(a["predicted_logical_bits"], b["predicted_logical_bits"]);
            assert_eq!(
                a["planning_pair_evaluations"],
                b["planning_possible_pair_evaluations"]
            );
            assert!(
                b["planning_pair_evaluations"].as_u64().unwrap()
                    < a["planning_pair_evaluations"].as_u64().unwrap()
            );
            for bank in a["banks"].as_array().unwrap() {
                let name = bank["file"].as_str().unwrap();
                assert_eq!(
                    fs::read(Path::new(&old).join(name)).unwrap(),
                    fs::read(Path::new(&new).join(name)).unwrap()
                );
            }
            let a = ok(&["query", "--index", &old, "--query", &query], "4");
            let b = ok(&["query", "--index", &new, "--query", &query], "4");
            assert_eq!(
                a.lines().skip(1).collect::<Vec<_>>(),
                b.lines().skip(1).collect::<Vec<_>>()
            );
        }
    }
}

#[test]
fn physical_column_order_preserves_partitions_scores_and_dataset_order() {
    let f = Fixture::new();
    let sequences: Vec<_> = (0..4).map(|id| dna(700, 777 + id)).collect();
    let paths: Vec<_> = sequences
        .iter()
        .enumerate()
        .map(|(id, seq)| f.write(&format!("{id}.fa"), format!(">r\n{seq}\n")))
        .collect();
    let origins = [0, 1, 0, 2, 1, 2, 0, 3, 0];
    let names: Vec<_> = (0..origins.len())
        .map(|id| format!("c{}", 20 - id))
        .collect();
    let inputs = f.write(
        "inputs.fof",
        origins
            .iter()
            .enumerate()
            .map(|(id, &origin)| format!("{}\t{}\n", names[id], paths[origin]))
            .collect::<String>(),
    );
    for planner in ["streaming", "all-pairs"] {
        let mut baseline = None;
        let mut baseline_partition = None;
        for order in ["width", "compression"] {
            let index = f.path(&format!("{planner}-{order}"));
            build(
                &inputs,
                &index,
                &[
                    "--merge-planner",
                    planner,
                    "--merge-column-order",
                    order,
                    "--merge-max-colors",
                    "2",
                ],
            );
            let info = ok(&["info", "--index", &index], "1");
            assert!(info.contains(&format!("merge_planner\t{planner}")));
            assert!(info.contains(&format!("merge_column_order\t{order}")));
            assert!(info.contains("merge_reverts_enabled\tfalse"));
            let manifest: serde_json::Value =
                serde_json::from_slice(&fs::read(Path::new(&index).join("manifest.json")).unwrap())
                    .unwrap();
            let mut partition = Vec::new();
            for bank in manifest["banks"].as_array().unwrap() {
                let mut previous_width = 0;
                for group in bank["groups"].as_array().unwrap() {
                    let width = group["width"].as_u64().unwrap();
                    assert!(width >= previous_width);
                    previous_width = width;
                    partition.push(
                        group["members"]
                            .as_array()
                            .unwrap()
                            .iter()
                            .map(|id| id.as_u64().unwrap())
                            .collect::<Vec<_>>(),
                    );
                }
            }
            partition.sort();
            if let Some(old) = &baseline_partition {
                assert_eq!(old, &partition);
            }
            baseline_partition = Some(partition);
            let mut outputs = Vec::new();
            for (origin, path) in paths.iter().enumerate() {
                let result = ok(&["query", "--index", &index, "--query", path], "4");
                // The first line includes a fresh elapsed wall time; compare
                // result rows, not this deliberately variable timing field.
                let rows = result.lines().skip(1).collect::<Vec<_>>().join("\n");
                let single = ok(&["query", "--index", &index, "--query", path], "1");
                assert_eq!(rows, single.lines().skip(1).collect::<Vec<_>>().join("\n"));
                for (id, &owner) in origins.iter().enumerate() {
                    if owner == origin {
                        assert_eq!(matches(&result, &names[id]), 670);
                    }
                }
                assert_eq!(
                    result
                        .lines()
                        .skip(2)
                        .map(|line| line.split('\t').next().unwrap())
                        .collect::<Vec<_>>(),
                    names
                );
                outputs.push(rows);
            }
            if let Some(old) = &baseline {
                assert_eq!(old, &outputs);
            }
            baseline = Some(outputs);
        }
    }
}

#[test]
fn slice_padding_does_not_revert_a_planned_merge() {
    let f = Fixture::new();
    // No 31-mers: two 8-bit singleton payloads cost exactly as much as a
    // 10-bit merged payload once the latter is padded to two bytes per row.
    let reference = f.write("short.fa", ">a\nACGT\n");
    let inputs = f.write("inputs.fof", format!("a\t{reference}\nb\t{reference}\n"));
    for planner in ["all-pairs", "streaming"] {
        let index = f.path(planner);
        build(
            &inputs,
            &index,
            &["--fingerprint-bits", "8", "--merge-planner", planner],
        );
        let info = ok(&["info", "--index", &index], "1");
        assert_eq!(stat(&info, "groups"), 1);
        assert_eq!(stat(&info, "planned_groups"), 1);
        assert_eq!(stat(&info, "rejected_groups"), 0);
        assert_eq!(
            stat(&info, "stack_bytes") * 8,
            stat(&info, "singleton_logical_bits")
        );
        assert!(stat(&info, "padding_bits") > 0);
    }
}

#[test]
fn non_default_hash_counts_and_feature_modes_use_matching_bank_addresses() {
    let f = Fixture::new();
    let seq = dna(500, 83);
    let reference = f.write("a.fa", format!(">a\n{seq}\n"));
    let inputs = f.write("inputs.fof", format!("a\t{reference}\nb\t{reference}\n"));
    for hashes in [2, 3, 4, 5, 6, 7, 8, 32] {
        let index = f.path(&format!("hashes-{hashes}"));
        build(
            &inputs,
            &index,
            &["--hashes", &hashes.to_string(), "--slot-scale", "2"],
        );
        let result = ok(&["query", "--index", &index, "--query", &reference], "4");
        assert_eq!(matches(&result, "a"), 470);
        assert_eq!(matches(&result, "b"), 470);
    }
    for mode in ["minimizers", "modimizers"] {
        let index = f.path(mode);
        let mut settings = vec!["--index-mode", mode];
        if mode == "minimizers" {
            settings.extend(["--k", "21", "--window", "11"]);
        }
        build(&inputs, &index, &settings);
        let result = ok(&["query", "--index", &index, "--query", &reference], "4");
        let total: u64 = result
            .lines()
            .next()
            .unwrap()
            .split('\t')
            .find_map(|s| s.strip_prefix("total_features="))
            .unwrap()
            .parse()
            .unwrap();
        assert!(total > 0);
        assert_eq!(matches(&result, "a"), total);
        assert_eq!(matches(&result, "b"), total);
    }
}

#[test]
fn incompatible_options_and_corrupt_mappings_are_rejected() {
    let f = Fixture::new();
    let reference = f.write("a.fa", format!(">a\n{}\n", dna(100, 2)));
    let inputs = f.write("inputs.fof", format!("a\t{reference}\n"));
    for extra in [
        vec!["--technical-layout", "split-large"],
        vec!["--union-graph"],
        vec!["--merge-sketch-size", "1000"],
        vec!["--merge-planner", "unknown"],
        vec!["--merge-max-colors", "64"],
        vec!["--stack-compression", "lz4"],
    ] {
        let output = f.path("invalid");
        let mut args = vec![
            "build",
            "--fof",
            &inputs,
            "--output",
            &output,
            "--k",
            "31",
            "--power-of-two-merge",
            "--fingerprint-bits",
            "16",
        ];
        args.extend(extra);
        assert!(!call(&args, "1").status.success());
        assert!(!Path::new(&output).exists());
    }
    let output = f.path("valid");
    build(&inputs, &output, &[]);
    let manifest = Path::new(&output).join("manifest.json");
    let mut json: serde_json::Value =
        serde_json::from_slice(&fs::read(&manifest).unwrap()).unwrap();
    json["banks"][0]["groups"][0]["members"][0] = serde_json::json!(9999);
    fs::write(&manifest, serde_json::to_vec(&json).unwrap()).unwrap();
    assert!(!call(&["info", "--index", &output], "1").status.success());
}

#[test]
fn full_64_bit_membership_payload_survives_unaligned_slice_serialization() {
    let f = Fixture::new();
    let a = f.write("a.fa", format!(">a\n{}\n", dna(700, 811)));
    let b = f.write("b.fa", format!(">b\n{}\n", dna(700, 923)));
    let mut fof = String::new();
    for id in 0..48 {
        fof.push_str(&format!("a{id}\t{a}\n"));
    }
    for id in 0..2 {
        fof.push_str(&format!("b{id}\t{b}\n"));
    }
    let inputs = f.write("inputs.fof", fof);
    for planner in ["all-pairs", "streaming"] {
        let index = f.path(planner);
        build(&inputs, &index, &["--merge-planner", planner]);
        let manifest: serde_json::Value =
            serde_json::from_slice(&fs::read(Path::new(&index).join("manifest.json")).unwrap())
                .unwrap();
        let groups = manifest["banks"]
            .as_array()
            .unwrap()
            .iter()
            .flat_map(|b| b["groups"].as_array().unwrap())
            .collect::<Vec<_>>();
        assert!(groups
            .iter()
            .any(|g| g["width"] == 64 && g["offset"].as_u64().unwrap() % 8 != 0));
        let query = ok(&["query", "--index", &index, "--query", &a], "4");
        for id in 0..48 {
            assert_eq!(matches(&query, &format!("a{id}")), 670);
        }
        for id in 0..2 {
            assert_eq!(matches(&query, &format!("b{id}")), 0);
        }
    }
}

#[test]
fn all_cell_storage_modes_build_identical_filters_and_queries() {
    let f = Fixture::new();
    let shared = dna(1024, 17);
    let mut inputs = String::new();
    for id in 0..10 {
        let seq = if id < 6 {
            shared.clone()
        } else {
            dna([100, 300, 600, 1200][id - 6], id as u64)
        };
        let private = dna(45 + id * 7, id as u64 + 111);
        let path = f.write(
            &format!("d{id}.fa"),
            format!(">r\n{seq}\n>private\n{private}\n"),
        );
        inputs.push_str(&format!("d{id}\t{path}\n"));
    }
    let fof = f.write("inputs.fof", inputs);
    let query = f.write(
        "query.fa",
        format!(">present\n{shared}\n>absent\n{}\n", dna(300, 999)),
    );
    for planner in ["streaming", "streaming-unpruned", "all-pairs"] {
        for fp in ["1", "2", "4", "8", "16"] {
            for order in ["width", "compression"] {
                let mut outputs = Vec::new();
                let mut manifests = Vec::new();
                let mut results = Vec::new();
                for (mode, storage) in [
                    (None, "two-pass-temporary-disk-packed-tiled-packing"),
                    (Some("--low-disk"), "two-pass-memory-packed-no-scratch"),
                    (Some("--disk-backed"), "temporary-disk-u64-tiled-packing"),
                ] {
                    let path = f.path(&format!("{planner}-{fp}-{order}-{storage}"));
                    let mut args = vec![
                        "--fingerprint-bits",
                        fp,
                        "--merge-planner",
                        planner,
                        "--merge-column-order",
                        order,
                    ];
                    args.extend(mode);
                    build(&fof, &path, &args);
                    let mut json: serde_json::Value = serde_json::from_slice(
                        &fs::read(Path::new(&path).join("manifest.json")).unwrap(),
                    )
                    .unwrap();
                    assert_eq!(json["construction_cell_storage"], storage);
                    assert_eq!(
                        json["construction_abandoned_storage"],
                        if mode == Some("--low-disk") {
                            "memory"
                        } else {
                            "temporary-disk-sorted-runs"
                        }
                    );
                    assert!(json["merge_memory_budget_bytes"].as_u64().unwrap() > 0);
                    for key in [
                        "construction_cell_storage",
                        "construction_abandoned_storage",
                        "sketch_seconds",
                        "plan_seconds",
                        "build_seconds",
                        // Auto budget snapshots available RAM separately for
                        // each build; it is provenance, not filter content.
                        "merge_memory_budget_bytes",
                    ] {
                        json.as_object_mut().unwrap().remove(key);
                    }
                    manifests.push(json);
                    let output = ok(
                        &["query", "--index", &path, "--query", &query],
                        if mode.is_some() { "1" } else { "4" },
                    );
                    results.push(
                        output
                            .lines()
                            .filter(|l| !l.starts_with('#'))
                            .collect::<Vec<_>>()
                            .join("\n"),
                    );
                    outputs.push(path);
                }
                for other in 1..outputs.len() {
                    assert_eq!(manifests[0], manifests[other]);
                    assert_eq!(results[0], results[other]);
                    for entry in fs::read_dir(&outputs[0]).unwrap() {
                        let name = entry.unwrap().file_name();
                        if name == "manifest.json" {
                            continue;
                        }
                        assert_eq!(
                            fs::read(Path::new(&outputs[0]).join(&name)).unwrap(),
                            fs::read(Path::new(&outputs[other]).join(&name)).unwrap(),
                            "{name:?}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn all_cell_storage_modes_preserve_all_feature_extractors() {
    let f = Fixture::new();
    let sequence = dna(500, 42);
    let a = f.write("a.fa", format!(">a\n{sequence}\n"));
    let fof = f.write("inputs.fof", format!("a\t{a}\nb\t{a}\n"));
    let query = f.write("query.fa", format!(">r\n{sequence}\n"));
    for (mode, settings) in [
        (
            "minimizers",
            vec!["--k", "31", "--index-mode", "minimizers", "--window", "19"],
        ),
        (
            "findere",
            vec!["--k", "31", "--index-mode", "findere", "--findere-z", "10"],
        ),
        ("modimizers", vec!["--index-mode", "modimizers"]),
        ("hpc", vec!["--homopolymer-compression"]),
        ("spaced", vec!["--k", "8", "--seed-pattern", "##_##_##"]),
        (
            "randstrobes",
            vec![
                "--k",
                "15",
                "--randstrobemers",
                "--strobe-min",
                "15",
                "--strobe-max",
                "30",
            ],
        ),
        ("no-recovery", vec!["--drop-abandoned"]),
    ] {
        let staged = f.path(&format!("{mode}-staged"));
        build(&fof, &staged, &settings);
        let one = ok(&["query", "--index", &staged, "--query", &query], "4");
        for flag in ["--low-disk", "--disk-backed"] {
            let other = f.path(&format!("{mode}{flag}"));
            let mut other_settings = settings.clone();
            other_settings.push(flag);
            build(&fof, &other, &other_settings);
            let two = ok(&["query", "--index", &other, "--query", &query], "1");
            assert_eq!(
                one.lines().skip(1).collect::<Vec<_>>(),
                two.lines().skip(1).collect::<Vec<_>>(),
                "{mode} {flag}"
            );
            for entry in fs::read_dir(&staged).unwrap() {
                let name = entry.unwrap().file_name();
                if name == "manifest.json" {
                    continue;
                }
                assert_eq!(
                    fs::read(Path::new(&staged).join(&name)).unwrap(),
                    fs::read(Path::new(&other).join(&name)).unwrap(),
                    "{mode} {flag} {name:?}"
                );
            }
        }
    }
    for flags in [
        ["--plain", "--disk-backed"],
        ["--plain", "--low-disk"],
        ["--low-disk", "--disk-backed"],
    ] {
        let invalid = f.path("invalid");
        assert!(!call(
            &["build", flags[0], flags[1], "--fof", &fof, "--output", &invalid, "--k", "31"],
            "1"
        )
        .status
        .success());
        assert!(!Path::new(&invalid).exists());
    }
}

#[test]
fn minimum_subindex_documents_promotes_without_changing_document_membership() {
    let f = Fixture::new();
    let mut fof = String::new();
    let mut paths = Vec::new();
    for id in 0..12 {
        let length = if id < 6 { 100 } else { 1000 };
        let sequence = dna(length, id + 1234);
        let path = f.write(&format!("minimum-{id}.fa"), format!(">r\n{sequence}\n"));
        fof.push_str(&format!("d{id}\t{path}\n"));
        paths.push(path);
    }
    let fof = f.write("minimum.fof", fof);
    for (minimum, maximum, expected_banks) in [
        ("1", None, 2),
        ("8", None, 1),
        ("128", None, 1),
        ("1", Some("1"), 1),
        ("1", Some("3"), 2),
        ("8", Some("1"), 1),
    ] {
        let index = f.path(&format!("minimum-{minimum}-maximum-{maximum:?}"));
        let mut extra = vec!["--merge-disabled", "--min-subindex-documents", minimum];
        if let Some(maximum) = maximum {
            extra.extend(["--max-subindexes", maximum]);
        }
        build(&fof, &index, &extra);
        let m: serde_json::Value =
            serde_json::from_slice(&fs::read(Path::new(&index).join("manifest.json")).unwrap())
                .unwrap();
        let banks = m["banks"].as_array().unwrap();
        assert_eq!(
            m["max_subindexes"].as_u64(),
            maximum.map(|m| m.parse::<u64>().unwrap())
        );
        assert_eq!(banks.len(), expected_banks);
        for bank in banks {
            assert!(bank["slots"].as_u64().unwrap().is_power_of_two());
            let docs: usize = bank["groups"]
                .as_array()
                .unwrap()
                .iter()
                .map(|g| g["members"].as_array().unwrap().len())
                .sum();
            assert!(docs >= minimum.parse::<usize>().unwrap().min(12));
        }
        for (id, path) in paths.iter().enumerate() {
            let output = ok(&["query", "--index", &index, "--query", path], "4");
            assert_eq!(
                matches(&output, &format!("d{id}")),
                if id < 6 { 70 } else { 970 }
            );
        }
    }
    let invalid = f.path("minimum-invalid");
    assert!(!call(
        &[
            "build",
            "--fof",
            &fof,
            "--output",
            &invalid,
            "--min-subindex-documents",
            "0"
        ],
        "1"
    )
    .status
    .success());
    assert!(!Path::new(&invalid).exists());
    let invalid = f.path("maximum-invalid");
    assert!(!call(
        &[
            "build",
            "--fof",
            &fof,
            "--output",
            &invalid,
            "--max-subindexes",
            "0"
        ],
        "1"
    )
    .status
    .success());
    assert!(!Path::new(&invalid).exists());
}

#[test]
fn record_queries_match_aggregate_counts_and_preserve_order_across_workers() {
    let f = Fixture::new();
    let shared = dna(90, 41);
    let private = dna(70, 42);
    let a = f.write("record-a.fa", format!(">a\n{shared}\n{private}\n"));
    let b = f.write("record-b.fa", format!(">b\n{shared}\n"));
    let fof = f.write("records.fof", format!("a\t{a}\nb\t{b}\n"));
    let sequences = [
        shared[..31].to_owned(),
        dna(31, 2026),
        format!("{}NNN{}", shared.to_lowercase(), private),
        "NNNN".to_owned(),
        String::new(),
    ];
    let fasta: String = sequences
        .iter()
        .enumerate()
        .map(|(id, seq)| format!(">q{id} description\r\n{seq}\r\n"))
        .collect();
    let query = f.write(
        "records.fa.zst",
        zstd::bulk::compress(fasta.as_bytes(), 1).unwrap(),
    );
    for bits in ["1", "2", "4", "8", "16"] {
        let index = f.path(&format!("records-fp{bits}"));
        build(&fof, &index, &["--fingerprint-bits", bits]);
        for recovery in [true, false] {
            let mut reference = Vec::new();
            for (id, seq) in sequences.iter().enumerate() {
                let single = f.write("single-record.fa", format!(">q{id}\n{seq}\n"));
                let mut args = vec!["query", "--index", &index, "--query", &single];
                if !recovery {
                    args.push("--ignore-abandoned");
                }
                let aggregate = ok(&args, "1");
                let total: u64 = aggregate
                    .lines()
                    .next()
                    .unwrap()
                    .split('\t')
                    .find_map(|field| field.strip_prefix("total_features="))
                    .unwrap()
                    .parse()
                    .unwrap();
                reference.push((total, [matches(&aggregate, "a"), matches(&aggregate, "b")]));
            }
            for threshold in ["0", "0.5", "1"] {
                let min_ratio: f64 = threshold.parse().unwrap();
                let mut expected = "query_id\tvalid_kmers\thits\n".to_owned();
                for (id, (total, scores)) in reference.iter().enumerate() {
                    let hits: Vec<_> = scores
                        .iter()
                        .enumerate()
                        .filter_map(|(origin, &score)| {
                            let ratio = if *total == 0 {
                                0.0
                            } else {
                                score as f64 / *total as f64
                            };
                            (ratio >= min_ratio).then(|| format!("{origin}:{score}"))
                        })
                        .collect();
                    expected.push_str(&format!("q{id}\t{total}\t{}\n", hits.join(",")));
                }
                for (threads, record_flag, batch_flag, batch) in [
                    ("1", "--records", "--record-batch-size", "1"),
                    ("4", "--records", "--record-batch-size", "3"),
                    ("4", "--bench-records", "--bench-batch", "4096"),
                ] {
                    let mut args = vec![
                        "query",
                        "--index",
                        &index,
                        "--query",
                        &query,
                        "--min-ratio",
                        threshold,
                        record_flag,
                        batch_flag,
                        batch,
                    ];
                    if !recovery {
                        args.push("--ignore-abandoned");
                    }
                    assert_eq!(ok(&args, threads), expected, "fp{bits} recovery={recovery}");
                }
            }
        }
    }
}

#[test]
fn sequential_queries_match_full_loading_for_seed_modes_widths_and_thresholds() {
    let f = Fixture::new();
    let shared = dna(700, 310);
    let large = dna(6000, 311);
    let mut inputs = String::new();
    // Input order interleaves size classes; shared small documents can merge.
    for (id, sequence) in [&large, &shared, &shared, &large[..2000].to_owned(), &shared]
        .into_iter()
        .enumerate()
    {
        let path = f.write(&format!("stream-{id}.fa"), format!(">doc\n{sequence}\n"));
        inputs.push_str(&format!("doc{id}\t{path}\n"));
    }
    let fof = f.write("stream.fof", inputs);
    let negative = dna(600, 312);
    let reverse = String::from_utf8(zorindex::reverse_complement(shared.as_bytes())).unwrap();
    let reads = format!(">duplicate description\n{shared}\n>negative\n{negative}\n>duplicate\n{reverse}\n>empty\n\n>ambiguous\nNNNN\n");
    let query = f.write(
        "stream-reads.fa.zst",
        zstd::bulk::compress(reads.as_bytes(), 1).unwrap(),
    );
    let modes: &[(&str, &[&str])] = &[
        ("kmers", &[]),
        (
            "minimizers",
            &["--index-mode", "minimizers", "--k", "31", "--window", "19"],
        ),
        (
            "simd-minimizers",
            &[
                "--index-mode",
                "minimizers",
                "--k",
                "31",
                "--window",
                "19",
                "--minimizer-scheme",
                "simd",
            ],
        ),
        (
            "spaced",
            &["--seed-pattern", "#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#"],
        ),
        (
            "spaced-minimizers",
            &[
                "--index-mode",
                "minimizers",
                "--k",
                "31",
                "--window",
                "19",
                "--seed-pattern",
                "#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#",
            ],
        ),
        ("hpc", &["--homopolymer-compression"]),
        (
            "strobes",
            &[
                "--randstrobemers",
                "--strobe-min",
                "9",
                "--strobe-max",
                "80",
            ],
        ),
        (
            "minimizer-strobes",
            &[
                "--index-mode",
                "minimizers",
                "--k",
                "31",
                "--window",
                "19",
                "--randstrobemers",
                "--strobe-min",
                "9",
                "--strobe-max",
                "80",
            ],
        ),
        (
            "modimizers",
            &["--index-mode", "modimizers", "--modimizer-sampling", "7"],
        ),
        ("findere", &["--index-mode", "findere"]),
    ];
    for bits in ["1", "2", "4", "8", "16"] {
        for &(name, mode) in modes {
            let index = f.path(&format!("stream-{bits}-{name}"));
            // Keep genuinely different size classes while still exercising
            // shared colored filters; otherwise this tiny fixture can collapse
            // entirely into one subindex for some seed schemes.
            let mut extra = vec!["--fingerprint-bits", bits, "--merge-max-colors", "2"];
            extra.extend_from_slice(mode);
            build(&fof, &index, &extra);
            let manifest: serde_json::Value =
                serde_json::from_slice(&fs::read(Path::new(&index).join("manifest.json")).unwrap())
                    .unwrap();
            assert!(
                manifest["banks"].as_array().unwrap().len() >= 2,
                "fp{bits} {name}"
            );
            for ignore in [false, true] {
                let mut args = vec!["query", "--index", &index, "--query", &query];
                if ignore {
                    args.push("--ignore-abandoned");
                }
                let streamed = call(&args, "1");
                assert!(
                    streamed.status.success(),
                    "{}",
                    String::from_utf8_lossy(&streamed.stderr)
                );
                assert!(
                    String::from_utf8_lossy(&streamed.stderr).contains("max_loaded_subindexes=1")
                );
                if matches!(
                    name,
                    "kmers"
                        | "spaced"
                        | "hpc"
                        | "minimizers"
                        | "simd-minimizers"
                        | "spaced-minimizers"
                ) {
                    assert!(
                        String::from_utf8_lossy(&streamed.stderr).contains("lookup=record-pruned"),
                        "fp{bits} {name}: direct dense/minimizer query path is required"
                    );
                }
                let streamed = String::from_utf8(streamed.stdout).unwrap();
                let mut full = args.clone();
                full.extend(["--query-memory", "all"]);
                let full = ok(&full, "4");
                // Query statistics and every result are identical; elapsed time differs.
                assert_eq!(
                    streamed
                        .lines()
                        .next()
                        .unwrap()
                        .split("\telapsed_s=")
                        .next(),
                    full.lines().next().unwrap().split("\telapsed_s=").next(),
                    "fp{bits} {name} ignore={ignore}"
                );
                assert_eq!(
                    streamed.lines().skip(1).collect::<Vec<_>>(),
                    full.lines().skip(1).collect::<Vec<_>>(),
                    "fp{bits} {name} ignore={ignore}"
                );
                if name != "findere" {
                    for ratio in ["0", "0.8", "1"] {
                        let mut records = args.clone();
                        records.extend([
                            "--records",
                            "--record-batch-size",
                            "2",
                            "--min-ratio",
                            ratio,
                        ]);
                        let streamed = ok(&records, "4");
                        records.extend(["--query-memory", "all"]);
                        assert_eq!(
                            streamed,
                            ok(&records, "1"),
                            "fp{bits} {name} ratio={ratio} ignore={ignore}"
                        );
                        assert!(fs::read_dir(&f.0).unwrap().all(|entry| !entry
                            .unwrap()
                            .file_name()
                            .to_string_lossy()
                            .starts_with(".zor-query-")));
                    }
                }
            }
        }
    }
}

#[test]
fn sequential_record_queries_keep_results_in_ram_for_one_or_multiple_subindexes() {
    let f = Fixture::new();
    let a = f.write("a.fa", format!(">a\n{}\n", dna(600, 313)));
    let b = f.write("b.fa", format!(">b\n{}\n", dna(6000, 314)));
    let fof = f.write("input.fof", format!("a\t{a}\nb\t{b}\n"));
    let index = f.path("multiple");
    build(&fof, &index, &[]);
    let invalid = f.write("invalid.fa", "ACGT\n");
    let args = ["query", "--index", &index, "--query", &invalid, "--records"];
    assert!(!call(&args, "2").status.success());
    assert!(fs::read_dir(&f.0).unwrap().all(|entry| !entry
        .unwrap()
        .file_name()
        .to_string_lossy()
        .starts_with(".zor-query-")));
    // Invalid input must still fail cleanly without creating any result files.
    let failure = call(
        &["query", "--index", &index, "--query", f.0.to_str().unwrap()],
        "1",
    );
    assert!(!failure.status.success());
    let single = f.path("single");
    build(&fof, &single, &["--max-subindexes", "1"]);
    let nonexistent = f.path("not-a-directory");
    let query = ["query", "--index", &single, "--query", &a, "--records"];
    let actual = ok(&query, "2");
    let expected = ok(
        &query
            .into_iter()
            .chain(["--query-memory", "all"])
            .collect::<Vec<_>>(),
        "1",
    );
    assert_eq!(actual, expected);
    assert!(!Path::new(&nonexistent).exists());
}

#[test]
fn dense_direct_queries_preserve_long_ambiguous_and_repeated_feature_counts() {
    let f = Fixture::new();
    let shared = dna(2400, 7319);
    let other = dna(3200, 6928);
    let a = f.write("dense-a.fa", format!(">a\n{shared}\n"));
    let b = f.write("dense-b.fa", format!(">b\n{other}\n"));
    let fof = f.write("dense.fof", format!("a\t{a}\nb\t{b}\n"));
    let reverse = String::from_utf8(zorindex::reverse_complement(shared.as_bytes())).unwrap();
    let reads = format!(
        ">long\n{}\n>reverse\n{reverse}\n>mixed\n{}NNNN{}\n>short\nACGT\n>empty\n\n>ambiguous\nNNNN\n>negative\n{}\n",
        shared.repeat(16), shared.repeat(2), other.repeat(2), dna(1000, 31267)
    );
    let query = f.write(
        "dense-queries.fa.zst",
        zstd::bulk::compress(reads.as_bytes(), 1).unwrap(),
    );
    for bits in ["4", "8", "16"] {
        let index = f.path(&format!("dense-direct-{bits}"));
        build(
            &fof,
            &index,
            &[
                "--fingerprint-bits",
                bits,
                "--slot-count",
                "8192",
                "--max-subindexes",
                "1",
                "--max-subindex-size",
                "4KiB",
            ],
        );
        for ignore in [false, true] {
            for ratio in ["0", "0.8", "1"] {
                let mut args = vec![
                    "query",
                    "--index",
                    &index,
                    "--query",
                    &query,
                    "--records",
                    "--min-ratio",
                    ratio,
                    "--record-batch-size",
                    "2",
                ];
                if ignore {
                    args.push("--ignore-abandoned");
                }
                let direct = call(&args, "4");
                assert!(
                    direct.status.success(),
                    "{}",
                    String::from_utf8_lossy(&direct.stderr)
                );
                assert!(String::from_utf8_lossy(&direct.stderr).contains("lookup=record-pruned"));
                args.extend(["--query-memory", "all"]);
                let full = call(&args, "1");
                assert!(
                    full.status.success(),
                    "{}",
                    String::from_utf8_lossy(&full.stderr)
                );
                assert_eq!(
                    direct.stdout, full.stdout,
                    "fp{bits} ratio{ratio} ignore{ignore}"
                );
            }
        }
    }
}

#[test]
fn capped_subindexes_route_hashes_and_preserve_membership_across_seed_modes_and_widths() {
    let f = Fixture::new();
    let shared = dna(900, 781);
    let other = dna(1500, 9182);
    let mut entries = String::new();
    for (id, seq) in [&shared, &other, &shared].into_iter().enumerate() {
        let path = f.write(&format!("partition-{id}.fa"), format!(">d{id}\n{seq}\n"));
        entries.push_str(&format!("d{id}\t{path}\n"));
    }
    let fof = f.write("partition.fof", entries);
    let reads = format!(
        ">shared\n{shared}\n>negative\n{}\n>empty\n\n",
        dna(500, 91237)
    );
    let query = f.write(
        "partition-queries.fa.zst",
        zstd::bulk::compress(reads.as_bytes(), 1).unwrap(),
    );
    let positive = f.write("partition-positive.fa", format!(">positive\n{shared}\n"));
    let modes: &[(&str, &[&str])] = &[
        ("kmers", &[]),
        (
            "minimizers",
            &["--index-mode", "minimizers", "--k", "31", "--window", "19"],
        ),
        (
            "spaced",
            &["--seed-pattern", "#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#"],
        ),
        (
            "strobes",
            &[
                "--randstrobemers",
                "--strobe-min",
                "9",
                "--strobe-max",
                "80",
            ],
        ),
        (
            "minimizer-strobes",
            &[
                "--index-mode",
                "minimizers",
                "--k",
                "31",
                "--window",
                "19",
                "--randstrobemers",
                "--strobe-min",
                "9",
                "--strobe-max",
                "80",
            ],
        ),
        ("findere", &["--index-mode", "findere"]),
    ];
    for bits in ["1", "2", "4", "8", "16"] {
        for &(mode, options) in modes {
            let index = f.path(&format!("cap-{bits}-{mode}"));
            let mut extra = vec![
                "--fingerprint-bits",
                bits,
                "--slot-count",
                "8192",
                "--max-subindexes",
                "1",
                "--max-subindex-size",
                "4KiB",
            ];
            extra.extend_from_slice(options);
            build(&fof, &index, &extra);
            let m: serde_json::Value =
                serde_json::from_slice(&fs::read(Path::new(&index).join("manifest.json")).unwrap())
                    .unwrap();
            assert_eq!(m["version"], 9);
            assert_eq!(m["max_subindex_size"], 4096);
            assert_eq!(m["max_subindexes"], 1);
            let banks = m["banks"].as_array().unwrap();
            let count = banks.len();
            assert!(count > 1 && count.is_power_of_two());
            for (part, bank) in banks.iter().enumerate() {
                let slots = bank["slots"].as_u64().unwrap();
                assert!(slots.is_power_of_two());
                assert_eq!(slots * count as u64, 8192);
                assert_eq!(bank["routing"]["partition"], part);
                assert_eq!(bank["routing"]["partitions"], count);
                assert_eq!(bank["routing"]["source_slots"], 8192);
                let raw = zstd::decode_all(
                    fs::File::open(Path::new(&index).join(bank["file"].as_str().unwrap())).unwrap(),
                )
                .unwrap();
                let payload = slots * bank["row_bits"].as_u64().unwrap().div_ceil(8);
                assert_eq!(raw.len() as u64, 40 + payload);
                assert!(payload <= 4096);
            }
            for ignored in [false, true] {
                let mut args = vec!["query", "--index", &index, "--query", &query];
                if ignored {
                    args.push("--ignore-abandoned");
                }
                let a = call(&args, "4");
                assert!(a.status.success(), "{}", String::from_utf8_lossy(&a.stderr));
                let log = String::from_utf8_lossy(&a.stderr);
                assert!(log.contains("input_passes=1\tseed_extraction_passes=1"));
                assert!(log.contains(if ignored {
                    "abandoned_file_scans=0"
                } else {
                    "abandoned_file_scans=1"
                }));
                assert!(log.contains("result_storage=RAM"));
                args.extend(["--query-memory", "all"]);
                let a = String::from_utf8(a.stdout).unwrap();
                let b = ok(&args, "1");
                assert_eq!(
                    a.lines().skip(1).collect::<Vec<_>>(),
                    b.lines().skip(1).collect::<Vec<_>>(),
                    "fp{bits} {mode}"
                );
                if mode != "findere" {
                    for threshold in ["0", "0.8", "1"] {
                        let mut records = vec![
                            "query",
                            "--index",
                            &index,
                            "--query",
                            &query,
                            "--records",
                            "--min-ratio",
                            threshold,
                        ];
                        if ignored {
                            records.push("--ignore-abandoned");
                        }
                        let a = ok(&records, "4");
                        records.extend(["--query-memory", "all"]);
                        assert_eq!(
                            a,
                            ok(&records, "1"),
                            "fp{bits} {mode} threshold={threshold}"
                        );
                    }
                }
            }
            let result = ok(&["query", "--index", &index, "--query", &positive], "4");
            let total = result
                .lines()
                .next()
                .unwrap()
                .split('\t')
                .find_map(|s| s.strip_prefix("total_features="))
                .unwrap()
                .parse::<u64>()
                .unwrap();
            assert!(total > 0);
            assert_eq!(matches(&result, "d0"), total, "fp{bits} {mode}");
            assert_eq!(matches(&result, "d2"), total, "fp{bits} {mode}");
            if bits == "8" && mode == "kmers" {
                // A single seed needs exactly one partition. Delete every
                // other payload and verify the query still resolves fully.
                let one = f.write("single-seed.fa", format!(">one\n{}\n", &shared[..31]));
                let expected = ok(&["query", "--index", &index, "--query", &one], "1");
                let mut key = 0;
                zorindex::scan_canonical_kmers(shared[..31].as_bytes(), 31, &mut |k| key = k);
                let selected =
                    (zorindex::splitmix64(key ^ 69 ^ 0x34B6_1D8F_29A7_C053) as usize) & (count - 1);
                for (part, bank) in banks.iter().enumerate() {
                    if part != selected {
                        fs::remove_file(Path::new(&index).join(bank["file"].as_str().unwrap()))
                            .unwrap();
                    }
                }
                let actual = ok(&["query", "--index", &index, "--query", &one], "4");
                assert_eq!(
                    actual.lines().skip(1).collect::<Vec<_>>(),
                    expected.lines().skip(1).collect::<Vec<_>>()
                );
            }
        }
    }
}

#[test]
fn partition_caps_validate_and_all_construction_storage_modes_agree() {
    let f = Fixture::new();
    let path = f.write("cap.fa", format!(">ref\n{}\n", dna(1300, 187)));
    let fof = f.write("cap.fof", format!("a\t{path}\nb\t{path}\n"));
    let mut expected: Option<serde_json::Value> = None;
    let mut expected_directory = String::new();
    for (name, flags) in [
        ("default", vec![]),
        ("ram", vec!["--low-disk"]),
        ("disk", vec!["--disk-backed"]),
        ("dropped", vec!["--drop-abandoned"]),
    ] {
        let index = f.path(name);
        let mut args = vec![
            "--slot-count",
            "8192",
            "--max-subindex-size",
            "2KiB",
            "--fingerprint-bits",
            "8",
        ];
        args.extend(flags);
        build(&fof, &index, &args);
        let manifest: serde_json::Value =
            serde_json::from_slice(&fs::read(Path::new(&index).join("manifest.json")).unwrap())
                .unwrap();
        assert!(manifest["banks"].as_array().unwrap().len() > 1);
        if name != "dropped" {
            if let Some(expected) = &expected {
                assert_eq!(manifest["banks"], expected["banks"]);
                for bank in manifest["banks"].as_array().unwrap() {
                    let file = bank["file"].as_str().unwrap();
                    assert_eq!(
                        fs::read(Path::new(&index).join(file)).unwrap(),
                        fs::read(Path::new(&expected_directory).join(file)).unwrap()
                    );
                }
                assert_eq!(
                    fs::read(Path::new(&index).join("abandoned.zormap.zst")).unwrap(),
                    fs::read(Path::new(&expected_directory).join("abandoned.zormap.zst")).unwrap()
                );
            } else {
                expected = Some(manifest);
                expected_directory = index.clone();
            }
        }
        let args = ["query", "--index", &index, "--query", &path];
        let a = ok(&args, "4");
        let b = ok(
            &args
                .into_iter()
                .chain(["--query-memory", "all"])
                .collect::<Vec<_>>(),
            "1",
        );
        assert_eq!(
            a.lines().skip(1).collect::<Vec<_>>(),
            b.lines().skip(1).collect::<Vec<_>>()
        );
        assert!(!Path::new(&index).join("partitions").exists());
    }
    for cap in [
        "0",
        "3KiB",
        "1.5GiB",
        "junk",
        "18446744073709551615GiB",
        "1B",
    ] {
        let index = f.path(&format!("invalid-{cap}"));
        let result = call(
            &[
                "build",
                "--fof",
                &fof,
                "--output",
                &index,
                "--k",
                "31",
                "--max-subindex-size",
                cap,
            ],
            "1",
        );
        assert!(!result.status.success(), "cap {cap}");
        assert!(!Path::new(&index).exists());
    }
}

#[test]
fn malformed_hash_partition_metadata_is_rejected_before_reading_payloads() {
    let f = Fixture::new();
    let path = f.write("reference.fa", format!(">ref\n{}\n", dna(1300, 8716)));
    let fof = f.write("input.fof", format!("a\t{path}\nb\t{path}\n"));
    let index = f.path("index");
    build(
        &fof,
        &index,
        &["--slot-count", "8192", "--max-subindex-size", "4KiB"],
    );
    let manifest_path = Path::new(&index).join("manifest.json");
    let original: serde_json::Value =
        serde_json::from_slice(&fs::read(&manifest_path).unwrap()).unwrap();
    let info = ok(&["info", "--index", &index], "1");
    assert_eq!(stat(&info, "size_classes"), 1);
    assert!(stat(&info, "largest_subindex_bytes") <= 4096);
    for corruption in [
        "missing",
        "duplicate",
        "not-power-two",
        "source-slots",
        "member-order",
        "old-version",
        "over-cap",
        "bank-order",
    ] {
        let mut m = original.clone();
        match corruption {
            "missing" => {
                m["banks"].as_array_mut().unwrap().remove(1);
            }
            "duplicate" => m["banks"][1]["routing"]["partition"] = 0.into(),
            "not-power-two" => m["banks"][0]["routing"]["partitions"] = 3.into(),
            "source-slots" => m["banks"][0]["routing"]["source_slots"] = 4096.into(),
            "member-order" => m["banks"][1]["groups"][0]["members"]
                .as_array_mut()
                .unwrap()
                .reverse(),
            "old-version" => m["version"] = 8.into(),
            "over-cap" => m["max_subindex_size"] = 512.into(),
            "bank-order" => m["banks"].as_array_mut().unwrap().swap(0, 1),
            _ => unreachable!(),
        }
        fs::write(&manifest_path, serde_json::to_vec(&m).unwrap()).unwrap();
        let result = call(&["query", "--index", &index, "--query", &path], "1");
        assert!(!result.status.success(), "{corruption}");
    }
    fs::write(manifest_path, serde_json::to_vec(&original).unwrap()).unwrap();
    assert!(call(&["query", "--index", &index, "--query", &path], "1")
        .status
        .success());
}

#[cfg(unix)]
#[test]
fn sequential_record_queries_parse_non_replayable_input_once() {
    use std::io::Write;
    use std::process::Stdio;
    let f = Fixture::new();
    let a = f.write("a.fa", format!(">a\n{}\n", dna(500, 823)));
    let b = f.write("b.fa", format!(">b\n{}\n", dna(6000, 719)));
    let fof = f.write("inputs.fof", format!("a\t{a}\nb\t{b}\n"));
    let index = f.path("index");
    build(&fof, &index, &[]);
    let mut child = Command::new(env!("CARGO_BIN_EXE_zorindex"))
        .args([
            "--threads",
            "2",
            "query",
            "--index",
            &index,
            "--query",
            "/dev/stdin",
            "--records",
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(&fs::read(&a).unwrap())
        .unwrap();
    let result = child.wait_with_output().unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        String::from_utf8(result.stdout).unwrap(),
        ok(
            &[
                "query",
                "--index",
                &index,
                "--query",
                &a,
                "--records",
                "--query-memory",
                "all"
            ],
            "1"
        )
    );
}

#[test]
fn legacy_manifest_without_planner_fields_still_loads() {
    let f = Fixture::new();
    let reference = f.write("a.fa", format!(">a\n{}\n", dna(700, 123)));
    let inputs = f.write("inputs.fof", format!("a\t{reference}\n"));
    let index = f.path("legacy");
    build(&inputs, &index, &[]);
    let path = Path::new(&index).join("manifest.json");
    let mut manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
    for field in [
        "merge_planner",
        "planning_pair_evaluations",
        "merge_reverts_enabled",
        "construction_cell_storage",
        "column_order",
    ] {
        manifest.as_object_mut().unwrap().remove(field);
    }
    fs::write(path, serde_json::to_vec(&manifest).unwrap()).unwrap();
    let info = ok(&["info", "--index", &index], "1");
    assert!(info.contains("merge_planner\tall-pairs"));
    assert!(info.contains("merge_reverts_enabled\ttrue"));
    assert!(info.contains("construction_cell_storage\tin-memory-u64"));
    assert!(info.contains("merge_column_order\twidth"));
    let query = ok(&["query", "--index", &index, "--query", &reference], "4");
    assert_eq!(matches(&query, "a"), 670);
}

#[test]
fn repacking_subindexes_preserves_filter_payloads_and_both_query_modes() {
    use std::io::Read;
    let f = Fixture::new();
    let a = dna(5000, 9921);
    let b = dna(3000, 9922);
    let pa = f.write("a.fa", format!(">a\n{a}\n"));
    let pb = f.write("b.fa", format!(">b\n{b}\n"));
    let fof = f.write("input.fof", format!("a\t{pa}\nb\t{pb}\n"));
    let index = f.path("before");
    let output = f.path("after");
    ok(
        &[
            "build",
            "--fof",
            &fof,
            "--output",
            &index,
            "--index-mode",
            "minimizers",
            "--k",
            "31",
            "--window",
            "19",
            "--min-subindex-documents",
            "1",
            "--hashes",
            "4",
        ],
        "4",
    );
    ok(
        &["repack-subindexes", "--index", &index, "--output", &output],
        "4",
    );
    assert_eq!(
        fs::read(Path::new(&index).join("manifest.json")).unwrap(),
        fs::read(Path::new(&output).join("manifest.json")).unwrap()
    );
    let manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(Path::new(&index).join("manifest.json")).unwrap())
            .unwrap();
    for bank in manifest["banks"].as_array().unwrap() {
        let name = bank["file"].as_str().unwrap();
        let mut before = Vec::new();
        let mut after = Vec::new();
        zorindex::open_index_reader(&Path::new(&index).join(name))
            .unwrap()
            .read_to_end(&mut before)
            .unwrap();
        zorindex::open_index_reader(&Path::new(&output).join(name))
            .unwrap()
            .read_to_end(&mut after)
            .unwrap();
        assert_eq!(before, after);
    }
    let query = f.write(
        "queries.fa",
        format!(
            ">a\n{}\n>negative\n{}\n>empty\nNNN\n",
            &a[50..1050],
            dna(1000, 9221)
        ),
    );
    for ratio in ["0", "0.8", "1"] {
        let expected = ok(
            &[
                "query",
                "--index",
                &index,
                "--query",
                &query,
                "--records",
                "--min-ratio",
                ratio,
                "--query-memory",
                "all",
            ],
            "4",
        );
        for dir in [&index, &output] {
            for mode in ["all", "subindex"] {
                assert_eq!(
                    expected,
                    ok(
                        &[
                            "query",
                            "--index",
                            dir,
                            "--query",
                            &query,
                            "--records",
                            "--min-ratio",
                            ratio,
                            "--query-memory",
                            mode
                        ],
                        "4"
                    )
                );
            }
        }
    }
    assert!(!call(
        &["repack-subindexes", "--index", &index, "--output", &output],
        "4"
    )
    .status
    .success());
}
