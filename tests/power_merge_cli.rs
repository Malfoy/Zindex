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
    ] {
        if !extra.contains(&flag) {
            args.extend([flag, value]);
        }
    }
    args.extend(extra);
    ok(&args, "4");
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
        "most-deg2",
        "--tie-scan",
        "8",
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
            "49",
            "--index-mode",
            "minimizers",
            "--minimizer-size",
            "31",
            "--fingerprint-bits",
            "16",
        ];
        args.extend(settings);
        ok(&args, "4");
        let manifest: serde_json::Value =
            serde_json::from_slice(&fs::read(Path::new(&index).join("manifest.json")).unwrap())
                .unwrap();
        assert_eq!(manifest["merge_planner"], "streaming");
        assert_eq!(manifest["column_order"], "compression");
        assert_eq!(
            manifest["construction_cell_storage"],
            "two-pass-memory-packed-no-scratch"
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
    assert!(info.contains("construction_cell_storage\ttwo-pass-memory-packed-no-scratch"));
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
    assert_eq!(manifest["version"], 3);
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
    for hashes in [2, 3, 4, 5, 8, 32] {
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
        build(&inputs, &index, &["--index-mode", mode]);
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
fn two_pass_and_disk_backed_build_identical_filters_and_queries() {
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
        for fp in ["1", "8", "16"] {
            for order in ["width", "compression"] {
                let mut outputs = Vec::new();
                let mut manifests = Vec::new();
                let mut results = Vec::new();
                for disk in [false, true] {
                    let path = f.path(&format!("{planner}-{fp}-{order}-{disk}"));
                    let mut args = vec![
                        "--fingerprint-bits",
                        fp,
                        "--merge-planner",
                        planner,
                        "--merge-column-order",
                        order,
                    ];
                    if disk {
                        args.push("--disk-backed");
                    }
                    build(&fof, &path, &args);
                    let mut json: serde_json::Value = serde_json::from_slice(
                        &fs::read(Path::new(&path).join("manifest.json")).unwrap(),
                    )
                    .unwrap();
                    assert_eq!(
                        json["construction_cell_storage"],
                        if disk {
                            "temporary-disk-u64-tiled-packing"
                        } else {
                            "two-pass-memory-packed-no-scratch"
                        }
                    );
                    for key in [
                        "construction_cell_storage",
                        "sketch_seconds",
                        "plan_seconds",
                        "build_seconds",
                    ] {
                        json.as_object_mut().unwrap().remove(key);
                    }
                    manifests.push(json);
                    let output = ok(
                        &["query", "--index", &path, "--query", &query],
                        if disk { "1" } else { "4" },
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
                assert_eq!(manifests[0], manifests[1]);
                assert_eq!(results[0], results[1]);
                for entry in fs::read_dir(&outputs[0]).unwrap() {
                    let name = entry.unwrap().file_name();
                    if name == "manifest.json" {
                        continue;
                    }
                    assert_eq!(
                        fs::read(Path::new(&outputs[0]).join(&name)).unwrap(),
                        fs::read(Path::new(&outputs[1]).join(&name)).unwrap(),
                        "{name:?}"
                    );
                }
            }
        }
    }
}

#[test]
fn scratch_free_mode_preserves_all_feature_extractors() {
    let f = Fixture::new();
    let sequence = dna(500, 42);
    let a = f.write("a.fa", format!(">a\n{sequence}\n"));
    let fof = f.write("inputs.fof", format!("a\t{a}\nb\t{a}\n"));
    let query = f.write("query.fa", format!(">r\n{sequence}\n"));
    for (mode, settings) in [
        (
            "minimizers",
            vec![
                "--k",
                "49",
                "--index-mode",
                "minimizers",
                "--minimizer-size",
                "31",
            ],
        ),
        (
            "findere",
            vec!["--k", "31", "--index-mode", "findere", "--findere-z", "10"],
        ),
        ("modimizers", vec!["--index-mode", "modimizers"]),
        ("no-recovery", vec!["--drop-abandoned"]),
    ] {
        let memory = f.path(&format!("{mode}-memory"));
        let disk = f.path(&format!("{mode}-disk"));
        build(&fof, &memory, &settings);
        let mut disk_settings = settings.clone();
        disk_settings.push("--disk-backed");
        build(&fof, &disk, &disk_settings);
        let one = ok(&["query", "--index", &memory, "--query", &query], "4");
        let two = ok(&["query", "--index", &disk, "--query", &query], "1");
        assert_eq!(
            one.lines().skip(1).collect::<Vec<_>>(),
            two.lines().skip(1).collect::<Vec<_>>(),
            "{mode}"
        );
        for entry in fs::read_dir(&memory).unwrap() {
            let name = entry.unwrap().file_name();
            if name == "manifest.json" {
                continue;
            }
            assert_eq!(
                fs::read(Path::new(&memory).join(&name)).unwrap(),
                fs::read(Path::new(&disk).join(&name)).unwrap(),
                "{mode} {name:?}"
            );
        }
    }
    let invalid = f.path("invalid");
    assert!(!call(
        &[
            "build",
            "--plain",
            "--disk-backed",
            "--fof",
            &fof,
            "--output",
            &invalid,
            "--k",
            "31"
        ],
        "1"
    )
    .status
    .success());
    assert!(!Path::new(&invalid).exists());
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
