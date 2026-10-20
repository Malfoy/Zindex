use std::{
    fs,
    path::PathBuf,
    process::Command,
    sync::atomic::{AtomicU64, Ordering},
};
use zorindex::{scan_index_features, FeatureConfig, IndexMode, MinimizerScheme};

static NEXT: AtomicU64 = AtomicU64::new(0);
struct Fixture(PathBuf);
impl Fixture {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "zor-extraction-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn path(&self, name: &str) -> String {
        self.0.join(name).to_str().unwrap().into()
    }
    fn write(&self, name: &str, data: impl AsRef<[u8]>) -> String {
        let path = self.path(name);
        fs::write(&path, data).unwrap();
        path
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn ok(args: &[&str], threads: &str) -> String {
    let output = Command::new(env!("CARGO_BIN_EXE_zorindex"))
        .args(["--threads", threads])
        .args(args)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{args:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}
fn dna(n: usize) -> Vec<u8> {
    (0..n)
        .map(|i| b"ACGT"[(zorindex::splitmix64(i as u64) & 3) as usize])
        .collect()
}
fn expanded(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .enumerate()
        .flat_map(|(i, &b)| std::iter::repeat_n(b, 1 + i % 5))
        .collect()
}
fn fasta(seq: &[u8]) -> Vec<u8> {
    let mut result = b">seq\n".to_vec();
    result.extend_from_slice(seq);
    result.push(b'\n');
    result
}
fn score(output: &str) -> (u64, u64) {
    let line = output.lines().find(|line| line.starts_with("a\t")).unwrap();
    let fields: Vec<_> = line.split('\t').collect();
    let total = output
        .lines()
        .next()
        .unwrap()
        .split('\t')
        .find_map(|field| field.strip_prefix("total_features="))
        .unwrap()
        .parse()
        .unwrap();
    (fields[1].parse().unwrap(), total)
}

#[test]
fn persisted_schemes_and_homopolymer_setting_match_build_query_and_plain_append() {
    let f = Fixture::new();
    let original = dna(1800);
    let compressed = zorindex::homopolymer_compress(&original);
    let reference = f.write("expanded.fa", fasta(&expanded(&compressed)));
    let query = f.write("compressed.fa", fasta(&compressed));
    let rc: Vec<_> = compressed
        .iter()
        .rev()
        .map(|&b| b"TGCA"[zorindex::encode_base(b).unwrap() as usize])
        .collect();
    let reversed = f.write("rc.fa", fasta(&rc));
    let fof = f.write("inputs.fof", format!("a\t{reference}\n"));
    for plain in [false, true] {
        for mode in [
            IndexMode::Kmers,
            IndexMode::Minimizers,
            IndexMode::Modimizers,
            IndexMode::Findere,
        ] {
            for scheme in [MinimizerScheme::Simd, MinimizerScheme::OpenClosedMod] {
                let index = f.path(&format!(
                    "index-{plain}-{}-{}",
                    mode.label(),
                    scheme.label()
                ));
                let mut args = vec![
                    "build",
                    "--fof",
                    &fof,
                    "--output",
                    &index,
                    "--k",
                    "31",
                    "--index-mode",
                    mode.label(),
                    "--minimizer-size",
                    "21",
                    "--minimizer-scheme",
                    scheme.label(),
                    "--homopolymer-compression",
                    "--fingerprint-bits",
                    "16",
                ];
                // Plain filters can abandon keys without recovery. Leave room
                // for this test to assert extraction matching independently.
                if plain {
                    args.extend(["--plain", "--slot-scale", "8"]);
                }
                ok(&args, "4");
                let info = ok(&["info", "--index", &index], "1");
                assert!(info.contains(&format!("minimizer_scheme\t{}", scheme.label())));
                assert!(info.contains("homopolymer_compression\ttrue"));
                let config = FeatureConfig::new(mode, 21, 16)
                    .with_minimizer_scheme(scheme)
                    .with_homopolymer_compression(true);
                let mut expected = 0;
                if mode == IndexMode::Findere {
                    expected = (compressed.len() - 31 + 1) as u64;
                } else {
                    scan_index_features(&compressed, 31, 69, config, &mut |_| expected += 1);
                }
                for path in [&reference, &query, &reversed] {
                    for threads in ["1", "4"] {
                        let output = ok(&["query", "--index", &index, "--query", path], threads);
                        assert_eq!(
                            score(&output),
                            (expected, expected),
                            "{plain} {mode:?} {scheme:?}: {output}"
                        );
                    }
                }
                if plain {
                    let update = f.write("update.fof", format!("b\t{query}\n"));
                    ok(&["append", "--index", &index, "--fof", &update], "4");
                    let output = ok(&["query", "--index", &index, "--query", &reference], "4");
                    assert_eq!(score(&output), (expected, expected));
                    assert!(output
                        .lines()
                        .any(|line| line.starts_with(&format!("b\t{expected}\t1.00000000\t"))));
                }
            }
        }
    }
}

#[test]
fn default_is_open_closed_and_legacy_metadata_keeps_simd() {
    let f = Fixture::new();
    let reference = f.write("a.fa", fasta(&dna(1900)));
    let fof = f.write("inputs.fof", format!("a\t{reference}\n"));
    let index = f.path("default");
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
            "49",
            "--minimizer-size",
            "31",
        ],
        "4",
    );
    assert!(ok(&["info", "--index", &index], "1").contains("minimizer_scheme\topen-closed-mod"));
    for plain in [false, true] {
        let index = f.path(&format!("legacy-{plain}"));
        let mut args = vec![
            "build",
            "--fof",
            &fof,
            "--output",
            &index,
            "--index-mode",
            "minimizers",
            "--k",
            "49",
            "--minimizer-size",
            "31",
            "--minimizer-scheme",
            "simd",
        ];
        if plain {
            args.push("--plain");
        }
        ok(&args, "4");
        let expected = ok(&["query", "--index", &index, "--query", &reference], "1");
        if plain {
            let mut bytes = zstd::stream::decode_all(fs::File::open(&index).unwrap()).unwrap();
            assert_eq!(u32::from_le_bytes(bytes[8..12].try_into().unwrap()), 10);
            // v8 has extraction settings but no seed-pattern/randstrobe tags.
            bytes[8..12].copy_from_slice(&8u32.to_le_bytes());
            bytes.drain(36..38);
            fs::write(
                &index,
                zstd::stream::encode_all(bytes.as_slice(), 0).unwrap(),
            )
            .unwrap();
            let v8 = ok(&["query", "--index", &index, "--query", &reference], "1");
            assert_eq!(score(&expected), score(&v8));
            // v7 has no scheme or compression byte after the feature header.
            bytes[8..12].copy_from_slice(&7u32.to_le_bytes());
            bytes.drain(34..36);
            fs::write(
                &index,
                zstd::stream::encode_all(bytes.as_slice(), 0).unwrap(),
            )
            .unwrap();
        } else {
            let path = PathBuf::from(&index).join("manifest.json");
            let mut manifest: serde_json::Value =
                serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
            manifest["version"] = 3.into();
            manifest.as_object_mut().unwrap().remove("minimizer_scheme");
            manifest
                .as_object_mut()
                .unwrap()
                .remove("homopolymer_compression");
            fs::write(path, serde_json::to_vec(&manifest).unwrap()).unwrap();
        }
        let info = ok(&["info", "--index", &index], "1");
        assert!(info.contains("minimizer_scheme\tsimd"));
        assert!(info.contains("homopolymer_compression\tfalse"));
        let actual = ok(&["query", "--index", &index, "--query", &reference], "1");
        assert_eq!(score(&expected), score(&actual));
    }
}

#[test]
fn homopolymer_compression_preserves_breaks_and_feature_settings() {
    assert_eq!(
        zorindex::homopolymer_compress(b"AaAACccNNNTttAA"),
        b"ACNNNTA"
    );
    for mode in [
        IndexMode::Kmers,
        IndexMode::Minimizers,
        IndexMode::Modimizers,
        IndexMode::Findere,
    ] {
        for scheme in [MinimizerScheme::Simd, MinimizerScheme::OpenClosedMod] {
            let config = FeatureConfig::new(mode, 3, 1)
                .with_findere_z(2)
                .with_minimizer_scheme(scheme)
                .with_homopolymer_compression(true);
            let input = b"AAAACCCCGGGGTTTTAAAAANNNAAAACCCCGGGGTTTTAAAAA";
            let mut actual = Vec::new();
            let mut expected = Vec::new();
            scan_index_features(input, 5, 69, config, &mut |key| actual.push(key));
            for _ in 0..2 {
                scan_index_features(
                    b"ACGTA",
                    5,
                    69,
                    config.with_homopolymer_compression(false),
                    &mut |key| expected.push(key),
                );
            }
            assert_eq!(actual, expected, "{mode:?} {scheme:?}");
        }
    }
}

#[test]
fn anchor_settings_roundtrip_queries_reverse_complements_and_plain_updates() {
    use zorindex::RandstrobeConfig;
    let f = Fixture::new();
    let seq = dna(1800);
    let reference = f.write("anchors.fa", fasta(&seq));
    let reversed = f.write(
        "anchors-rc.fa",
        fasta(
            &seq.iter()
                .rev()
                .map(|&b| b"TGCA"[zorindex::encode_base(b).unwrap() as usize])
                .collect::<Vec<_>>(),
        ),
    );
    let fof = f.write("anchors.fof", format!("a\t{reference}\n"));
    for plain in [false, true] {
        for mode in [IndexMode::Kmers, IndexMode::Minimizers] {
            for scheme in [MinimizerScheme::Simd, MinimizerScheme::OpenClosedMod] {
                for choice in 0..3 {
                    for hpc in [false, true] {
                        let k = if mode == IndexMode::Kmers { 9 } else { 31 };
                        let k_text = k.to_string();
                        let index = f.path(&format!(
                            "anchors-{plain}-{mode:?}-{scheme:?}-{choice}-{hpc}"
                        ));
                        let masked = choice != 1;
                        let paired = choice != 0;
                        let pattern = "##__#__##";
                        let mut args = vec![
                            "build",
                            "--fof",
                            &fof,
                            "--output",
                            &index,
                            "--k",
                            &k_text,
                            "--index-mode",
                            mode.label(),
                            "--minimizer-size",
                            "9",
                            "--minimizer-scheme",
                            scheme.label(),
                            "--fingerprint-bits",
                            "16",
                        ];
                        if masked {
                            args.extend(["--seed-pattern", pattern]);
                        }
                        if paired {
                            args.extend([
                                "--randstrobemers",
                                "--strobe-min",
                                "9",
                                "--strobe-max",
                                "40",
                            ]);
                        }
                        if hpc {
                            args.push("--homopolymer-compression");
                        }
                        if plain {
                            args.extend(["--plain", "--slot-scale", "8"]);
                        }
                        ok(&args, "4");
                        let info = ok(&["info", "--index", &index], "1");
                        assert!(info.contains(&format!(
                            "seed_pattern\t{}",
                            if masked { pattern } else { "none" }
                        )));
                        assert!(info.contains(&format!("randstrobemers\t{paired}")));
                        if paired {
                            assert!(info.contains("strobe_min\t9\nstrobe_max\t40"));
                        }
                        let config = FeatureConfig::new(mode, 9, 16)
                            .with_minimizer_scheme(scheme)
                            .with_homopolymer_compression(hpc)
                            .with_seed_pattern(masked.then(|| pattern.parse().unwrap()))
                            .with_randstrobemers(paired.then_some(RandstrobeConfig {
                                min_distance: 9,
                                max_distance: 40,
                            }));
                        let mut expected = 0;
                        scan_index_features(&seq, k, 69, config, &mut |_| expected += 1);
                        assert!(expected > 0);
                        for query in [&reference, &reversed] {
                            for threads in ["1", "4"] {
                                assert_eq!(
                                    score(&ok(
                                        &["query", "--index", &index, "--query", query],
                                        threads
                                    )),
                                    (expected, expected)
                                );
                            }
                        }
                        if plain {
                            let update = f.write("anchors-update.fof", format!("b\t{reversed}\n"));
                            ok(&["append", "--index", &index, "--fof", &update], "4");
                            let result =
                                ok(&["query", "--index", &index, "--query", &reference], "4");
                            assert_eq!(score(&result), (expected, expected));
                            assert!(result
                                .lines()
                                .any(|line| line
                                    .starts_with(&format!("b\t{expected}\t1.00000000\t"))));
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn seed_pattern_ignores_only_masked_bases_and_is_reused_from_index() {
    let f = Fixture::new();
    let reference = f.write("masked.fa", ">a\nACATACCG\n");
    let query = f.write("mutated.fa", ">a\nACGTAGCG\n");
    let fof = f.write("masked.fof", format!("a\t{reference}\n"));
    for plain in [false, true] {
        let index = f.path(&format!("masked-{plain}"));
        let mut args = vec![
            "build",
            "--fof",
            &fof,
            "--output",
            &index,
            "--k",
            "8",
            "--seed-pattern",
            "##_##_##",
            "--fingerprint-bits",
            "16",
        ];
        if plain {
            args.extend(["--plain", "--slot-scale", "8"]);
        }
        ok(&args, "1");
        assert_eq!(
            score(&ok(&["query", "--index", &index, "--query", &query], "4")),
            (1, 1)
        );
    }
}

#[test]
fn invalid_anchor_arguments_fail_before_construction() {
    let f = Fixture::new();
    let reference = f.write("invalid.fa", fasta(&dna(500)));
    let fof = f.write("invalid.fof", format!("a\t{reference}\n"));
    let cases: &[(&[&str], &str)] = &[
        (&["--k", "8", "--seed-pattern", "##_##_#"], "symmetric"),
        (&["--k", "8", "--seed-pattern", "#_#"], "must match k"),
        (
            &[
                "--index-mode",
                "minimizers",
                "--k",
                "31",
                "--minimizer-size",
                "9",
                "--seed-pattern",
                "#_#",
            ],
            "minimizer size m",
        ),
        (&["--seed-pattern", "#x#"], "only '#' and '_'"),
        (&["--seed-pattern", "___"], "at least one '#'"),
        (
            &["--randstrobemers", "--strobe-min", "0"],
            "greater than zero",
        ),
        (
            &[
                "--randstrobemers",
                "--strobe-min",
                "20",
                "--strobe-max",
                "10",
            ],
            "must not exceed",
        ),
        (&["--strobe-min", "5"], "--randstrobemers"),
        (
            &["--index-mode", "findere", "--randstrobemers"],
            "kmers or minimizers",
        ),
        (
            &[
                "--index-mode",
                "modimizers",
                "--seed-pattern",
                "###############################",
            ],
            "kmers or minimizers",
        ),
    ];
    for (i, &(extra, message)) in cases.iter().enumerate() {
        let index = f.path(&format!("invalid-{i}"));
        let mut command = Command::new(env!("CARGO_BIN_EXE_zorindex"));
        command.args(["build", "--fof", &fof, "--output", &index]);
        if !extra.contains(&"--k") {
            command.args(["--k", "31"]);
        }
        let output = command.args(extra).output().unwrap();
        assert!(!output.status.success(), "{extra:?}");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(stderr.contains(message), "{extra:?}: {stderr}");
        assert!(!PathBuf::from(index).exists());
    }
}

#[test]
fn downstream_strobes_keep_strands_separate_and_reject_custom_legacy_indexes() {
    let f = Fixture::new();
    let seq = dna(1800);
    let rc = zorindex::reverse_complement(&seq);
    let reference = f.write("reference.fa", fasta(&seq));
    let mixed = f.write(
        "mixed.fa",
        format!(
            ">forward\n{}\n>reverse\n{}\n",
            String::from_utf8_lossy(&seq),
            String::from_utf8_lossy(&rc)
        ),
    );
    let reverse_fastq = f.write(
        "reverse.fastq.zst",
        zstd::stream::encode_all(
            format!(
                "@reverse\n{}\n+\n{}\n",
                String::from_utf8_lossy(&rc),
                "I".repeat(rc.len())
            )
            .as_bytes(),
            0,
        )
        .unwrap(),
    );
    let fof = f.write("input.fof", format!("a\t{reference}\n"));
    for plain in [false, true] {
        let index = f.path(&format!("strands-{plain}"));
        let mut args = vec![
            "build",
            "--fof",
            &fof,
            "--output",
            &index,
            "--k",
            "11",
            "--randstrobemers",
            "--strobe-min",
            "15",
            "--strobe-max",
            "50",
            "--fingerprint-bits",
            "16",
        ];
        if plain {
            args.extend(["--plain", "--slot-scale", "8"]);
        }
        ok(&args, "4");
        let total = (seq.len() - 11 + 1) as u64;
        assert_eq!(
            score(&ok(&["query", "--index", &index, "--query", &mixed], "4")),
            (2 * total, 2 * total)
        );
        assert_eq!(
            score(&ok(
                &["query", "--index", &index, "--query", &reverse_fastq],
                "4"
            )),
            (total, total)
        );
        if !plain {
            for threads in ["1", "4"] {
                let records = ok(
                    &[
                        "query",
                        "--index",
                        &index,
                        "--query",
                        &mixed,
                        "--records",
                        "--min-ratio",
                        "1",
                    ],
                    threads,
                );
                assert_eq!(
                    records.lines().skip(1).collect::<Vec<_>>(),
                    [
                        format!("forward\t{total}\t0:{total}"),
                        format!("reverse\t{total}\t0:{total}"),
                    ]
                );
            }
            let path = PathBuf::from(&index).join("manifest.json");
            let mut manifest: serde_json::Value =
                serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
            assert_eq!(manifest["version"], 6);
            assert_eq!(manifest["colors"][0]["occurrences"], total);
            assert_eq!(manifest["randstrobe_scheme"], zorindex::RANDSTROBE_SCHEME);
            manifest["version"] = 5.into();
            manifest
                .as_object_mut()
                .unwrap()
                .remove("randstrobe_scheme");
            fs::write(path, serde_json::to_vec(&manifest).unwrap()).unwrap();
        } else {
            let mut bytes = zstd::stream::decode_all(fs::File::open(&index).unwrap()).unwrap();
            bytes[8..12].copy_from_slice(&9u32.to_le_bytes());
            fs::write(
                &index,
                zstd::stream::encode_all(bytes.as_slice(), 0).unwrap(),
            )
            .unwrap();
        }
        let output = Command::new(env!("CARGO_BIN_EXE_zorindex"))
            .args(["query", "--index", &index, "--query", &mixed])
            .output()
            .unwrap();
        assert!(!output.status.success());
        assert!(String::from_utf8_lossy(&output.stderr).contains("rebuild"));
    }
}

#[test]
fn randstrobe_distance_defaults_follow_strobe_length() {
    let f = Fixture::new();
    let reference = f.write("defaults.fa", fasta(&dna(600)));
    let fof = f.write("defaults.fof", format!("a\t{reference}\n"));
    for mode in ["kmers", "minimizers"] {
        let index = f.path(mode);
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
                mode,
                "--minimizer-size",
                "9",
                "--randstrobemers",
            ],
            "4",
        );
        let expected = if mode == "kmers" {
            "strobe_min\t31\nstrobe_max\t93"
        } else {
            "strobe_min\t9\nstrobe_max\t27"
        };
        assert!(ok(&["info", "--index", &index], "1").contains(expected));
    }
}

#[test]
fn invalid_persisted_anchor_settings_are_rejected() {
    let f = Fixture::new();
    let reference = f.write("corrupt.fa", fasta(&dna(500)));
    let fof = f.write("corrupt.fof", format!("a\t{reference}\n"));
    for plain in [false, true] {
        let index = f.path(&format!("corrupt-{plain}"));
        let mut args = vec![
            "build",
            "--fof",
            &fof,
            "--output",
            &index,
            "--k",
            "9",
            "--seed-pattern",
            "##__#__##",
            "--randstrobemers",
        ];
        if plain {
            args.extend(["--plain", "--slot-scale", "8"]);
        }
        ok(&args, "1");
        if plain {
            let mut bytes = zstd::stream::decode_all(fs::File::open(&index).unwrap()).unwrap();
            // Pattern width follows the scheme and compression bytes.
            bytes[36] = 32;
            fs::write(
                &index,
                zstd::stream::encode_all(bytes.as_slice(), 0).unwrap(),
            )
            .unwrap();
        } else {
            let path = PathBuf::from(&index).join("manifest.json");
            let mut manifest: serde_json::Value =
                serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
            manifest["randstrobemers"]["min_distance"] = 0.into();
            fs::write(path, serde_json::to_vec(&manifest).unwrap()).unwrap();
        }
        let result = Command::new(env!("CARGO_BIN_EXE_zorindex"))
            .args(["info", "--index", &index])
            .output()
            .unwrap();
        assert!(!result.status.success());
        let error = String::from_utf8_lossy(&result.stderr);
        assert!(
            error.contains(if plain {
                "invalid seed pattern width"
            } else {
                "minimum distance must be greater than zero"
            }),
            "{error}"
        );
    }
}

#[test]
#[ignore = "requires ZOR_OCM_BEFORE and ZOR_OCM_GENOMES_FOF for before/after CLI timings"]
fn real_collection_before_after_benchmark() {
    let before =
        std::env::var("ZOR_OCM_BEFORE").expect("path to saved pre-optimization release binary");
    let fof = std::env::var("ZOR_OCM_GENOMES_FOF").expect("collection FOF");
    let inputs = zorindex::read_fof(std::path::Path::new(&fof)).unwrap();
    let f = Fixture::new();
    let normalized = inputs
        .iter()
        .map(|input| {
            format!(
                "{}\t{}\n",
                input.name,
                fs::canonicalize(&input.path).unwrap().display()
            )
        })
        .collect::<String>();
    let fof = f.write("normalized.fof", normalized);
    let mut manifests = Vec::new();
    let mut builds = Vec::new();
    let binaries = [before.as_str(), env!("CARGO_BIN_EXE_zorindex")];
    let indexes = [f.path("before"), f.path("after")];
    for (binary, index) in binaries.into_iter().zip(&indexes) {
        let start = std::time::Instant::now();
        let output = Command::new(binary)
            .args([
                "--threads",
                "4",
                "build",
                "--fof",
                &fof,
                "--output",
                index,
                "--k",
                "49",
                "--minimizer-size",
                "31",
                "--index-mode",
                "minimizers",
                "--fingerprint-bits",
                "16",
            ])
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        builds.push(start.elapsed().as_secs_f64());
        let manifest: serde_json::Value =
            serde_json::from_slice(&fs::read(PathBuf::from(index).join("manifest.json")).unwrap())
                .unwrap();
        manifests.push(manifest);
    }
    // Timings are the only manifest differences. The complete stored filter
    // bytes, origin map, dataset metadata and merge decisions must be unchanged.
    for field in [
        "colors",
        "banks",
        "planned_groups",
        "predicted_logical_bits",
        "singleton_logical_bits",
        "planning_pair_evaluations",
    ] {
        assert_eq!(manifests[0][field], manifests[1][field], "{field}");
    }
    for bank in manifests[0]["banks"].as_array().unwrap() {
        let file = bank["file"].as_str().unwrap();
        assert_eq!(
            fs::read(PathBuf::from(&indexes[0]).join(file)).unwrap(),
            fs::read(PathBuf::from(&indexes[1]).join(file)).unwrap(),
            "bank {file}"
        );
    }
    let map = "abandoned.zormap.zst";
    assert_eq!(
        fs::read(PathBuf::from(&indexes[0]).join(map)).unwrap(),
        fs::read(PathBuf::from(&indexes[1]).join(map)).unwrap()
    );
    eprintln!(
        "collection={} build_before_s={:.6} build_after_s={:.6} stored_bank_bytes_identical=true",
        inputs.len(),
        builds[0],
        builds[1]
    );
    let mut queries = vec![("genome", fs::canonicalize(&inputs[1].path).unwrap())];
    if let Ok(path) = std::env::var("ZOR_OCM_QUERY") {
        queries.push(("short-queries", PathBuf::from(path)));
    }
    for (label, query) in queries {
        for threads in ["1", "4"] {
            let mut times: [Vec<f64>; 2] = std::array::from_fn(|_| Vec::new());
            let mut expected = None;
            for repeat in 0..7 {
                for step in 0..2 {
                    let id = (step + repeat) % 2;
                    let start = std::time::Instant::now();
                    let result = Command::new(binaries[id])
                        .args([
                            "--threads",
                            threads,
                            "query",
                            "--index",
                            &indexes[id],
                            "--query",
                            query.to_str().unwrap(),
                        ])
                        .output()
                        .unwrap();
                    assert!(
                        result.status.success(),
                        "{}",
                        String::from_utf8_lossy(&result.stderr)
                    );
                    let elapsed = start.elapsed().as_secs_f64();
                    let text = String::from_utf8(result.stdout).unwrap();
                    let first = text.lines().next().unwrap();
                    // Header differs only in elapsed time; counts/probes match.
                    let header = first.split("\telapsed_s=").next().unwrap().to_string();
                    let scores = text.lines().skip(1).collect::<Vec<_>>().join("\n");
                    assert_eq!(
                        expected.get_or_insert((header.clone(), scores.clone())),
                        &(header, scores)
                    );
                    times[id].push(elapsed);
                }
            }
            for (id, mut values) in times.into_iter().enumerate() {
                values.sort_by(f64::total_cmp);
                eprintln!(
                    "query={label} threads={threads} path={} median_wall_s={:.6} min_wall_s={:.6}",
                    if id == 0 { "before" } else { "after" },
                    values[3],
                    values[0]
                );
            }
            eprintln!("verified_query_header={}", expected.unwrap().0);
        }
    }
}
