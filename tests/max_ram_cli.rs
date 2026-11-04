//! Fixed total RSS policy: units, validation, measured stop and unchanged data.
use std::{
    fs,
    path::{Path, PathBuf},
    process::{Command, Output},
    sync::atomic::{AtomicU64, Ordering},
};
static NEXT: AtomicU64 = AtomicU64::new(0);
struct Fixture(PathBuf);
impl Fixture {
    fn new(n: usize) -> Self {
        let path = std::env::temp_dir().join(format!(
            "zor-max-ram-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        let sequence: String = (0..n)
            .map(|i| b"ACGT"[(zorindex::splitmix64(i as u64) & 3) as usize] as char)
            .collect();
        fs::write(path.join("genome.fa"), format!(">g\n{sequence}\n")).unwrap();
        fs::write(
            path.join("query.fa"),
            format!(">q\n{}\n", &sequence[..1000.min(n)]),
        )
        .unwrap();
        fs::write(
            path.join("fof.tsv"),
            format!("g\t{}\n", path.join("genome.fa").display()),
        )
        .unwrap();
        Self(path)
    }
    fn build(&self, name: &str, extra: &[&str]) -> Output {
        Command::new(env!("CARGO_BIN_EXE_zorindex"))
            .args([
                "--threads",
                &std::thread::available_parallelism()
                    .unwrap()
                    .get()
                    .to_string(),
                "build",
                "--fof",
                self.0.join("fof.tsv").to_str().unwrap(),
                "--output",
                self.0.join(name).to_str().unwrap(),
                "--k",
                "31",
                "--index-mode",
                "minimizers",
                "--window",
                "19",
                "--minimizer-scheme",
                "open-closed-mod",
                "--hashes",
                "4",
                "--fingerprint-bits",
                "8",
                "--min-subindex-documents",
                "1",
                "--skip-oversized-documents=false",
            ])
            .args(extra)
            .output()
            .unwrap()
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn files(path: &Path, base: &Path, out: &mut std::collections::BTreeMap<PathBuf, Vec<u8>>) {
    for entry in fs::read_dir(path).unwrap() {
        let path = entry.unwrap().path();
        if path.is_dir() {
            files(&path, base, out);
        } else if path.file_name().unwrap() != "manifest.json" {
            out.insert(
                path.strip_prefix(base).unwrap().to_owned(),
                fs::read(&path).unwrap(),
            );
        }
    }
}
#[test]
fn invalid_total_budget_and_conflicting_budgets_are_rejected() {
    let f = Fixture::new(2000);
    for value in ["0", "-1", "NaN", "inf", "1e100", "0.0000000001"] {
        let out = f.build("invalid", &["--max-ram-gb", value]);
        assert!(!out.status.success(), "accepted {value}");
        assert!(!f.0.join("invalid").exists());
    }
    let out = f.build(
        "invalid",
        &["--max-ram-gb", "500", "--build-memory-gib", "1"],
    );
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("cannot be used with"));
}
#[test]
fn measured_rss_above_small_total_fails_without_publishing() {
    let f = Fixture::new(2000);
    let out = f.build("capped", &["--max-ram-gb", "0.000001"]);
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    let evidence = stderr
        .lines()
        .find_map(|line| line.strip_prefix("memory-limit-exceeded\t"))
        .expect("actual RSS evidence");
    let evidence: serde_json::Value = serde_json::from_str(evidence).unwrap();
    assert_eq!(evidence["kind"], "actual_rss_limit");
    assert_eq!(evidence["limit_bytes"], 1000);
    assert!(evidence["rss_bytes"].as_u64().unwrap() > 1000);
    assert!(!f.0.join("capped").exists());
}
#[test]
fn watchdog_stops_actual_growth_during_native_construction() {
    let f = Fixture::new(10_000_000);
    let out = f.build("growing", &["--max-ram-gb", "0.02"]);
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("sketching"),
        "must reach native work before exceeding cap: {stderr}"
    );
    assert_eq!(
        out.status.code(),
        Some(86),
        "background watchdog must stop growth: {stderr}"
    );
    let evidence = stderr
        .lines()
        .find_map(|line| line.strip_prefix("memory-limit-exceeded\t"))
        .expect("actual RSS evidence");
    let evidence: serde_json::Value = serde_json::from_str(evidence).unwrap();
    assert_eq!(evidence["kind"], "actual_rss_limit");
    assert_eq!(evidence["limit_bytes"], 20_000_000);
    assert!(evidence["rss_bytes"].as_u64().unwrap() > 20_000_000);
    assert!(!f.0.join("growing").exists());
}
#[test]
fn fixed_budget_keeps_native_filter_bytes_and_queries_identical() {
    let f = Fixture::new(50_000);
    let baseline = f.build("baseline", &[]);
    assert!(
        baseline.status.success(),
        "{}",
        String::from_utf8_lossy(&baseline.stderr)
    );
    let fixed = f.build("fixed", &["--max-ram-gb", "500"]);
    assert!(
        fixed.status.success(),
        "{}",
        String::from_utf8_lossy(&fixed.stderr)
    );
    assert!(String::from_utf8_lossy(&fixed.stderr).contains("max_ram_bytes=500000000000"));
    let (mut a, mut b) = (Default::default(), Default::default());
    files(&f.0.join("baseline"), &f.0.join("baseline"), &mut a);
    files(&f.0.join("fixed"), &f.0.join("fixed"), &mut b);
    assert_eq!(a, b);
    let query = |name| {
        Command::new(env!("CARGO_BIN_EXE_zorindex"))
            .args([
                "--threads",
                &std::thread::available_parallelism()
                    .unwrap()
                    .get()
                    .to_string(),
                "query",
                "--index",
                f.0.join(name).to_str().unwrap(),
                "--query",
                f.0.join("query.fa").to_str().unwrap(),
                "--records",
                "--output-format",
                "compact",
                "--min-ratio",
                "1",
            ])
            .output()
            .unwrap()
    };
    let (a, b) = (query("baseline"), query("fixed"));
    assert!(a.status.success() && b.status.success());
    assert_eq!(a.stdout, b.stdout);
    assert!(String::from_utf8_lossy(&a.stdout).contains("0:"));
}
