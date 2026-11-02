//! Memory-guard fixtures exercise exclusions, retained membership and failures.
use std::{
    fs,
    io::Write,
    path::PathBuf,
    process::{Command, Output},
    sync::atomic::{AtomicU64, Ordering},
};
static NEXT: AtomicU64 = AtomicU64::new(0);
struct Fixture(PathBuf);
impl Fixture {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "zor-memory-guard-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&p).unwrap();
        Self(p)
    }
    fn path(&self, name: &str) -> String {
        self.0.join(name).to_string_lossy().into_owned()
    }
    fn write(&self, name: &str, text: impl AsRef<[u8]>) -> String {
        let p = self.path(name);
        fs::write(&p, text).unwrap();
        p
    }
    fn call(&self, fof: &str, name: &str, extra: &[&str]) -> Output {
        let threads = std::thread::available_parallelism()
            .unwrap()
            .get()
            .to_string();
        Command::new(env!("CARGO_BIN_EXE_zorindex"))
            .args([
                "--threads",
                &threads,
                "build",
                "--fof",
                fof,
                "--output",
                &self.path(name),
                "--k",
                "31",
                "--hashes",
                "8",
                "--fingerprint-bits",
                "8",
                "--min-subindex-documents",
                "1",
                "--max-document-memory-gib",
                "0.01",
            ])
            .args(extra)
            .output()
            .unwrap()
    }
    fn report(&self, name: &str) -> serde_json::Value {
        serde_json::from_slice(&fs::read(self.path(name)).unwrap()).unwrap()
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn dna(n: usize) -> String {
    (0..n)
        .map(|i| b"ACGT"[(zorindex::splitmix64(i as u64) & 3) as usize] as char)
        .collect()
}

#[test]
fn oversized_document_is_skipped_without_losing_included_document_identity() {
    let f = Fixture::new();
    let seq = dna(3000);
    let small = f.write("small.fa", format!(">small\n{seq}\n"));
    let large = f.write("large.fa", format!(">large\n{}\n", dna(600_000)));
    let fof = f.write(
        "inputs.tsv",
        format!("before\t{small}\noversized\t{large}\nafter\t{small}\n"),
    );
    let output = f.call(&fof, "index", &[]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let manifest = f.report("index/manifest.json");
    assert_eq!(manifest["input_document_count"], 3);
    assert_eq!(manifest["colors"].as_array().unwrap().len(), 2);
    assert_eq!(manifest["source_document_ids"], serde_json::json!([0, 2]));
    assert_eq!(manifest["skipped_documents"][0]["source_document_id"], 1);
    assert_eq!(manifest["skipped_documents"][0]["name"], "oversized");
    assert_eq!(
        f.report("index.skipped-documents.json")["status"],
        "indexed"
    );
    assert_eq!(
        f.report("index/skipped-documents.json")["indexed_documents"],
        2
    );
    let query = f.write("query.fa", format!(">q0\n{}\n", &seq[..1000]));
    let out = Command::new(env!("CARGO_BIN_EXE_zorindex"))
        .args([
            "query",
            "--index",
            &f.path("index"),
            "--query",
            &query,
            "--records",
            "--min-ratio",
            "1",
        ])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let text = String::from_utf8(out.stdout).unwrap();
    assert!(text.contains("0:970"), "{text}");
    assert!(text.contains("1:970"), "{text}");
}

#[test]
fn all_skipped_leaves_a_report_without_publishing_an_index() {
    let f = Fixture::new();
    let file = f.write("large.fa", format!(">large\n{}\n", dna(600_000)));
    let fof = f.write("fof.txt", format!("large\t{file}\n"));
    let out = f.call(&fof, "index", &[]);
    assert!(!out.status.success());
    assert!(!f.0.join("index").exists());
    assert_eq!(
        f.report("index.skipped-documents.json")["status"],
        "all_skipped"
    );
    assert!(!fs::read_dir(&f.0).unwrap().any(|p| p
        .unwrap()
        .file_name()
        .to_string_lossy()
        .starts_with(".zor-pow2-building-")));
}

#[test]
fn compressed_size_underestimate_is_caught_after_exact_extraction() {
    let f = Fixture::new();
    let mut repeated = String::new();
    for _ in 0..2 {
        repeated.push_str(&dna(100_000));
    }
    let mut encoder = deko::write::AnyEncoder::new(
        Vec::new(),
        deko::Format::Gz,
        deko::write::Compression::Default,
    )
    .unwrap();
    encoder
        .write_all(format!(">large\n{repeated}\n").as_bytes())
        .unwrap();
    let file = f.write("large.fa.gz", encoder.finish().unwrap());
    let small = f.write("small.fa", format!(">small\n{}\n", dna(1000)));
    let fof = f.write("fof.txt", format!("large\t{file}\nsmall\t{small}\n"));
    let out = f.call(&fof, "index", &[]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let report = f.report("index.skipped-documents.json");
    assert_eq!(report["skipped_documents"][0]["stage"], "singleton-graph");
}

#[test]
fn compressed_expansion_is_bounded_and_not_accepted_as_truncated_input() {
    let f = Fixture::new();
    let mut encoder = deko::write::AnyEncoder::new(
        Vec::new(),
        deko::Format::Gz,
        deko::write::Compression::Default,
    )
    .unwrap();
    encoder
        .write_all(format!(">large\n{}\n", "A".repeat(4_000_000)).as_bytes())
        .unwrap();
    let file = f.write("large.fa.gz", encoder.finish().unwrap());
    let fof = f.write("fof.txt", format!("large\t{file}\n"));
    let out = f.call(&fof, "index", &[]);
    assert!(!out.status.success());
    assert_eq!(
        f.report("index.skipped-documents.json")["skipped_documents"][0]["stage"],
        "decompressed-input"
    );
    assert!(!f.0.join("index").exists());
}

#[test]
fn corruption_is_an_error_and_not_a_memory_skip() {
    let f = Fixture::new();
    let file = f.write("bad.fa.gz", b"\x1f\x8bgarbage");
    let fof = f.write("fof.txt", format!("bad\t{file}\n"));
    let out = f.call(&fof, "index", &[]);
    assert!(!out.status.success());
    assert!(!f.0.join("index.skipped-documents.json").exists());
}

#[test]
fn oversized_merge_is_rejected_during_planning_and_all_documents_still_query() {
    let f = Fixture::new();
    let seq = dna(50_000);
    let input = f.write("shared.fa", format!(">shared\n{seq}\n"));
    let fof = f.write(
        "four.tsv",
        (0..4)
            .map(|i| format!("doc{i}\t{input}\n"))
            .collect::<String>(),
    );
    let threads = std::thread::available_parallelism()
        .unwrap()
        .get()
        .to_string();
    for mode in ["streaming", "streaming-unpruned", "all-pairs"] {
        let name = format!("index-{mode}");
        let out = Command::new(env!("CARGO_BIN_EXE_zorindex"))
            .args([
                "--threads",
                &threads,
                "build",
                "--fof",
                &fof,
                "--output",
                &f.path(&name),
                "--k",
                "31",
                "--hashes",
                "8",
                "--fingerprint-bits",
                "8",
                "--min-subindex-documents",
                "1",
                "--merge-plan-margin",
                "0",
                "--max-document-memory-gib",
                "0.01953125",
                "--merge-planner",
                mode,
            ])
            .output()
            .unwrap();
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        let manifest = f.report(&format!("{name}/manifest.json"));
        assert_eq!(manifest["colors"].as_array().unwrap().len(), 4);
        assert_eq!(manifest["planned_groups"], 2);
        assert_eq!(manifest["memory_guard_split_groups"], 0);
        assert_eq!(manifest["merge_memory_budget_bytes"], 20 << 20);
        let query = f.write("query.fa", format!(">q0\n{}\n", &seq[..1000]));
        let out = Command::new(env!("CARGO_BIN_EXE_zorindex"))
            .args([
                "--threads",
                &threads,
                "query",
                "--index",
                &f.path(&name),
                "--query",
                &query,
                "--records",
                "--min-ratio",
                "1",
            ])
            .output()
            .unwrap();
        assert!(out.status.success());
        let text = String::from_utf8(out.stdout).unwrap();
        for id in 0..4 {
            assert!(text.contains(&format!("{id}:970")), "{mode}: {text}");
        }
    }
}

#[test]
fn file_size_hint_does_not_reject_sparse_minimizer_input() {
    let f = Fixture::new();
    let input = f.write("genome.fa", format!(">large\n{}\n", dna(600_000)));
    // File*24 exceeds10MiB, but OCM extraction plus its exact graph fit.
    let fof = f.write("fof.tsv", format!("large\t{input}\n"));
    let threads = std::thread::available_parallelism()
        .unwrap()
        .get()
        .to_string();
    let out = Command::new(env!("CARGO_BIN_EXE_zorindex"))
        .args([
            "--threads",
            &threads,
            "build",
            "--fof",
            &fof,
            "--output",
            &f.path("index"),
            "--k",
            "49",
            "--index-mode",
            "minimizers",
            "--minimizer-size",
            "31",
            "--minimizer-scheme",
            "open-closed-mod",
            "--hashes",
            "4",
            "--fingerprint-bits",
            "8",
            "--min-subindex-documents",
            "1",
            "--max-document-memory-gib",
            "0.01",
        ])
        .env("ZORINDEX_MEMORY_TRACE", "1")
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let manifest = f.report("index/manifest.json");
    assert_eq!(manifest["colors"].as_array().unwrap().len(), 1);
    let trace = String::from_utf8(out.stderr).unwrap();
    assert!(trace.contains("document-extracted"));
    assert!(trace.contains("graph-finished"));
}

#[test]
fn admission_uses_resident_rss_instead_of_capacity_bound() {
    let f = Fixture::new();
    let seq = dna(150_000);
    let input = f.write("genome.fa", format!(">g\n{seq}\n"));
    let fof = f.write("fof.tsv", format!("g\t{input}\n"));
    let threads = std::thread::available_parallelism()
        .unwrap()
        .get()
        .to_string();
    // 18MiB fits the resident prediction plus its10% margin, whereas the old
    // all-capacities estimate would reject this exact construction.
    let output = Command::new(env!("CARGO_BIN_EXE_zorindex"))
        .args([
            "--threads",
            &threads,
            "build",
            "--fof",
            &fof,
            "--output",
            &f.path("index"),
            "--k",
            "31",
            "--hashes",
            "4",
            "--fingerprint-bits",
            "8",
            "--min-subindex-documents",
            "1",
            "--max-document-memory-gib",
            "0.017578125",
        ])
        .env("ZORINDEX_MEMORY_TRACE", "1")
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let manifest = f.report("index/manifest.json");
    assert_eq!(manifest["colors"].as_array().unwrap().len(), 1);
    let stderr = String::from_utf8(output.stderr).unwrap();
    let events: Vec<serde_json::Value> = stderr
        .lines()
        .filter_map(|line| {
            line.strip_prefix("memory-trace\t")
                .map(|s| serde_json::from_str(s).unwrap())
        })
        .collect();
    let event = events
        .iter()
        .find(|v| v["stage"] == "document-extracted")
        .unwrap();
    let d = &event["detail"];
    assert_eq!(d["prediction_model"], "compact-resident-v2");
    assert!(d["admission_workspace_bytes"].as_u64().unwrap() < 18 << 20);
    assert!(d["conservative_capacity_bound_bytes"].as_u64().unwrap() > 18 << 20);
    assert!(events.iter().any(|v| v["stage"] == "graph-finished"));
}
