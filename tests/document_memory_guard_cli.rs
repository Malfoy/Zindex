//! The single memory budget must fail whole builds rather than exclude datasets.
use std::{
    fs,
    path::PathBuf,
    process::{Command, Output},
    sync::atomic::{AtomicU64, Ordering},
};
static NEXT: AtomicU64 = AtomicU64::new(0);
struct Fixture(PathBuf);
impl Fixture {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "zindex-memory-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn write(&self, name: &str, data: impl AsRef<[u8]>) {
        fs::write(self.0.join(name), data).unwrap();
    }
    fn run(&self, extra: &[&str]) -> Output {
        Command::new(env!("CARGO_BIN_EXE_zorindex"))
            .current_dir(&self.0)
            .args([
                "-t",
                "2",
                "build",
                "-f",
                "inputs.fof",
                "-o",
                "index",
                "-k",
                "31",
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
fn dna(n: usize) -> String {
    (0..n)
        .map(|i| b"ACGT"[(zorindex::splitmix64(i as u64) & 3) as usize] as char)
        .collect()
}
#[test]
fn every_planner_requires_all_datasets_within_the_single_memory_limit() {
    let f = Fixture::new();
    f.write("small.fa", format!(">small\n{}\n", dna(1000)));
    f.write("large.fa", format!(">large\n{}\n", dna(2_000_000)));
    f.write("inputs.fof", "small\tsmall.fa\nlarge\tlarge.fa\n");
    for planner in ["streaming", "full", "naive"] {
        let result = f.run(&["--merge-planner", planner, "--max-ram-gb", "0.06"]);
        assert!(!result.status.success());
        assert!(!f.0.join("index").exists());
        assert!(!f.0.join("index.skipped-documents.json").exists());
    }
}
#[test]
fn compressed_expansion_is_a_failed_build_not_a_partial_collection() {
    let f = Fixture::new();
    let compressed =
        zstd::bulk::compress(format!(">large\n{}\n", "A".repeat(8_000_000)).as_bytes(), 1).unwrap();
    f.write("large.fa.zst", compressed);
    f.write("inputs.fof", "large\tlarge.fa.zst\n");
    let result = f.run(&["--max-ram-gb", "0.06"]);
    assert!(!result.status.success());
    assert!(!f.0.join("index").exists());
}
#[test]
fn corruption_is_never_treated_as_an_excluded_dataset() {
    let f = Fixture::new();
    f.write("bad.fa.zst", [0x28, 0xb5, 0x2f, 0xfd, 1, 2, 3]);
    f.write("inputs.fof", "bad\tbad.fa.zst\n");
    let result = f.run(&["--max-ram-gb", "1"]);
    assert!(!result.status.success());
    assert!(!f.0.join("index").exists());
    assert!(!String::from_utf8_lossy(&result.stderr).contains("status=incomplete"));
}
#[test]
fn sparse_minimizer_input_keeps_every_dataset() {
    let f = Fixture::new();
    f.write("source.fa", format!(">source\n{}\n", dna(600_000)));
    f.write("inputs.fof", "source\tsource.fa\n");
    let result = Command::new(env!("CARGO_BIN_EXE_zorindex"))
        .current_dir(&f.0)
        .args([
            "-t",
            "2",
            "build",
            "-f",
            "inputs.fof",
            "-o",
            "index",
            "-k",
            "49",
            "-i",
            "minimizers",
            "-m",
            "31",
            "-g",
            "0.5",
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(String::from_utf8_lossy(&result.stderr)
        .contains("status=complete\trequested=1\tindexed=1\tskipped=0"));
}
