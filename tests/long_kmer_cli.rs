use std::{
    fs,
    path::PathBuf,
    process::Command,
    sync::atomic::{AtomicU64, Ordering},
};

static NEXT: AtomicU64 = AtomicU64::new(0);
struct Fixture(PathBuf);
impl Fixture {
    fn new() -> Self {
        let directory = std::env::temp_dir().join(format!(
            "zindex-long-k-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&directory).unwrap();
        Self(directory)
    }
    fn run(&self, args: &[&str]) -> std::process::Output {
        let result = Command::new(env!("CARGO_BIN_EXE_zorindex"))
            .current_dir(&self.0)
            .args(["-t", "2"])
            .args(args)
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{args:?}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        result
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
fn long_kmers_preserve_source_matches_across_query_paths_storage_and_fusion() {
    let f = Fixture::new();
    let reference = dna(1800);
    fs::write(f.0.join("source.fa"), format!(">reference\n{reference}\n")).unwrap();
    fs::write(f.0.join("sources.fof"), "source\tsource.fa\n").unwrap();
    fs::write(
        f.0.join("queries.fa"),
        format!(
            ">forward\n{reference}\n>reverse\n{}\n>short\nACGT\n",
            String::from_utf8(zorindex::reverse_complement(reference.as_bytes())).unwrap()
        ),
    )
    .unwrap();
    for (k, storage) in [
        (32, "disk"),
        (33, "cached"),
        (49, "memory"),
        (64, "disk"),
        (65, "cached"),
        (255, "memory"),
        (256, "disk"),
        (511, "cached"),
        (1025, "memory"),
    ] {
        let index = format!("index-{k}");
        f.run(&[
            "build",
            "-f",
            "sources.fof",
            "-k",
            &k.to_string(),
            "-o",
            &index,
            "--build-storage",
            storage,
            "-b",
            "16",
        ]);
        for memory in ["subindex", "all"] {
            let result = f.run(&[
                "query",
                "-i",
                &index,
                "-q",
                "queries.fa",
                "-r",
                "-a",
                "1",
                "-M",
                memory,
            ]);
            let result = String::from_utf8(result.stdout).unwrap();
            let rows: Vec<_> = result
                .lines()
                .skip(1)
                .map(|line| line.split('\t').collect::<Vec<_>>())
                .collect();
            let expected = (reference.len() - k + 1).to_string();
            assert_eq!(rows[0][3], expected);
            assert_eq!(rows[0][4], expected);
            assert_eq!(rows[1][3], expected);
            assert_eq!(rows[1][4], expected);
            assert_eq!(rows[2][1], ".");
        }
        f.run(&["fusion", "-i", &index, &index, "-o", "fused"]);
        let result = f.run(&["query", "-i", "fused", "-q", "queries.fa", "-r", "-a", "1"]);
        assert_eq!(
            String::from_utf8(result.stdout)
                .unwrap()
                .lines()
                .filter(|line| line.starts_with("forward\t"))
                .count(),
            2
        );
        fs::remove_dir_all(f.0.join("fused")).unwrap();
    }
}

#[test]
fn omitted_output_creates_distinct_timestamped_build_and_fusion_directories() {
    let f = Fixture::new();
    fs::write(f.0.join("source.fa"), format!(">reference\n{}\n", dna(300))).unwrap();
    fs::write(f.0.join("sources.fof"), "source\tsource.fa\n").unwrap();
    for _ in 0..2 {
        let output = f.run(&["build", "-f", "sources.fof", "-k", "65"]);
        assert!(String::from_utf8_lossy(&output.stderr).contains("Output index: myindex-"));
    }
    let directories = || {
        let mut paths: Vec<_> = fs::read_dir(&f.0)
            .unwrap()
            .map(|entry| entry.unwrap().path())
            .filter(|p| {
                p.is_dir()
                    && p.file_name()
                        .unwrap()
                        .to_string_lossy()
                        .starts_with("myindex-")
            })
            .collect();
        paths.sort();
        paths
    };
    let inputs = directories();
    assert_eq!(inputs.len(), 2);
    f.run(&[
        "fusion",
        "-i",
        inputs[0].to_str().unwrap(),
        inputs[1].to_str().unwrap(),
    ]);
    assert_eq!(directories().len(), 3);
}
