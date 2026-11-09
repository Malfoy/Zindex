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
        let path = std::env::temp_dir().join(format!(
            "zor-one-bit-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn path(&self, name: &str) -> String {
        self.0.join(name).to_str().unwrap().to_owned()
    }
    fn write(&self, name: &str, value: &str) -> String {
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
fn ok(args: &[&str], threads: &str) -> String {
    let format =
        if args.first() == Some(&"query") && args.iter().any(|a| ["--records", "-r"].contains(a)) {
            vec!["--output-format", "compact"]
        } else {
            vec![]
        };
    let output = Command::new(env!("CARGO_BIN_EXE_zorindex"))
        .args(["--threads", threads])
        .args(args)
        .args(format)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{args:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}
fn stat(output: &str, field: &str) -> u64 {
    output
        .lines()
        .find_map(|l| {
            l.strip_prefix(&format!("{field}\t"))
                .and_then(|s| s.parse().ok())
        })
        .unwrap()
}
fn rows(output: &str) -> Vec<&str> {
    output.lines().filter(|l| !l.starts_with('#')).collect()
}
fn scores(output: &str) -> Vec<u64> {
    rows(output)
        .iter()
        .skip(1)
        .map(|l| l.split('\t').nth(1).unwrap().parse().unwrap())
        .collect()
}
fn dna(size: usize, seed: u64) -> String {
    (0..size)
        .map(|i| b"ACGT"[(zorindex::splitmix64(i as u64 ^ seed) & 3) as usize] as char)
        .collect()
}

#[test]
fn packed_power_banks_query_one_bit_for_all_planners_and_word_tails() {
    let f = Fixture::new();
    let seq = dna(1054, 13);
    let reference = f.write("ref.fa", &format!(">ref\n{seq}\n"));
    for count in [1usize, 9, 65] {
        let inputs = f.write(
            "inputs.fof",
            &(0..count)
                .map(|i| format!("d{i}\t{reference}\n"))
                .collect::<String>(),
        );
        for planner in ["streaming", "full"] {
            let index = f.path(&format!("{planner}-{count}"));
            ok(
                &[
                    "build",
                    "--fof",
                    &inputs,
                    "--output",
                    &index,
                    "--k",
                    "31",
                    "--fingerprint-bits",
                    "1",
                    "--merge-planner",
                    planner,
                ],
                "4",
            );
            let info = ok(&["info", "--index", &index], "1");
            assert_eq!(stat(&info, "fingerprint_bits"), 1);
            assert_eq!(stat(&info, "groups"), count as u64); // Origin bits make these merges unprofitable.
            assert_eq!(stat(&info, "logical_bits"), count as u64 * 1024);
            assert_eq!(stat(&info, "stack_bytes"), count.div_ceil(8) as u64 * 1024);
            assert!(stat(&info, "abandoned_keys") > 0);
            let one = ok(&["query", "--index", &index, "--query", &reference], "1");
            let many = ok(&["query", "--index", &index, "--query", &reference], "4");
            assert_eq!(rows(&one), rows(&many));
            assert_eq!(scores(&one), vec![1024; count]); // Abandoned-key recovery stays additive.
        }
    }
}

#[test]
fn packed_minimizer_and_findere_queries_keep_present_features() {
    let f = Fixture::new();
    let reference = f.write("ref.fa", &format!(">ref\n{}\n", dna(800, 987)));
    let inputs = f.write("inputs.fof", &format!("a\t{reference}\n"));
    for bits in ["1", "2", "4"] {
        for mode in ["minimizers", "findere"] {
            let index = f.path(&format!("{bits}-{mode}-merged"));
            let mut args = vec![
                "build",
                "--fof",
                &inputs,
                "--output",
                &index,
                "--k",
                "49",
                "--fingerprint-bits",
                bits,
                "--index-mode",
                mode,
            ];
            if mode == "findere" {
                args.extend(["--findere-z", "20"]);
            } else {
                args.extend(["--minimizer-size", "21"]);
            }

            ok(&args, "4");
            let one = ok(&["query", "--index", &index, "--query", &reference], "1");
            let four = ok(&["query", "--index", &index, "--query", &reference], "4");
            assert_eq!(rows(&one), rows(&four));
            let record = rows(&one)[1];
            assert!(record.split('\t').nth(1).unwrap().parse::<u64>().unwrap() > 0);
            assert_eq!(record.split('\t').nth(2), Some("1.00000000"));
        }
    }
}
