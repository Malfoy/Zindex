use std::{
    fs,
    path::PathBuf,
    process::{Command, Output},
};

struct Fixture(PathBuf);
impl Fixture {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!("zindex-user-output-{}", std::process::id()));
        fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn path(&self, name: &str) -> String {
        self.0.join(name).to_str().unwrap().to_owned()
    }
    fn write(&self, name: &str, text: impl AsRef<[u8]>) -> String {
        let path = self.path(name);
        fs::write(&path, text).unwrap();
        path
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn call(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_zorindex"))
        .args(["-t", "2"])
        .args(args)
        .output()
        .unwrap()
}
fn ok(args: &[&str]) -> Output {
    let result = call(args);
    assert!(
        result.status.success(),
        "{args:?}: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    result
}

#[test]
fn named_rows_match_compact_counts_and_preserve_unmatched_queries_in_every_query_path() {
    let f = Fixture::new();
    let sequence: String = (0..2000)
        .map(|i| b"ACGT"[(zorindex::splitmix64(i) & 3) as usize] as char)
        .collect();
    let source = f.write("source.fa", format!(">s\n{sequence}\n"));
    let fof = f.write(
        "samples.fof",
        format!("sample one\t{source}\nsample one\t{source}\n"),
    );
    let queries = f.write(
        "queries.fa",
        format!(
            ">source\n{}\n>negative\n{}\n>short\nA\n>ambiguous\nNNNNNN\n",
            &sequence[..500],
            "A".repeat(150)
        ),
    );
    for (mode, settings) in [
        ("kmers", vec![]),
        ("minimizers", vec!["-i", "minimizers", "-m", "31"]),
        ("strobes", vec!["--randstrobemers"]),
        (
            "sampled-strobes",
            vec!["-i", "minimizers", "-m", "31", "--randstrobemers"],
        ),
    ] {
        let index = f.path(mode);
        let mut args = vec![
            "build",
            "-f",
            &fof,
            "-o",
            &index,
            "-k",
            if settings.contains(&"minimizers") {
                "49"
            } else {
                "31"
            },
            "-b",
            "16",
            "-D",
            "1",
        ];
        args.extend(settings);
        let result = ok(&args);
        assert!(String::from_utf8_lossy(&result.stderr)
            .contains("collection\tstatus=complete\trequested=2\tindexed=2\tskipped=0"));
        let info = String::from_utf8(ok(&["info", "-i", &index]).stdout).unwrap();
        assert!(info.contains("collection_status\tcomplete\nindexed_documents\t2\n"));
        for memory in ["all", "subindex"] {
            for threshold in ["0", "0.8"] {
                let args = [
                    "query", "-i", &index, "-q", &queries, "-r", "-M", memory, "-a", threshold,
                    "-B", "1",
                ];
                let named = String::from_utf8(ok(&args).stdout).unwrap();
                assert!(named.starts_with(
                    "query_id\tdataset_id\tdataset_name\tmatches\tqueried_seeds\tratio\n"
                ));
                let compact_args: Vec<_> = args.into_iter().chain(["-F", "compact"]).collect();
                let compact = String::from_utf8(ok(&compact_args).stdout).unwrap();
                let mut named_rows = named.lines().skip(1);
                for line in compact.lines().skip(1) {
                    let fields: Vec<_> = line.split('\t').collect();
                    let total: u64 = fields[1].parse().unwrap();
                    let hits: Vec<_> = fields[2].split(',').filter(|s| !s.is_empty()).collect();
                    if hits.is_empty() {
                        assert_eq!(
                            named_rows.next().unwrap(),
                            format!("{}\t.\t.\t0\t{total}\t0.00000000", fields[0])
                        );
                    } else {
                        for hit in hits {
                            let (id, count) = hit.split_once(':').unwrap();
                            let count: u64 = count.parse().unwrap();
                            let ratio = if total == 0 {
                                0.0
                            } else {
                                count as f64 / total as f64
                            };
                            assert_eq!(
                                named_rows.next().unwrap(),
                                format!(
                                    "{}\t{id}\tsample one\t{count}\t{total}\t{ratio:.8}",
                                    fields[0]
                                )
                            );
                        }
                    }
                }
                assert!(named_rows.next().is_none());
            }
        }
        let aggregate = ok(&["query", "-i", &index, "-q", &queries]);
        assert!(aggregate
            .stdout
            .starts_with(b"color\tmatches\tratio\tunique_kmers\n"));
        assert!(String::from_utf8_lossy(&aggregate.stderr).contains("#query\t"));
    }
}
