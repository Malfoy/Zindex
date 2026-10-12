use deko::{
    write::{AnyEncoder, Compression},
    Format,
};
use std::{fs, io::Write, path::PathBuf, process::Command};

#[test]
fn compressed_build_and_query_preserve_memberships() {
    let root = std::env::temp_dir().join(format!("zor-compression-test-{}", std::process::id()));
    fs::create_dir_all(&root).unwrap();
    let sequence: String = (0..1800)
        .map(|i| b"ACGT"[(zorindex::splitmix64(i ^ 42) & 3) as usize] as char)
        .collect();
    let fasta = format!(
        ">first\n{}\n>second\n{}\n",
        &sequence[..1000],
        &sequence[1000..]
    );
    let threads = std::thread::available_parallelism()
        .unwrap()
        .get()
        .to_string();
    let run = |binary: &str, args: Vec<String>| {
        let output = Command::new(binary)
            .args(["--threads", &threads])
            .args(&args)
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{args:?}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        String::from_utf8(output.stdout)
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with('#'))
            .collect::<Vec<_>>()
            .join("\n")
    };
    for binary in [
        env!("CARGO_BIN_EXE_bloomindex"),
        env!("CARGO_BIN_EXE_zorindex"),
    ] {
        for mode in ["kmers", "findere"] {
            let mut expected = None;
            for (suffix, format) in [
                ("fasta", Format::Verbatim),
                ("fasta.gz", Format::Gz),
                ("fasta.zst", Format::Zstd),
                ("fasta.xz", Format::Xz),
                ("fasta.bz2", Format::Bz),
                ("fasta.bz", Format::Bz),
            ] {
                let input = root.join(format!("input.{suffix}"));
                let mut encoder =
                    AnyEncoder::new(Vec::new(), format, Compression::Default).unwrap();
                encoder.write_all(fasta.as_bytes()).unwrap();
                fs::write(&input, encoder.finish().unwrap()).unwrap();
                let fof = root.join("inputs.tsv");
                fs::write(&fof, format!("reference\t{}\n", input.display())).unwrap();
                let index = root.join(format!(
                    "{}-{mode}-{suffix}.index",
                    PathBuf::from(binary).file_name().unwrap().to_str().unwrap()
                ));
                let mut build = vec![
                    "build".into(),
                    "--fof".into(),
                    fof.display().to_string(),
                    "--output".into(),
                    index.display().to_string(),
                    "--k".into(),
                    "31".into(),
                    "--index-mode".into(),
                    mode.into(),
                    "--hashes".into(),
                    "8".into(),
                ];
                if mode == "findere" {
                    build.extend(["--findere-z".into(), "2".into()]);
                }
                if binary.ends_with("zorindex") {
                    build.extend(["--fingerprint-bits".into(), "8".into()]);
                } else {
                    build.extend(["--false-positive-rate".into(), "0.00390625".into()]);
                }
                run(binary, build);
                let query = vec![
                    "query".into(),
                    "--index".into(),
                    index.display().to_string(),
                    "--query".into(),
                    input.display().to_string(),
                    "--min-ratio".into(),
                    "0".into(),
                ];
                let actual = run(binary, query);
                if let Some(ref expected) = expected {
                    assert_eq!(expected, &actual, "{binary} {mode} {suffix}");
                } else {
                    expected = Some(actual);
                }
            }
        }
    }
    fs::remove_dir_all(root).unwrap();
}
