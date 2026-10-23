use anyhow::{ensure, Context, Result};
use clap::Parser;
use minimizer_fasta::{generate, Config, MinimizerScheme, Stats};
use rayon::prelude::*;
use serde::Serialize;
use std::{
    fs::{self, File, OpenOptions},
    io::{self, BufRead, BufReader, BufWriter, Read, Write},
    path::{Path, PathBuf},
    time::Instant,
};

#[derive(Debug, Parser)]
#[command(
    version,
    about = "Parallel streaming canonical SIMD/open-closed-mod minimizer FASTA generator"
)]
struct Args {
    /// Input FASTA, optionally gzip/zstd/xz/bzip2 compressed; '-' reads stdin.
    #[arg(conflicts_with = "fof")]
    input: Option<PathBuf>,
    /// Output FASTA; .gz enables gzip; '-' writes stdout. Never overwrites.
    #[arg(conflicts_with = "output_dir")]
    output: Option<PathBuf>,
    /// Collection: one path or name<TAB>path per line; relative paths use cwd.
    #[arg(long, requires = "output_dir", conflicts_with = "input")]
    fof: Option<PathBuf>,
    /// Collection output directory; numbered FASTAs plus minimizers.fof.
    #[arg(long, requires = "fof", conflicts_with = "output")]
    output_dir: Option<PathBuf>,
    /// Length of each selected seed (window span is k + window - 1).
    #[arg(short, long, default_value_t = 31)]
    k: usize,
    /// Number of candidate k-mers in each minimizer window.
    #[arg(short, long, alias = "w", default_value_t = 19)]
    window: usize,
    /// Selection scheme; simd preserves legacy output, open-closed-mod matches Zindex.
    #[arg(long, value_enum, default_value = "simd")]
    minimizer_scheme: MinimizerScheme,
    /// Open-closed-mod hash seed; ignored by the legacy SIMD selector.
    #[arg(long, default_value_t = 69)]
    seed: u64,
    /// Rayon worker count; defaults to all CPUs available to this process.
    #[arg(short, long)]
    threads: Option<usize>,
    /// Compress output with gzip (also inferred from .gz/.gz.partial).
    #[arg(long)]
    gzip: bool,
    /// Gzip compression level, 0..=9. Independent chunks compress in parallel.
    #[arg(long, default_value_t = 1)]
    gzip_level: u32,
    /// Owned seed-start positions per work chunk; context is added automatically.
    #[arg(long, default_value_t = 262_144)]
    chunk_bases: usize,
}

fn output_is_gzip(path: &Path, force: bool) -> bool {
    force
        || path
            .to_str()
            .is_some_and(|p| p.ends_with(".gz") || p.ends_with(".gz.partial"))
}

fn open_input(path: &Path) -> Result<Box<dyn Read + Send>> {
    let input: Box<dyn Read + Send> = if path == Path::new("-") {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(path).with_context(|| format!("opening {}", path.display()))?)
    };
    // generate() performs magic-based decoding, including concatenated members.
    Ok(Box::new(BufReader::with_capacity(1024 * 1024, input)))
}

fn create_output(path: &Path) -> Result<Box<dyn Write + Send>> {
    let output: Box<dyn Write + Send> = if path == Path::new("-") {
        Box::new(io::stdout())
    } else {
        Box::new(
            OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(path)
                .with_context(|| {
                    format!(
                        "creating {} (outputs are never overwritten)",
                        path.display()
                    )
                })?,
        )
    };
    Ok(Box::new(BufWriter::with_capacity(1024 * 1024, output)))
}

fn convert(input: &Path, output: &Path, config: Config) -> Result<Stats> {
    let reader = open_input(input)?;
    let writer = create_output(output)?;
    generate(reader, writer, config).with_context(|| format!("converting {}", input.display()))
}

fn collection(args: &Args, config: Config, threads: usize) -> Result<(Stats, usize, PathBuf)> {
    let directory = args.output_dir.as_ref().unwrap();
    let mut files = Vec::new();
    for line in BufReader::new(File::open(args.fof.as_ref().unwrap())?).lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let (name, input) = if let Some((name, path)) = line.split_once('\t') {
            ensure!(
                !name.is_empty() && !path.is_empty() && !path.contains('\t'),
                "expected name<TAB>path in FOF"
            );
            (name.to_owned(), PathBuf::from(path))
        } else {
            (format!("{:08}", files.len()), PathBuf::from(line))
        };
        ensure!(input.is_file(), "input is not a file: {}", input.display());
        let suffix = if args.gzip {
            ".minimizers.fa.gz"
        } else {
            ".minimizers.fa"
        };
        let output = directory.join(format!("{:08}{suffix}", files.len()));
        ensure!(
            !output.exists(),
            "output already exists: {}",
            output.display()
        );
        files.push((name, input, output));
    }
    ensure!(!files.is_empty(), "empty FOF");
    let manifest = directory.join("minimizers.fof");
    ensure!(
        !manifest.exists(),
        "output manifest already exists: {}",
        manifest.display()
    );
    fs::create_dir_all(directory)?;
    // Avoid threads^2 in-flight chunk buffers when many files run concurrently.
    // Nested Rayon work uses the same pool, never a pool per genome.
    let config = Config {
        batch_chunks: threads.div_ceil(files.len()).max(1),
        ..config
    };
    let results: Vec<Stats> = files
        .par_iter()
        .map(|(_, input, output)| convert(input, output, config))
        .collect::<Result<_>>()?;
    let mut total = Stats::default();
    for stats in results {
        total += stats;
    }
    let mut out = BufWriter::new(
        OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(&manifest)?,
    );
    for (name, _, output) in &files {
        // Absolute paths make the generated FOF independent of the next command's cwd.
        writeln!(out, "{name}\t{}", fs::canonicalize(output)?.display())?;
    }
    out.flush()?;
    Ok((total, files.len(), manifest))
}

#[derive(Serialize)]
struct Summary {
    k: usize,
    w: usize,
    minimizer_scheme: MinimizerScheme,
    seed: u64,
    canonical: bool,
    #[serde(flatten)]
    stats: Stats,
    files: usize,
    threads: usize,
    simd: &'static str,
    wall_seconds: f64,
    short_run_policy: &'static str,
    #[serde(skip_serializing_if = "Option::is_none")]
    output_fof: Option<PathBuf>,
}

fn main() -> Result<()> {
    let args = Args::parse();
    ensure!(
        args.fof.is_some() || (args.input.is_some() && args.output.is_some()),
        "provide INPUT OUTPUT or --fof FILE --output-dir DIR"
    );
    let threads = args
        .threads
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, usize::from));
    ensure!(threads > 0, "threads must be positive");
    ensure!(args.gzip_level <= 9, "gzip level must be in 0..=9");
    let gzip = args
        .output
        .as_deref()
        .map_or(args.gzip, |p| output_is_gzip(p, args.gzip));
    let config = Config {
        k: args.k,
        w: args.window,
        minimizer_scheme: args.minimizer_scheme,
        seed: args.seed,
        chunk_bases: args.chunk_bases,
        batch_chunks: threads,
        gzip_level: gzip.then_some(args.gzip_level),
    };
    config.validate()?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()?;
    let start = Instant::now();
    let (stats, files, output_fof) = pool.install(|| -> Result<_> {
        if args.fof.is_some() {
            let (stats, files, manifest) = collection(&args, config, threads)?;
            Ok((stats, files, Some(manifest)))
        } else {
            Ok((
                convert(
                    args.input.as_ref().unwrap(),
                    args.output.as_ref().unwrap(),
                    config,
                )?,
                1,
                None,
            ))
        }
    })?;
    let ocm_avx2 = {
        #[cfg(target_arch = "x86_64")]
        {
            std::is_x86_feature_detected!("avx2")
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            false
        }
    };
    let simd = if args.minimizer_scheme == MinimizerScheme::OpenClosedMod {
        if ocm_avx2 {
            "AVX2"
        } else {
            "scalar"
        }
    } else if cfg!(target_feature = "avx2") {
        "AVX2"
    } else if cfg!(target_feature = "neon") {
        "NEON"
    } else {
        "scalar"
    };
    let summary = Summary {
        k: args.k,
        w: args.window,
        minimizer_scheme: args.minimizer_scheme,
        seed: args.seed,
        canonical: true,
        stats,
        files,
        threads,
        simd,
        wall_seconds: start.elapsed().as_secs_f64(),
        short_run_policy: "no minimizer without a full k+w-1 base window",
        output_fof,
    };
    if args.output.as_deref() == Some(Path::new("-")) {
        eprintln!("{}", serde_json::to_string(&summary)?);
    } else {
        println!("{}", serde_json::to_string(&summary)?);
    }
    Ok(())
}
