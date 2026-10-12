use std::collections::HashSet;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

use anyhow::{ensure, Context, Result};
use clap::{Args, Parser, Subcommand};
use rayon::prelude::*;
use zorindex::{
    create_index_writer, finish_index_writer, mixsplit, open_index_reader,
    parse_features_from_file, read_fof, splitmix64, ColorInput, FeatureConfig, IndexMode,
    DEFAULT_MINIMIZER_SIZE, DEFAULT_MODIMIZER_SAMPLING, MAX_HASHES,
};

const MAGIC: &[u8; 8] = b"BLMIDX1\0";
const FORMAT_VERSION: u32 = 2;
const DEFAULT_FALSE_POSITIVE_RATE: f64 = 1.0 / 256.0;
const SPLIT_SEED_XOR: u64 = 0xA54D_BA5E_1F2C_8D93;

#[derive(Parser, Debug)]
#[command(
    name = "bloomindex",
    version,
    about = "Approximate colored k-mer index using interleaved per-color Bloom filters"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,

    /// Number of rayon worker threads. Defaults to all available hardware threads.
    #[arg(short = 't', long, global = true)]
    threads: Option<usize>,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Build an index from a file-of-files.
    Build(BuildArgs),
    /// Append colors from a file-of-files to an existing index.
    Append(AppendArgs),
    /// Query all indexed features from a FASTA/FASTQ file.
    Query(QueryArgs),
    /// Print index metadata.
    Info(IndexArg),
}

#[derive(Args, Debug)]
struct BuildArgs {
    /// File containing one input FASTA/FASTQ path per line, or name<TAB>path.
    #[arg(long)]
    fof: PathBuf,

    /// Query k-mer length. Findere stores shorter (k-z)-mers of length 1..=31.
    #[arg(short = 'k', long)]
    k: usize,

    /// Output index path.
    #[arg(short, long)]
    output: PathBuf,

    /// Number of Bloom hash locations per feature.
    #[arg(long, default_value_t = 4)]
    hashes: usize,

    /// Shared hash seed.
    #[arg(long, default_value_t = 69)]
    seed: u64,

    /// Feature extraction mode stored in the index.
    #[arg(long, value_enum, default_value = "kmers")]
    index_mode: IndexMode,

    /// Findere extension: index (k-z)-mers and require z+1 consecutive hits per dataset.
    #[arg(long, alias = "z", default_value_t = zorindex::findere::DEFAULT_Z)]
    findere_z: usize,

    /// Canonical minimizer length when --index-mode minimizers is used.
    #[arg(long, default_value_t = DEFAULT_MINIMIZER_SIZE)]
    minimizer_size: usize,

    /// Keep roughly 1/N canonical k-mers when --index-mode modimizers is used.
    #[arg(long, default_value_t = DEFAULT_MODIMIZER_SAMPLING)]
    modimizer_sampling: u64,

    /// Target Bloom false positive rate used to size each color filter.
    #[arg(long, default_value_t = DEFAULT_FALSE_POSITIVE_RATE)]
    false_positive_rate: f64,

    /// Override Bloom bits per largest-color feature.
    #[arg(long)]
    bits_per_item: Option<f64>,

    /// Override the common Bloom bit count for every color.
    #[arg(long)]
    bit_count: Option<usize>,

    /// Internal color layout: flat keeps one technical bin per dataset, split-large partitions outliers.
    #[arg(long, value_enum, default_value = "flat")]
    technical_layout: TechnicalLayoutMode,

    /// Explicit maximum number of unique features per technical bin in split-large mode.
    #[arg(long)]
    target_technical_features: Option<usize>,

    /// Maximum number of technical bins created for one dataset in split-large mode.
    #[arg(long, default_value_t = 64)]
    max_split_parts: usize,

    /// Maximum number of datasets parsed/built concurrently. Defaults to all Rayon workers.
    #[arg(long)]
    build_jobs: Option<usize>,
}

#[derive(Args, Debug)]
struct AppendArgs {
    /// Existing index path to extend.
    #[arg(short, long)]
    index: PathBuf,

    /// File containing new FASTA/FASTQ inputs: one path per line, or name<TAB>path.
    #[arg(long)]
    fof: PathBuf,

    /// Output path. Defaults to overwriting --index after it is fully loaded.
    #[arg(short, long)]
    output: Option<PathBuf>,
}

#[derive(Args, Debug)]
struct QueryArgs {
    /// Index path.
    #[arg(short, long)]
    index: PathBuf,

    /// Query FASTA/FASTQ path.
    #[arg(short, long)]
    query: PathBuf,

    /// Only report colors with hit ratio at least this value.
    #[arg(long, default_value_t = 0.0)]
    min_ratio: f64,
}

#[derive(Args, Debug)]
struct IndexArg {
    /// Index path.
    #[arg(short, long)]
    index: PathBuf,
}

#[derive(Clone, Debug)]
struct ColorMeta {
    name: String,
    path: String,
    occurrences: u64,
    unique_kmers: u64,
}

#[derive(Clone, Debug)]
struct TechnicalBinMeta {
    color_id: usize,
    part_id: usize,
    part_count: usize,
    unique_kmers: u64,
}

#[derive(Clone, Debug)]
struct ColorSummary {
    input: ColorInput,
    occurrences: u64,
    unique_kmers: usize,
}

#[derive(Clone, Debug)]
struct ColorBuildPlan {
    input: ColorInput,
    color_id: usize,
    occurrences: u64,
    unique_kmers: usize,
    first_technical_bin: usize,
    part_count: usize,
}

#[derive(Clone, Debug)]
struct TechnicalLayoutPlan {
    mode: TechnicalLayoutMode,
    target_features: usize,
    max_part_features: usize,
    max_part_count: usize,
    part_counts: Vec<usize>,
}

#[derive(Clone, Copy, Debug)]
struct AppendBloomSpec {
    first_technical_bin: usize,
    part_count: usize,
}

#[derive(Clone, Debug)]
struct ParsedColor {
    input: ColorInput,
    kmers: Vec<u64>,
    occurrences: u64,
}

#[derive(Clone, Debug)]
struct BloomIndex {
    k: usize,
    hashes: usize,
    seed: u64,
    target_fp_rate: f64,
    feature_config: FeatureConfig,
    bit_count: usize,
    technical_layout: TechnicalLayoutMode,
    split_target_features: u64,
    colors: Vec<ColorMeta>,
    technical_bins: Vec<TechnicalBinMeta>,
    bitstack: Vec<u64>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
enum TechnicalLayoutMode {
    /// One internal technical bin per output color.
    Flat,
    /// Split large output colors across several internal technical bins.
    #[value(name = "split-large")]
    SplitLarge,
}

impl TechnicalLayoutMode {
    fn as_u8(self) -> u8 {
        match self {
            Self::Flat => 0,
            Self::SplitLarge => 1,
        }
    }

    fn from_u8(value: u8) -> Result<Self> {
        match value {
            0 => Ok(Self::Flat),
            1 => Ok(Self::SplitLarge),
            _ => anyhow::bail!("unsupported technical layout mode {}", value),
        }
    }

    fn label(self) -> &'static str {
        match self {
            Self::Flat => "flat",
            Self::SplitLarge => "split-large",
        }
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let threads = match cli.threads {
        Some(threads) => {
            ensure!(threads > 0, "--threads must be greater than 0");
            threads
        }
        None => std::thread::available_parallelism()
            .map(|threads| threads.get())
            .unwrap_or(1),
    };
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global();

    match cli.command {
        Command::Build(args) => build_command(args),
        Command::Append(args) => append_command(args),
        Command::Query(args) => query_command(args),
        Command::Info(args) => info_command(args),
    }
}

fn build_command(args: BuildArgs) -> Result<()> {
    validate_build_args(&args)?;
    let start = Instant::now();
    let build_jobs = bounded_build_jobs(args.build_jobs)?;
    let feature_config = FeatureConfig::new(
        args.index_mode,
        args.minimizer_size,
        args.modimizer_sampling,
    )
    .with_findere_z(args.findere_z);
    let inputs = read_fof(&args.fof).with_context(|| format!("reading {}", args.fof.display()))?;
    ensure!(
        !inputs.is_empty(),
        "{} contains no inputs",
        args.fof.display()
    );

    eprintln!(
        "counting {} color files with k={} and index_mode={} using build_jobs={}",
        inputs.len(),
        args.k,
        feature_config.mode.label(),
        build_jobs
    );
    let mut summaries = summarize_inputs(&inputs, args.k, args.seed, feature_config, build_jobs)?;
    summaries.sort_by(|left, right| left.input.name.cmp(&right.input.name));

    let technical_plan = choose_bloom_technical_layout(&summaries, &args)?;
    let corrected_fp_rate =
        split_corrected_fp_rate(args.false_positive_rate, technical_plan.max_part_count);
    let bit_count = choose_bit_count(
        technical_plan.max_part_features,
        args.hashes,
        corrected_fp_rate,
        args.bits_per_item,
        args.bit_count,
    )?;
    let (color_plans, mut colors, mut technical_bins) =
        make_color_build_plan(&summaries, &technical_plan);
    let technical_bin_count = technical_bins.len();
    let color_words = color_words(technical_bin_count);
    let stack_words = bit_count
        .checked_mul(color_words)
        .context("Bloom stack size overflow")?;

    eprintln!(
        "building {} Bloom technical bins from {} colors: technical_layout={}, target_features={}, max_part_features={}, max_part_count={}, bit_count={}, color_words={}, hashes={}, target_fp_rate={:.8}, corrected_part_fp_rate={:.8}, stack_bytes={}",
        technical_bin_count,
        colors.len(),
        technical_plan.mode.label(),
        technical_plan.target_features,
        technical_plan.max_part_features,
        technical_plan.max_part_count,
        bit_count,
        color_words,
        args.hashes,
        args.false_positive_rate,
        corrected_fp_rate,
        stack_words * 8
    );
    let stack = make_atomic_stack(stack_words);
    for chunk in color_plans.chunks(build_jobs) {
        let results = chunk
            .par_iter()
            .map(|plan| {
                let parsed = parse_color(&plan.input, args.k, args.seed, feature_config)?;
                ensure!(
                    parsed.kmers.len() == plan.unique_kmers,
                    "{} changed while building: counted {} unique features but parsed {}",
                    plan.input.name,
                    plan.unique_kmers,
                    parsed.kmers.len()
                );
                ensure!(
                    parsed.occurrences == plan.occurrences,
                    "{} changed while building: counted {} features but parsed {}",
                    plan.input.name,
                    plan.occurrences,
                    parsed.occurrences
                );
                let partitions =
                    partition_keys(&parsed.kmers, plan.part_count, args.seed ^ SPLIT_SEED_XOR);
                for (part_id, keys) in partitions.iter().enumerate() {
                    set_color_bits(
                        &stack,
                        bit_count,
                        color_words,
                        plan.first_technical_bin + part_id,
                        args.hashes,
                        args.seed,
                        keys,
                    );
                }
                Ok((plan.color_id, parsed.occurrences, partitions))
            })
            .collect::<Result<Vec<_>>>()?;
        for (color_id, occurrences, partitions) in results {
            let mut unique = 0u64;
            for (part_id, keys) in partitions.iter().enumerate() {
                let technical_id = color_plans[color_id].first_technical_bin + part_id;
                technical_bins[technical_id].unique_kmers = keys.len() as u64;
                unique += keys.len() as u64;
            }
            colors[color_id].occurrences = occurrences;
            colors[color_id].unique_kmers = unique;
        }
    }
    let bitstack = stack
        .into_iter()
        .map(AtomicU64::into_inner)
        .collect::<Vec<_>>();
    let total_occurrences = colors.iter().map(|c| c.occurrences).sum::<u64>();
    let total_unique_by_color = colors.iter().map(|c| c.unique_kmers).sum::<u64>();

    let index = BloomIndex {
        k: args.k,
        hashes: args.hashes,
        seed: args.seed,
        target_fp_rate: args.false_positive_rate,
        feature_config,
        bit_count,
        technical_layout: technical_plan.mode,
        split_target_features: technical_plan.target_features as u64,
        colors,
        technical_bins,
        bitstack,
    };
    index
        .save(&args.output)
        .with_context(|| format!("writing {}", args.output.display()))?;

    eprintln!(
        "built {} in {:.3}s: colors={}, technical_bins={}, occurrences={}, unique_by_color={}, bit_count={}, stack_bytes={}",
        args.output.display(),
        start.elapsed().as_secs_f64(),
        index.colors.len(),
        index.technical_bins.len(),
        total_occurrences,
        total_unique_by_color,
        bit_count,
        index.bitstack.len() * 8
    );
    Ok(())
}

fn append_command(args: AppendArgs) -> Result<()> {
    let start = Instant::now();
    let mut index = BloomIndex::load(&args.index)
        .with_context(|| format!("loading {}", args.index.display()))?;
    let inputs = read_fof(&args.fof).with_context(|| format!("reading {}", args.fof.display()))?;
    ensure!(
        !inputs.is_empty(),
        "{} contains no inputs",
        args.fof.display()
    );

    eprintln!(
        "parsing {} new color files with k={} and index_mode={}",
        inputs.len(),
        index.k,
        index.feature_config.mode.label()
    );
    let mut parsed = inputs
        .par_iter()
        .map(|input| parse_color(input, index.k, index.seed, index.feature_config))
        .collect::<Result<Vec<_>>>()?;
    parsed.sort_by(|left, right| left.input.name.cmp(&right.input.name));

    let mut names = index
        .colors
        .iter()
        .map(|color| color.name.as_str())
        .collect::<HashSet<_>>();
    for color in &parsed {
        ensure!(
            names.insert(color.input.name.as_str()),
            "duplicate color name {}",
            color.input.name
        );
    }
    drop(names);

    let old_color_count = index.colors.len();
    let old_technical_bin_count = index.technical_bins.len();
    let append_specs = parsed
        .iter()
        .scan(old_technical_bin_count, |next_technical_bin, color| {
            let part_count = appended_part_count(
                index.technical_layout,
                index.split_target_features,
                color.kmers.len(),
            );
            let spec = AppendBloomSpec {
                first_technical_bin: *next_technical_bin,
                part_count,
            };
            *next_technical_bin += part_count;
            Some(spec)
        })
        .collect::<Vec<_>>();
    let appended_technical_bins = append_specs
        .iter()
        .map(|spec| spec.part_count)
        .sum::<usize>();
    let new_technical_bin_count = old_technical_bin_count
        .checked_add(appended_technical_bins)
        .context("technical bin count overflow")?;
    let old_color_words = color_words(old_technical_bin_count);
    let new_color_words = color_words(new_technical_bin_count);
    let mut new_stack = if new_color_words == old_color_words {
        index.bitstack
    } else {
        eprintln!(
            "expanding packed color words from {} to {}",
            old_color_words, new_color_words
        );
        expand_color_words(
            &index.bitstack,
            index.bit_count,
            old_color_words,
            new_color_words,
        )?
    };
    let stack_words = index
        .bit_count
        .checked_mul(new_color_words)
        .context("Bloom stack size overflow")?;
    ensure!(
        new_stack.len() == stack_words,
        "internal stack expansion produced wrong length"
    );

    let stack = new_stack
        .drain(..)
        .map(AtomicU64::new)
        .collect::<Vec<AtomicU64>>();
    let append_sizes = parsed
        .par_iter()
        .zip(append_specs.par_iter())
        .map(|(color, spec)| {
            let partitions =
                partition_keys(&color.kmers, spec.part_count, index.seed ^ SPLIT_SEED_XOR);
            for (part_id, keys) in partitions.iter().enumerate() {
                set_color_bits(
                    &stack,
                    index.bit_count,
                    new_color_words,
                    spec.first_technical_bin + part_id,
                    index.hashes,
                    index.seed,
                    keys,
                );
            }
            partitions
                .iter()
                .map(|keys| keys.len() as u64)
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    index.bitstack = stack
        .into_iter()
        .map(AtomicU64::into_inner)
        .collect::<Vec<_>>();
    for (offset, color) in parsed.iter().enumerate() {
        let color_id = old_color_count + offset;
        let part_sizes = &append_sizes[offset];
        let unique_kmers = part_sizes.iter().sum::<u64>();
        index.colors.push(ColorMeta {
            name: color.input.name.clone(),
            path: color.input.path.display().to_string(),
            occurrences: color.occurrences,
            unique_kmers,
        });
        let part_count = part_sizes.len().max(1);
        for (part_id, &unique_kmers) in part_sizes.iter().enumerate() {
            index.technical_bins.push(TechnicalBinMeta {
                color_id,
                part_id,
                part_count,
                unique_kmers,
            });
        }
    }
    if index.technical_layout == TechnicalLayoutMode::Flat {
        let max_appended = parsed
            .iter()
            .map(|color| color.kmers.len() as u64)
            .max()
            .unwrap_or(0);
        index.split_target_features = index.split_target_features.max(max_appended);
    }

    for color_id in old_color_count..index.colors.len() {
        let estimate = index.estimated_user_fp_rate(color_id);
        if estimate > index.target_fp_rate * 1.25 {
            eprintln!(
                "warning: appended color {} has estimated_fp_rate={:.8}, above stored target {:.8}",
                index.colors[color_id].name, estimate, index.target_fp_rate
            );
        }
    }

    let output = args.output.as_deref().unwrap_or(&args.index);
    index
        .save(output)
        .with_context(|| format!("writing {}", output.display()))?;

    let appended_unique = parsed.iter().map(|c| c.kmers.len() as u64).sum::<u64>();
    eprintln!(
        "appended {} colors ({} technical bins) to {} in {:.3}s: appended_unique={}, colors={}, stack_bytes={}",
        parsed.len(),
        appended_technical_bins,
        output.display(),
        start.elapsed().as_secs_f64(),
        appended_unique,
        index.colors.len(),
        index.bitstack.len() * 8
    );
    Ok(())
}

fn query_command(args: QueryArgs) -> Result<()> {
    ensure!(
        (0.0..=1.0).contains(&args.min_ratio),
        "--min-ratio must be in [0, 1]"
    );
    let index = BloomIndex::load(&args.index)
        .with_context(|| format!("loading {}", args.index.display()))?;
    let start = Instant::now();
    let mut scores = vec![0u64; index.colors.len()];
    let mut scratch = BloomQueryScratch::default();
    let mut total_features = 0u64;
    let feature_queries;

    if index.feature_config.mode == IndexMode::Findere {
        let result = zorindex::findere::query_file(
            &args.query,
            index.k,
            index.feature_config.findere_z,
            index.colors.len(),
            BloomQueryScratch::default,
            |scratch, key, hits| {
                index.query_kmer_into_with_scratch(key, hits, scratch);
                true
            },
        )
        .with_context(|| format!("querying {}", args.query.display()))?;
        total_features = result.total_kmers;
        feature_queries = result.smer_queries;
        scores = result.scores;
    } else {
        parse_features_from_file(
            &args.query,
            index.k,
            index.seed,
            index.feature_config,
            |kmer| {
                total_features += 1;
                index.query_kmer_into_with_scratch(kmer, &mut scores, &mut scratch);
            },
        )
        .with_context(|| format!("querying {}", args.query.display()))?;
        feature_queries = total_features;
    }

    println!(
        "#query\t{}\tk={}\tindex_mode={}\ttotal_features={}\tbit_probes={}\telapsed_s={:.3}",
        args.query.display(),
        index.k,
        index.feature_config.mode.label(),
        total_features,
        feature_queries.saturating_mul(index.hashes as u64),
        start.elapsed().as_secs_f64()
    );
    println!("color\tmatches\tratio\tunique_kmers\testimated_fp_rate");
    for (color_id, (color, &matches)) in index.colors.iter().zip(scores.iter()).enumerate() {
        let ratio = if total_features == 0 {
            0.0
        } else {
            matches as f64 / total_features as f64
        };
        if ratio >= args.min_ratio {
            println!(
                "{}\t{}\t{:.8}\t{}\t{:.8}",
                color.name,
                matches,
                ratio,
                color.unique_kmers,
                index.estimated_user_fp_rate(color_id)
            );
        }
    }
    Ok(())
}

fn info_command(args: IndexArg) -> Result<()> {
    let index = BloomIndex::load(&args.index)
        .with_context(|| format!("loading {}", args.index.display()))?;
    let total_unique = index.colors.iter().map(|c| c.unique_kmers).sum::<u64>();
    let max_unique = index
        .colors
        .iter()
        .map(|c| c.unique_kmers)
        .max()
        .unwrap_or(0);
    println!("index\t{}", args.index.display());
    println!("k\t{}", index.k);
    println!("hashes\t{}", index.hashes);
    println!("seed\t{}", index.seed);
    println!("index_mode\t{}", index.feature_config.mode.label());
    if index.feature_config.mode == IndexMode::Findere {
        println!("findere_z\t{}", index.feature_config.findere_z);
        println!("indexed_k\t{}", index.k - index.feature_config.findere_z);
    }
    println!("minimizer_size\t{}", index.feature_config.minimizer_size);
    println!(
        "modimizer_sampling\t{}",
        index.feature_config.modimizer_sampling
    );
    println!("colors\t{}", index.colors.len());
    println!("technical_layout\t{}", index.technical_layout.label());
    println!("split_target_features\t{}", index.split_target_features);
    println!("technical_bins\t{}", index.technical_bins.len());
    println!("bits_per_technical_bin\t{}", index.bit_count);
    println!("color_words\t{}", color_words(index.technical_bins.len()));
    println!("stack_words\t{}", index.bitstack.len());
    println!("stack_bytes\t{}", index.bitstack.len() * 8);
    println!("unique_by_color\t{}", total_unique);
    println!("max_unique_color\t{}", max_unique);
    println!("target_fp_rate\t{:.8}", index.target_fp_rate);
    let max_estimated_fp = (0..index.colors.len())
        .map(|color_id| index.estimated_user_fp_rate(color_id))
        .fold(0.0, f64::max);
    println!("max_color_estimated_fp_rate\t{:.8}", max_estimated_fp);
    println!("color\tunique_kmers\toccurrences\testimated_fp_rate\tpath");
    for (color_id, color) in index.colors.iter().enumerate() {
        println!(
            "{}\t{}\t{}\t{:.8}\t{}",
            color.name,
            color.unique_kmers,
            color.occurrences,
            index.estimated_user_fp_rate(color_id),
            color.path
        );
    }
    println!("technical_bin\tcolor\tpart\tparts\tunique_kmers");
    for (technical_id, bin) in index.technical_bins.iter().enumerate() {
        let color_name = index
            .colors
            .get(bin.color_id)
            .map(|color| color.name.as_str())
            .unwrap_or("<invalid>");
        println!(
            "{}\t{}\t{}\t{}\t{}",
            technical_id, color_name, bin.part_id, bin.part_count, bin.unique_kmers
        );
    }
    Ok(())
}

fn validate_build_args(args: &BuildArgs) -> Result<()> {
    if args.index_mode != IndexMode::Findere {
        ensure!((1..=31).contains(&args.k), "k must be in 1..=31");
    }
    ensure!(
        (1..=MAX_HASHES).contains(&args.hashes),
        "--hashes must be in 1..={}",
        MAX_HASHES
    );
    validate_feature_config(
        args.k,
        FeatureConfig::new(
            args.index_mode,
            args.minimizer_size,
            args.modimizer_sampling,
        )
        .with_findere_z(args.findere_z),
    )?;
    ensure!(
        args.false_positive_rate.is_finite()
            && args.false_positive_rate > 0.0
            && args.false_positive_rate < 1.0,
        "--false-positive-rate must be finite and in (0, 1)"
    );
    if let Some(bits_per_item) = args.bits_per_item {
        ensure!(
            bits_per_item.is_finite() && bits_per_item > 0.0,
            "--bits-per-item must be finite and > 0"
        );
    }
    if let Some(bit_count) = args.bit_count {
        ensure!(bit_count > 0, "--bit-count must be greater than 0");
    }
    ensure!(
        args.max_split_parts > 0,
        "--max-split-parts must be greater than 0"
    );
    if let Some(target) = args.target_technical_features {
        ensure!(
            target > 0,
            "--target-technical-features must be greater than 0"
        );
    }
    if let Some(build_jobs) = args.build_jobs {
        ensure!(build_jobs > 0, "--build-jobs must be greater than 0");
    }
    Ok(())
}

fn validate_feature_config(k: usize, config: FeatureConfig) -> Result<()> {
    ensure!(
        (1..=31).contains(&config.minimizer_size),
        "--minimizer-size must be in 1..=31"
    );
    ensure!(
        config.modimizer_sampling > 0,
        "--modimizer-sampling must be greater than 0"
    );
    match config.mode {
        IndexMode::Findere => {
            zorindex::findere::indexed_length(k, config.findere_z)?;
        }
        IndexMode::Kmers => {}
        IndexMode::Minimizers => {
            ensure!(
                k % 2 == 1,
                "minimizer indexes require odd k with canonical SIMD minimizers"
            );
            ensure!(
                (1..=k).contains(&config.minimizer_size),
                "--minimizer-size must be in 1..=k for minimizer indexes"
            );
        }
        IndexMode::Modimizers => {}
    }
    Ok(())
}

fn choose_bit_count(
    max_unique: usize,
    hashes: usize,
    target_fp_rate: f64,
    bits_per_item: Option<f64>,
    bit_count: Option<usize>,
) -> Result<usize> {
    if let Some(bit_count) = bit_count {
        return Ok(align_to_64(bit_count));
    }
    let items = max_unique.max(1);
    let bits_per_item =
        bits_per_item.unwrap_or_else(|| bits_per_item_for_fp(target_fp_rate, hashes));
    let bits = ((items as f64) * bits_per_item).ceil() as usize;
    Ok(align_to_64(bits.max(1)))
}

fn bits_per_item_for_fp(fp_rate: f64, hashes: usize) -> f64 {
    let root = fp_rate.powf(1.0 / hashes as f64);
    -(hashes as f64) / (1.0 - root).ln()
}

fn estimated_fp_rate(unique: u64, bit_count: usize, hashes: usize) -> f64 {
    if unique == 0 || bit_count == 0 || hashes == 0 {
        return 0.0;
    }
    let exponent = -((hashes as f64) * (unique as f64)) / (bit_count as f64);
    (1.0 - exponent.exp()).powi(hashes as i32)
}

fn split_corrected_fp_rate(target_fp_rate: f64, split_parts: usize) -> f64 {
    if split_parts <= 1 {
        target_fp_rate
    } else {
        1.0 - (1.0 - target_fp_rate).powf(1.0 / split_parts as f64)
    }
}

fn align_to_64(value: usize) -> usize {
    value.div_ceil(64) * 64
}

fn color_words(color_count: usize) -> usize {
    color_count.div_ceil(64).max(1)
}

fn make_atomic_stack(words: usize) -> Vec<AtomicU64> {
    let mut stack = Vec::with_capacity(words);
    stack.resize_with(words, || AtomicU64::new(0));
    stack
}

fn expand_color_words(
    old_stack: &[u64],
    bit_count: usize,
    old_color_words: usize,
    new_color_words: usize,
) -> Result<Vec<u64>> {
    ensure!(
        old_stack.len() == bit_count * old_color_words,
        "old Bloom stack length does not match dimensions"
    );
    ensure!(
        new_color_words >= old_color_words,
        "new color word count cannot shrink"
    );
    let mut new_stack = vec![0u64; bit_count * new_color_words];
    new_stack
        .par_chunks_mut(new_color_words)
        .enumerate()
        .for_each(|(bit, row)| {
            let old_start = bit * old_color_words;
            row[..old_color_words]
                .copy_from_slice(&old_stack[old_start..old_start + old_color_words]);
        });
    Ok(new_stack)
}

fn set_color_bits(
    stack: &[AtomicU64],
    bit_count: usize,
    color_words: usize,
    color_id: usize,
    hashes: usize,
    seed: u64,
    keys: &[u64],
) {
    let color_word = color_id / 64;
    let color_mask = 1u64 << (color_id % 64);
    let mut positions = [0usize; MAX_HASHES];
    for &key in keys {
        fill_bloom_positions(key, bit_count, hashes, seed, &mut positions[..hashes]);
        for &position in &positions[..hashes] {
            stack[position * color_words + color_word].fetch_or(color_mask, Ordering::Relaxed);
        }
    }
}

#[inline(always)]
fn fill_bloom_positions(key: u64, bit_count: usize, hashes: usize, seed: u64, out: &mut [usize]) {
    debug_assert!(out.len() >= hashes);
    let hash1 = mixsplit(key, seed);
    let hash2 = splitmix64(hash1);
    let bit_count = bit_count as u64;
    let mut hash = hash1;
    for slot in out.iter_mut().take(hashes) {
        *slot = reduce_u64(hash, bit_count);
        hash = hash.wrapping_add(hash2);
    }
}

#[inline(always)]
fn reduce_u64(hash: u64, modulus: u64) -> usize {
    (((hash as u128) * (modulus as u128)) >> 64) as usize
}

fn parse_color(
    input: &ColorInput,
    k: usize,
    seed: u64,
    feature_config: FeatureConfig,
) -> Result<ParsedColor> {
    let start = Instant::now();
    let mut kmers = Vec::new();
    let mut occurrences = 0u64;
    let stats = parse_features_from_file(&input.path, k, seed, feature_config, |kmer| {
        kmers.push(kmer);
        occurrences += 1;
    })?;
    kmers.sort_unstable();
    kmers.dedup();

    eprintln!(
        "parsed {}\tfeatures={}\tunique={}\tbases={}\tchunks={}\telapsed_s={:.3}",
        input.name,
        occurrences,
        kmers.len(),
        stats.bases,
        stats.chunks,
        start.elapsed().as_secs_f64()
    );

    Ok(ParsedColor {
        input: input.clone(),
        kmers,
        occurrences,
    })
}

fn bounded_build_jobs(requested: Option<usize>) -> Result<usize> {
    let jobs = requested.unwrap_or_else(|| rayon::current_num_threads().max(1));
    ensure!(jobs > 0, "--build-jobs must be greater than 0");
    Ok(jobs)
}

fn summarize_inputs(
    inputs: &[ColorInput],
    k: usize,
    seed: u64,
    feature_config: FeatureConfig,
    build_jobs: usize,
) -> Result<Vec<ColorSummary>> {
    let mut summaries = Vec::with_capacity(inputs.len());
    for chunk in inputs.chunks(build_jobs) {
        let mut chunk_summaries = chunk
            .par_iter()
            .map(|input| summarize_color(input, k, seed, feature_config))
            .collect::<Result<Vec<_>>>()?;
        summaries.append(&mut chunk_summaries);
    }
    Ok(summaries)
}

fn summarize_color(
    input: &ColorInput,
    k: usize,
    seed: u64,
    feature_config: FeatureConfig,
) -> Result<ColorSummary> {
    let start = Instant::now();
    let mut kmers = Vec::new();
    let mut occurrences = 0u64;
    let stats = parse_features_from_file(&input.path, k, seed, feature_config, |kmer| {
        kmers.push(kmer);
        occurrences += 1;
    })?;
    kmers.sort_unstable();
    kmers.dedup();
    eprintln!(
        "counted {}\tfeatures={}\tunique={}\tbases={}\tchunks={}\telapsed_s={:.3}",
        input.name,
        occurrences,
        kmers.len(),
        stats.bases,
        stats.chunks,
        start.elapsed().as_secs_f64()
    );
    Ok(ColorSummary {
        input: input.clone(),
        occurrences,
        unique_kmers: kmers.len(),
    })
}

fn choose_bloom_technical_layout(
    summaries: &[ColorSummary],
    args: &BuildArgs,
) -> Result<TechnicalLayoutPlan> {
    let counts = summaries
        .iter()
        .map(|summary| summary.unique_kmers)
        .collect::<Vec<_>>();
    let max_unique = counts.iter().copied().max().unwrap_or(0).max(1);
    if args.technical_layout == TechnicalLayoutMode::Flat {
        return Ok(TechnicalLayoutPlan {
            mode: TechnicalLayoutMode::Flat,
            target_features: max_unique,
            max_part_features: max_unique,
            max_part_count: 1,
            part_counts: vec![1; summaries.len()],
        });
    }

    let candidates = if let Some(target) = args.target_technical_features {
        vec![target.max(1)]
    } else {
        technical_target_candidates(&counts)
    };
    let color_count = summaries.len().max(1);
    let mut best: Option<(f64, usize, usize, usize, Vec<usize>)> = None;
    for target in candidates {
        let (part_counts, max_part_features) =
            split_part_counts(&counts, target, args.max_split_parts);
        let max_part_count = part_counts.iter().copied().max().unwrap_or(1);
        let technical_bin_count = part_counts.iter().sum::<usize>().max(1);
        let corrected_fp_rate = split_corrected_fp_rate(args.false_positive_rate, max_part_count);
        let bit_count = choose_bit_count(
            max_part_features,
            args.hashes,
            corrected_fp_rate,
            args.bits_per_item,
            args.bit_count,
        )?;
        let estimated_words = bit_count
            .checked_mul(color_words(technical_bin_count))
            .context("technical layout size overflow")?;
        let technical_bin_penalty = technical_bin_count as f64 / color_count as f64;
        let estimated_cost = estimated_words as f64 * technical_bin_penalty;
        let replace = best
            .as_ref()
            .map(|(best_cost, _, _, _, _)| estimated_cost < *best_cost)
            .unwrap_or(true);
        if replace {
            best = Some((
                estimated_cost,
                target,
                max_part_features,
                max_part_count,
                part_counts,
            ));
        }
    }
    let Some((_, target_features, max_part_features, max_part_count, part_counts)) = best else {
        anyhow::bail!("no technical layout candidates were generated");
    };
    Ok(TechnicalLayoutPlan {
        mode: args.technical_layout,
        target_features,
        max_part_features,
        max_part_count,
        part_counts,
    })
}

fn technical_target_candidates(counts: &[usize]) -> Vec<usize> {
    let mut nonzero = counts
        .iter()
        .copied()
        .filter(|&count| count > 0)
        .collect::<Vec<_>>();
    if nonzero.is_empty() {
        return vec![1];
    }
    nonzero.sort_unstable();
    let max_unique = *nonzero.last().unwrap();
    let mut candidates = vec![
        max_unique,
        percentile_sorted(&nonzero, 99),
        percentile_sorted(&nonzero, 95),
        percentile_sorted(&nonzero, 90),
        percentile_sorted(&nonzero, 75),
        percentile_sorted(&nonzero, 50),
    ];
    candidates.retain(|&target| target > 0);
    candidates.sort_unstable();
    candidates.dedup();
    candidates
}

fn percentile_sorted(sorted: &[usize], percentile: usize) -> usize {
    debug_assert!(!sorted.is_empty());
    let index = ((sorted.len() - 1) * percentile).div_ceil(100);
    sorted[index.min(sorted.len() - 1)].max(1)
}

fn split_part_counts(
    counts: &[usize],
    target_features: usize,
    max_split_parts: usize,
) -> (Vec<usize>, usize) {
    let target = target_features.max(1);
    let max_parts = max_split_parts.max(1);
    let mut max_part_features = 1usize;
    let part_counts = counts
        .iter()
        .map(|&count| {
            let parts = count.div_ceil(target).max(1).min(max_parts);
            max_part_features = max_part_features.max(count.div_ceil(parts));
            parts
        })
        .collect::<Vec<_>>();
    (part_counts, max_part_features)
}

fn make_color_build_plan(
    summaries: &[ColorSummary],
    technical_plan: &TechnicalLayoutPlan,
) -> (Vec<ColorBuildPlan>, Vec<ColorMeta>, Vec<TechnicalBinMeta>) {
    let mut color_plans = Vec::with_capacity(summaries.len());
    let mut colors = Vec::with_capacity(summaries.len());
    let mut technical_bins = Vec::new();
    for (color_id, summary) in summaries.iter().enumerate() {
        let part_count = technical_plan.part_counts[color_id].max(1);
        let first_technical_bin = technical_bins.len();
        colors.push(ColorMeta {
            name: summary.input.name.clone(),
            path: summary.input.path.display().to_string(),
            occurrences: summary.occurrences,
            unique_kmers: summary.unique_kmers as u64,
        });
        color_plans.push(ColorBuildPlan {
            input: summary.input.clone(),
            color_id,
            occurrences: summary.occurrences,
            unique_kmers: summary.unique_kmers,
            first_technical_bin,
            part_count,
        });
        for part_id in 0..part_count {
            technical_bins.push(TechnicalBinMeta {
                color_id,
                part_id,
                part_count,
                unique_kmers: 0,
            });
        }
    }
    (color_plans, colors, technical_bins)
}

fn flat_technical_bins(colors: &[ColorMeta]) -> Vec<TechnicalBinMeta> {
    colors
        .iter()
        .enumerate()
        .map(|(color_id, color)| TechnicalBinMeta {
            color_id,
            part_id: 0,
            part_count: 1,
            unique_kmers: color.unique_kmers,
        })
        .collect()
}

fn partition_keys(keys: &[u64], part_count: usize, seed: u64) -> Vec<Vec<u64>> {
    if part_count <= 1 {
        return vec![keys.to_vec()];
    }
    let mut partitions = (0..part_count).map(|_| Vec::new()).collect::<Vec<_>>();
    for &key in keys {
        let part = reduce_u64(mixsplit(key, seed), part_count as u64);
        partitions[part].push(key);
    }
    partitions
}

fn appended_part_count(
    technical_layout: TechnicalLayoutMode,
    split_target_features: u64,
    unique_kmers: usize,
) -> usize {
    if technical_layout == TechnicalLayoutMode::SplitLarge {
        let target = (split_target_features as usize).max(1);
        unique_kmers.div_ceil(target).max(1)
    } else {
        1
    }
}

#[derive(Default)]
struct BloomQueryScratch {
    seen_users: Vec<u32>,
    epoch: u32,
}

impl BloomQueryScratch {
    fn next_epoch(&mut self, user_count: usize) -> u32 {
        if self.seen_users.len() != user_count {
            self.seen_users.resize(user_count, 0);
            self.epoch = 0;
        }
        self.epoch = self.epoch.wrapping_add(1);
        if self.epoch == 0 {
            self.seen_users.fill(0);
            self.epoch = 1;
        }
        self.epoch
    }
}

impl BloomIndex {
    fn technical_words(&self) -> usize {
        color_words(self.technical_bins.len())
    }

    fn is_flat_technical_layout(&self) -> bool {
        self.technical_bins.len() == self.colors.len()
            && self
                .technical_bins
                .iter()
                .enumerate()
                .all(|(idx, bin)| bin.color_id == idx && bin.part_id == 0 && bin.part_count == 1)
    }

    fn estimated_user_fp_rate(&self, color_id: usize) -> f64 {
        let mut no_fp = 1.0;
        for bin in self
            .technical_bins
            .iter()
            .filter(|bin| bin.color_id == color_id)
        {
            let part_fp = estimated_fp_rate(bin.unique_kmers, self.bit_count, self.hashes);
            no_fp *= 1.0 - part_fp;
        }
        1.0 - no_fp
    }

    #[cfg(test)]
    fn query_kmer_into(&self, key: u64, scores: &mut [u64]) {
        let mut scratch = BloomQueryScratch::default();
        self.query_kmer_into_with_scratch(key, scores, &mut scratch);
    }

    fn query_kmer_into_with_scratch(
        &self,
        key: u64,
        scores: &mut [u64],
        scratch: &mut BloomQueryScratch,
    ) {
        debug_assert_eq!(scores.len(), self.colors.len());
        let color_words = self.technical_words();
        let mut positions = [0usize; MAX_HASHES];
        fill_bloom_positions(
            key,
            self.bit_count,
            self.hashes,
            self.seed,
            &mut positions[..self.hashes],
        );
        if !self.is_flat_technical_layout() {
            let epoch = scratch.next_epoch(self.colors.len());
            for word in 0..color_words {
                let mut mask = !0u64;
                for &position in &positions[..self.hashes] {
                    mask &= self.bitstack[position * color_words + word];
                    if mask == 0 {
                        break;
                    }
                }
                while mask != 0 {
                    let bit = mask.trailing_zeros() as usize;
                    let technical_id = word * 64 + bit;
                    if technical_id < self.technical_bins.len() {
                        let color_id = self.technical_bins[technical_id].color_id;
                        if scratch.seen_users[color_id] != epoch {
                            scratch.seen_users[color_id] = epoch;
                            scores[color_id] += 1;
                        }
                    }
                    mask &= mask - 1;
                }
            }
            return;
        }
        for word in 0..color_words {
            let mut mask = !0u64;
            for &position in &positions[..self.hashes] {
                mask &= self.bitstack[position * color_words + word];
                if mask == 0 {
                    break;
                }
            }
            while mask != 0 {
                let bit = mask.trailing_zeros() as usize;
                let color = word * 64 + bit;
                if color < scores.len() {
                    scores[color] += 1;
                }
                mask &= mask - 1;
            }
        }
    }

    fn save(&self, path: &Path) -> Result<()> {
        let mut writer = create_index_writer(path)?;
        writer.write_all(MAGIC)?;
        write_u32(&mut writer, FORMAT_VERSION)?;
        write_u8(&mut writer, self.k as u8)?;
        write_u8(&mut writer, self.hashes as u8)?;
        write_u16(&mut writer, 0)?;
        write_u64(&mut writer, self.seed)?;
        write_u8(&mut writer, self.feature_config.mode.as_u8())?;
        write_u8(&mut writer, self.feature_config.minimizer_size as u8)?;
        write_u64(&mut writer, self.feature_config.modimizer_sampling)?;
        // Mode-specific extension: legacy modes keep their original byte layout.
        // Old readers reject mode 3 rather than misinterpreting a findere index.
        if self.feature_config.mode == IndexMode::Findere {
            write_u8(&mut writer, self.feature_config.findere_z as u8)?;
        }
        write_f64(&mut writer, self.target_fp_rate)?;
        write_u64(&mut writer, self.bit_count as u64)?;
        write_u64(&mut writer, self.colors.len() as u64)?;
        for color in &self.colors {
            write_string(&mut writer, &color.name)?;
            write_string(&mut writer, &color.path)?;
            write_u64(&mut writer, color.occurrences)?;
            write_u64(&mut writer, color.unique_kmers)?;
        }
        write_u8(&mut writer, self.technical_layout.as_u8())?;
        write_u64(&mut writer, self.split_target_features)?;
        write_u64(&mut writer, self.technical_bins.len() as u64)?;
        for bin in &self.technical_bins {
            write_u64(&mut writer, bin.color_id as u64)?;
            write_u64(&mut writer, bin.part_id as u64)?;
            write_u64(&mut writer, bin.part_count as u64)?;
            write_u64(&mut writer, bin.unique_kmers)?;
        }
        write_u64(&mut writer, self.bitstack.len() as u64)?;
        write_u64_slice(&mut writer, &self.bitstack)?;
        finish_index_writer(writer)
    }

    fn load(path: &Path) -> Result<Self> {
        let mut reader = open_index_reader(path)?;
        let mut magic = [0u8; 8];
        reader.read_exact(&mut magic)?;
        ensure!(&magic == MAGIC, "not a BLOOMINDEX file");
        let version = read_u32(&mut reader)?;
        ensure!(
            version == 1 || version == FORMAT_VERSION,
            "unsupported Bloom index version {}",
            version
        );
        let k = read_u8(&mut reader)? as usize;
        let hashes = read_u8(&mut reader)? as usize;
        let _reserved = read_u16(&mut reader)?;
        let seed = read_u64(&mut reader)?;
        let mode = IndexMode::from_u8(read_u8(&mut reader)?)?;
        let minimizer_size = read_u8(&mut reader)? as usize;
        let modimizer_sampling = read_u64(&mut reader)?;
        let findere_z = if mode == IndexMode::Findere {
            read_u8(&mut reader)? as usize
        } else {
            zorindex::findere::DEFAULT_Z
        };
        let target_fp_rate = read_f64(&mut reader)?;
        let bit_count = read_u64(&mut reader)? as usize;
        let color_count = read_u64(&mut reader)? as usize;
        if mode != IndexMode::Findere {
            ensure!((1..=31).contains(&k), "invalid k in index");
        }
        ensure!(
            (1..=MAX_HASHES).contains(&hashes),
            "invalid Bloom hash count in index"
        );
        let feature_config = FeatureConfig {
            mode,
            minimizer_size,
            modimizer_sampling,
            findere_z,
        };
        validate_feature_config(k, feature_config)?;
        ensure!(bit_count > 0, "invalid zero bit count in index");

        let mut colors = Vec::with_capacity(color_count);
        for _ in 0..color_count {
            colors.push(ColorMeta {
                name: read_string(&mut reader)?,
                path: read_string(&mut reader)?,
                occurrences: read_u64(&mut reader)?,
                unique_kmers: read_u64(&mut reader)?,
            });
        }
        let (technical_layout, split_target_features, technical_bins) = if version >= 2 {
            let technical_layout = TechnicalLayoutMode::from_u8(read_u8(&mut reader)?)?;
            let split_target_features = read_u64(&mut reader)?;
            let technical_count = read_u64(&mut reader)? as usize;
            ensure!(technical_count > 0, "corrupt index: no technical bins");
            let mut technical_bins = Vec::with_capacity(technical_count);
            for _ in 0..technical_count {
                let color_id = read_u64(&mut reader)? as usize;
                let part_id = read_u64(&mut reader)? as usize;
                let part_count = read_u64(&mut reader)? as usize;
                let unique_kmers = read_u64(&mut reader)?;
                ensure!(
                    color_id < color_count,
                    "corrupt index: technical bin color id out of bounds"
                );
                ensure!(part_count > 0, "corrupt index: zero split part count");
                ensure!(
                    part_id < part_count,
                    "corrupt index: split part id out of bounds"
                );
                technical_bins.push(TechnicalBinMeta {
                    color_id,
                    part_id,
                    part_count,
                    unique_kmers,
                });
            }
            (technical_layout, split_target_features, technical_bins)
        } else {
            (
                TechnicalLayoutMode::Flat,
                colors
                    .iter()
                    .map(|color| color.unique_kmers)
                    .max()
                    .unwrap_or(0),
                flat_technical_bins(&colors),
            )
        };
        let bitstack_len = read_u64(&mut reader)? as usize;
        let expected_len = bit_count
            .checked_mul(color_words(technical_bins.len()))
            .context("index dimensions overflow")?;
        ensure!(
            bitstack_len == expected_len,
            "corrupt index: bitstack length {} != expected {}",
            bitstack_len,
            expected_len
        );
        let bitstack = read_u64_vec(&mut reader, bitstack_len)?;

        Ok(Self {
            k,
            hashes,
            seed,
            target_fp_rate,
            feature_config,
            bit_count,
            technical_layout,
            split_target_features,
            colors,
            technical_bins,
            bitstack,
        })
    }
}

fn write_u8<W: Write>(writer: &mut W, value: u8) -> Result<()> {
    writer.write_all(&[value])?;
    Ok(())
}

fn write_u16<W: Write>(writer: &mut W, value: u16) -> Result<()> {
    writer.write_all(&value.to_le_bytes())?;
    Ok(())
}

fn write_u32<W: Write>(writer: &mut W, value: u32) -> Result<()> {
    writer.write_all(&value.to_le_bytes())?;
    Ok(())
}

fn write_u64<W: Write>(writer: &mut W, value: u64) -> Result<()> {
    writer.write_all(&value.to_le_bytes())?;
    Ok(())
}

fn write_f64<W: Write>(writer: &mut W, value: f64) -> Result<()> {
    write_u64(writer, value.to_bits())
}

fn write_string<W: Write>(writer: &mut W, value: &str) -> Result<()> {
    ensure!(
        value.len() <= u32::MAX as usize,
        "string too large for index format"
    );
    write_u32(writer, value.len() as u32)?;
    writer.write_all(value.as_bytes())?;
    Ok(())
}

fn write_u64_slice<W: Write>(writer: &mut W, values: &[u64]) -> Result<()> {
    #[cfg(target_endian = "little")]
    {
        let bytes = unsafe {
            std::slice::from_raw_parts(values.as_ptr() as *const u8, std::mem::size_of_val(values))
        };
        writer.write_all(bytes)?;
    }
    #[cfg(not(target_endian = "little"))]
    {
        for &value in values {
            write_u64(writer, value)?;
        }
    }
    Ok(())
}

fn read_u8<R: Read>(reader: &mut R) -> Result<u8> {
    let mut bytes = [0u8; 1];
    reader.read_exact(&mut bytes)?;
    Ok(bytes[0])
}

fn read_u16<R: Read>(reader: &mut R) -> Result<u16> {
    let mut bytes = [0u8; 2];
    reader.read_exact(&mut bytes)?;
    Ok(u16::from_le_bytes(bytes))
}

fn read_u32<R: Read>(reader: &mut R) -> Result<u32> {
    let mut bytes = [0u8; 4];
    reader.read_exact(&mut bytes)?;
    Ok(u32::from_le_bytes(bytes))
}

fn read_u64<R: Read>(reader: &mut R) -> Result<u64> {
    let mut bytes = [0u8; 8];
    reader.read_exact(&mut bytes)?;
    Ok(u64::from_le_bytes(bytes))
}

fn read_f64<R: Read>(reader: &mut R) -> Result<f64> {
    Ok(f64::from_bits(read_u64(reader)?))
}

fn read_string<R: Read>(reader: &mut R) -> Result<String> {
    let len = read_u32(reader)? as usize;
    let mut bytes = vec![0u8; len];
    reader.read_exact(&mut bytes)?;
    Ok(String::from_utf8(bytes)?)
}

fn read_u64_vec<R: Read>(reader: &mut R, len: usize) -> Result<Vec<u64>> {
    let byte_len = len.checked_mul(8).context("u64 vector size overflow")?;
    let mut values = vec![0u64; len];
    #[cfg(target_endian = "little")]
    {
        let bytes =
            unsafe { std::slice::from_raw_parts_mut(values.as_mut_ptr() as *mut u8, byte_len) };
        reader.read_exact(bytes)?;
    }
    #[cfg(not(target_endian = "little"))]
    {
        for value in &mut values {
            *value = read_u64(reader)?;
        }
    }
    Ok(values)
}

#[cfg(test)]
mod tests {
    use super::*;
    use zorindex::{scan_canonical_kmers, ZSTD_FRAME_MAGIC};

    fn temp_index_path(name: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "bloomindex_{name}_{}_{}.blmidx",
            std::process::id(),
            splitmix64(name.len() as u64)
        ))
    }

    fn temp_test_dir(name: &str) -> PathBuf {
        let path = std::env::temp_dir().join(format!(
            "bloomindex_{name}_{}_{}",
            std::process::id(),
            splitmix64(name.len() as u64)
        ));
        let _ = std::fs::remove_dir_all(&path);
        std::fs::create_dir_all(&path).unwrap();
        path
    }

    fn write_fasta(path: &Path, name: &str, seq: &str) {
        std::fs::write(path, format!(">{name}\n{seq}\n")).unwrap();
    }

    fn write_fof(path: &Path, entries: &[(&str, &Path)]) {
        let mut text = String::new();
        for (name, fasta) in entries {
            text.push_str(name);
            text.push('\t');
            text.push_str(&fasta.display().to_string());
            text.push('\n');
        }
        std::fs::write(path, text).unwrap();
    }

    #[test]
    fn bit_count_for_one_over_256_with_four_hashes_is_reasonable() {
        let bits = bits_per_item_for_fp(1.0 / 256.0, 4);
        assert!(bits > 13.8 && bits < 14.0);
    }

    #[test]
    fn bloom_query_scores_inserted_colors() {
        let colors = vec![
            vec![1u64, 2, 3, 4],
            vec![4u64, 5, 6, 7],
            vec![9u64, 10, 11, 12],
        ];
        let bit_count = 4096;
        let color_words = color_words(colors.len());
        let stack = make_atomic_stack(bit_count * color_words);
        for (color_id, keys) in colors.iter().enumerate() {
            set_color_bits(&stack, bit_count, color_words, color_id, 4, 99, keys);
        }
        let bitstack = stack
            .into_iter()
            .map(AtomicU64::into_inner)
            .collect::<Vec<_>>();
        let color_metas = (0..3)
            .map(|idx| ColorMeta {
                name: format!("c{idx}"),
                path: String::new(),
                occurrences: 0,
                unique_kmers: colors[idx].len() as u64,
            })
            .collect::<Vec<_>>();
        let technical_bins = flat_technical_bins(&color_metas);
        let index = BloomIndex {
            k: 3,
            hashes: 4,
            seed: 99,
            target_fp_rate: DEFAULT_FALSE_POSITIVE_RATE,
            feature_config: FeatureConfig {
                mode: IndexMode::Kmers,
                findere_z: zorindex::findere::DEFAULT_Z,
                minimizer_size: DEFAULT_MINIMIZER_SIZE,
                modimizer_sampling: DEFAULT_MODIMIZER_SAMPLING,
            },
            bit_count,
            technical_layout: TechnicalLayoutMode::Flat,
            split_target_features: 0,
            colors: color_metas,
            technical_bins,
            bitstack,
        };
        let mut scores = vec![0u64; 3];
        index.query_kmer_into(4, &mut scores);
        assert_eq!(scores[0], 1);
        assert_eq!(scores[1], 1);
        assert_eq!(scores[2], 0);
    }

    #[test]
    fn split_bloom_technical_bins_count_once_per_user_color() {
        let bit_count = 4096;
        let technical_words = color_words(2);
        let stack = make_atomic_stack(bit_count * technical_words);
        set_color_bits(&stack, bit_count, technical_words, 0, 4, 99, &[4]);
        set_color_bits(&stack, bit_count, technical_words, 1, 4, 99, &[4]);
        let bitstack = stack
            .into_iter()
            .map(AtomicU64::into_inner)
            .collect::<Vec<_>>();
        let colors = vec![ColorMeta {
            name: "split".into(),
            path: "split.fa".into(),
            occurrences: 0,
            unique_kmers: 2,
        }];
        let technical_bins = vec![
            TechnicalBinMeta {
                color_id: 0,
                part_id: 0,
                part_count: 2,
                unique_kmers: 1,
            },
            TechnicalBinMeta {
                color_id: 0,
                part_id: 1,
                part_count: 2,
                unique_kmers: 1,
            },
        ];
        let index = BloomIndex {
            k: 3,
            hashes: 4,
            seed: 99,
            target_fp_rate: DEFAULT_FALSE_POSITIVE_RATE,
            feature_config: FeatureConfig::legacy_kmers(),
            bit_count,
            technical_layout: TechnicalLayoutMode::SplitLarge,
            split_target_features: 1,
            colors,
            technical_bins,
            bitstack,
        };
        let mut scores = vec![0u64; 1];
        index.query_kmer_into(4, &mut scores);
        assert_eq!(scores, [1]);
    }

    #[test]
    fn append_splits_large_color_in_split_large_index() {
        let dir = temp_test_dir("append_split_bloom");
        let a = dir.join("a.fa");
        let b = dir.join("b.fa");
        let build_fof = dir.join("build.fof");
        let append_fof = dir.join("append.fof");
        let index_path = dir.join("base.blmidx");
        let appended_path = dir.join("appended.blmidx");
        write_fasta(&a, "a", "ACGTACGTACGTACGT");
        write_fasta(&b, "b", "AAAACCCCCGGGGGTTTTTACGTACGT");
        write_fof(&build_fof, &[("a", &a)]);
        write_fof(&append_fof, &[("b", &b)]);

        build_command(BuildArgs {
            fof: build_fof,
            k: 5,
            output: index_path.clone(),
            hashes: 4,
            seed: 99,
            index_mode: IndexMode::Kmers,
            findere_z: zorindex::findere::DEFAULT_Z,
            minimizer_size: DEFAULT_MINIMIZER_SIZE,
            modimizer_sampling: DEFAULT_MODIMIZER_SAMPLING,
            false_positive_rate: DEFAULT_FALSE_POSITIVE_RATE,
            bits_per_item: None,
            bit_count: Some(4096),
            technical_layout: TechnicalLayoutMode::SplitLarge,
            target_technical_features: Some(2),
            max_split_parts: 64,
            build_jobs: Some(1),
        })
        .unwrap();
        append_command(AppendArgs {
            index: index_path,
            fof: append_fof,
            output: Some(appended_path.clone()),
        })
        .unwrap();

        let loaded = BloomIndex::load(&appended_path).unwrap();
        assert_eq!(loaded.colors.len(), 2);
        assert_eq!(loaded.technical_layout, TechnicalLayoutMode::SplitLarge);
        let appended_bins = loaded
            .technical_bins
            .iter()
            .filter(|bin| bin.color_id == 1)
            .count();
        assert!(appended_bins > 1);
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn canonical_encoding_matches_reverse_complement() {
        let mut left = Vec::new();
        let mut right = Vec::new();
        scan_canonical_kmers(b"ACGTAC", 4, &mut |k| left.push(k));
        scan_canonical_kmers(b"GTACGT", 4, &mut |k| right.push(k));
        left.sort_unstable();
        right.sort_unstable();
        assert_eq!(left, right);
    }

    #[test]
    fn saved_bloom_index_is_zstd_serialized_and_loadable() {
        let path = temp_index_path("saved_bloom");
        let bit_count = 128;
        let colors = vec![
            ColorMeta {
                name: "a".into(),
                path: "a.fa".into(),
                occurrences: 7,
                unique_kmers: 5,
            },
            ColorMeta {
                name: "b".into(),
                path: "b.fa".into(),
                occurrences: 9,
                unique_kmers: 6,
            },
        ];
        let bitstack = vec![0x55aa_55aa_55aa_55aau64; bit_count * color_words(colors.len())];
        let technical_bins = flat_technical_bins(&colors);
        let index = BloomIndex {
            k: 5,
            hashes: 3,
            seed: 99,
            target_fp_rate: DEFAULT_FALSE_POSITIVE_RATE,
            feature_config: FeatureConfig::legacy_kmers(),
            bit_count,
            technical_layout: TechnicalLayoutMode::Flat,
            split_target_features: 0,
            colors,
            technical_bins,
            bitstack,
        };
        index.save(&path).unwrap();

        let mut file = std::fs::File::open(&path).unwrap();
        let mut magic = [0u8; 4];
        file.read_exact(&mut magic).unwrap();
        assert_eq!(magic, ZSTD_FRAME_MAGIC);

        let loaded = BloomIndex::load(&path).unwrap();
        assert_eq!(loaded.k, index.k);
        assert_eq!(loaded.hashes, index.hashes);
        assert_eq!(loaded.seed, index.seed);
        assert_eq!(loaded.bit_count, index.bit_count);
        assert_eq!(loaded.colors.len(), index.colors.len());
        assert_eq!(loaded.bitstack, index.bitstack);
        let _ = std::fs::remove_file(path);
    }
}
