use std::cmp::{self, Reverse};
use std::collections::{BinaryHeap, HashSet, VecDeque};
use std::fs::File;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{ensure, Context, Result};
use bitpacking::{BitPacker, BitPacker8x};
use clap::{Args, Parser, Subcommand, ValueEnum};
use deko::read::AnyDecoder;
use fastpfor::AnyLenCodec;
use helicase::input::{FromMmap, FromSlice};
use helicase::{FastxParser, HelicaseParser};
use rayon::prelude::*;
use zorindex::power_merge::{
    DegreeQueue, GraphIndex, IndexedDegreeHeap, PeelGraph, TypedPeelGraph,
};
use zorindex::{
    create_index_writer, finish_index_writer, is_compressed_path, mixsplit, open_index_reader,
    parse_features_from_file, read_fof, scan_index_features, splitmix64, ColorInput, FastxStats,
    FeatureConfig, IndexMode, DEFAULT_MINIMIZER_SIZE, DEFAULT_MODIMIZER_SAMPLING, MAX_HASHES,
    PARSER_CONFIG,
};

const MAGIC: &[u8; 8] = b"ZORIDX1\0";
const FORMAT_VERSION: u32 = 7;
const FLAG_UNION_GRAPH: u8 = 1;
const MAX_SEGMENT_LENGTH_LOG: u32 = 12;
const DEFAULT_STACK_BLOCK_ROWS: usize = 4096;
const QUERY_BATCH_FEATURES: usize = 262_144;
const FINGERPRINT_SEED: u64 = 0x8E2D_4F6A_9B13_C5D7;
const UNION_SEED_XOR: u64 = 0xC6BC_2796_92B5_CC83;
const SPLIT_SEED_XOR: u64 = 0xA54D_BA5E_1F2C_8D93;

#[derive(Parser, Debug)]
#[command(
    name = "zorindex",
    version,
    about = "Approximate colored k-mer index using stacked pure ZOR filters"
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
    /// Rewrite an index with a different stack compression mode.
    Repack(RepackArgs),
    /// Query all k-mers from a FASTA/FASTQ file against the color stack.
    Query(QueryArgs),
    /// Print index metadata.
    Info(IndexArg),
}

#[derive(Args, Debug)]
struct BuildArgs {
    #[command(flatten)]
    merge: zorindex::power_merge::MergeArgs,

    /// Use the legacy single-stack index instead of the default merged power-of-two banks.
    #[arg(long, conflicts_with_all = ["power_of_two_merge", "drop_abandoned", "disabled"])]
    plain: bool,

    /// File containing one input FASTA/FASTQ path per line, or name<TAB>path.
    #[arg(long)]
    fof: PathBuf,

    /// K-mer window length. Full k-mer/modimizer modes are limited to 1..=31;
    /// minimizer/findere modes may use larger windows because only shorter features are encoded.
    #[arg(short = 'k', long)]
    k: usize,

    /// Output index path.
    #[arg(short, long)]
    output: PathBuf,

    /// Number of hash locations per ZOR equation.
    #[arg(long, default_value_t = 4)]
    hashes: usize,

    /// Fingerprint width in bits. Supported values: 8 or 16.
    #[arg(long, default_value_t = 8)]
    fingerprint_bits: u8,

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

    /// Scale applied to the largest color cardinality to size each color layer.
    #[arg(long, default_value_t = 1.0)]
    slot_scale: f64,

    /// Override the common slot count for all color layers.
    #[arg(long)]
    slot_count: Option<usize>,

    /// Override the ZOR segment length. Must be a non-zero power of two.
    #[arg(long)]
    segment_length: Option<usize>,

    /// Cycle-breaking heuristic used when the pure ZOR peel reaches a core.
    #[arg(long, value_enum, default_value = "most-deg2")]
    cycle_break: CycleBreakHeuristic,

    /// Number of tied min-degree cells to scan for heuristic cycle breaking.
    #[arg(long, default_value_t = 8)]
    tie_scan: usize,

    /// Add a separate union pure ZOR filter for graph/topology negative checks.
    #[arg(long, default_value_t = false)]
    union_graph: bool,

    /// Compression scheme for the slot-major color stack.
    #[arg(long, value_enum, default_value = "none")]
    stack_compression: StackCompression,

    /// Deprecated: slice codecs always compress one slot row per slice.
    #[arg(long, default_value_t = DEFAULT_STACK_BLOCK_ROWS)]
    stack_block_rows: usize,

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

    /// Cycle-breaking heuristic used for newly appended color layers.
    #[arg(long, value_enum, default_value = "most-deg2")]
    cycle_break: CycleBreakHeuristic,

    /// Number of tied min-degree cells to scan for heuristic cycle breaking.
    #[arg(long, default_value_t = 8)]
    tie_scan: usize,

    /// Override stack compression for the output index.
    #[arg(long, value_enum)]
    stack_compression: Option<StackCompression>,

    /// Deprecated: slice codecs always compress one slot row per slice.
    #[arg(long)]
    stack_block_rows: Option<usize>,
}

#[derive(Args, Debug)]
struct RepackArgs {
    /// Existing index path.
    #[arg(short, long)]
    index: PathBuf,

    /// Output index path.
    #[arg(short, long)]
    output: PathBuf,

    /// Compression scheme for the slot-major color stack.
    #[arg(long, value_enum)]
    stack_compression: StackCompression,

    /// Deprecated: slice codecs always compress one slot row per slice.
    #[arg(long, default_value_t = DEFAULT_STACK_BLOCK_ROWS)]
    stack_block_rows: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum CycleBreakHeuristic {
    /// Keep the first active key in a min-degree cell.
    NoHeuristic,
    /// Keep the key that leaves the smallest degree mass abandoned.
    Lightest,
    /// Keep the key that leaves the largest degree mass abandoned.
    Heaviest,
    /// Prefer abandoned sets concentrated in low-degree cells.
    MostDeg2,
    /// Keep the key that minimizes the abandoned maximum degree.
    MinMaxDegree,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum StackCompression {
    /// Store the slot-major color stack as plain bytes.
    None,
    /// Legacy StreamVByte 1/2/3/4-byte u32 coding. Larger than plain for fingerprint rows.
    #[value(skip)]
    Svb32,
    /// Compress each slot row with SIMD StreamVByte 0/1/2/4-byte u32 coding.
    #[value(name = "svb32-0124", alias = "svb32-sparse")]
    Svb32Sparse,
    /// Legacy SIMD-BP128-style bit-packing. Larger than plain for fingerprint rows.
    #[value(skip)]
    Bitpack,
    /// Compress each slot row with LZ4 block compression.
    Lz4,
    /// Compress each slot row with liblz4 block compression.
    #[value(name = "lz4-lib")]
    Lz4Lib,
    /// Compress each slot row with Snappy raw block compression.
    Snappy,
    /// Compress each slot row with pure-Rust FastPFOR-256 plus variable-byte tail.
    Fastpfor256,
    /// Compress each slot row with liblz4 HC at the minimum HC level.
    #[value(name = "lz4-hc")]
    Lz4Hc,
    /// Group fingerprints into u32 words before FastPFOR-256 coding.
    #[value(name = "fastpfor-pack")]
    FastpforPack,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
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

impl StackCompression {
    fn as_u8(self) -> u8 {
        match self {
            Self::None => 0,
            Self::Svb32 => 1,
            Self::Svb32Sparse => 2,
            Self::Bitpack => 3,
            Self::Lz4 => 5,
            Self::Snappy => 6,
            Self::Fastpfor256 => 7,
            Self::Lz4Lib => 8,
            Self::Lz4Hc => 9,
            Self::FastpforPack => 10,
        }
    }

    fn from_u8(value: u8) -> Result<Self> {
        match value {
            0 => Ok(Self::None),
            1 => Ok(Self::Svb32),
            2 => Ok(Self::Svb32Sparse),
            3 => Ok(Self::Bitpack),
            5 => Ok(Self::Lz4),
            6 => Ok(Self::Snappy),
            7 => Ok(Self::Fastpfor256),
            8 => Ok(Self::Lz4Lib),
            9 => Ok(Self::Lz4Hc),
            10 => Ok(Self::FastpforPack),
            _ => anyhow::bail!("unsupported stack compression {}", value),
        }
    }

    fn label(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::Svb32 => "svb32",
            Self::Svb32Sparse => "svb32-0124",
            Self::Bitpack => "bitpack",
            Self::Lz4 => "lz4",
            Self::Lz4Lib => "lz4-lib",
            Self::Snappy => "snappy",
            Self::Fastpfor256 => "fastpfor256",
            Self::Lz4Hc => "lz4-hc",
            Self::FastpforPack => "fastpfor-pack",
        }
    }
}

impl Default for CycleBreakHeuristic {
    fn default() -> Self {
        Self::NoHeuristic
    }
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

    /// Skip the global abandoned-key recovery map for merged indexes.
    #[arg(long)]
    ignore_abandoned: bool,

    /// Do not use the optional union graph to skip color-stack probes.
    #[arg(long, default_value_t = false)]
    ignore_union: bool,

    /// Decode a compressed stack to the plain SIMD layout before scanning queries.
    #[arg(long, default_value_t = false)]
    decode_stack: bool,
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
    abandoned: u64,
}

#[derive(Clone, Debug)]
struct TechnicalBinMeta {
    color_id: usize,
    part_id: usize,
    part_count: usize,
    unique_kmers: u64,
}

#[derive(Clone, Debug)]
struct ParsedColor {
    input: ColorInput,
    kmers: Vec<u64>,
    occurrences: u64,
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
    part_counts: Vec<usize>,
}

#[derive(Debug)]
struct ZorPartBuild {
    technical_bin_id: usize,
    unique_kmers: usize,
    abandoned: u64,
    cells: Vec<u8>,
}

#[derive(Debug)]
struct ZorColorBuild {
    color_id: usize,
    occurrences: u64,
    unique_kmers: usize,
    abandoned: u64,
    parts: Vec<ZorPartBuild>,
}

#[derive(Debug)]
struct AppendedZorPart {
    unique_kmers: usize,
    abandoned: u64,
    cells: Vec<u8>,
}

#[derive(Debug)]
struct AppendedZorColor {
    color_offset: usize,
    input: ColorInput,
    occurrences: u64,
    unique_kmers: usize,
    abandoned: u64,
    parts: Vec<AppendedZorPart>,
}

#[derive(Clone, Debug)]
struct PureZorBuild {
    cells: Vec<u8>,
    abandoned: u64,
    abandoned_keys: Vec<u64>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Layout {
    segment_length: usize,
    segment_length_mask: usize,
    segment_count: usize,
    segment_count_length: usize,
    array_length: usize,
}

#[derive(Clone, Copy, Debug)]
struct PureZorConfig {
    layout: Layout,
    hashes: usize,
    seed: u64,
    fingerprint_bits: u8,
    cycle_break: CycleBreakHeuristic,
    tie_scan: usize,
}

#[derive(Clone, Debug)]
struct UnionGraph {
    layout: Layout,
    unique_kmers: u64,
    abandoned: u64,
    cells: Vec<u8>,
}

#[derive(Debug)]
struct CompressedStack {
    codec: StackCompression,
    row_bytes: usize,
    fingerprint_bytes: usize,
    rows: usize,
    uncompressed_len: usize,
    offsets: Vec<u64>,
    data: Vec<u8>,
}

#[derive(Debug)]
enum StackStorage {
    Plain(Vec<u8>),
    Compressed(CompressedStack),
}

#[derive(Debug)]
struct ZorIndex {
    k: usize,
    hashes: usize,
    fingerprint_bits: u8,
    seed: u64,
    feature_config: FeatureConfig,
    layout: Layout,
    technical_layout: TechnicalLayoutMode,
    split_target_features: u64,
    colors: Vec<ColorMeta>,
    technical_bins: Vec<TechnicalBinMeta>,
    stack: StackStorage,
    union_graph: Option<UnionGraph>,
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
        Command::Repack(args) => repack_command(args),
        Command::Query(args) => query_command(args),
        Command::Info(args) => info_command(args),
    }
}

fn build_command(args: BuildArgs) -> Result<()> {
    validate_build_args(&args)?;
    if args.merge.power_of_two_merge && !args.plain {
        ensure!(
            args.technical_layout == TechnicalLayoutMode::Flat,
            "power-of-two merging cannot be combined with --technical-layout split-large; use --plain for the legacy layout"
        );
        ensure!(
            !args.union_graph,
            "power-of-two banks do not yet support --union-graph; use --plain for the legacy layout"
        );
        ensure!(
            args.stack_compression == StackCompression::None,
            "power-of-two banks currently support plain slices only (files are zstd compressed); use --plain for legacy stack codecs"
        );
        ensure!(
            args.segment_length.is_none(),
            "power-of-two banks select their segment length automatically; use --plain for manual segment lengths"
        );
        let config = zorindex::power_merge::BuildConfig {
            backend: zorindex::power_merge::Backend::Zor,
            k: args.k,
            hashes: args.hashes,
            seed: args.seed,
            fingerprint_bits: args.fingerprint_bits as usize,
            feature_config: FeatureConfig::new(
                args.index_mode,
                args.minimizer_size,
                args.modimizer_sampling,
            )
            .with_findere_z(args.findere_z),
            slot_scale: args.slot_scale,
            minimum_slots: args.slot_count,
            merge: args.merge.clone(),
            build_jobs: bounded_build_jobs(args.build_jobs)?,
        };
        return zorindex::power_merge::build(&args.fof, &args.output, config, |keys, slots| {
            let segment_length = (slots / args.hashes)
                .max(1)
                .min(1 << MAX_SEGMENT_LENGTH_LOG);
            let segment_length = 1 << segment_length.ilog2();
            let layout = calculate_layout_with_segment_length(slots, args.hashes, segment_length)?;
            ensure!(layout.array_length == slots, "invalid power-of-two layout");
            peel_zor_graph(
                keys,
                PureZorConfig {
                    layout,
                    hashes: args.hashes,
                    seed: args.seed,
                    fingerprint_bits: args.fingerprint_bits,
                    cycle_break: args.cycle_break,
                    tie_scan: args.tie_scan,
                },
            )
        });
    }
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

    let technical_plan = choose_zor_technical_layout(&summaries, &args)?;
    let target_slots = match args.slot_count {
        Some(slots) => slots,
        None => scaled_slots(technical_plan.max_part_features, args.slot_scale),
    }
    .max(args.hashes)
    .max(1);
    let layout = calculate_layout(target_slots, args.hashes, args.segment_length)?;
    let build_config = PureZorConfig {
        layout,
        hashes: args.hashes,
        seed: args.seed,
        fingerprint_bits: args.fingerprint_bits,
        cycle_break: args.cycle_break,
        tie_scan: args.tie_scan,
    };
    let fingerprint_bytes = fingerprint_bytes(args.fingerprint_bits)?;
    let (color_plans, mut colors, mut technical_bins) =
        make_color_build_plan(&summaries, &technical_plan);
    let color_count = colors.len();
    let technical_bin_count = technical_bins.len();
    let row_bytes = technical_bin_count
        .checked_mul(fingerprint_bytes)
        .context("row width overflow")?;
    let total_bytes = layout
        .array_length
        .checked_mul(row_bytes)
        .context("interlaced size overflow")?;
    let mut interlaced = vec![0u8; total_bytes];

    eprintln!(
        "building {} pure ZOR technical bins from {} colors: technical_layout={}, target_features={}, max_part_features={}, target_slots={}, array_slots={}, segment_length={}, segment_count={}, hashes={}, fingerprint_bits={}, cycle_break={:?}, tie_scan={}",
        technical_bin_count,
        color_count,
        technical_plan.mode.label(),
        technical_plan.target_features,
        technical_plan.max_part_features,
        target_slots,
        layout.array_length,
        layout.segment_length,
        layout.segment_count,
        args.hashes,
        args.fingerprint_bits,
        args.cycle_break,
        args.tie_scan
    );

    for chunk in color_plans.chunks(build_jobs) {
        let results = chunk
            .par_iter()
            .map(|plan| {
                build_zor_color_parts(plan, args.k, args.seed, feature_config, build_config)
            })
            .collect::<Result<Vec<_>>>()?;
        for result in results {
            for part in result.parts {
                technical_bins[part.technical_bin_id].unique_kmers = part.unique_kmers as u64;
                write_layer_into_interlaced(
                    &mut interlaced,
                    layout,
                    technical_bin_count,
                    fingerprint_bytes,
                    part.technical_bin_id,
                    &part.cells,
                )?;
            }
            colors[result.color_id].occurrences = result.occurrences;
            colors[result.color_id].unique_kmers = result.unique_kmers as u64;
            colors[result.color_id].abandoned = result.abandoned;
        }
    }

    let stack = StackStorage::from_plain(
        interlaced,
        args.stack_compression,
        row_bytes,
        fingerprint_bytes,
        layout.array_length,
    )?;

    let union_graph = if args.union_graph {
        Some(build_union_graph_streaming(
            &color_plans,
            &args,
            feature_config,
            build_config,
        )?)
    } else {
        None
    };

    let total_occurrences = colors.iter().map(|c| c.occurrences).sum::<u64>();
    let total_unique_by_color = colors.iter().map(|c| c.unique_kmers).sum::<u64>();
    let total_abandoned = colors.iter().map(|c| c.abandoned).sum::<u64>();

    let index = ZorIndex {
        k: args.k,
        hashes: args.hashes,
        fingerprint_bits: args.fingerprint_bits,
        seed: args.seed,
        feature_config,
        layout,
        technical_layout: technical_plan.mode,
        split_target_features: technical_plan.target_features as u64,
        colors,
        technical_bins,
        stack,
        union_graph,
    };
    index
        .save(&args.output)
        .with_context(|| format!("writing {}", args.output.display()))?;

    eprintln!(
        "built {} in {:.3}s: colors={}, technical_bins={}, occurrences={}, unique_by_color={}, abandoned={} ({:.4}%), stack_compression={}, stack_stored_bytes={}",
        args.output.display(),
        start.elapsed().as_secs_f64(),
        color_count,
        technical_bin_count,
        total_occurrences,
        total_unique_by_color,
        total_abandoned,
        percent(total_abandoned, total_unique_by_color),
        index.stack.compression().label(),
        index.stack.stored_bytes()
    );
    Ok(())
}

fn append_command(args: AppendArgs) -> Result<()> {
    ensure!(
        !args.index.is_dir(),
        "append is not supported for power-of-two banks; rebuild from the complete FOF"
    );
    ensure!(args.tie_scan > 0, "--tie-scan must be greater than 0");
    let start = Instant::now();
    let mut index =
        ZorIndex::load(&args.index).with_context(|| format!("loading {}", args.index.display()))?;
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

    let build_config = PureZorConfig {
        layout: index.layout,
        hashes: index.hashes,
        seed: index.seed,
        fingerprint_bits: index.fingerprint_bits,
        cycle_break: args.cycle_break,
        tie_scan: args.tie_scan,
    };
    let fingerprint_bytes = fingerprint_bytes(index.fingerprint_bits)?;
    let output_compression = args
        .stack_compression
        .unwrap_or_else(|| index.stack.compression());
    if let Some(stack_block_rows) = args.stack_block_rows {
        ensure!(
            stack_block_rows > 0,
            "--stack-block-rows must be greater than 0"
        );
    }
    let append_results = parsed
        .par_iter()
        .enumerate()
        .map(|(offset, color)| build_appended_zor_color(&index, color, offset, build_config))
        .collect::<Result<Vec<_>>>()?;
    let appended_technical_bins = append_results
        .iter()
        .map(|result| result.parts.len())
        .sum::<usize>();

    eprintln!(
        "building {} appended pure ZOR technical bins from {} colors: slots={}, hashes={}, fingerprint_bits={}, cycle_break={:?}, tie_scan={}",
        appended_technical_bins,
        parsed.len(),
        index.layout.array_length,
        index.hashes,
        index.fingerprint_bits,
        args.cycle_break,
        args.tie_scan
    );

    let old_color_count = index.colors.len();
    let old_technical_bin_count = index.technical_bins.len();
    let mut builds = Vec::with_capacity(appended_technical_bins);
    let mut new_colors = Vec::with_capacity(append_results.len());
    let mut new_technical_bins = Vec::with_capacity(appended_technical_bins);
    for result in append_results {
        let color_id = old_color_count + result.color_offset;
        let part_count = result.parts.len().max(1);
        new_colors.push(ColorMeta {
            name: result.input.name.clone(),
            path: result.input.path.display().to_string(),
            occurrences: result.occurrences,
            unique_kmers: result.unique_kmers as u64,
            abandoned: result.abandoned,
        });
        for (part_id, part) in result.parts.into_iter().enumerate() {
            new_technical_bins.push(TechnicalBinMeta {
                color_id,
                part_id,
                part_count,
                unique_kmers: part.unique_kmers as u64,
            });
            builds.push(PureZorBuild {
                cells: part.cells,
                abandoned: part.abandoned,
                abandoned_keys: Vec::new(),
            });
        }
    }

    let existing_stack = index.stack.to_plain_vec()?;
    let interlaced = append_interlaced(
        &existing_stack,
        index.layout,
        old_technical_bin_count,
        fingerprint_bytes,
        &builds,
    )?;
    index.colors.extend(new_colors);
    index.technical_bins.extend(new_technical_bins);
    if index.technical_layout == TechnicalLayoutMode::Flat {
        let max_appended = parsed
            .iter()
            .map(|color| color.kmers.len() as u64)
            .max()
            .unwrap_or(0);
        index.split_target_features = index.split_target_features.max(max_appended);
    }
    let row_bytes = index
        .technical_bins
        .len()
        .checked_mul(fingerprint_bytes)
        .context("row width overflow")?;
    index.stack = StackStorage::from_plain(
        interlaced,
        output_compression,
        row_bytes,
        fingerprint_bytes,
        index.layout.array_length,
    )?;

    if index.union_graph.is_some() {
        eprintln!(
            "dropping existing union graph: append cannot update it exactly because old k-mers are not stored"
        );
        index.union_graph = None;
    }

    let output = args.output.as_deref().unwrap_or(&args.index);
    index
        .save(output)
        .with_context(|| format!("writing {}", output.display()))?;

    let appended_unique = parsed.iter().map(|c| c.kmers.len() as u64).sum::<u64>();
    let appended_abandoned = builds.iter().map(|b| b.abandoned).sum::<u64>();
    eprintln!(
        "appended {} colors ({} technical bins) to {} in {:.3}s: appended_unique={}, appended_abandoned={} ({:.4}%), stack_compression={}, stack_stored_bytes={}",
        parsed.len(),
        appended_technical_bins,
        output.display(),
        start.elapsed().as_secs_f64(),
        appended_unique,
        appended_abandoned,
        percent(appended_abandoned, appended_unique),
        index.stack.compression().label(),
        index.stack.stored_bytes()
    );
    Ok(())
}

fn repack_command(args: RepackArgs) -> Result<()> {
    ensure!(
        !args.index.is_dir(),
        "repack is not supported for power-of-two banks"
    );
    ensure!(
        args.stack_block_rows > 0,
        "--stack-block-rows must be greater than 0"
    );
    let start = Instant::now();
    let mut index =
        ZorIndex::load(&args.index).with_context(|| format!("loading {}", args.index.display()))?;
    let fingerprint_bytes = fingerprint_bytes(index.fingerprint_bits)?;
    let row_bytes = index
        .technical_bins
        .len()
        .checked_mul(fingerprint_bytes)
        .context("row width overflow")?;
    let original_codec = index.stack.compression();
    let original_stored = index.stack.stored_bytes();
    let plain = index.stack.to_plain_vec()?;
    index.stack = StackStorage::from_plain(
        plain,
        args.stack_compression,
        row_bytes,
        fingerprint_bytes,
        index.layout.array_length,
    )?;
    let new_stored = index.stack.stored_bytes();
    index
        .save(&args.output)
        .with_context(|| format!("writing {}", args.output.display()))?;
    eprintln!(
        "repacked {} -> {} in {:.3}s: {} {} bytes -> {} {} bytes",
        args.index.display(),
        args.output.display(),
        start.elapsed().as_secs_f64(),
        original_codec.label(),
        original_stored,
        index.stack.compression().label(),
        new_stored
    );
    Ok(())
}

fn query_command(args: QueryArgs) -> Result<()> {
    if args.index.is_dir() {
        return zorindex::power_merge::query_command_with_recovery(
            &args.index,
            &args.query,
            args.min_ratio,
            !args.ignore_abandoned,
        );
    }
    ensure!(
        (0.0..=1.0).contains(&args.min_ratio),
        "--min-ratio must be in [0, 1]"
    );
    let mut index =
        ZorIndex::load(&args.index).with_context(|| format!("loading {}", args.index.display()))?;
    if args.decode_stack && !matches!(index.stack, StackStorage::Plain(_)) {
        let decode_start = Instant::now();
        let plain = index.stack.to_plain_vec()?;
        let original_codec = index.stack.compression();
        let original_stored = index.stack.stored_bytes();
        index.stack = StackStorage::Plain(plain);
        eprintln!(
            "decoded stack before query: {} {} bytes -> plain {} bytes in {:.3}s",
            original_codec.label(),
            original_stored,
            index.stack.len_uncompressed(),
            decode_start.elapsed().as_secs_f64()
        );
    }
    let start = Instant::now();
    let query_result = query_features_from_file_parallel(&index, &args.query, args.ignore_union)
        .with_context(|| format!("querying {}", args.query.display()))?;
    let scores = query_result.accum.scores;
    let total_features = query_result.accum.total_features;
    let stack_queries = query_result.accum.stack_queries;

    println!(
        "#query\t{}\tk={}\tindex_mode={}\ttotal_features={}\tstack_probes={}\telapsed_s={:.3}",
        args.query.display(),
        index.k,
        index.feature_config.mode.label(),
        total_features,
        stack_queries,
        start.elapsed().as_secs_f64()
    );
    println!("color\tmatches\tratio\tunique_kmers\tabandoned");
    for (color, &matches) in index.colors.iter().zip(scores.iter()) {
        let ratio = if total_features == 0 {
            0.0
        } else {
            matches as f64 / total_features as f64
        };
        if ratio >= args.min_ratio {
            println!(
                "{}\t{}\t{:.8}\t{}\t{}",
                color.name, matches, ratio, color.unique_kmers, color.abandoned
            );
        }
    }
    Ok(())
}

fn info_command(args: IndexArg) -> Result<()> {
    if args.index.is_dir() {
        return zorindex::power_merge::info_command(&args.index);
    }
    let index =
        ZorIndex::load(&args.index).with_context(|| format!("loading {}", args.index.display()))?;
    let total_unique = index.colors.iter().map(|c| c.unique_kmers).sum::<u64>();
    let total_abandoned = index.colors.iter().map(|c| c.abandoned).sum::<u64>();
    println!("index\t{}", args.index.display());
    println!("k\t{}", index.k);
    println!("hashes\t{}", index.hashes);
    println!("fingerprint_bits\t{}", index.fingerprint_bits);
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
    println!("slots_per_technical_bin\t{}", index.layout.array_length);
    println!("segment_length\t{}", index.layout.segment_length);
    println!("segment_count\t{}", index.layout.segment_count);
    println!("stack_bytes\t{}", index.stack.len_uncompressed());
    println!("stack_stored_bytes\t{}", index.stack.stored_bytes());
    println!("stack_compression\t{}", index.stack.compression().label());
    println!("stack_block_rows\t{}", index.stack.block_rows());
    println!("unique_by_color\t{}", total_unique);
    println!(
        "abandoned_by_color\t{}\t{:.4}%",
        total_abandoned,
        percent(total_abandoned, total_unique)
    );
    if let Some(union) = &index.union_graph {
        println!("union_slots\t{}", union.layout.array_length);
        println!("union_segment_length\t{}", union.layout.segment_length);
        println!("union_segment_count\t{}", union.layout.segment_count);
        println!("union_unique\t{}", union.unique_kmers);
        println!(
            "union_abandoned\t{}\t{:.4}%",
            union.abandoned,
            percent(union.abandoned, union.unique_kmers)
        );
    } else {
        println!("union_graph\tabsent");
    }
    println!("color\tunique_kmers\toccurrences\tabandoned\tpath");
    for color in &index.colors {
        println!(
            "{}\t{}\t{}\t{}\t{}",
            color.name, color.unique_kmers, color.occurrences, color.abandoned, color.path
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
    ensure!(
        (1..=u8::MAX as usize).contains(&args.k),
        "k must be in 1..=255"
    );
    ensure!(
        (2..=32).contains(&args.hashes),
        "--hashes must be in 2..=32"
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
    fingerprint_bytes(args.fingerprint_bits)
        .with_context(|| "--fingerprint-bits must be 8 or 16")?;
    ensure!(
        args.slot_scale.is_finite() && args.slot_scale > 0.0,
        "--slot-scale must be finite and > 0"
    );
    ensure!(args.tie_scan > 0, "--tie-scan must be greater than 0");
    ensure!(
        args.stack_block_rows > 0,
        "--stack-block-rows must be greater than 0"
    );
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
    if let Some(segment_length) = args.segment_length {
        ensure!(
            segment_length.is_power_of_two(),
            "--segment-length must be a non-zero power of two"
        );
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
        IndexMode::Kmers => {
            ensure!(
                (1..=31).contains(&k),
                "k must be in 1..=31 for full k-mer indexes"
            );
        }
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
        IndexMode::Modimizers => {
            ensure!(
                (1..=31).contains(&k),
                "k must be in 1..=31 for modimizer indexes"
            );
        }
    }
    Ok(())
}

fn fingerprint_bytes(fingerprint_bits: u8) -> Result<usize> {
    match fingerprint_bits {
        8 => Ok(1),
        16 => Ok(2),
        _ => anyhow::bail!("unsupported fingerprint width {}", fingerprint_bits),
    }
}

impl StackStorage {
    fn from_plain(
        data: Vec<u8>,
        codec: StackCompression,
        row_bytes: usize,
        fingerprint_bytes: usize,
        rows: usize,
    ) -> Result<Self> {
        ensure!(row_bytes > 0, "stack row width must be greater than 0");
        ensure!(
            fingerprint_bytes == 1 || fingerprint_bytes == 2,
            "fingerprint bytes must be 1 or 2"
        );
        ensure!(
            row_bytes % fingerprint_bytes == 0,
            "row width must be a multiple of fingerprint bytes"
        );
        let expected_len = rows
            .checked_mul(row_bytes)
            .context("stack dimensions overflow")?;
        ensure!(
            data.len() == expected_len,
            "stack length {} does not match {} rows of {} bytes",
            data.len(),
            rows,
            row_bytes
        );
        if codec == StackCompression::None {
            return Ok(Self::Plain(data));
        }

        let mut offsets = Vec::with_capacity(rows + 1);
        let mut compressed_data = Vec::new();
        offsets.push(0);
        for row in data.chunks_exact(row_bytes) {
            compress_stack_row(codec, row, fingerprint_bytes, &mut compressed_data)?;
            offsets.push(compressed_data.len() as u64);
        }
        Ok(Self::Compressed(CompressedStack {
            codec,
            row_bytes,
            fingerprint_bytes,
            rows,
            uncompressed_len: data.len(),
            offsets,
            data: compressed_data,
        }))
    }

    fn compression(&self) -> StackCompression {
        match self {
            Self::Plain(_) => StackCompression::None,
            Self::Compressed(stack) => stack.codec,
        }
    }

    fn block_rows(&self) -> usize {
        match self {
            Self::Plain(_) => 0,
            Self::Compressed(_) => 1,
        }
    }

    fn len_uncompressed(&self) -> usize {
        match self {
            Self::Plain(data) => data.len(),
            Self::Compressed(stack) => stack.uncompressed_len,
        }
    }

    fn stored_bytes(&self) -> usize {
        match self {
            Self::Plain(data) => data.len(),
            Self::Compressed(stack) => {
                stack.data.len() + stack.offsets.len() * std::mem::size_of::<u64>()
            }
        }
    }

    fn as_plain(&self) -> Option<&[u8]> {
        match self {
            Self::Plain(data) => Some(data),
            Self::Compressed(_) => None,
        }
    }

    fn to_plain_vec(&self) -> Result<Vec<u8>> {
        match self {
            Self::Plain(data) => Ok(data.clone()),
            Self::Compressed(stack) => {
                let mut out = vec![0u8; stack.uncompressed_len];
                out.par_chunks_mut(stack.row_bytes)
                    .enumerate()
                    .try_for_each(|(row, dst)| {
                        let mut scratch = RowCodecScratch::default();
                        stack.decode_row_into(row, dst, &mut scratch)
                    })?;
                Ok(out)
            }
        }
    }

    fn xor_row_into(
        &self,
        row: usize,
        acc: &mut [u8],
        scratch: &mut RowCodecScratch,
    ) -> Result<()> {
        match self {
            Self::Plain(_) => anyhow::bail!("plain stack rows are read through the SIMD fast path"),
            Self::Compressed(stack) => stack.xor_row_into(row, acc, scratch),
        }
    }
}

impl CompressedStack {
    fn row_data(&self, row: usize) -> Result<&[u8]> {
        ensure!(row < self.rows, "stack row {} out of bounds", row);
        let start = self.offsets[row] as usize;
        let end = self.offsets[row + 1] as usize;
        ensure!(
            start <= end && end <= self.data.len(),
            "corrupt index: compressed row offsets out of bounds"
        );
        Ok(&self.data[start..end])
    }

    fn decode_row_into(
        &self,
        row: usize,
        out: &mut [u8],
        scratch: &mut RowCodecScratch,
    ) -> Result<()> {
        ensure!(
            out.len() == self.row_bytes,
            "row decode output has length {} != {}",
            out.len(),
            self.row_bytes
        );
        decompress_stack_row(
            self.codec,
            self.row_data(row)?,
            self.fingerprint_bytes,
            out,
            scratch,
        )
    }

    fn xor_row_into(
        &self,
        row: usize,
        acc: &mut [u8],
        scratch: &mut RowCodecScratch,
    ) -> Result<()> {
        ensure!(
            acc.len() == self.row_bytes,
            "row xor accumulator has length {} != {}",
            acc.len(),
            self.row_bytes
        );
        xor_compressed_stack_row(
            self.codec,
            self.row_data(row)?,
            self.fingerprint_bytes,
            acc,
            scratch,
        )
    }
}

#[derive(Default)]
struct RowCodecScratch {
    values: Vec<u32>,
    words: Vec<u32>,
    bytes: Vec<u8>,
    fastpfor256: fastpfor::FastPFor256,
}

fn compress_stack_row(
    codec: StackCompression,
    row: &[u8],
    fingerprint_bytes: usize,
    out: &mut Vec<u8>,
) -> Result<()> {
    match codec {
        StackCompression::None => out.extend_from_slice(row),
        StackCompression::Svb32 => {
            let values = row_to_u32_values(row, fingerprint_bytes)?;
            svb::u32::U32Classic.encode_into(&values, out);
        }
        StackCompression::Svb32Sparse => {
            let values = row_to_u32_values(row, fingerprint_bytes)?;
            svb::u32::U32Variant0124.encode_into(&values, out);
        }
        StackCompression::Bitpack => compress_bitpack_row(row, fingerprint_bytes, out)?,
        StackCompression::Lz4 => out.extend_from_slice(&lz4_flex::compress_prepend_size(row)),
        StackCompression::Snappy => {
            let mut encoder = snap::raw::Encoder::new();
            let compressed = encoder
                .compress_vec(row)
                .context("snappy row encode failed")?;
            out.extend_from_slice(&compressed);
        }
        StackCompression::Fastpfor256 => compress_fastpfor256_row(row, fingerprint_bytes, out)?,
        StackCompression::Lz4Lib => {
            let mut compressed = vec![0u8; lzzzz::lz4::max_compressed_size(row.len())];
            let written = lzzzz::lz4::compress(row, &mut compressed, lzzzz::lz4::ACC_LEVEL_DEFAULT)
                .context("liblz4 row encode failed")?;
            out.extend_from_slice(&compressed[..written]);
        }
        StackCompression::Lz4Hc => {
            let mut compressed = vec![0u8; lzzzz::lz4::max_compressed_size(row.len())];
            let written = lzzzz::lz4_hc::compress(row, &mut compressed, lzzzz::lz4_hc::CLEVEL_MIN)
                .context("lz4-hc row encode failed")?;
            out.extend_from_slice(&compressed[..written]);
        }
        StackCompression::FastpforPack => compress_fastpfor_pack_row(row, fingerprint_bytes, out)?,
    }
    Ok(())
}

fn decompress_stack_row(
    codec: StackCompression,
    data: &[u8],
    fingerprint_bytes: usize,
    out: &mut [u8],
    scratch: &mut RowCodecScratch,
) -> Result<()> {
    match codec {
        StackCompression::None => {
            ensure!(
                data.len() == out.len(),
                "plain row length {} != expected {}",
                data.len(),
                out.len()
            );
            out.copy_from_slice(data);
        }
        StackCompression::Svb32 => {
            decode_svb32_row(data, fingerprint_bytes, out, false, scratch)?;
        }
        StackCompression::Svb32Sparse => {
            decode_svb32_row(data, fingerprint_bytes, out, true, scratch)?;
        }
        StackCompression::Bitpack => decompress_bitpack_row(data, fingerprint_bytes, out, scratch)?,
        StackCompression::Lz4 => decompress_lz4_row(data, out)?,
        StackCompression::Snappy => decompress_snappy_row(data, out)?,
        StackCompression::Fastpfor256 => {
            decompress_fastpfor256_row(data, fingerprint_bytes, out, scratch)?
        }
        StackCompression::Lz4Lib => decompress_lz4_lib_row(data, out)?,
        StackCompression::Lz4Hc => decompress_lz4_lib_row(data, out)?,
        StackCompression::FastpforPack => {
            decompress_fastpfor_pack_row(data, fingerprint_bytes, out, scratch)?
        }
    }
    Ok(())
}

fn xor_compressed_stack_row(
    codec: StackCompression,
    data: &[u8],
    fingerprint_bytes: usize,
    acc: &mut [u8],
    scratch: &mut RowCodecScratch,
) -> Result<()> {
    match codec {
        StackCompression::None => xor_bytes_into(acc, data),
        StackCompression::Svb32 => {
            decode_svb32_values(data, acc.len() / fingerprint_bytes, false, scratch)?;
            xor_u32_values_into_row(&scratch.values, fingerprint_bytes, acc)?;
        }
        StackCompression::Svb32Sparse => {
            decode_svb32_values(data, acc.len() / fingerprint_bytes, true, scratch)?;
            xor_u32_values_into_row(&scratch.values, fingerprint_bytes, acc)?;
        }
        StackCompression::Bitpack => {
            decode_bitpack_values(data, acc.len() / fingerprint_bytes, scratch)?;
            xor_u32_values_into_row(
                &scratch.values[..acc.len() / fingerprint_bytes],
                fingerprint_bytes,
                acc,
            )?;
        }
        StackCompression::Lz4 => {
            scratch.bytes.resize(acc.len(), 0);
            decompress_lz4_row(data, &mut scratch.bytes)?;
            xor_bytes_into(acc, &scratch.bytes);
        }
        StackCompression::Snappy => {
            scratch.bytes.resize(acc.len(), 0);
            decompress_snappy_row(data, &mut scratch.bytes)?;
            xor_bytes_into(acc, &scratch.bytes);
        }
        StackCompression::Fastpfor256 => {
            decode_fastpfor256_values(data, acc.len() / fingerprint_bytes, scratch)?;
            xor_u32_values_into_row(
                &scratch.values[..acc.len() / fingerprint_bytes],
                fingerprint_bytes,
                acc,
            )?;
        }
        StackCompression::Lz4Lib => {
            scratch.bytes.resize(acc.len(), 0);
            decompress_lz4_lib_row(data, &mut scratch.bytes)?;
            xor_bytes_into(acc, &scratch.bytes);
        }
        StackCompression::Lz4Hc => {
            scratch.bytes.resize(acc.len(), 0);
            decompress_lz4_lib_row(data, &mut scratch.bytes)?;
            xor_bytes_into(acc, &scratch.bytes);
        }
        StackCompression::FastpforPack => {
            decode_fastpfor_pack_values(data, acc.len(), fingerprint_bytes, scratch)?;
            xor_packed_u32_values_into_row(&scratch.values, fingerprint_bytes, acc)?;
        }
    }
    Ok(())
}

fn row_to_u32_values(row: &[u8], fingerprint_bytes: usize) -> Result<Vec<u32>> {
    match fingerprint_bytes {
        1 => Ok(row.iter().map(|&value| value as u32).collect()),
        2 => {
            ensure!(row.len() % 2 == 0, "u16 row has odd byte length");
            Ok(row
                .chunks_exact(2)
                .map(|bytes| u16::from_le_bytes([bytes[0], bytes[1]]) as u32)
                .collect())
        }
        _ => anyhow::bail!("unsupported fingerprint byte width {}", fingerprint_bytes),
    }
}

fn write_u32_values_to_row(values: &[u32], fingerprint_bytes: usize, out: &mut [u8]) -> Result<()> {
    match fingerprint_bytes {
        1 => {
            ensure!(values.len() == out.len(), "decoded value count mismatch");
            for (dst, &value) in out.iter_mut().zip(values.iter()) {
                ensure!(value <= u8::MAX as u32, "decoded u8 value out of range");
                *dst = value as u8;
            }
        }
        2 => {
            ensure!(
                values.len() * 2 == out.len(),
                "decoded value count mismatch"
            );
            for (idx, &value) in values.iter().enumerate() {
                ensure!(value <= u16::MAX as u32, "decoded u16 value out of range");
                out[idx * 2..idx * 2 + 2].copy_from_slice(&(value as u16).to_le_bytes());
            }
        }
        _ => anyhow::bail!("unsupported fingerprint byte width {}", fingerprint_bytes),
    }
    Ok(())
}

fn decode_svb32_row(
    data: &[u8],
    fingerprint_bytes: usize,
    out: &mut [u8],
    sparse: bool,
    scratch: &mut RowCodecScratch,
) -> Result<()> {
    let value_count = out.len() / fingerprint_bytes;
    decode_svb32_values(data, value_count, sparse, scratch)?;
    write_u32_values_to_row(&scratch.values, fingerprint_bytes, out)
}

fn decode_svb32_values(
    data: &[u8],
    value_count: usize,
    sparse: bool,
    scratch: &mut RowCodecScratch,
) -> Result<()> {
    scratch.values.clear();
    if sparse {
        svb::u32::U32Variant0124
            .decode_into(data, value_count, &mut scratch.values)
            .context("svb32-sparse row decode failed")?;
    } else {
        svb::u32::U32Classic
            .decode_into(data, value_count, &mut scratch.values)
            .context("svb32 row decode failed")?;
    }
    ensure!(
        scratch.values.len() == value_count,
        "decoded row value count {} != expected {}",
        scratch.values.len(),
        value_count
    );
    Ok(())
}

fn xor_u32_values_into_row(values: &[u32], fingerprint_bytes: usize, acc: &mut [u8]) -> Result<()> {
    match fingerprint_bytes {
        1 => {
            ensure!(values.len() == acc.len(), "decoded value count mismatch");
            for (dst, &value) in acc.iter_mut().zip(values.iter()) {
                ensure!(value <= u8::MAX as u32, "decoded u8 value out of range");
                *dst ^= value as u8;
            }
        }
        2 => {
            ensure!(
                values.len() * 2 == acc.len(),
                "decoded value count mismatch"
            );
            for (idx, &value) in values.iter().enumerate() {
                ensure!(value <= u16::MAX as u32, "decoded u16 value out of range");
                let bytes = (value as u16).to_le_bytes();
                acc[idx * 2] ^= bytes[0];
                acc[idx * 2 + 1] ^= bytes[1];
            }
        }
        _ => anyhow::bail!("unsupported fingerprint byte width {}", fingerprint_bytes),
    }
    Ok(())
}

fn compress_bitpack_row(row: &[u8], fingerprint_bytes: usize, out: &mut Vec<u8>) -> Result<()> {
    let mut values = row_to_u32_values(row, fingerprint_bytes)?;
    let bitpacker = BitPacker8x::new();
    let block_len = BitPacker8x::BLOCK_LEN;
    let padded_len = values.len().div_ceil(block_len) * block_len;
    values.resize(padded_len, 0);
    let mut compressed = vec![0u8; block_len * 4];
    for block in values.chunks_exact(block_len) {
        let num_bits = bitpacker.num_bits(block);
        out.push(num_bits);
        let written = bitpacker.compress(block, &mut compressed, num_bits);
        out.extend_from_slice(&compressed[..written]);
    }
    Ok(())
}

fn decompress_bitpack_row(
    data: &[u8],
    fingerprint_bytes: usize,
    out: &mut [u8],
    scratch: &mut RowCodecScratch,
) -> Result<()> {
    let value_count = out.len() / fingerprint_bytes;
    decode_bitpack_values(data, value_count, scratch)?;
    write_u32_values_to_row(&scratch.values[..value_count], fingerprint_bytes, out)
}

fn decode_bitpack_values(
    data: &[u8],
    value_count: usize,
    scratch: &mut RowCodecScratch,
) -> Result<()> {
    let bitpacker = BitPacker8x::new();
    let block_len = BitPacker8x::BLOCK_LEN;
    let padded_len = value_count.div_ceil(block_len) * block_len;
    scratch.values.resize(padded_len, 0);
    let mut input_offset = 0usize;
    for block in scratch.values.chunks_exact_mut(block_len) {
        ensure!(input_offset < data.len(), "truncated bitpack row");
        let num_bits = data[input_offset];
        input_offset += 1;
        ensure!(num_bits <= 32, "invalid bitpack width {}", num_bits);
        let byte_len = BitPacker8x::compressed_block_size(num_bits);
        ensure!(
            input_offset + byte_len <= data.len(),
            "truncated bitpack row block"
        );
        bitpacker.decompress(
            &data[input_offset..input_offset + byte_len],
            block,
            num_bits,
        );
        input_offset += byte_len;
    }
    ensure!(
        input_offset == data.len(),
        "bitpack row has {} trailing bytes",
        data.len() - input_offset
    );
    Ok(())
}

fn compress_fastpfor256_row(row: &[u8], fingerprint_bytes: usize, out: &mut Vec<u8>) -> Result<()> {
    let values = row_to_u32_values(row, fingerprint_bytes)?;
    let mut words = Vec::new();
    let mut codec = fastpfor::FastPFor256::default();
    codec
        .encode(&values, &mut words)
        .context("fastpfor256 row encode failed")?;
    write_u32_words_le(&words, out);
    Ok(())
}

fn decompress_fastpfor256_row(
    data: &[u8],
    fingerprint_bytes: usize,
    out: &mut [u8],
    scratch: &mut RowCodecScratch,
) -> Result<()> {
    let value_count = out.len() / fingerprint_bytes;
    decode_fastpfor256_values(data, value_count, scratch)?;
    write_u32_values_to_row(&scratch.values[..value_count], fingerprint_bytes, out)
}

fn decode_fastpfor256_values(
    data: &[u8],
    value_count: usize,
    scratch: &mut RowCodecScratch,
) -> Result<()> {
    read_u32_words_le_into(data, &mut scratch.words)?;
    scratch.values.clear();
    scratch
        .fastpfor256
        .decode(
            &scratch.words,
            &mut scratch.values,
            Some(
                value_count
                    .try_into()
                    .context("fastpfor256 value count exceeds u32")?,
            ),
        )
        .context("fastpfor256 row decode failed")?;
    ensure!(
        scratch.values.len() == value_count,
        "decoded fastpfor256 value count {} != expected {}",
        scratch.values.len(),
        value_count
    );
    Ok(())
}

fn compress_fastpfor_pack_row(
    row: &[u8],
    fingerprint_bytes: usize,
    out: &mut Vec<u8>,
) -> Result<()> {
    let values = row_to_packed_u32_values(row, fingerprint_bytes)?;
    let mut words = Vec::new();
    let mut codec = fastpfor::FastPFor256::default();
    codec
        .encode(&values, &mut words)
        .context("fastpfor-pack row encode failed")?;
    write_u32_words_le(&words, out);
    Ok(())
}

fn decompress_fastpfor_pack_row(
    data: &[u8],
    fingerprint_bytes: usize,
    out: &mut [u8],
    scratch: &mut RowCodecScratch,
) -> Result<()> {
    decode_fastpfor_pack_values(data, out.len(), fingerprint_bytes, scratch)?;
    write_packed_u32_values_to_row(&scratch.values, fingerprint_bytes, out)
}

fn decode_fastpfor_pack_values(
    data: &[u8],
    row_bytes: usize,
    fingerprint_bytes: usize,
    scratch: &mut RowCodecScratch,
) -> Result<()> {
    let value_count = packed_u32_value_count(row_bytes, fingerprint_bytes)?;
    read_u32_words_le_into(data, &mut scratch.words)?;
    scratch.values.clear();
    scratch
        .fastpfor256
        .decode(
            &scratch.words,
            &mut scratch.values,
            Some(
                value_count
                    .try_into()
                    .context("fastpfor-pack value count exceeds u32")?,
            ),
        )
        .context("fastpfor-pack row decode failed")?;
    ensure!(
        scratch.values.len() == value_count,
        "decoded fastpfor-pack value count {} != expected {}",
        scratch.values.len(),
        value_count
    );
    Ok(())
}

fn row_to_packed_u32_values(row: &[u8], fingerprint_bytes: usize) -> Result<Vec<u32>> {
    ensure!(
        fingerprint_bytes == 1 || fingerprint_bytes == 2,
        "unsupported fingerprint byte width {}",
        fingerprint_bytes
    );
    ensure!(
        row.len() % fingerprint_bytes == 0,
        "row byte length is not a multiple of fingerprint width"
    );
    let mut values = Vec::with_capacity(packed_u32_value_count(row.len(), fingerprint_bytes)?);
    for chunk in row.chunks(4) {
        let mut bytes = [0u8; 4];
        bytes[..chunk.len()].copy_from_slice(chunk);
        values.push(u32::from_le_bytes(bytes));
    }
    Ok(values)
}

fn packed_u32_value_count(row_bytes: usize, fingerprint_bytes: usize) -> Result<usize> {
    ensure!(
        fingerprint_bytes == 1 || fingerprint_bytes == 2,
        "unsupported fingerprint byte width {}",
        fingerprint_bytes
    );
    ensure!(
        row_bytes % fingerprint_bytes == 0,
        "row byte length is not a multiple of fingerprint width"
    );
    Ok(row_bytes.div_ceil(4))
}

fn write_packed_u32_values_to_row(
    values: &[u32],
    fingerprint_bytes: usize,
    out: &mut [u8],
) -> Result<()> {
    ensure!(
        values.len() == packed_u32_value_count(out.len(), fingerprint_bytes)?,
        "decoded packed value count mismatch"
    );
    for (word_idx, &value) in values.iter().enumerate() {
        let bytes = value.to_le_bytes();
        let offset = word_idx * 4;
        let end = (offset + 4).min(out.len());
        out[offset..end].copy_from_slice(&bytes[..end - offset]);
    }
    Ok(())
}

fn xor_packed_u32_values_into_row(
    values: &[u32],
    fingerprint_bytes: usize,
    acc: &mut [u8],
) -> Result<()> {
    ensure!(
        values.len() == packed_u32_value_count(acc.len(), fingerprint_bytes)?,
        "decoded packed value count mismatch"
    );
    for (word_idx, &value) in values.iter().enumerate() {
        let bytes = value.to_le_bytes();
        let offset = word_idx * 4;
        let end = (offset + 4).min(acc.len());
        for byte in 0..end - offset {
            acc[offset + byte] ^= bytes[byte];
        }
    }
    Ok(())
}

fn decompress_lz4_row(data: &[u8], out: &mut [u8]) -> Result<()> {
    ensure!(data.len() >= 4, "truncated lz4 row header");
    let expected = read_u32_le(data, 0) as usize;
    ensure!(
        expected == out.len(),
        "lz4 row length {} != expected {}",
        expected,
        out.len()
    );
    let written = lz4_flex::decompress_into(&data[4..], out).context("lz4 row decode failed")?;
    ensure!(
        written == out.len(),
        "lz4 row wrote {} bytes != expected {}",
        written,
        out.len()
    );
    Ok(())
}

fn decompress_lz4_lib_row(data: &[u8], out: &mut [u8]) -> Result<()> {
    let written = lzzzz::lz4::decompress(data, out).context("liblz4 row decode failed")?;
    ensure!(
        written == out.len(),
        "liblz4 row wrote {} bytes != expected {}",
        written,
        out.len()
    );
    Ok(())
}

fn decompress_snappy_row(data: &[u8], out: &mut [u8]) -> Result<()> {
    let mut decoder = snap::raw::Decoder::new();
    let written = decoder
        .decompress(data, out)
        .context("snappy row decode failed")?;
    ensure!(
        written == out.len(),
        "snappy row wrote {} bytes != expected {}",
        written,
        out.len()
    );
    Ok(())
}

fn write_u32_words_le(words: &[u32], out: &mut Vec<u8>) {
    out.reserve(words.len() * 4);
    for &word in words {
        out.extend_from_slice(&word.to_le_bytes());
    }
}

fn read_u32_words_le_into(data: &[u8], out: &mut Vec<u32>) -> Result<()> {
    ensure!(
        data.len() % 4 == 0,
        "u32 codec row byte length {} is not divisible by 4",
        data.len()
    );
    out.clear();
    out.reserve(data.len() / 4);
    for chunk in data.chunks_exact(4) {
        out.push(u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]));
    }
    Ok(())
}

fn decompress_v5_stack_block(codec: u8, data: &[u8], uncompressed_len: usize) -> Result<Vec<u8>> {
    match codec {
        1 => lz4_flex::decompress_size_prepended(data).context("v5 lz4 decompression failed"),
        2 => {
            let mut decoder = snap::raw::Decoder::new();
            decoder
                .decompress_vec(data)
                .context("v5 snappy decompression failed")
        }
        3 => zstd::bulk::decompress(data, uncompressed_len).context("v5 zstd decompression failed"),
        _ => anyhow::bail!("unsupported v5 stack compression {}", codec),
    }
}

fn append_interlaced(
    existing: &[u8],
    layout: Layout,
    old_color_count: usize,
    fingerprint_bytes: usize,
    builds: &[PureZorBuild],
) -> Result<Vec<u8>> {
    let new_color_count = builds.len();
    let total_color_count = old_color_count
        .checked_add(new_color_count)
        .context("color count overflow")?;
    let old_row_bytes = old_color_count
        .checked_mul(fingerprint_bytes)
        .context("row width overflow")?;
    let new_row_bytes = total_color_count
        .checked_mul(fingerprint_bytes)
        .context("row width overflow")?;
    ensure!(
        existing.len() == layout.array_length * old_row_bytes,
        "existing interlaced data length does not match layout"
    );
    for build in builds {
        ensure!(
            build.cells.len() == layout.array_length * fingerprint_bytes,
            "new build cell length does not match layout"
        );
    }

    let mut interlaced = vec![0u8; layout.array_length * new_row_bytes];
    interlaced
        .par_chunks_mut(new_row_bytes)
        .enumerate()
        .for_each(|(slot, row)| {
            let old_src = slot * old_row_bytes;
            row[..old_row_bytes].copy_from_slice(&existing[old_src..old_src + old_row_bytes]);
            for (new_color_id, build) in builds.iter().enumerate() {
                let src = slot * fingerprint_bytes;
                let dst = (old_color_count + new_color_id) * fingerprint_bytes;
                row[dst..dst + fingerprint_bytes]
                    .copy_from_slice(&build.cells[src..src + fingerprint_bytes]);
            }
        });

    Ok(interlaced)
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

fn choose_zor_technical_layout(
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
            part_counts: vec![1; summaries.len()],
        });
    }

    let candidates = if let Some(target) = args.target_technical_features {
        vec![target.max(1)]
    } else {
        technical_target_candidates(&counts)
    };
    let fingerprint_bytes = fingerprint_bytes(args.fingerprint_bits)?;
    let color_count = summaries.len().max(1);
    let mut best: Option<(f64, usize, usize, Vec<usize>)> = None;
    for target in candidates {
        let (part_counts, max_part_features) =
            split_part_counts(&counts, target, args.max_split_parts);
        let technical_bin_count = part_counts.iter().sum::<usize>().max(1);
        let target_slots = match args.slot_count {
            Some(slots) => slots,
            None => scaled_slots(max_part_features, args.slot_scale),
        }
        .max(args.hashes)
        .max(1);
        let layout = calculate_layout(target_slots, args.hashes, args.segment_length)?;
        let estimated_bytes = layout
            .array_length
            .checked_mul(technical_bin_count)
            .and_then(|bytes| bytes.checked_mul(fingerprint_bytes))
            .context("technical layout size overflow")?;
        let technical_bin_penalty = technical_bin_count as f64 / color_count as f64;
        let estimated_cost = estimated_bytes as f64 * technical_bin_penalty;
        let replace = best
            .as_ref()
            .map(|(best_cost, _, _, _)| estimated_cost < *best_cost)
            .unwrap_or(true);
        if replace {
            best = Some((estimated_cost, target, max_part_features, part_counts));
        }
    }
    let Some((_, target_features, max_part_features, part_counts)) = best else {
        anyhow::bail!("no technical layout candidates were generated");
    };
    Ok(TechnicalLayoutPlan {
        mode: args.technical_layout,
        target_features,
        max_part_features,
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
            abandoned: 0,
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

fn build_zor_color_parts(
    plan: &ColorBuildPlan,
    k: usize,
    seed: u64,
    feature_config: FeatureConfig,
    build_config: PureZorConfig,
) -> Result<ZorColorBuild> {
    let parsed = parse_color(&plan.input, k, seed, feature_config)?;
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
    let partitions = partition_keys(&parsed.kmers, plan.part_count, seed ^ SPLIT_SEED_XOR);
    let (parts, abandoned) =
        build_zor_part_sequence(partitions, plan.first_technical_bin, build_config)?;
    Ok(ZorColorBuild {
        color_id: plan.color_id,
        occurrences: parsed.occurrences,
        unique_kmers: parsed.kmers.len(),
        abandoned,
        parts,
    })
}

fn build_appended_zor_color(
    index: &ZorIndex,
    color: &ParsedColor,
    color_offset: usize,
    build_config: PureZorConfig,
) -> Result<AppendedZorColor> {
    let part_count = appended_part_count(
        index.technical_layout,
        index.split_target_features,
        color.kmers.len(),
    );
    let partitions = partition_keys(&color.kmers, part_count, index.seed ^ SPLIT_SEED_XOR);
    let (parts, abandoned) = build_zor_part_sequence(partitions, 0, build_config)?;
    let parts = parts
        .into_iter()
        .map(|part| AppendedZorPart {
            unique_kmers: part.unique_kmers,
            abandoned: part.abandoned,
            cells: part.cells,
        })
        .collect();
    Ok(AppendedZorColor {
        color_offset,
        input: color.input.clone(),
        occurrences: color.occurrences,
        unique_kmers: color.kmers.len(),
        abandoned,
        parts,
    })
}

fn build_zor_part_sequence(
    partitions: Vec<Vec<u64>>,
    first_technical_bin: usize,
    build_config: PureZorConfig,
) -> Result<(Vec<ZorPartBuild>, u64)> {
    let mut parts = Vec::with_capacity(partitions.len());
    let mut carry = Vec::new();
    let mut final_abandoned = 0u64;

    for (part_id, mut keys) in partitions.into_iter().enumerate() {
        if !carry.is_empty() {
            keys.extend_from_slice(&carry);
            keys.sort_unstable();
            keys.dedup();
        }
        let build = build_pure_zor(&keys, build_config)?;
        final_abandoned = build.abandoned;
        carry = build.abandoned_keys;
        parts.push(ZorPartBuild {
            technical_bin_id: first_technical_bin + part_id,
            unique_kmers: keys.len(),
            abandoned: build.abandoned,
            cells: build.cells,
        });
    }

    Ok((parts, final_abandoned))
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

#[inline(always)]
fn reduce_u64(hash: u64, modulus: u64) -> usize {
    (((hash as u128) * (modulus as u128)) >> 64) as usize
}

fn write_layer_into_interlaced(
    interlaced: &mut [u8],
    layout: Layout,
    technical_bin_count: usize,
    fingerprint_bytes: usize,
    technical_bin_id: usize,
    cells: &[u8],
) -> Result<()> {
    ensure!(
        cells.len() == layout.array_length * fingerprint_bytes,
        "build cell length does not match layout"
    );
    ensure!(
        technical_bin_id < technical_bin_count,
        "technical bin id out of bounds"
    );
    let row_bytes = technical_bin_count
        .checked_mul(fingerprint_bytes)
        .context("row width overflow")?;
    ensure!(
        interlaced.len() == layout.array_length * row_bytes,
        "interlaced stack length does not match layout"
    );
    for slot in 0..layout.array_length {
        let src = slot * fingerprint_bytes;
        let dst = slot * row_bytes + technical_bin_id * fingerprint_bytes;
        interlaced[dst..dst + fingerprint_bytes]
            .copy_from_slice(&cells[src..src + fingerprint_bytes]);
    }
    Ok(())
}

fn build_union_graph_streaming(
    color_plans: &[ColorBuildPlan],
    args: &BuildArgs,
    feature_config: FeatureConfig,
    build_config: PureZorConfig,
) -> Result<UnionGraph> {
    eprintln!("building optional union graph filter");
    let mut union = Vec::new();
    for plan in color_plans {
        let parsed = parse_color(&plan.input, args.k, args.seed, feature_config)?;
        union.extend_from_slice(&parsed.kmers);
    }
    union.sort_unstable();
    union.dedup();
    let union_target_slots = scaled_slots(union.len(), args.slot_scale)
        .max(args.hashes)
        .max(1);
    let union_layout = calculate_layout(union_target_slots, args.hashes, args.segment_length)?;
    let build = build_pure_zor(
        &union,
        PureZorConfig {
            layout: union_layout,
            seed: args.seed ^ UNION_SEED_XOR,
            ..build_config
        },
    )?;
    Ok(UnionGraph {
        layout: union_layout,
        unique_kmers: union.len() as u64,
        abandoned: build.abandoned,
        cells: build.cells,
    })
}

#[derive(Debug)]
struct QueryRunResult {
    _stats: FastxStats,
    accum: QueryAccum,
}

#[derive(Debug)]
struct QueryAccum {
    scores: Vec<u64>,
    total_features: u64,
    stack_queries: u64,
}

impl QueryAccum {
    fn new(color_count: usize) -> Self {
        Self {
            scores: vec![0; color_count],
            total_features: 0,
            stack_queries: 0,
        }
    }

    fn merge(&mut self, other: QueryAccum) {
        self.total_features += other.total_features;
        self.stack_queries += other.stack_queries;
        for (left, right) in self.scores.iter_mut().zip(other.scores) {
            *left += right;
        }
    }
}

fn query_features_from_file_parallel(
    index: &ZorIndex,
    path: &Path,
    ignore_union: bool,
) -> Result<QueryRunResult> {
    if index.feature_config.mode == IndexMode::Findere {
        let result = zorindex::findere::query_file(
            path,
            index.k,
            index.feature_config.findere_z,
            index.colors.len(),
            QueryScratch::default,
            |scratch, key, hits| {
                if !ignore_union {
                    if let Some(union) = &index.union_graph {
                        if !contains_raw(
                            key,
                            &union.cells,
                            union.layout,
                            index.hashes,
                            index.seed ^ UNION_SEED_XOR,
                            index.fingerprint_bits,
                        ) {
                            return false;
                        }
                    }
                }
                index.query_kmer_into_with_scratch(key, hits, scratch);
                true
            },
        )?;
        return Ok(QueryRunResult {
            _stats: result.stats,
            accum: QueryAccum {
                scores: result.scores,
                total_features: result.total_kmers,
                stack_queries: result.backend_queries,
            },
        });
    }
    if is_compressed_path(path) {
        let file =
            File::open(path).with_context(|| format!("unable to open {}", path.display()))?;
        let mut decoder = AnyDecoder::new(file);
        let mut data = Vec::new();
        decoder
            .read_to_end(&mut data)
            .with_context(|| format!("decompressing {}", path.display()))?;
        ensure!(!data.is_empty(), "{} is empty", path.display());
        let parser = FastxParser::<PARSER_CONFIG>::from_slice(&data)
            .with_context(|| format!("parsing {}", path.display()))?;
        return query_features_from_parser_parallel(parser, index, ignore_union);
    }

    let parser = FastxParser::<PARSER_CONFIG>::from_file_mmap(path)
        .with_context(|| format!("unable to mmap {}", path.display()))?;
    query_features_from_parser_parallel(parser, index, ignore_union)
}

fn query_features_from_parser_parallel(
    mut parser: FastxParser<'_, PARSER_CONFIG>,
    index: &ZorIndex,
    ignore_union: bool,
) -> Result<QueryRunResult> {
    let mut stats = FastxStats::default();
    let mut accum = QueryAccum::new(index.colors.len());
    let mut batch = Vec::with_capacity(QUERY_BATCH_FEATURES);

    while parser.next().is_some() {
        let seq = parser.get_dna_string();
        stats.bases += seq.len() as u64;
        stats.chunks += 1;
        scan_index_features(
            seq,
            index.k,
            index.seed,
            index.feature_config,
            &mut |feature| {
                batch.push(feature);
                if batch.len() >= QUERY_BATCH_FEATURES {
                    accum.merge(query_feature_batch(index, &batch, ignore_union));
                    batch.clear();
                }
            },
        );
    }

    if !batch.is_empty() {
        accum.merge(query_feature_batch(index, &batch, ignore_union));
    }

    Ok(QueryRunResult {
        _stats: stats,
        accum,
    })
}

fn query_feature_batch(index: &ZorIndex, features: &[u64], ignore_union: bool) -> QueryAccum {
    if features.is_empty() {
        return QueryAccum::new(index.colors.len());
    }

    let worker_count = rayon::current_num_threads().max(1);
    let chunk_len = (features.len() / (worker_count * 8)).max(1);
    features
        .par_chunks(chunk_len)
        .map(|features| query_feature_chunk(index, features, ignore_union))
        .reduce(
            || QueryAccum::new(index.colors.len()),
            |mut left, right| {
                left.merge(right);
                left
            },
        )
}

fn query_feature_chunk(index: &ZorIndex, features: &[u64], ignore_union: bool) -> QueryAccum {
    let mut accum = QueryAccum::new(index.colors.len());
    let mut scratch = QueryScratch::default();

    for &kmer in features {
        accum.total_features += 1;
        if !ignore_union {
            if let Some(union) = &index.union_graph {
                if !contains_raw(
                    kmer,
                    &union.cells,
                    union.layout,
                    index.hashes,
                    index.seed ^ UNION_SEED_XOR,
                    index.fingerprint_bits,
                ) {
                    continue;
                }
            }
        }
        accum.stack_queries += 1;
        index.query_kmer_into_with_scratch(kmer, &mut accum.scores, &mut scratch);
    }

    accum
}

fn scaled_slots(unique: usize, scale: f64) -> usize {
    if unique == 0 {
        1
    } else {
        ((unique as f64) * scale).ceil() as usize
    }
}

fn calculate_layout(
    target_slots: usize,
    hashes: usize,
    segment_length_override: Option<usize>,
) -> Result<Layout> {
    ensure!((2..=MAX_HASHES).contains(&hashes), "invalid hash count");
    ensure!(target_slots > 0, "target slot count must be positive");

    let segment_length = if let Some(segment_length) = segment_length_override {
        ensure!(
            segment_length.is_power_of_two(),
            "segment length must be a non-zero power of two"
        );
        segment_length
    } else {
        let mut segment_length = segment_length_for(hashes, target_slots).max(1);
        let max_segment_length = 1usize << MAX_SEGMENT_LENGTH_LOG;
        if segment_length > max_segment_length {
            segment_length = max_segment_length;
        }
        while segment_length > target_slots {
            segment_length >>= 1;
            if segment_length == 0 {
                segment_length = 1;
                break;
            }
        }
        segment_length
    };

    calculate_layout_with_segment_length(target_slots, hashes, segment_length)
}

fn calculate_layout_with_segment_length(
    target_slots: usize,
    hashes: usize,
    segment_length: usize,
) -> Result<Layout> {
    ensure!(
        segment_length > 0 && segment_length.is_power_of_two(),
        "segment length must be a non-zero power of two"
    );
    let mut total_segments = target_slots.div_ceil(segment_length);
    if total_segments < hashes {
        total_segments = hashes;
    }
    let segment_count = total_segments.saturating_sub(hashes - 1).max(1);
    let total_segments_with_overlap = segment_count
        .checked_add(hashes - 1)
        .context("layout segment count overflow")?;
    let array_length = segment_length
        .checked_mul(total_segments_with_overlap)
        .context("layout array length overflow")?;
    let segment_count_length = segment_length
        .checked_mul(segment_count)
        .context("layout segment count length overflow")?;

    Ok(Layout {
        segment_length,
        segment_length_mask: segment_length - 1,
        segment_count,
        segment_count_length,
        array_length,
    })
}

fn segment_length_for(hashes: usize, key_count: usize) -> usize {
    let size = cmp::max(key_count, 1) as f64;
    let log_size = size.ln();
    let (base, offset) = if hashes <= 3 {
        (3.33_f64, 2.25_f64)
    } else {
        (2.91_f64, -0.5_f64)
    };
    let shift = (log_size / base.ln() + offset).floor() as i32;
    let clamped = shift.clamp(1, MAX_SEGMENT_LENGTH_LOG as i32);
    1usize << clamped
}

fn validate_layout(layout: Layout, hashes: usize) -> Result<()> {
    ensure!(
        layout.segment_length > 0 && layout.segment_length.is_power_of_two(),
        "invalid layout segment length"
    );
    ensure!(
        layout.segment_length_mask == layout.segment_length - 1,
        "invalid layout segment mask"
    );
    ensure!(layout.segment_count > 0, "invalid layout segment count");
    ensure!(
        layout.segment_count_length == layout.segment_length * layout.segment_count,
        "invalid layout segment count length"
    );
    ensure!(
        layout.array_length == layout.segment_length * (layout.segment_count + hashes - 1),
        "invalid layout array length"
    );
    Ok(())
}

fn peel_zor_graph(keys: &[u64], config: PureZorConfig) -> Result<zorindex::power_merge::PeelGraph> {
    if compact_graph_eligible(keys.len(), config.layout.array_length, config.hashes) {
        Ok(PeelGraph::Compact(peel_zor_graph_impl::<
            u32,
            IndexedDegreeHeap,
        >(keys, config)?))
    } else {
        Ok(PeelGraph::Wide(peel_zor_graph_impl::<
            usize,
            BinaryHeap<(Reverse<u32>, usize)>,
        >(keys, config)?))
    }
}

fn compact_graph_eligible(edges: usize, slots: usize, hashes: usize) -> bool {
    edges <= u32::MAX as usize
        && slots <= u32::MAX as usize
        && edges
            .checked_mul(hashes)
            .is_some_and(|n| n <= u32::MAX as usize)
}

fn peel_zor_graph_impl<I: GraphIndex, Q: DegreeQueue>(
    keys: &[u64],
    config: PureZorConfig,
) -> Result<TypedPeelGraph<I>> {
    ensure!(config.hashes >= 2, "hash count must be at least 2");
    ensure!(
        config.hashes <= MAX_HASHES,
        "hash count must be at most {}",
        MAX_HASHES
    );
    ensure!(config.tie_scan > 0, "tie scan must be greater than 0");
    fingerprint_bytes(config.fingerprint_bits)?;
    validate_layout(config.layout, config.hashes)?;
    if keys.is_empty() {
        return Ok(TypedPeelGraph {
            positions: Vec::new(),
            order: Vec::new(),
            dropped: Vec::new(),
        });
    }

    let edge_count = keys.len();
    let incidence_count = edge_count
        .checked_mul(config.hashes)
        .context("too many key/hash incidences")?;
    let mut edge_slots = vec![I::from_usize(0); incidence_count];
    let mut degrees = vec![0u32; config.layout.array_length];

    for (edge, &key) in keys.iter().enumerate() {
        let mut locations = [0usize; MAX_HASHES];
        fill_positions(
            key,
            config.layout,
            config.hashes,
            config.seed,
            &mut locations[..config.hashes],
        );
        for (stored, &slot) in edge_slots[edge * config.hashes..][..config.hashes]
            .iter_mut()
            .zip(&locations)
        {
            *stored = I::from_usize(slot);
            degrees[slot] += 1;
        }
    }

    let mut offsets = vec![I::from_usize(0); config.layout.array_length + 1];
    for slot in 0..config.layout.array_length {
        offsets[slot + 1] = I::from_usize(offsets[slot].as_usize() + degrees[slot] as usize);
    }
    let mut cursor = offsets[..config.layout.array_length].to_vec();
    let mut adjacency = vec![I::from_usize(0); incidence_count];
    for edge in 0..edge_count {
        for &slot in &edge_slots[edge * config.hashes..edge * config.hashes + config.hashes] {
            let slot = slot.as_usize();
            let pos = cursor[slot].as_usize();
            adjacency[pos] = I::from_usize(edge);
            cursor[slot] = I::from_usize(pos + 1);
        }
    }

    let mut queue = VecDeque::new();
    let mut multi_heap = Q::new(config.layout.array_length);
    for (slot, &degree) in degrees.iter().enumerate() {
        match degree {
            1 => queue.push_back(I::from_usize(slot)),
            degree if degree > 1 => multi_heap.push((Reverse(degree), slot)),
            _ => {}
        }
    }

    let mut active = vec![true; edge_count];
    let mut remaining = edge_count;
    let mut abandoned_edges = Vec::new();
    let mut peel_order = Vec::with_capacity(edge_count);

    while remaining > 0 {
        let mut progress = false;

        while let Some(slot) = queue.pop_front() {
            let slot = slot.as_usize();
            if degrees[slot] != 1 {
                continue;
            }
            let Some(edge) = find_one_active_edge(slot, &offsets, &adjacency, &active) else {
                degrees[slot] = 0;
                continue;
            };
            if !active[edge] {
                continue;
            }
            progress = true;
            active[edge] = false;
            remaining -= 1;
            peel_order.push((I::from_usize(edge), I::from_usize(slot)));
            decrement_edge(
                edge,
                config.hashes,
                &edge_slots,
                &mut degrees,
                &mut queue,
                &mut multi_heap,
            );
        }

        if remaining == 0 {
            break;
        }

        if progress {
            continue;
        }

        let candidate = choose_cycle_break_candidate(
            &mut degrees,
            &offsets,
            &adjacency,
            &active,
            &edge_slots,
            config.hashes,
            config.cycle_break,
            config.tie_scan,
            &mut multi_heap,
        );

        let Some((slot, keep)) = candidate else {
            let Some((edge, _)) = active.iter().enumerate().find(|(_, active)| **active) else {
                break;
            };
            active[edge] = false;
            remaining -= 1;
            abandoned_edges.push(I::from_usize(edge));
            decrement_edge(
                edge,
                config.hashes,
                &edge_slots,
                &mut degrees,
                &mut queue,
                &mut multi_heap,
            );
            continue;
        };

        let mut to_abandon = Vec::new();
        for edge in adjacency[offsets[slot].as_usize()..offsets[slot + 1].as_usize()]
            .iter()
            .map(|&edge| edge.as_usize())
        {
            if edge != keep && active[edge] {
                to_abandon.push(edge);
            }
        }

        if to_abandon.is_empty() {
            if active[keep] {
                degrees[slot] = 1;
                queue.push_back(I::from_usize(slot));
            }
            continue;
        }

        for edge in to_abandon {
            active[edge] = false;
            remaining -= 1;
            abandoned_edges.push(I::from_usize(edge));
            decrement_edge(
                edge,
                config.hashes,
                &edge_slots,
                &mut degrees,
                &mut queue,
                &mut multi_heap,
            );
        }
        if degrees[slot] == 1 {
            queue.push_back(I::from_usize(slot));
        }
    }

    for (edge, is_active) in active.iter_mut().enumerate() {
        if *is_active {
            *is_active = false;
            abandoned_edges.push(I::from_usize(edge));
        }
    }

    Ok(TypedPeelGraph {
        positions: edge_slots,
        order: peel_order,
        dropped: abandoned_edges,
    })
}

fn build_pure_zor(keys: &[u64], config: PureZorConfig) -> Result<PureZorBuild> {
    match peel_zor_graph(keys, config)? {
        PeelGraph::Compact(graph) => finish_pure_zor(keys, config, graph),
        PeelGraph::Wide(graph) => finish_pure_zor(keys, config, graph),
    }
}

fn finish_pure_zor<I: GraphIndex>(
    keys: &[u64],
    config: PureZorConfig,
    graph: TypedPeelGraph<I>,
) -> Result<PureZorBuild> {
    let TypedPeelGraph {
        positions: edge_slots,
        order: peel_order,
        dropped: mut abandoned_edges,
    } = graph;
    let cells = assign_cells(
        keys,
        &edge_slots,
        &peel_order,
        config.layout,
        config.hashes,
        config.seed,
        config.fingerprint_bits,
    )?;

    abandoned_edges.sort_unstable();
    abandoned_edges.dedup();
    let abandoned_keys = abandoned_edges
        .iter()
        .filter_map(|&edge| {
            let key = keys[edge.as_usize()];
            if contains_raw(
                key,
                &cells,
                config.layout,
                config.hashes,
                config.seed,
                config.fingerprint_bits,
            ) {
                None
            } else {
                Some(key)
            }
        })
        .collect::<Vec<_>>();
    let abandoned = abandoned_keys.len() as u64;

    Ok(PureZorBuild {
        cells,
        abandoned,
        abandoned_keys,
    })
}

fn assign_cells<I: GraphIndex>(
    keys: &[u64],
    edge_slots: &[I],
    peel_order: &[(I, I)],
    layout: Layout,
    hashes: usize,
    seed: u64,
    fingerprint_bits: u8,
) -> Result<Vec<u8>> {
    match fingerprint_bits {
        8 => {
            let mut cells = vec![0u8; layout.array_length];
            for &(edge, pivot) in peel_order.iter().rev() {
                let edge = edge.as_usize();
                let pivot = pivot.as_usize();
                let mut value = fingerprint8(keys[edge], seed);
                for &slot in &edge_slots[edge * hashes..edge * hashes + hashes] {
                    let slot = slot.as_usize();
                    if slot != pivot {
                        value ^= cells[slot];
                    }
                }
                cells[pivot] = value;
            }
            Ok(cells)
        }
        16 => {
            let mut cells = vec![0u8; layout.array_length * 2];
            for &(edge, pivot) in peel_order.iter().rev() {
                let edge = edge.as_usize();
                let pivot = pivot.as_usize();
                let mut value = fingerprint16(keys[edge], seed);
                for &slot in &edge_slots[edge * hashes..edge * hashes + hashes] {
                    let slot = slot.as_usize();
                    if slot != pivot {
                        value ^= read_cell_u16(&cells, slot);
                    }
                }
                write_cell_u16(&mut cells, pivot, value);
            }
            Ok(cells)
        }
        _ => anyhow::bail!("unsupported fingerprint width {}", fingerprint_bits),
    }
}

#[inline(always)]
fn decrement_edge<I: GraphIndex, Q: DegreeQueue>(
    edge: usize,
    hashes: usize,
    edge_slots: &[I],
    degrees: &mut [u32],
    queue: &mut VecDeque<I>,
    multi_heap: &mut Q,
) {
    for &slot in &edge_slots[edge * hashes..edge * hashes + hashes] {
        let slot = slot.as_usize();
        if degrees[slot] > 0 {
            degrees[slot] -= 1;
            if degrees[slot] == 1 {
                queue.push_back(I::from_usize(slot));
            }
            multi_heap.update(slot, degrees[slot]);
        }
    }
}

fn find_one_active_edge<I: GraphIndex>(
    slot: usize,
    offsets: &[I],
    adjacency: &[I],
    active: &[bool],
) -> Option<usize> {
    adjacency[offsets[slot].as_usize()..offsets[slot + 1].as_usize()]
        .iter()
        .map(|&edge| edge.as_usize())
        .find(|&edge| active[edge])
}

#[derive(Clone, Copy)]
struct KeyStats {
    sum_degrees: u64,
    max_degree: u32,
    deg2_count: u32,
    degrees: [u32; MAX_HASHES],
    len: usize,
}

#[derive(Clone, Copy)]
struct AbandonStats {
    sum_degrees: u64,
    max_degree: u32,
    deg2_count: u32,
}

fn choose_cycle_break_candidate<I: GraphIndex, Q: DegreeQueue>(
    degrees: &mut [u32],
    offsets: &[I],
    adjacency: &[I],
    active: &[bool],
    edge_slots: &[I],
    hashes: usize,
    cycle_break: CycleBreakHeuristic,
    tie_scan: usize,
    multi_heap: &mut Q,
) -> Option<(usize, usize)> {
    if cycle_break == CycleBreakHeuristic::NoHeuristic {
        loop {
            let (Reverse(recorded_degree), cell) = multi_heap.pop()?;
            let current_degree = degrees[cell];
            if current_degree <= 1 || current_degree != recorded_degree {
                continue;
            }
            let Some(keep_edge) = find_one_active_edge(cell, offsets, adjacency, active) else {
                degrees[cell] = 0;
                continue;
            };
            return Some((cell, keep_edge));
        }
    }

    loop {
        let (Reverse(recorded_degree), cell) = multi_heap.pop()?;
        let current_degree = degrees[cell];
        if current_degree <= 1 || current_degree != recorded_degree {
            continue;
        }

        let mut best_cell = cell;
        let mut best_edge = None;
        let mut best_stats = AbandonStats {
            sum_degrees: 0,
            max_degree: 0,
            deg2_count: 0,
        };
        let mut best_degrees = Vec::new();
        let mut scanned_cells = Vec::new();
        let mut scanned = 0usize;

        if let Some((edge, stats, degrees_vec)) = best_key_for_cell(
            cell,
            active,
            degrees,
            offsets,
            adjacency,
            edge_slots,
            hashes,
            cycle_break,
        ) {
            best_edge = Some(edge);
            best_stats = stats;
            best_degrees = degrees_vec;
            scanned += 1;
        } else {
            degrees[cell] = 0;
        }

        while scanned < tie_scan {
            let Some((Reverse(next_degree), _)) = multi_heap.peek() else {
                break;
            };
            if next_degree != recorded_degree {
                break;
            }
            let (_, other_cell) = multi_heap.pop().unwrap();
            let current_degree = degrees[other_cell];
            if current_degree <= 1 || current_degree != recorded_degree {
                continue;
            }
            if let Some((edge, stats, degrees_vec)) = best_key_for_cell(
                other_cell,
                active,
                degrees,
                offsets,
                adjacency,
                edge_slots,
                hashes,
                cycle_break,
            ) {
                if best_edge.is_some() {
                    if better_abandon(
                        cycle_break,
                        &stats,
                        &degrees_vec,
                        &best_stats,
                        &best_degrees,
                    ) {
                        scanned_cells.push((Reverse(recorded_degree), best_cell));
                        best_cell = other_cell;
                        best_edge = Some(edge);
                        best_stats = stats;
                        best_degrees = degrees_vec;
                    } else {
                        scanned_cells.push((Reverse(recorded_degree), other_cell));
                    }
                } else {
                    best_cell = other_cell;
                    best_edge = Some(edge);
                    best_stats = stats;
                    best_degrees = degrees_vec;
                }
                scanned += 1;
            } else {
                degrees[other_cell] = 0;
            }
        }

        for entry in scanned_cells {
            multi_heap.push(entry);
        }

        if let Some(best_edge) = best_edge {
            return Some((best_cell, best_edge));
        }
    }
}

fn best_key_for_cell<I: GraphIndex>(
    cell: usize,
    active: &[bool],
    degrees: &[u32],
    offsets: &[I],
    adjacency: &[I],
    edge_slots: &[I],
    hashes: usize,
    cycle_break: CycleBreakHeuristic,
) -> Option<(usize, AbandonStats, Vec<u32>)> {
    let mut active_edges = Vec::new();
    let mut total_sum = 0u64;
    let mut total_deg2 = 0u32;
    let mut max_degree = 0u32;
    let mut max_count = 0u32;
    let mut second_max = 0u32;

    for edge in adjacency[offsets[cell].as_usize()..offsets[cell + 1].as_usize()]
        .iter()
        .map(|&edge| edge.as_usize())
    {
        if !active[edge] {
            continue;
        }
        let stats = key_stats(edge, degrees, edge_slots, hashes);
        total_sum += stats.sum_degrees;
        total_deg2 += stats.deg2_count;
        if stats.max_degree > max_degree {
            second_max = max_degree;
            max_degree = stats.max_degree;
            max_count = 1;
        } else if stats.max_degree == max_degree {
            max_count += 1;
        } else if stats.max_degree > second_max {
            second_max = stats.max_degree;
        }
        active_edges.push((edge, stats));
    }

    if active_edges.is_empty() {
        return None;
    }

    let mut best_edge = None;
    let mut best_stats = AbandonStats {
        sum_degrees: 0,
        max_degree: 0,
        deg2_count: 0,
    };
    let mut best_degrees = Vec::new();

    for (edge, stats) in active_edges.iter().copied() {
        let abandon_max = if stats.max_degree == max_degree && max_count == 1 {
            second_max
        } else {
            max_degree
        };
        let abandon_stats = AbandonStats {
            sum_degrees: total_sum - stats.sum_degrees,
            max_degree: abandon_max,
            deg2_count: total_deg2 - stats.deg2_count,
        };

        let mut abandon_degrees = Vec::new();
        if cycle_break == CycleBreakHeuristic::MostDeg2 {
            abandon_degrees.reserve((active_edges.len().saturating_sub(1)) * hashes);
            for (other_edge, other_stats) in &active_edges {
                if *other_edge != edge {
                    abandon_degrees.extend_from_slice(&other_stats.degrees[..other_stats.len]);
                }
            }
            abandon_degrees.sort_unstable();
        }

        if best_edge.is_none()
            || better_abandon(
                cycle_break,
                &abandon_stats,
                &abandon_degrees,
                &best_stats,
                &best_degrees,
            )
        {
            best_edge = Some(edge);
            best_stats = abandon_stats;
            best_degrees = abandon_degrees;
        }
    }

    best_edge.map(|edge| (edge, best_stats, best_degrees))
}

fn key_stats<I: GraphIndex>(
    edge: usize,
    degrees: &[u32],
    edge_slots: &[I],
    hashes: usize,
) -> KeyStats {
    let mut sum_degrees = 0u64;
    let mut max_degree = 0u32;
    let mut deg2_count = 0u32;
    let mut degrees_buf = [0u32; MAX_HASHES];

    for (idx, &slot) in edge_slots[edge * hashes..edge * hashes + hashes]
        .iter()
        .enumerate()
    {
        let degree = degrees[slot.as_usize()];
        sum_degrees += degree as u64;
        max_degree = max_degree.max(degree);
        if degree == 2 {
            deg2_count += 1;
        }
        degrees_buf[idx] = degree;
    }
    degrees_buf[..hashes].sort_unstable();

    KeyStats {
        sum_degrees,
        max_degree,
        deg2_count,
        degrees: degrees_buf,
        len: hashes,
    }
}

fn better_abandon(
    cycle_break: CycleBreakHeuristic,
    candidate: &AbandonStats,
    candidate_degrees: &[u32],
    best: &AbandonStats,
    best_degrees: &[u32],
) -> bool {
    match cycle_break {
        CycleBreakHeuristic::NoHeuristic => false,
        CycleBreakHeuristic::Lightest => {
            candidate.sum_degrees < best.sum_degrees
                || (candidate.sum_degrees == best.sum_degrees
                    && (candidate.deg2_count > best.deg2_count
                        || (candidate.deg2_count == best.deg2_count
                            && candidate.max_degree < best.max_degree)))
        }
        CycleBreakHeuristic::Heaviest => {
            candidate.sum_degrees > best.sum_degrees
                || (candidate.sum_degrees == best.sum_degrees
                    && (candidate.deg2_count < best.deg2_count
                        || (candidate.deg2_count == best.deg2_count
                            && candidate.max_degree > best.max_degree)))
        }
        CycleBreakHeuristic::MostDeg2 => {
            let len = candidate_degrees.len().min(best_degrees.len());
            for idx in 0..len {
                let candidate_degree = candidate_degrees[idx];
                let best_degree = best_degrees[idx];
                if candidate_degree != best_degree {
                    return candidate_degree < best_degree;
                }
            }
            if candidate.sum_degrees != best.sum_degrees {
                candidate.sum_degrees < best.sum_degrees
            } else {
                candidate.max_degree < best.max_degree
            }
        }
        CycleBreakHeuristic::MinMaxDegree => {
            candidate.max_degree < best.max_degree
                || (candidate.max_degree == best.max_degree
                    && (candidate.sum_degrees < best.sum_degrees
                        || (candidate.sum_degrees == best.sum_degrees
                            && candidate.deg2_count > best.deg2_count)))
        }
    }
}

fn contains_raw(
    key: u64,
    cells: &[u8],
    layout: Layout,
    hashes: usize,
    seed: u64,
    fingerprint_bits: u8,
) -> bool {
    if cells.is_empty() {
        return false;
    }
    let mut slots = [0usize; MAX_HASHES];
    fill_positions(key, layout, hashes, seed, &mut slots[..hashes]);
    match fingerprint_bits {
        8 => {
            let mut value = 0u8;
            for &slot in &slots[..hashes] {
                value ^= cells[slot];
            }
            value == fingerprint8(key, seed)
        }
        16 => {
            let mut value = 0u16;
            for &slot in &slots[..hashes] {
                value ^= read_cell_u16(cells, slot);
            }
            value == fingerprint16(key, seed)
        }
        _ => false,
    }
}

#[inline(always)]
fn fill_positions(key: u64, layout: Layout, hashes: usize, seed: u64, out: &mut [usize]) {
    debug_assert!(out.len() >= hashes);
    if hashes == 0 || layout.array_length == 0 {
        return;
    }
    let hash1 = mixsplit(key, seed);
    let hash2 = splitmix64(hash1);
    let segment_count = layout.segment_count as u64;
    let base_segment = (((hash1 as u128) * (segment_count as u128)) >> 64) as u64;
    let segment_length = layout.segment_length as u64;
    let base_offset = base_segment * segment_length;
    let mask = layout.segment_length_mask as u64;

    let mut hash = hash1;
    for (idx, slot) in out.iter_mut().take(hashes).enumerate() {
        let variation = (hash ^ (hash >> 33)) & mask;
        *slot = (base_offset + (idx as u64) * segment_length + variation) as usize;
        hash = hash.wrapping_add(hash2);
    }
}

#[inline(always)]
fn fingerprint8(key: u64, seed: u64) -> u8 {
    splitmix64(key ^ seed ^ FINGERPRINT_SEED) as u8
}

#[inline(always)]
fn fingerprint16(key: u64, seed: u64) -> u16 {
    splitmix64(key ^ seed ^ FINGERPRINT_SEED) as u16
}

#[inline(always)]
fn repeated_byte_u64(byte: u8) -> u64 {
    (byte as u64) * 0x0101_0101_0101_0101
}

#[inline(always)]
fn repeated_u16_u64(value: u16) -> u64 {
    (value as u64) * 0x0001_0001_0001_0001
}

#[inline(always)]
fn read_u64_le(bytes: &[u8], offset: usize) -> u64 {
    debug_assert!(offset + 8 <= bytes.len());
    // SAFETY: the debug assertion documents the required bounds, and `[u8; 8]`
    // has alignment 1. `read_unaligned` keeps this valid for any row offset.
    let array = unsafe { std::ptr::read_unaligned(bytes.as_ptr().add(offset) as *const [u8; 8]) };
    u64::from_le_bytes(array)
}

#[inline(always)]
fn read_u32_le(bytes: &[u8], offset: usize) -> u32 {
    debug_assert!(offset + 4 <= bytes.len());
    let array = unsafe { std::ptr::read_unaligned(bytes.as_ptr().add(offset) as *const [u8; 4]) };
    u32::from_le_bytes(array)
}

#[inline(always)]
fn read_u16_le(bytes: &[u8], offset: usize) -> u16 {
    debug_assert!(offset + 2 <= bytes.len());
    u16::from_le_bytes([bytes[offset], bytes[offset + 1]])
}

#[inline(always)]
fn read_cell_u16(cells: &[u8], slot: usize) -> u16 {
    read_u16_le(cells, slot * 2)
}

#[inline(always)]
fn write_cell_u16(cells: &mut [u8], slot: usize, value: u16) {
    let offset = slot * 2;
    cells[offset..offset + 2].copy_from_slice(&value.to_le_bytes());
}

#[inline(always)]
fn increment_matching_bytes(diff: u64, scores: &mut [u64], offset: usize) {
    for lane in 0..8 {
        if ((diff >> (lane * 8)) & 0xFF) == 0 {
            scores[offset + lane] += 1;
        }
    }
}

#[inline(always)]
fn increment_matching_u16_lanes(diff: u64, scores: &mut [u64], offset: usize) {
    for lane in 0..4 {
        if ((diff >> (lane * 16)) & 0xFFFF) == 0 {
            scores[offset + lane] += 1;
        }
    }
}

#[derive(Default)]
struct QueryScratch {
    acc: Vec<u8>,
    row_codec: RowCodecScratch,
    seen_users: Vec<u32>,
    epoch: u32,
}

impl QueryScratch {
    fn ensure(&mut self, row_bytes: usize) {
        if self.acc.len() != row_bytes {
            self.acc.resize(row_bytes, 0);
        }
    }

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

#[inline(always)]
fn xor_bytes_into(acc: &mut [u8], row: &[u8]) {
    debug_assert_eq!(acc.len(), row.len());
    let mut offset = 0usize;
    while offset + 8 <= acc.len() {
        let value = read_u64_le(acc, offset) ^ read_u64_le(row, offset);
        acc[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
        offset += 8;
    }
    while offset < acc.len() {
        acc[offset] ^= row[offset];
        offset += 1;
    }
}

impl ZorIndex {
    fn technical_bin_count(&self) -> usize {
        self.technical_bins.len()
    }

    fn is_flat_technical_layout(&self) -> bool {
        self.technical_bins.len() == self.colors.len()
            && self
                .technical_bins
                .iter()
                .enumerate()
                .all(|(idx, bin)| bin.color_id == idx && bin.part_id == 0 && bin.part_count == 1)
    }

    #[cfg(test)]
    fn query_kmer_into(&self, key: u64, scores: &mut [u64]) {
        let mut scratch = QueryScratch::default();
        self.query_kmer_into_with_scratch(key, scores, &mut scratch);
    }

    fn query_kmer_into_with_scratch(
        &self,
        key: u64,
        scores: &mut [u64],
        scratch: &mut QueryScratch,
    ) {
        debug_assert_eq!(scores.len(), self.colors.len());
        if !self.is_flat_technical_layout() {
            match self.fingerprint_bits {
                8 => self.query_kmer_into_mapped_8(key, scores, scratch),
                16 => self.query_kmer_into_mapped_16(key, scores, scratch),
                _ => unreachable!("invalid fingerprint width in loaded index"),
            }
            return;
        }
        match self.fingerprint_bits {
            8 => self.query_kmer_into_8(key, scores, scratch),
            16 => self.query_kmer_into_16(key, scores, scratch),
            _ => unreachable!("invalid fingerprint width in loaded index"),
        }
    }

    fn query_kmer_into_mapped_8(&self, key: u64, scores: &mut [u64], scratch: &mut QueryScratch) {
        let mut slots = [0usize; MAX_HASHES];
        fill_positions(
            key,
            self.layout,
            self.hashes,
            self.seed,
            &mut slots[..self.hashes],
        );
        let target = fingerprint8(key, self.seed);
        let row_bytes = self.technical_bin_count();
        scratch.ensure(row_bytes);
        scratch.acc.fill(0);
        if let Some(interlaced) = self.stack.as_plain() {
            for &slot in &slots[..self.hashes] {
                let start = slot * row_bytes;
                xor_bytes_into(&mut scratch.acc, &interlaced[start..start + row_bytes]);
            }
        } else {
            for &slot in &slots[..self.hashes] {
                self.stack
                    .xor_row_into(slot, &mut scratch.acc, &mut scratch.row_codec)
                    .expect("failed to decode compressed stack slice");
            }
        }
        self.score_mapped_row_8(target, scores, scratch);
    }

    fn query_kmer_into_mapped_16(&self, key: u64, scores: &mut [u64], scratch: &mut QueryScratch) {
        let mut slots = [0usize; MAX_HASHES];
        fill_positions(
            key,
            self.layout,
            self.hashes,
            self.seed,
            &mut slots[..self.hashes],
        );
        let target = fingerprint16(key, self.seed);
        let row_bytes = self.technical_bin_count() * 2;
        scratch.ensure(row_bytes);
        scratch.acc.fill(0);
        if let Some(interlaced) = self.stack.as_plain() {
            for &slot in &slots[..self.hashes] {
                let start = slot * row_bytes;
                xor_bytes_into(&mut scratch.acc, &interlaced[start..start + row_bytes]);
            }
        } else {
            for &slot in &slots[..self.hashes] {
                self.stack
                    .xor_row_into(slot, &mut scratch.acc, &mut scratch.row_codec)
                    .expect("failed to decode compressed stack slice");
            }
        }
        self.score_mapped_row_16(target, scores, scratch);
    }

    fn score_mapped_row_8(&self, target: u8, scores: &mut [u64], scratch: &mut QueryScratch) {
        let epoch = scratch.next_epoch(self.colors.len());
        for (technical_id, &value) in scratch.acc.iter().enumerate() {
            if value != target {
                continue;
            }
            let color_id = self.technical_bins[technical_id].color_id;
            if scratch.seen_users[color_id] != epoch {
                scratch.seen_users[color_id] = epoch;
                scores[color_id] += 1;
            }
        }
    }

    fn score_mapped_row_16(&self, target: u16, scores: &mut [u64], scratch: &mut QueryScratch) {
        let epoch = scratch.next_epoch(self.colors.len());
        for technical_id in 0..self.technical_bin_count() {
            if read_u16_le(&scratch.acc, technical_id * 2) != target {
                continue;
            }
            let color_id = self.technical_bins[technical_id].color_id;
            if scratch.seen_users[color_id] != epoch {
                scratch.seen_users[color_id] = epoch;
                scores[color_id] += 1;
            }
        }
    }

    fn query_kmer_into_8(&self, key: u64, scores: &mut [u64], scratch: &mut QueryScratch) {
        let mut slots = [0usize; MAX_HASHES];
        fill_positions(
            key,
            self.layout,
            self.hashes,
            self.seed,
            &mut slots[..self.hashes],
        );
        let target = fingerprint8(key, self.seed);
        let slots = &slots[..self.hashes];
        let mut offset = 0usize;

        if let Some(interlaced) = self.stack.as_plain() {
            #[cfg(target_arch = "x86_64")]
            {
                // SAFETY: x86_64 guarantees SSE2. The helper only reads complete
                // 16-byte chunks from rows of length `color_count`.
                offset =
                    unsafe { self.query_kmer_into_sse2(interlaced, slots, target, scores, offset) };
            }

            offset = self.query_kmer_into_word_chunks(interlaced, slots, target, scores, offset);
            self.query_kmer_into_scalar_tail(interlaced, slots, target, scores, offset);
        } else {
            self.query_kmer_into_compressed_8(slots, target, scores, scratch);
        }
    }

    fn query_kmer_into_16(&self, key: u64, scores: &mut [u64], scratch: &mut QueryScratch) {
        let mut slots = [0usize; MAX_HASHES];
        fill_positions(
            key,
            self.layout,
            self.hashes,
            self.seed,
            &mut slots[..self.hashes],
        );
        let target = fingerprint16(key, self.seed);
        let slots = &slots[..self.hashes];
        let mut offset = 0usize;

        if let Some(interlaced) = self.stack.as_plain() {
            #[cfg(target_arch = "x86_64")]
            {
                // SAFETY: x86_64 guarantees SSE2. The helper only reads complete
                // 16-byte chunks, i.e. 8 little-endian u16 color cells.
                offset = unsafe {
                    self.query_kmer_into_sse2_16(interlaced, slots, target, scores, offset)
                };
            }

            offset =
                self.query_kmer_into_u16_word_chunks(interlaced, slots, target, scores, offset);
            self.query_kmer_into_u16_scalar_tail(interlaced, slots, target, scores, offset);
        } else {
            self.query_kmer_into_compressed_16(slots, target, scores, scratch);
        }
    }

    fn query_kmer_into_compressed_8(
        &self,
        slots: &[usize],
        target: u8,
        scores: &mut [u64],
        scratch: &mut QueryScratch,
    ) {
        let row_bytes = self.colors.len();
        scratch.ensure(row_bytes);
        scratch.acc.fill(0);
        for &slot in slots {
            self.stack
                .xor_row_into(slot, &mut scratch.acc, &mut scratch.row_codec)
                .expect("failed to decode compressed stack slice");
        }
        self.score_decoded_row_8(&scratch.acc, target, scores);
    }

    fn query_kmer_into_compressed_16(
        &self,
        slots: &[usize],
        target: u16,
        scores: &mut [u64],
        scratch: &mut QueryScratch,
    ) {
        let row_bytes = self.colors.len() * 2;
        scratch.ensure(row_bytes);
        scratch.acc.fill(0);
        for &slot in slots {
            self.stack
                .xor_row_into(slot, &mut scratch.acc, &mut scratch.row_codec)
                .expect("failed to decode compressed stack slice");
        }
        self.score_decoded_row_16(&scratch.acc, target, scores);
    }

    fn score_decoded_row_8(&self, row: &[u8], target: u8, scores: &mut [u64]) {
        let mut offset = 0usize;
        #[cfg(target_arch = "x86_64")]
        {
            // SAFETY: x86_64 guarantees SSE2. The helper only reads complete
            // 16-byte chunks from the decoded accumulator row.
            offset = unsafe { self.score_decoded_row_sse2(row, target, scores, offset) };
        }
        let color_count = self.colors.len();
        let target_word = repeated_byte_u64(target);
        while offset + 8 <= color_count {
            increment_matching_bytes(read_u64_le(row, offset) ^ target_word, scores, offset);
            offset += 8;
        }
        for color in offset..color_count {
            if row[color] == target {
                scores[color] += 1;
            }
        }
    }

    fn score_decoded_row_16(&self, row: &[u8], target: u16, scores: &mut [u64]) {
        let mut offset = 0usize;
        #[cfg(target_arch = "x86_64")]
        {
            // SAFETY: x86_64 guarantees SSE2. The helper only reads complete
            // 16-byte chunks, i.e. 8 little-endian u16 color cells.
            offset = unsafe { self.score_decoded_row_sse2_16(row, target, scores, offset) };
        }
        let color_count = self.colors.len();
        let target_word = repeated_u16_u64(target);
        while offset + 4 <= color_count {
            increment_matching_u16_lanes(
                read_u64_le(row, offset * 2) ^ target_word,
                scores,
                offset,
            );
            offset += 4;
        }
        for color in offset..color_count {
            if read_u16_le(row, color * 2) == target {
                scores[color] += 1;
            }
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "sse2")]
    unsafe fn score_decoded_row_sse2(
        &self,
        row: &[u8],
        target: u8,
        scores: &mut [u64],
        mut offset: usize,
    ) -> usize {
        use std::arch::x86_64::{
            __m128i, _mm_cmpeq_epi8, _mm_loadu_si128, _mm_movemask_epi8, _mm_set1_epi8,
        };

        let color_count = self.colors.len();
        let end = color_count & !15;
        let target_vec = _mm_set1_epi8(target as i8);
        while offset < end {
            let ptr = row.as_ptr().add(offset);
            let value = _mm_loadu_si128(ptr as *const __m128i);
            let matches = _mm_cmpeq_epi8(value, target_vec);
            let mut mask = _mm_movemask_epi8(matches) as u32;
            while mask != 0 {
                let lane = mask.trailing_zeros() as usize;
                scores[offset + lane] += 1;
                mask &= mask - 1;
            }
            offset += 16;
        }
        offset
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "sse2")]
    unsafe fn score_decoded_row_sse2_16(
        &self,
        row: &[u8],
        target: u16,
        scores: &mut [u64],
        mut offset: usize,
    ) -> usize {
        use std::arch::x86_64::{
            __m128i, _mm_cmpeq_epi16, _mm_loadu_si128, _mm_movemask_epi8, _mm_set1_epi16,
        };

        let color_count = self.colors.len();
        let end = color_count & !7;
        let target_vec = _mm_set1_epi16(target as i16);
        while offset < end {
            let ptr = row.as_ptr().add(offset * 2);
            let value = _mm_loadu_si128(ptr as *const __m128i);
            let matches = _mm_cmpeq_epi16(value, target_vec);
            let mask = _mm_movemask_epi8(matches) as u32;
            for lane in 0..8 {
                let lane_bits = 0b11 << (lane * 2);
                if mask & lane_bits == lane_bits {
                    scores[offset + lane] += 1;
                }
            }
            offset += 8;
        }
        offset
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "sse2")]
    unsafe fn query_kmer_into_sse2(
        &self,
        interlaced: &[u8],
        slots: &[usize],
        target: u8,
        scores: &mut [u64],
        mut offset: usize,
    ) -> usize {
        use std::arch::x86_64::{
            __m128i, _mm_cmpeq_epi8, _mm_loadu_si128, _mm_movemask_epi8, _mm_set1_epi8,
            _mm_setzero_si128, _mm_xor_si128,
        };

        let color_count = self.colors.len();
        let end = color_count & !15;
        let target_vec = _mm_set1_epi8(target as i8);

        while offset < end {
            let mut value = _mm_setzero_si128();
            for &slot in slots {
                let ptr = interlaced.as_ptr().add(slot * color_count + offset);
                let row = _mm_loadu_si128(ptr as *const __m128i);
                value = _mm_xor_si128(value, row);
            }
            let matches = _mm_cmpeq_epi8(value, target_vec);
            let mut mask = _mm_movemask_epi8(matches) as u32;
            while mask != 0 {
                let lane = mask.trailing_zeros() as usize;
                scores[offset + lane] += 1;
                mask &= mask - 1;
            }
            offset += 16;
        }

        offset
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "sse2")]
    unsafe fn query_kmer_into_sse2_16(
        &self,
        interlaced: &[u8],
        slots: &[usize],
        target: u16,
        scores: &mut [u64],
        mut offset: usize,
    ) -> usize {
        use std::arch::x86_64::{
            __m128i, _mm_cmpeq_epi16, _mm_loadu_si128, _mm_movemask_epi8, _mm_set1_epi16,
            _mm_setzero_si128, _mm_xor_si128,
        };

        let color_count = self.colors.len();
        let row_bytes = color_count * 2;
        let end = color_count & !7;
        let target_vec = _mm_set1_epi16(target as i16);

        while offset < end {
            let mut value = _mm_setzero_si128();
            for &slot in slots {
                let ptr = interlaced.as_ptr().add(slot * row_bytes + offset * 2);
                let row = _mm_loadu_si128(ptr as *const __m128i);
                value = _mm_xor_si128(value, row);
            }
            let matches = _mm_cmpeq_epi16(value, target_vec);
            let mask = _mm_movemask_epi8(matches) as u32;
            for lane in 0..8 {
                let lane_bits = 0b11 << (lane * 2);
                if mask & lane_bits == lane_bits {
                    scores[offset + lane] += 1;
                }
            }
            offset += 8;
        }

        offset
    }

    fn query_kmer_into_word_chunks(
        &self,
        interlaced: &[u8],
        slots: &[usize],
        target: u8,
        scores: &mut [u64],
        mut offset: usize,
    ) -> usize {
        let color_count = self.colors.len();
        let target_word = repeated_byte_u64(target);
        while offset + 8 <= color_count {
            let mut value = 0u64;
            for &slot in slots {
                value ^= read_u64_le(interlaced, slot * color_count + offset);
            }
            increment_matching_bytes(value ^ target_word, scores, offset);
            offset += 8;
        }
        offset
    }

    fn query_kmer_into_scalar_tail(
        &self,
        interlaced: &[u8],
        slots: &[usize],
        target: u8,
        scores: &mut [u64],
        offset: usize,
    ) {
        let color_count = self.colors.len();
        for color in offset..color_count {
            let mut value = 0u8;
            for &slot in slots {
                value ^= interlaced[slot * color_count + color];
            }
            if value == target {
                scores[color] += 1;
            }
        }
    }

    fn query_kmer_into_u16_word_chunks(
        &self,
        interlaced: &[u8],
        slots: &[usize],
        target: u16,
        scores: &mut [u64],
        mut offset: usize,
    ) -> usize {
        let color_count = self.colors.len();
        let row_bytes = color_count * 2;
        let target_word = repeated_u16_u64(target);
        while offset + 4 <= color_count {
            let mut value = 0u64;
            for &slot in slots {
                value ^= read_u64_le(interlaced, slot * row_bytes + offset * 2);
            }
            increment_matching_u16_lanes(value ^ target_word, scores, offset);
            offset += 4;
        }
        offset
    }

    fn query_kmer_into_u16_scalar_tail(
        &self,
        interlaced: &[u8],
        slots: &[usize],
        target: u16,
        scores: &mut [u64],
        offset: usize,
    ) {
        let color_count = self.colors.len();
        let row_bytes = color_count * 2;
        for color in offset..color_count {
            let mut value = 0u16;
            for &slot in slots {
                value ^= read_u16_le(interlaced, slot * row_bytes + color * 2);
            }
            if value == target {
                scores[color] += 1;
            }
        }
    }

    fn save(&self, path: &Path) -> Result<()> {
        let mut writer = create_index_writer(path)?;
        writer.write_all(MAGIC)?;
        write_u32(&mut writer, FORMAT_VERSION)?;
        write_u8(&mut writer, self.k as u8)?;
        write_u8(&mut writer, self.hashes as u8)?;
        write_u8(&mut writer, self.fingerprint_bits)?;
        write_u8(
            &mut writer,
            if self.union_graph.is_some() {
                FLAG_UNION_GRAPH
            } else {
                0
            },
        )?;
        write_u64(&mut writer, self.seed)?;
        write_u8(&mut writer, self.feature_config.mode.as_u8())?;
        write_u8(&mut writer, self.feature_config.minimizer_size as u8)?;
        write_u64(&mut writer, self.feature_config.modimizer_sampling)?;
        // Mode-specific extension: legacy modes keep their original byte layout.
        // Old readers reject mode 3 rather than misinterpreting a findere index.
        if self.feature_config.mode == IndexMode::Findere {
            write_u8(&mut writer, self.feature_config.findere_z as u8)?;
        }
        write_layout(&mut writer, self.layout)?;
        write_u64(&mut writer, self.colors.len() as u64)?;
        for color in &self.colors {
            write_string(&mut writer, &color.name)?;
            write_string(&mut writer, &color.path)?;
            write_u64(&mut writer, color.occurrences)?;
            write_u64(&mut writer, color.unique_kmers)?;
            write_u64(&mut writer, color.abandoned)?;
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
        write_stack(&mut writer, &self.stack)?;
        if let Some(union) = &self.union_graph {
            write_layout(&mut writer, union.layout)?;
            write_u64(&mut writer, union.unique_kmers)?;
            write_u64(&mut writer, union.abandoned)?;
            write_u64(&mut writer, union.cells.len() as u64)?;
            writer.write_all(&union.cells)?;
        }
        finish_index_writer(writer)
    }

    fn load(path: &Path) -> Result<Self> {
        let mut reader = open_index_reader(path)?;
        let mut magic = [0u8; 8];
        reader.read_exact(&mut magic)?;
        ensure!(&magic == MAGIC, "not a ZORINDEX file");
        let version = read_u32(&mut reader)?;
        ensure!(
            version == 2
                || version == 3
                || version == 4
                || version == 5
                || version == FORMAT_VERSION,
            "unsupported index version {}",
            version
        );
        let k = read_u8(&mut reader)? as usize;
        let hashes = read_u8(&mut reader)? as usize;
        let fingerprint_bits = read_u8(&mut reader)?;
        let fingerprint_bytes = fingerprint_bytes(fingerprint_bits)?;
        let flags = read_u8(&mut reader)?;
        let seed = read_u64(&mut reader)?;
        let feature_config = if version >= 4 {
            let mode = IndexMode::from_u8(read_u8(&mut reader)?)?;
            let minimizer_size = read_u8(&mut reader)? as usize;
            let modimizer_sampling = read_u64(&mut reader)?;
            let findere_z = if mode == IndexMode::Findere {
                read_u8(&mut reader)? as usize
            } else {
                zorindex::findere::DEFAULT_Z
            };
            FeatureConfig {
                mode,
                minimizer_size,
                modimizer_sampling,
                findere_z,
            }
        } else {
            FeatureConfig::legacy_kmers()
        };
        let layout = read_layout(&mut reader)?;
        let color_count = read_u64(&mut reader)? as usize;
        ensure!(k > 0, "invalid k in index");
        ensure!((2..=32).contains(&hashes), "invalid hash count in index");
        validate_feature_config(k, feature_config)?;
        validate_layout(layout, hashes)?;

        let mut colors = Vec::with_capacity(color_count);
        for _ in 0..color_count {
            colors.push(ColorMeta {
                name: read_string(&mut reader)?,
                path: read_string(&mut reader)?,
                occurrences: read_u64(&mut reader)?,
                unique_kmers: read_u64(&mut reader)?,
                abandoned: read_u64(&mut reader)?,
            });
        }
        let (technical_layout, split_target_features, technical_bins) = if version >= 7 {
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
        let stack = read_stack(
            &mut reader,
            version,
            layout,
            technical_bins.len(),
            fingerprint_bytes,
        )?;

        let union_graph = if flags & FLAG_UNION_GRAPH != 0 {
            let union_layout = read_layout(&mut reader)?;
            validate_layout(union_layout, hashes)?;
            let unique_kmers = read_u64(&mut reader)?;
            let abandoned = read_u64(&mut reader)?;
            let cell_len = read_u64(&mut reader)? as usize;
            ensure!(
                cell_len == union_layout.array_length * fingerprint_bytes,
                "corrupt index: union cells length does not match slots"
            );
            let mut cells = vec![0u8; cell_len];
            reader.read_exact(&mut cells)?;
            Some(UnionGraph {
                layout: union_layout,
                unique_kmers,
                abandoned,
                cells,
            })
        } else {
            None
        };

        Ok(Self {
            k,
            hashes,
            fingerprint_bits,
            seed,
            feature_config,
            layout,
            technical_layout,
            split_target_features,
            colors,
            technical_bins,
            stack,
            union_graph,
        })
    }
}

fn write_stack<W: Write>(writer: &mut W, stack: &StackStorage) -> Result<()> {
    write_u8(writer, stack.compression().as_u8())?;
    write_u64(writer, stack.len_uncompressed() as u64)?;
    match stack {
        StackStorage::Plain(data) => {
            write_u64(writer, data.len() as u64)?;
            writer.write_all(data)?;
        }
        StackStorage::Compressed(compressed) => {
            write_u64(writer, compressed.row_bytes as u64)?;
            write_u64(writer, compressed.fingerprint_bytes as u64)?;
            write_u64(writer, compressed.rows as u64)?;
            write_u64(writer, compressed.offsets.len() as u64)?;
            write_u64(writer, compressed.data.len() as u64)?;
            for &offset in &compressed.offsets {
                write_u64(writer, offset)?;
            }
            writer.write_all(&compressed.data)?;
        }
    }
    Ok(())
}

fn read_stack<R: Read>(
    reader: &mut R,
    version: u32,
    layout: Layout,
    color_count: usize,
    fingerprint_bytes: usize,
) -> Result<StackStorage> {
    let row_bytes = color_count
        .checked_mul(fingerprint_bytes)
        .context("index dimensions overflow")?;
    let expected_len = layout
        .array_length
        .checked_mul(row_bytes)
        .context("index dimensions overflow")?;

    if version < 5 {
        let interlaced_len = read_u64(reader)? as usize;
        ensure!(
            interlaced_len == expected_len,
            "corrupt index: interlaced length {} != expected {}",
            interlaced_len,
            expected_len
        );
        let mut interlaced = vec![0u8; interlaced_len];
        reader.read_exact(&mut interlaced)?;
        return Ok(StackStorage::Plain(interlaced));
    }
    if version == 5 {
        return read_v5_stack(reader, layout, row_bytes, expected_len);
    }

    let codec = StackCompression::from_u8(read_u8(reader)?)?;
    let uncompressed_len = read_u64(reader)? as usize;
    ensure!(
        uncompressed_len == expected_len,
        "corrupt index: stack length {} != expected {}",
        uncompressed_len,
        expected_len
    );
    if codec == StackCompression::None {
        let stored_len = read_u64(reader)? as usize;
        ensure!(
            stored_len == expected_len,
            "corrupt index: plain stack length {} != expected {}",
            stored_len,
            expected_len
        );
        let mut data = vec![0u8; stored_len];
        reader.read_exact(&mut data)?;
        return Ok(StackStorage::Plain(data));
    }

    let stored_row_bytes = read_u64(reader)? as usize;
    let stored_fingerprint_bytes = read_u64(reader)? as usize;
    let rows = read_u64(reader)? as usize;
    let offset_count = read_u64(reader)? as usize;
    let data_len = read_u64(reader)? as usize;
    ensure!(
        stored_row_bytes == row_bytes,
        "corrupt index: stack row bytes {} != expected {}",
        stored_row_bytes,
        row_bytes
    );
    ensure!(
        stored_fingerprint_bytes == fingerprint_bytes,
        "corrupt index: stack fingerprint bytes {} != expected {}",
        stored_fingerprint_bytes,
        fingerprint_bytes
    );
    ensure!(
        rows == layout.array_length,
        "corrupt index: stack rows {} != expected {}",
        rows,
        layout.array_length
    );
    ensure!(
        offset_count == rows + 1,
        "corrupt index: stack offset count {} != expected {}",
        offset_count,
        rows + 1
    );

    let mut offsets = Vec::with_capacity(offset_count);
    let mut previous = 0u64;
    for idx in 0..offset_count {
        let offset = read_u64(reader)?;
        ensure!(
            offset >= previous,
            "corrupt index: stack offset {} is not monotonic",
            idx
        );
        ensure!(
            offset <= data_len as u64,
            "corrupt index: stack offset {} exceeds data length",
            idx
        );
        previous = offset;
        offsets.push(offset);
    }
    ensure!(
        offsets.first() == Some(&0) && offsets.last() == Some(&(data_len as u64)),
        "corrupt index: stack offsets do not cover compressed data"
    );
    let mut data = vec![0u8; data_len];
    reader.read_exact(&mut data)?;
    Ok(StackStorage::Compressed(CompressedStack {
        codec,
        row_bytes,
        fingerprint_bytes,
        rows,
        uncompressed_len: expected_len,
        offsets,
        data,
    }))
}

fn read_v5_stack<R: Read>(
    reader: &mut R,
    layout: Layout,
    row_bytes: usize,
    expected_len: usize,
) -> Result<StackStorage> {
    let codec = read_u8(reader)?;
    let uncompressed_len = read_u64(reader)? as usize;
    ensure!(
        uncompressed_len == expected_len,
        "corrupt v5 index: stack length {} != expected {}",
        uncompressed_len,
        expected_len
    );
    if codec == 0 {
        let stored_len = read_u64(reader)? as usize;
        ensure!(
            stored_len == expected_len,
            "corrupt v5 index: plain stack length {} != expected {}",
            stored_len,
            expected_len
        );
        let mut data = vec![0u8; stored_len];
        reader.read_exact(&mut data)?;
        return Ok(StackStorage::Plain(data));
    }

    let block_rows = read_u64(reader)? as usize;
    let stored_row_bytes = read_u64(reader)? as usize;
    let rows = read_u64(reader)? as usize;
    let block_count = read_u64(reader)? as usize;
    ensure!(block_rows > 0, "corrupt v5 index: zero stack block rows");
    ensure!(
        stored_row_bytes == row_bytes,
        "corrupt v5 index: stack row bytes {} != expected {}",
        stored_row_bytes,
        row_bytes
    );
    ensure!(
        rows == layout.array_length,
        "corrupt v5 index: stack rows {} != expected {}",
        rows,
        layout.array_length
    );
    ensure!(
        block_count == rows.div_ceil(block_rows),
        "corrupt v5 index: stack block count does not match row/block layout"
    );
    let mut out = Vec::with_capacity(expected_len);
    for block_id in 0..block_count {
        let block_uncompressed_len = read_u64(reader)? as usize;
        let stored_len = read_u64(reader)? as usize;
        let remaining_rows = rows.saturating_sub(block_id * block_rows);
        let expected_rows = remaining_rows.min(block_rows);
        let expected_block_len = expected_rows
            .checked_mul(row_bytes)
            .context("v5 stack block dimensions overflow")?;
        ensure!(
            block_uncompressed_len == expected_block_len,
            "corrupt v5 index: stack block length {} != expected {}",
            block_uncompressed_len,
            expected_block_len
        );
        let mut data = vec![0u8; stored_len];
        reader.read_exact(&mut data)?;
        let decoded = decompress_v5_stack_block(codec, &data, block_uncompressed_len)?;
        ensure!(
            decoded.len() == block_uncompressed_len,
            "corrupt v5 index: decoded block length {} != expected {}",
            decoded.len(),
            block_uncompressed_len
        );
        out.extend_from_slice(&decoded);
    }
    ensure!(
        out.len() == expected_len,
        "corrupt v5 index: decoded stack length {} != expected {}",
        out.len(),
        expected_len
    );
    Ok(StackStorage::Plain(out))
}

fn write_u8<W: Write>(writer: &mut W, value: u8) -> Result<()> {
    writer.write_all(&[value])?;
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

fn write_layout<W: Write>(writer: &mut W, layout: Layout) -> Result<()> {
    write_u64(writer, layout.segment_length as u64)?;
    write_u64(writer, layout.segment_count as u64)?;
    write_u64(writer, layout.segment_count_length as u64)?;
    write_u64(writer, layout.array_length as u64)?;
    Ok(())
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

fn read_u8<R: Read>(reader: &mut R) -> Result<u8> {
    let mut bytes = [0u8; 1];
    reader.read_exact(&mut bytes)?;
    Ok(bytes[0])
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

fn read_layout<R: Read>(reader: &mut R) -> Result<Layout> {
    let segment_length = read_u64(reader)? as usize;
    let segment_count = read_u64(reader)? as usize;
    let segment_count_length = read_u64(reader)? as usize;
    let array_length = read_u64(reader)? as usize;
    ensure!(segment_length > 0, "invalid zero segment length in index");
    Ok(Layout {
        segment_length,
        segment_length_mask: segment_length - 1,
        segment_count,
        segment_count_length,
        array_length,
    })
}

fn read_string<R: Read>(reader: &mut R) -> Result<String> {
    let len = read_u32(reader)? as usize;
    let mut bytes = vec![0u8; len];
    reader.read_exact(&mut bytes)?;
    Ok(String::from_utf8(bytes)?)
}

fn percent(part: u64, total: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        (part as f64 * 100.0) / total as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use zorindex::{scan_canonical_kmers, ZSTD_FRAME_MAGIC};

    #[test]
    fn cli_defaults_choose_optimized_merged_layout_and_allow_plain_opt_out() {
        let base = [
            "zorindex",
            "build",
            "--fof",
            "samples.fof",
            "--output",
            "samples.zindex",
            "--k",
            "49",
        ];
        let Command::Build(defaults) = Cli::try_parse_from(base).unwrap().command else {
            panic!()
        };
        assert!(defaults.merge.power_of_two_merge);
        assert!(!defaults.plain);
        assert_eq!(
            defaults.merge.planner,
            zorindex::power_merge::MergePlanner::Streaming
        );
        assert_eq!(
            defaults.merge.column_order,
            zorindex::power_merge::ColumnOrder::Compression
        );
        assert_eq!(defaults.slot_scale, 1.0);
        assert_eq!(defaults.cycle_break, CycleBreakHeuristic::MostDeg2);
        assert_eq!(defaults.tie_scan, 8);

        let Command::Build(plain) = Cli::try_parse_from(base.into_iter().chain(["--plain"]))
            .unwrap()
            .command
        else {
            panic!()
        };
        assert!(plain.plain);
        assert!(
            Cli::try_parse_from(base.into_iter().chain(["--plain", "--power-of-two-merge"]))
                .is_err()
        );
        assert!(
            Cli::try_parse_from(base.into_iter().chain(["--plain", "--drop-abandoned"])).is_err()
        );
        assert!(
            Cli::try_parse_from(base.into_iter().chain(["--plain", "--merge-disabled"])).is_err()
        );
        let Command::Build(dropped) =
            Cli::try_parse_from(base.into_iter().chain(["--drop-abandoned"]))
                .unwrap()
                .command
        else {
            panic!()
        };
        assert!(dropped.merge.drop_abandoned);

        let Command::Build(custom) = Cli::try_parse_from(base.into_iter().chain([
            "--merge-planner",
            "all-pairs",
            "--merge-column-order",
            "width",
            "--cycle-break",
            "no-heuristic",
            "--tie-scan",
            "1",
        ]))
        .unwrap()
        .command
        else {
            panic!()
        };
        assert_eq!(
            custom.merge.planner,
            zorindex::power_merge::MergePlanner::AllPairs
        );
        assert_eq!(
            custom.merge.column_order,
            zorindex::power_merge::ColumnOrder::Width
        );
        assert_eq!(custom.cycle_break, CycleBreakHeuristic::NoHeuristic);
        assert_eq!(custom.tie_scan, 1);

        let Command::Append(append) = Cli::try_parse_from([
            "zorindex",
            "append",
            "--index",
            "old.zoridx",
            "--fof",
            "more.fof",
        ])
        .unwrap()
        .command
        else {
            panic!()
        };
        assert_eq!(append.cycle_break, CycleBreakHeuristic::MostDeg2);
        assert_eq!(append.tie_scan, 8);
    }

    #[test]
    fn compact_graph_preserves_lazy_heap_decisions() {
        for n in [0usize, 17, 129, 513] {
            let keys: Vec<_> = (0..n).map(|i| splitmix64(i as u64)).collect();
            for hashes in [2, 3, 4, 8] {
                for heuristic in [
                    CycleBreakHeuristic::NoHeuristic,
                    CycleBreakHeuristic::Lightest,
                    CycleBreakHeuristic::Heaviest,
                    CycleBreakHeuristic::MostDeg2,
                    CycleBreakHeuristic::MinMaxDegree,
                ] {
                    for tie_scan in [1, 4, 8] {
                        let config = PureZorConfig {
                            layout: calculate_layout(n.max(16).next_power_of_two(), hashes, None)
                                .unwrap(),
                            hashes,
                            seed: 69,
                            fingerprint_bits: 16,
                            cycle_break: heuristic,
                            tie_scan,
                        };
                        let old = peel_zor_graph_impl::<usize, BinaryHeap<(Reverse<u32>, usize)>>(
                            &keys, config,
                        )
                        .unwrap();
                        let new =
                            peel_zor_graph_impl::<u32, IndexedDegreeHeap>(&keys, config).unwrap();
                        assert_eq!(
                            old.positions,
                            new.positions
                                .iter()
                                .map(|&i| i as usize)
                                .collect::<Vec<_>>()
                        );
                        assert_eq!(
                            old.order,
                            new.order
                                .iter()
                                .map(|&(a, b)| (a as usize, b as usize))
                                .collect::<Vec<_>>()
                        );
                        assert_eq!(
                            old.dropped,
                            new.dropped.iter().map(|&i| i as usize).collect::<Vec<_>>()
                        );
                        for bits in [8, 16] {
                            assert_eq!(
                                assign_cells(
                                    &keys,
                                    &old.positions,
                                    &old.order,
                                    config.layout,
                                    hashes,
                                    config.seed,
                                    bits
                                )
                                .unwrap(),
                                assign_cells(
                                    &keys,
                                    &new.positions,
                                    &new.order,
                                    config.layout,
                                    hashes,
                                    config.seed,
                                    bits
                                )
                                .unwrap()
                            );
                        }
                    }
                }
            }
        }
        assert!(!compact_graph_eligible(u32::MAX as usize, 16, 4));
        assert!(!compact_graph_eligible(16, u32::MAX as usize + 1, 4));
    }

    fn temp_index_path(name: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "zorindex_{name}_{}_{}.zoridx",
            std::process::id(),
            splitmix64(name.len() as u64)
        ))
    }

    fn temp_test_dir(name: &str) -> PathBuf {
        let path = std::env::temp_dir().join(format!(
            "zorindex_{name}_{}_{}",
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
    fn modimizer_sampling_one_matches_all_kmers() {
        let seq = b"ACGTACGT";
        let mut all = Vec::new();
        let mut sampled = Vec::new();
        scan_canonical_kmers(seq, 4, &mut |k| all.push(k));
        scan_index_features(
            seq,
            4,
            123,
            FeatureConfig {
                mode: IndexMode::Modimizers,
                findere_z: zorindex::findere::DEFAULT_Z,
                minimizer_size: DEFAULT_MINIMIZER_SIZE,
                modimizer_sampling: 1,
            },
            &mut |k| sampled.push(k),
        );
        assert_eq!(sampled, all);
    }

    #[test]
    fn minimizer_mode_emits_canonical_subkmers() {
        let seq = b"ACGTACGT";
        let k = 5;
        let m = 3;
        let mut minimizers = Vec::new();
        scan_index_features(
            seq,
            k,
            123,
            FeatureConfig {
                mode: IndexMode::Minimizers,
                findere_z: zorindex::findere::DEFAULT_Z,
                minimizer_size: m,
                modimizer_sampling: DEFAULT_MODIMIZER_SAMPLING,
            },
            &mut |k| minimizers.push(k),
        );
        assert!(!minimizers.is_empty());

        let mut subkmers = Vec::new();
        scan_canonical_kmers(seq, m, &mut |k| subkmers.push(k));
        for minimizer in minimizers {
            assert!(subkmers.contains(&minimizer));
        }
    }

    #[test]
    fn pure_zor_contains_most_inserted_keys() {
        let keys = (0..10_000u64)
            .map(|x| splitmix64(x) & ((1u64 << 40) - 1))
            .collect::<Vec<_>>();
        let layout = calculate_layout(keys.len(), 4, None).unwrap();
        let build = build_pure_zor(
            &keys,
            PureZorConfig {
                layout,
                hashes: 4,
                seed: 123,
                fingerprint_bits: 8,
                cycle_break: CycleBreakHeuristic::MostDeg2,
                tie_scan: 1,
            },
        )
        .unwrap();
        let hits = keys
            .iter()
            .filter(|&&key| contains_raw(key, &build.cells, layout, 4, 123, 8))
            .count();
        assert!(hits >= keys.len() - build.abandoned as usize);
    }

    #[test]
    fn pure_zor_16_contains_most_inserted_keys() {
        let keys = (0..10_000u64)
            .map(|x| splitmix64(x) & ((1u64 << 40) - 1))
            .collect::<Vec<_>>();
        let layout = calculate_layout(keys.len(), 4, None).unwrap();
        let build = build_pure_zor(
            &keys,
            PureZorConfig {
                layout,
                hashes: 4,
                seed: 123,
                fingerprint_bits: 16,
                cycle_break: CycleBreakHeuristic::MostDeg2,
                tie_scan: 1,
            },
        )
        .unwrap();
        let hits = keys
            .iter()
            .filter(|&&key| contains_raw(key, &build.cells, layout, 4, 123, 16))
            .count();
        assert!(hits >= keys.len() - build.abandoned as usize);
    }

    #[test]
    fn interlaced_query_scores_colors() {
        let c0 = vec![1, 2, 3, 4, 5];
        let c1 = vec![4, 5, 6, 7, 8];
        let layout = calculate_layout(16, 3, None).unwrap();
        let config = PureZorConfig {
            layout,
            hashes: 3,
            seed: 9,
            fingerprint_bits: 8,
            cycle_break: CycleBreakHeuristic::MostDeg2,
            tie_scan: 1,
        };
        let b0 = build_pure_zor(&c0, config).unwrap();
        let b1 = build_pure_zor(&c1, config).unwrap();
        let mut interlaced = vec![0; layout.array_length * 2];
        for slot in 0..layout.array_length {
            interlaced[slot * 2] = b0.cells[slot];
            interlaced[slot * 2 + 1] = b1.cells[slot];
        }
        let colors = vec![
            ColorMeta {
                name: "a".into(),
                path: "a.fa".into(),
                occurrences: 0,
                unique_kmers: c0.len() as u64,
                abandoned: b0.abandoned,
            },
            ColorMeta {
                name: "b".into(),
                path: "b.fa".into(),
                occurrences: 0,
                unique_kmers: c1.len() as u64,
                abandoned: b1.abandoned,
            },
        ];
        let technical_bins = flat_technical_bins(&colors);
        let index = ZorIndex {
            k: 3,
            hashes: 3,
            fingerprint_bits: 8,
            seed: 9,
            feature_config: FeatureConfig::legacy_kmers(),
            layout,
            technical_layout: TechnicalLayoutMode::Flat,
            split_target_features: 0,
            colors,
            technical_bins,
            stack: StackStorage::Plain(interlaced),
            union_graph: None,
        };
        let mut scores = [0, 0];
        index.query_kmer_into(4, &mut scores);
        assert_eq!(scores, [1, 1]);
    }

    #[test]
    fn split_technical_bins_count_once_per_user_color() {
        let layout = calculate_layout(32, 3, None).unwrap();
        let seed = 91;
        let key = (1..1000)
            .find(|&candidate| fingerprint8(candidate, seed) != 0)
            .unwrap();
        let target = fingerprint8(key, seed);
        let mut slots = [0usize; MAX_HASHES];
        fill_positions(key, layout, 3, seed, &mut slots[..3]);
        let mut interlaced = vec![0u8; layout.array_length * 2];
        interlaced[slots[0] * 2] = target;
        interlaced[slots[0] * 2 + 1] = target;
        let colors = vec![ColorMeta {
            name: "split".into(),
            path: "split.fa".into(),
            occurrences: 0,
            unique_kmers: 2,
            abandoned: 0,
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
        let index = ZorIndex {
            k: 5,
            hashes: 3,
            fingerprint_bits: 8,
            seed,
            feature_config: FeatureConfig::legacy_kmers(),
            layout,
            technical_layout: TechnicalLayoutMode::SplitLarge,
            split_target_features: 1,
            colors,
            technical_bins,
            stack: StackStorage::Plain(interlaced),
            union_graph: None,
        };
        let mut scores = [0u64];
        index.query_kmer_into(key, &mut scores);
        assert_eq!(scores, [1]);
    }

    #[test]
    fn split_build_carries_abandoned_keys_to_next_part() {
        let layout = calculate_layout(32, 3, None).unwrap();
        let config = PureZorConfig {
            layout,
            hashes: 3,
            seed: 123,
            fingerprint_bits: 16,
            cycle_break: CycleBreakHeuristic::NoHeuristic,
            tie_scan: 1,
        };
        let mut first = (0..200).map(splitmix64).collect::<Vec<_>>();
        let mut second = (200..240).map(splitmix64).collect::<Vec<_>>();
        first.sort_unstable();
        first.dedup();
        second.sort_unstable();
        second.dedup();

        let first_build = build_pure_zor(&first, config).unwrap();
        assert!(!first_build.abandoned_keys.is_empty());

        let mut expected_second_input = second.clone();
        expected_second_input.extend_from_slice(&first_build.abandoned_keys);
        expected_second_input.sort_unstable();
        expected_second_input.dedup();

        let (parts, final_abandoned) =
            build_zor_part_sequence(vec![first, second], 0, config).unwrap();
        assert_eq!(parts.len(), 2);
        assert_eq!(parts[0].abandoned, first_build.abandoned);
        assert_eq!(parts[1].unique_kmers, expected_second_input.len());
        assert_eq!(final_abandoned, parts[1].abandoned);
    }

    #[test]
    fn append_splits_large_color_in_split_large_index() {
        let dir = temp_test_dir("append_split_zor");
        let a = dir.join("a.fa");
        let b = dir.join("b.fa");
        let build_fof = dir.join("build.fof");
        let append_fof = dir.join("append.fof");
        let index_path = dir.join("base.zoridx");
        let appended_path = dir.join("appended.zoridx");
        write_fasta(&a, "a", "ACGTACGTACGTACGT");
        write_fasta(&b, "b", "AAAACCCCCGGGGGTTTTTACGTACGT");
        write_fof(&build_fof, &[("a", &a)]);
        write_fof(&append_fof, &[("b", &b)]);

        build_command(BuildArgs {
            merge: zorindex::power_merge::MergeArgs::default(),
            plain: true,
            fof: build_fof,
            k: 5,
            output: index_path.clone(),
            hashes: 3,
            fingerprint_bits: 8,
            seed: 99,
            index_mode: IndexMode::Kmers,
            findere_z: zorindex::findere::DEFAULT_Z,
            minimizer_size: DEFAULT_MINIMIZER_SIZE,
            modimizer_sampling: DEFAULT_MODIMIZER_SAMPLING,
            slot_scale: 4.0,
            slot_count: None,
            segment_length: None,
            cycle_break: CycleBreakHeuristic::MostDeg2,
            tie_scan: 1,
            union_graph: false,
            stack_compression: StackCompression::None,
            stack_block_rows: DEFAULT_STACK_BLOCK_ROWS,
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
            cycle_break: CycleBreakHeuristic::MostDeg2,
            tie_scan: 1,
            stack_compression: Some(StackCompression::None),
            stack_block_rows: None,
        })
        .unwrap();

        let loaded = ZorIndex::load(&appended_path).unwrap();
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
    fn chunked_query_matches_scalar_reference() {
        let color_count = 41;
        let hashes = 4;
        let layout = calculate_layout(512, hashes, None).unwrap();
        let mut interlaced = vec![0u8; layout.array_length * color_count];
        let mut state = 1u64;
        for byte in &mut interlaced {
            state = splitmix64(state);
            *byte = state as u8;
        }
        let colors = (0..color_count)
            .map(|idx| ColorMeta {
                name: format!("c{idx}"),
                path: String::new(),
                occurrences: 0,
                unique_kmers: 0,
                abandoned: 0,
            })
            .collect::<Vec<_>>();
        let technical_bins = flat_technical_bins(&colors);
        let index = ZorIndex {
            k: 7,
            hashes,
            fingerprint_bits: 8,
            seed: 99,
            feature_config: FeatureConfig::legacy_kmers(),
            layout,
            technical_layout: TechnicalLayoutMode::Flat,
            split_target_features: 0,
            colors,
            technical_bins,
            stack: StackStorage::Plain(interlaced),
            union_graph: None,
        };

        for key in (0..128).map(|idx| splitmix64(idx * 17)) {
            let mut optimized = vec![0u64; color_count];
            let mut scalar = vec![0u64; color_count];
            index.query_kmer_into(key, &mut optimized);
            query_kmer_into_scalar_reference(&index, key, &mut scalar);
            assert_eq!(optimized, scalar);
        }
    }

    #[test]
    fn chunked_query_16_matches_scalar_reference() {
        let color_count = 37;
        let hashes = 4;
        let layout = calculate_layout(512, hashes, None).unwrap();
        let mut interlaced = vec![0u8; layout.array_length * color_count * 2];
        let mut state = 1u64;
        for slot in 0..layout.array_length {
            for color in 0..color_count {
                state = splitmix64(state);
                let value = state as u16;
                let offset = (slot * color_count + color) * 2;
                interlaced[offset..offset + 2].copy_from_slice(&value.to_le_bytes());
            }
        }
        let colors = (0..color_count)
            .map(|idx| ColorMeta {
                name: format!("c{idx}"),
                path: String::new(),
                occurrences: 0,
                unique_kmers: 0,
                abandoned: 0,
            })
            .collect::<Vec<_>>();
        let technical_bins = flat_technical_bins(&colors);
        let index = ZorIndex {
            k: 7,
            hashes,
            fingerprint_bits: 16,
            seed: 99,
            feature_config: FeatureConfig::legacy_kmers(),
            layout,
            technical_layout: TechnicalLayoutMode::Flat,
            split_target_features: 0,
            colors,
            technical_bins,
            stack: StackStorage::Plain(interlaced),
            union_graph: None,
        };

        for key in (0..128).map(|idx| splitmix64(idx * 17)) {
            let mut optimized = vec![0u64; color_count];
            let mut scalar = vec![0u64; color_count];
            index.query_kmer_into(key, &mut optimized);
            query_kmer_into_scalar_reference(&index, key, &mut scalar);
            assert_eq!(optimized, scalar);
        }
    }

    #[test]
    fn compressed_query_matches_scalar_reference() {
        let color_count = 23;
        let hashes = 4;
        let layout = calculate_layout(512, hashes, None).unwrap();
        let mut interlaced = vec![0u8; layout.array_length * color_count];
        let mut state = 7u64;
        for byte in &mut interlaced {
            state = splitmix64(state);
            *byte = state as u8;
        }
        let colors = (0..color_count)
            .map(|idx| ColorMeta {
                name: format!("c{idx}"),
                path: String::new(),
                occurrences: 0,
                unique_kmers: 0,
                abandoned: 0,
            })
            .collect::<Vec<_>>();
        let technical_bins = flat_technical_bins(&colors);

        for codec in [
            StackCompression::Svb32Sparse,
            StackCompression::Lz4,
            StackCompression::Lz4Lib,
            StackCompression::Lz4Hc,
            StackCompression::Snappy,
            StackCompression::Fastpfor256,
            StackCompression::FastpforPack,
        ] {
            let index = ZorIndex {
                k: 7,
                hashes,
                fingerprint_bits: 8,
                seed: 99,
                feature_config: FeatureConfig::legacy_kmers(),
                layout,
                technical_layout: TechnicalLayoutMode::Flat,
                split_target_features: 0,
                colors: colors.clone(),
                technical_bins: technical_bins.clone(),
                stack: StackStorage::from_plain(
                    interlaced.clone(),
                    codec,
                    color_count,
                    1,
                    layout.array_length,
                )
                .unwrap(),
                union_graph: None,
            };

            for key in (0..128).map(|idx| splitmix64(idx * 19)) {
                let mut optimized = vec![0u64; color_count];
                let mut scalar = vec![0u64; color_count];
                index.query_kmer_into(key, &mut optimized);
                query_kmer_into_scalar_reference(&index, key, &mut scalar);
                assert_eq!(optimized, scalar, "codec={codec:?}");
            }
        }
    }

    #[test]
    fn compressed_query_16_matches_scalar_reference() {
        let color_count = 19;
        let hashes = 4;
        let layout = calculate_layout(512, hashes, None).unwrap();
        let mut interlaced = vec![0u8; layout.array_length * color_count * 2];
        let mut state = 11u64;
        for slot in 0..layout.array_length {
            for color in 0..color_count {
                state = splitmix64(state);
                let value = state as u16;
                let offset = (slot * color_count + color) * 2;
                interlaced[offset..offset + 2].copy_from_slice(&value.to_le_bytes());
            }
        }
        let colors = (0..color_count)
            .map(|idx| ColorMeta {
                name: format!("c{idx}"),
                path: String::new(),
                occurrences: 0,
                unique_kmers: 0,
                abandoned: 0,
            })
            .collect::<Vec<_>>();
        let technical_bins = flat_technical_bins(&colors);

        for codec in [
            StackCompression::Svb32Sparse,
            StackCompression::Lz4,
            StackCompression::Lz4Lib,
            StackCompression::Lz4Hc,
            StackCompression::Snappy,
            StackCompression::Fastpfor256,
            StackCompression::FastpforPack,
        ] {
            let index = ZorIndex {
                k: 7,
                hashes,
                fingerprint_bits: 16,
                seed: 99,
                feature_config: FeatureConfig::legacy_kmers(),
                layout,
                technical_layout: TechnicalLayoutMode::Flat,
                split_target_features: 0,
                colors: colors.clone(),
                technical_bins: technical_bins.clone(),
                stack: StackStorage::from_plain(
                    interlaced.clone(),
                    codec,
                    color_count * 2,
                    2,
                    layout.array_length,
                )
                .unwrap(),
                union_graph: None,
            };

            for key in (0..128).map(|idx| splitmix64(idx * 23)) {
                let mut optimized = vec![0u64; color_count];
                let mut scalar = vec![0u64; color_count];
                index.query_kmer_into(key, &mut optimized);
                query_kmer_into_scalar_reference(&index, key, &mut scalar);
                assert_eq!(optimized, scalar, "codec={codec:?}");
            }
        }
    }

    #[test]
    fn saved_zor_index_is_zstd_serialized_and_loadable() {
        let path = temp_index_path("saved_zor");
        let layout = calculate_layout(32, 3, None).unwrap();
        let colors = vec![
            ColorMeta {
                name: "a".into(),
                path: "a.fa".into(),
                occurrences: 7,
                unique_kmers: 5,
                abandoned: 0,
            },
            ColorMeta {
                name: "b".into(),
                path: "b.fa".into(),
                occurrences: 9,
                unique_kmers: 6,
                abandoned: 1,
            },
        ];
        let stack = StackStorage::Plain(vec![0x5au8; layout.array_length * colors.len()]);
        let technical_bins = flat_technical_bins(&colors);
        let index = ZorIndex {
            k: 5,
            hashes: 3,
            fingerprint_bits: 8,
            seed: 99,
            feature_config: FeatureConfig::legacy_kmers(),
            layout,
            technical_layout: TechnicalLayoutMode::Flat,
            split_target_features: 0,
            colors,
            technical_bins,
            stack,
            union_graph: None,
        };
        index.save(&path).unwrap();

        let mut file = File::open(&path).unwrap();
        let mut magic = [0u8; 4];
        file.read_exact(&mut magic).unwrap();
        assert_eq!(magic, ZSTD_FRAME_MAGIC);

        let loaded = ZorIndex::load(&path).unwrap();
        assert_eq!(loaded.k, index.k);
        assert_eq!(loaded.hashes, index.hashes);
        assert_eq!(loaded.fingerprint_bits, index.fingerprint_bits);
        assert_eq!(loaded.seed, index.seed);
        assert_eq!(loaded.layout, index.layout);
        assert_eq!(loaded.colors.len(), index.colors.len());
        assert_eq!(
            loaded.stack.to_plain_vec().unwrap(),
            index.stack.to_plain_vec().unwrap()
        );
        let _ = std::fs::remove_file(path);
    }

    fn query_kmer_into_scalar_reference(index: &ZorIndex, key: u64, scores: &mut [u64]) {
        let mut slots = [0usize; MAX_HASHES];
        fill_positions(
            key,
            index.layout,
            index.hashes,
            index.seed,
            &mut slots[..index.hashes],
        );
        let color_count = index.colors.len();
        let interlaced = index.stack.to_plain_vec().unwrap();
        match index.fingerprint_bits {
            8 => {
                let target = fingerprint8(key, index.seed);
                for color in 0..color_count {
                    let mut value = 0u8;
                    for &slot in &slots[..index.hashes] {
                        value ^= interlaced[slot * color_count + color];
                    }
                    if value == target {
                        scores[color] += 1;
                    }
                }
            }
            16 => {
                let target = fingerprint16(key, index.seed);
                let row_bytes = color_count * 2;
                for color in 0..color_count {
                    let mut value = 0u16;
                    for &slot in &slots[..index.hashes] {
                        value ^= read_u16_le(&interlaced, slot * row_bytes + color * 2);
                    }
                    if value == target {
                        scores[color] += 1;
                    }
                }
            }
            _ => unreachable!(),
        }
    }
}
