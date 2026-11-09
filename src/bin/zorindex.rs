use std::cmp::Reverse;
use std::collections::{BinaryHeap, VecDeque};
use std::ffi::{OsStr, OsString};
use std::path::PathBuf;

use anyhow::{ensure, Context, Result};
use clap::{Args, CommandFactory, FromArgMatches, Parser, Subcommand, ValueEnum};
use zorindex::power_merge::{
    DegreeQueue, GraphIndex, IndexedDegreeHeap, PeelGraph, TypedPeelGraph,
};
use zorindex::{
    mixsplit, splitmix64, FeatureConfig, IndexMode, MinimizerScheme, RandstrobeConfig, SpacedSeed,
    DEFAULT_MINIMIZER_SIZE, DEFAULT_MODIMIZER_SAMPLING, MAX_HASHES,
};

#[derive(Parser, Debug)]
#[command(name = "zorindex", version = env!("ZINDEX_BUILD_VERSION"),
    about = "Approximate sequence-collection queries using interleaved ZOR filters")]
struct Cli {
    #[command(subcommand)]
    command: Command,
    /// Worker threads. Defaults to available CPUs.
    #[arg(
        short = 't',
        long,
        global = true,
        help_heading = "Important options",
        display_order = 100
    )]
    threads: Option<usize>,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Build an index from a file-of-files.
    Build(BuildArgs),
    /// Query FASTA/FASTQ using the index's stored seed settings.
    Query(QueryArgs),
    /// Print index metadata.
    Info(IndexArg),
    /// Fuse two or more compatible indexes without rereading their source datasets.
    Fusion(FusionArgs),
}

#[derive(Args, Debug)]
struct FusionArgs {
    /// Input index directories, in output dataset order (at least two in total).
    #[arg(short = 'i', long, required = true, num_args = 1..)]
    index: Vec<PathBuf>,
    /// Output directory. Defaults to a new myindex-DATE-TIME directory.
    #[arg(short = 'o', long)]
    output: Option<PathBuf>,
}

#[derive(Args, Debug)]
struct BuildArgs {
    #[command(flatten)]
    merge: zorindex::power_merge::MergeArgs,
    /// One FASTA/FASTQ path per line, optionally preceded by a name and a tab.
    #[arg(short = 'f', long)]
    fof: PathBuf,
    /// Positive k-mer length, or window span in bases for minimizers.
    #[arg(short = 'k', long)]
    k: usize,
    /// Selected minimizer length in bases (1..=31). Required for minimizers.
    #[arg(short = 'm', long, required_if_eq("index_mode", "minimizers"))]
    minimizer_size: Option<usize>,
    /// Output directory. Defaults to a new myindex-DATE-TIME directory.
    #[arg(short = 'o', long)]
    output: Option<PathBuf>,
    /// Hash locations per ZOR equation, from 2 to 32.
    ///
    /// More hash functions can make the index smaller, but slow queries.
    #[arg(short = 'H', long, default_value_t = 4)]
    hashes: usize,
    /// Fingerprint width in bits: 1, 2, 4, 8 or 16.
    #[arg(short = 'b', long, default_value_t = 8)]
    fingerprint_bits: u8,
    /// Shared hash seed, stored for querying.
    ///
    /// Select a reproducible hash realization for extraction, sketches and
    /// filter construction. Different values can change selection, layout and
    /// abandoned equations. Queries use the stored value automatically.
    #[arg(short = 's', long, default_value_t = 69)]
    seed: u64,
    /// Feature extraction mode.
    #[arg(short = 'i', long, value_enum, default_value = "kmers")]
    index_mode: IndexMode,
    /// Index (k-z)-mers and require z+1 consecutive matches in findere mode.
    #[arg(short = 'z', long, default_value_t = zorindex::findere::DEFAULT_Z)]
    findere_z: usize,
    /// Minimizer selector, stored for querying.
    #[arg(short = 'M', long, value_enum, default_value = "open-closed-mod")]
    minimizer_scheme: MinimizerScheme,
    /// Collapse homopolymer runs before extracting seeds.
    #[arg(short = 'c', long)]
    homopolymer_compression: bool,
    /// Symmetric #/_ pattern of width k, or m in minimizer mode.
    #[arg(short = 'p', long)]
    seed_pattern: Option<SpacedSeed>,
    /// Link two k-mers or selected minimizers into canonical randstrobemers.
    #[arg(short = 'r', long)]
    randstrobemers: bool,
    /// Minimum strobe-start distance. Defaults to the strobe length (k or m).
    #[arg(short = 'l', long, requires = "randstrobemers")]
    strobe_min: Option<usize>,
    /// Maximum strobe-start distance. Defaults to three times the strobe length.
    #[arg(short = 'u', long, requires = "randstrobemers")]
    strobe_max: Option<usize>,
    /// Keep roughly one in N canonical k-mers in modimizer mode.
    #[arg(short = 'd', long, default_value_t = DEFAULT_MODIMIZER_SAMPLING)]
    modimizer_sampling: u64,
    /// Positive finite scale for filter capacity before power-of-two rounding.
    ///
    /// Trade filter cells against recovery work by scaling estimated feature
    /// counts. Larger values provide more cells. Capacity changes in power-of-two
    /// steps, and segment geometry still imposes a minimum capacity.
    #[arg(short = 'S', long, default_value_t = 1.0)]
    slot_scale: f64,
    /// Power-of-two hash segment length, 1..=4096 cells, default 64.
    ///
    /// Tune address locality and construction graph geometry. Each hash position
    /// occupies a separate segment. Changing this changes the filter layout,
    /// which is stored for queries and must match when fusing indexes.
    #[arg(short = 'L', long)]
    segment_length: Option<usize>,
    /// Strategy used when ZOR graph peeling reaches a cyclic core.
    /// Select which equation remains active when other equations must be abandoned.
    /// Recovery retains abandoned keys unless explicitly disabled.
    ///
    /// The default takes the first active equation in a minimum-degree cell.
    /// Other strategies score candidates using incident-cell degrees, spending
    /// extra construction work to change which equations enter recovery.
    /// This is a construction experiment, not a query-accuracy setting.
    #[arg(short = 'B', long, value_enum, default_value = "no-heuristic")]
    cycle_break: CycleBreakHeuristic,
    /// Positive number of tied cells examined by cycle-breaking heuristics, default 1.
    ///
    /// With a non-default cycle-breaking strategy, compare this many cells tied
    /// at the minimum active degree. Larger scans do more scoring work.
    /// The no-heuristic strategy does not use this setting.
    #[arg(short = 'T', long, default_value_t = 1)]
    tie_scan: usize,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
enum CycleBreakHeuristic {
    /// Keep the first active key in a min-degree cell.
    #[default]
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

#[derive(Args, Debug)]
struct QueryArgs {
    /// Index directory.
    #[arg(short = 'i', long)]
    index: PathBuf,
    /// Query FASTA/FASTQ path.
    #[arg(short = 'q', long)]
    query: PathBuf,
    /// Minimum matched fraction to report.
    #[arg(short = 'a', long, default_value_t = 0.0)]
    min_ratio: f64,
    /// Score each FASTA sequence or FASTQ read separately.
    ///
    /// Without this flag, counts are combined across the whole query file.
    /// Findere currently supports only combined file counts.
    #[arg(short = 'r', long)]
    records: bool,
    /// Per-record output: named TSV rows (default) or compact dataset_id:count hits.
    #[arg(
        short = 'F',
        long,
        value_enum,
        default_value = "tsv",
        requires = "records"
    )]
    output_format: zorindex::power_merge::RecordFormat,
    /// Positive number of records per parallel reporting batch, default 4096.
    ///
    /// Bound records retained per reporting batch and tune parallel work size.
    /// Smaller batches can reduce record-buffer memory, while larger batches
    /// can amortize scheduling and output overhead. Match semantics are unchanged.
    #[arg(short = 'B', long, default_value_t = 4096)]
    record_batch_size: usize,
    /// Load one subindex at a time or all subindexes together.
    #[arg(short = 'M', long, value_enum, default_value = "subindex")]
    query_memory: zorindex::power_merge::QueryMemory,
}

#[derive(Args, Debug)]
struct IndexArg {
    /// Index directory.
    #[arg(short = 'i', long)]
    index: PathBuf,
}

const ADVANCED_BUILD_OPTIONS: &[&str] = &[
    "max_ram_gb",
    "planner",
    "build_storage",
    "min_subindex_documents",
    "max_subindexes",
    "sketch_size",
    "plan_margin",
    "disabled",
    "seed",
    "slot_scale",
    "segment_length",
    "drop_abandoned",
    "seed_pattern",
    "randstrobemers",
    "strobe_min",
    "strobe_max",
    "homopolymer_compression",
    "modimizer_sampling",
    "findere_z",
    "cycle_break",
    "tie_scan",
];
const ADVANCED_QUERY_OPTIONS: &[&str] = &["query_memory", "record_batch_size", "output_format"];
const ADVANCED_HELP_HINT: &str = "-ah, --advanced-help  Show regular and advanced options.";

fn cli_command(advanced: bool) -> clap::Command {
    let mut command = Cli::command().after_help(ADVANCED_HELP_HINT).arg(
        clap::Arg::new("advanced-help")
            .long("advanced-help")
            .global(true)
            .hide(true)
            .action(clap::ArgAction::HelpLong),
    );
    for subcommand in command.get_subcommands_mut() {
        let advanced_options = match subcommand.get_name() {
            "build" => ADVANCED_BUILD_OPTIONS,
            "query" => ADVANCED_QUERY_OPTIONS,
            _ => &[],
        };
        let mut configured = subcommand.clone().after_help(ADVANCED_HELP_HINT);
        let mandatory: &[&str] = match subcommand.get_name() {
            "build" => &["fof", "k"],
            "query" => &["index", "query"],
            "info" => &["index"],
            "fusion" => &["index"],
            _ => &[],
        };
        for (order, &id) in mandatory.iter().enumerate() {
            configured = configured.mut_arg(id, |arg| {
                arg.help_heading("Mandatory options").display_order(order)
            });
        }
        let important: &[&str] = match subcommand.get_name() {
            "build" => &[
                "output",
                "hashes",
                "fingerprint_bits",
                "index_mode",
                "minimizer_size",
                "minimizer_scheme",
                "max_subindex_size",
            ],
            "query" => &["min_ratio", "records"],
            "fusion" => &["output"],
            _ => &[],
        };
        for (order, &id) in important.iter().enumerate() {
            configured = configured.mut_arg(id, |arg| {
                arg.help_heading("Important options")
                    .display_order(order + 101)
            });
        }
        for (order, &id) in advanced_options.iter().enumerate() {
            configured = configured.mut_arg(id, |arg| {
                arg.hide(!advanced)
                    .help_heading("Advanced options")
                    .display_order(order + 1000)
            });
        }
        *subcommand = configured;
    }
    command
}

fn parse_cli_from<I, T>(args: I) -> std::result::Result<Cli, clap::Error>
where
    I: IntoIterator<Item = T>,
    T: Into<OsString>,
{
    let mut args: Vec<OsString> = args.into_iter().map(Into::into).collect();
    let mut advanced = false;
    for arg in args.iter_mut().skip(1) {
        if arg == OsStr::new("--") {
            break;
        }
        if arg == OsStr::new("-ah") || arg == OsStr::new("--advanced-help") {
            advanced = true;
            *arg = OsString::from("--advanced-help");
        }
    }
    let matches = cli_command(advanced).try_get_matches_from(args)?;
    Cli::from_arg_matches(&matches)
}

#[derive(Clone, Copy, Debug)]
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

fn main() -> Result<()> {
    let cli = parse_cli_from(std::env::args_os()).unwrap_or_else(|error| error.exit());
    let threads = cli.threads.unwrap_or_else(|| {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    });
    ensure!(threads > 0, "--threads must be greater than 0");
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global();
    match cli.command {
        Command::Build(args) => build_command(args),
        Command::Query(args) => query_command(args),
        Command::Fusion(args) => {
            zorindex::power_merge::fuse_many_indexes(&args.index, &output_directory(args.output)?)
        }
        Command::Info(args) => {
            ensure!(
                args.index.is_dir(),
                "--index must name a merged index directory"
            );
            zorindex::power_merge::info_command(&args.index)
        }
    }
}

fn build_command(args: BuildArgs) -> Result<()> {
    validate_build_args(&args)?;
    let output = output_directory(args.output.clone())?;
    let segment_length = args
        .segment_length
        .unwrap_or(zorindex::power_merge::DEFAULT_SEGMENT_LENGTH);
    let config = zorindex::power_merge::BuildConfig {
        backend: zorindex::power_merge::Backend::Zor,
        k: extraction_span_from_args(&args)?,
        hashes: args.hashes,
        seed: args.seed,
        fingerprint_bits: args.fingerprint_bits as usize,
        feature_config: feature_config_from_args(&args),
        slot_scale: args.slot_scale,
        segment_length: Some(segment_length),
        merge: args.merge.clone(),
    };
    zorindex::power_merge::build(&args.fof, &output, config, |keys, slots| {
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
    })
}

fn output_directory(path: Option<PathBuf>) -> Result<PathBuf> {
    if let Some(path) = path {
        return Ok(path);
    }
    let now = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH)?;
    let days = now.as_secs() / 86400;
    // Gregorian civil date from days since 1970-01-01.
    let z = days as i64 + 719468;
    let era = z / 146097;
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let mut year = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = mp + if mp < 10 { 3 } else { -9 };
    year += i64::from(month <= 2);
    let seconds = now.as_secs() % 86400;
    let hour = seconds / 3600;
    let minute = seconds / 60 % 60;
    let second = seconds % 60;
    let stem = format!(
        "myindex-{year:04}{month:02}{day:02}-{hour:02}{minute:02}{second:02}-{:09}",
        now.subsec_nanos()
    );
    let mut path = PathBuf::from(&stem);
    let mut collision = 0u64;
    while path.exists() {
        collision += 1;
        path = PathBuf::from(format!("{stem}-{collision}"));
    }
    eprintln!("Output index: {}", path.display());
    Ok(path)
}

fn query_command(args: QueryArgs) -> Result<()> {
    ensure!(
        args.index.is_dir(),
        "--index must name a merged index directory"
    );
    if args.records {
        zorindex::power_merge::query_records_command_with_memory(
            &args.index,
            &args.query,
            args.min_ratio,
            true,
            args.record_batch_size,
            args.query_memory,
            args.output_format,
        )
    } else {
        zorindex::power_merge::query_command_with_memory(
            &args.index,
            &args.query,
            args.min_ratio,
            true,
            args.query_memory,
        )
    }
}

fn extraction_span_from_args(args: &BuildArgs) -> Result<usize> {
    if args.index_mode == IndexMode::Minimizers {
        args.minimizer_size
            .context("minimizer indexes require --minimizer-size")?;
        ensure!(
            args.k % 2 == 1,
            "canonical minimizers require an odd --k window span"
        );
    } else {
        ensure!(
            args.minimizer_size.is_none(),
            "--minimizer-size requires --index-mode minimizers"
        );
    }
    Ok(args.k)
}

fn feature_config_from_args(args: &BuildArgs) -> FeatureConfig {
    let minimizer_size = args.minimizer_size.unwrap_or(DEFAULT_MINIMIZER_SIZE);
    let strobe_length = if args.index_mode == IndexMode::Minimizers {
        minimizer_size
    } else {
        args.k
    };
    FeatureConfig::new(args.index_mode, minimizer_size, args.modimizer_sampling)
        .with_findere_z(args.findere_z)
        .with_minimizer_scheme(args.minimizer_scheme)
        .with_homopolymer_compression(args.homopolymer_compression)
        .with_seed_pattern(args.seed_pattern)
        .with_randstrobemers(args.randstrobemers.then_some(RandstrobeConfig {
            min_distance: args.strobe_min.unwrap_or(strobe_length),
            max_distance: args.strobe_max.unwrap_or(strobe_length.saturating_mul(3)),
        }))
}

fn validate_build_args(args: &BuildArgs) -> Result<()> {
    ensure!(args.k > 0, "--k must be positive");
    ensure!(
        (2..=32).contains(&args.hashes),
        "--hashes must be in 2..=32"
    );
    validate_feature_config(
        extraction_span_from_args(args)?,
        feature_config_from_args(args),
    )?;
    fingerprint_bytes(args.fingerprint_bits)
        .with_context(|| "--fingerprint-bits must be 1, 2, 4, 8 or 16")?;
    ensure!(
        args.slot_scale.is_finite() && args.slot_scale > 0.0,
        "--slot-scale must be finite and > 0"
    );
    ensure!(args.tie_scan > 0, "--tie-scan must be greater than 0");
    if let Some(segment_length) = args.segment_length {
        ensure!(
            segment_length.is_power_of_two(),
            "--segment-length must be a non-zero power of two"
        );
    }
    Ok(())
}

fn validate_feature_config(k: usize, config: FeatureConfig) -> Result<()> {
    config.validate_anchors(k)?;
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
            ensure!(k > 0, "--k must be positive");
        }
        IndexMode::Minimizers => {
            ensure!(
                k > 0 && k <= 255 && k % 2 == 1,
                "canonical minimizer indexes require an odd window span in 1..=255"
            );
            ensure!(
                (1..=k).contains(&config.minimizer_size),
                "minimizer seed length must be in 1..=window span"
            );
        }
        IndexMode::Modimizers => {
            ensure!(k > 0, "--k must be positive");
        }
    }
    Ok(())
}

fn fingerprint_bytes(fingerprint_bits: u8) -> Result<usize> {
    match fingerprint_bits {
        // Zero denotes sub-byte packed cells, never byte-per-fingerprint.
        1 | 2 | 4 => Ok(0),
        8 => Ok(1),
        16 => Ok(2),
        _ => anyhow::bail!("unsupported fingerprint width {}", fingerprint_bits),
    }
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

    zorindex::power_merge::trace_memory(
        "graph-finished",
        serde_json::json!({
            "distinct":edge_count,"slots":config.layout.array_length,"hashes":config.hashes,
            "index_bytes":std::mem::size_of::<I>(),
            "incidence_capacity_bytes":(edge_slots.capacity()+adjacency.capacity())*std::mem::size_of::<I>(),
            "degree_capacity_bytes":degrees.capacity()*4,
            "offset_cursor_capacity_bytes":(offsets.capacity()+cursor.capacity())*std::mem::size_of::<I>(),
            "active_capacity_bytes":active.capacity(),
            "peel_order_capacity_bytes":peel_order.capacity()*std::mem::size_of::<(I,I)>(),
            "abandoned_capacity_bytes":abandoned_edges.capacity()*std::mem::size_of::<I>(),
            "degree_one_queue_capacity_bytes":queue.capacity()*std::mem::size_of::<I>(),
            "degree_heap_capacity_bytes":multi_heap.allocated_bytes(),
            "peeled":peel_order.len(),"abandoned":abandoned_edges.len()
        }),
    );
    Ok(TypedPeelGraph {
        positions: edge_slots,
        order: peel_order,
        dropped: abandoned_edges,
    })
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

#[inline(always)]
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn minimizer_cli_uses_window_span_and_minimizer_size() {
        for (span, size) in [(49, 31), (49, 30), (31, 31), (255, 1)] {
            let Command::Build(args) = Cli::try_parse_from([
                "zorindex",
                "build",
                "--fof",
                "inputs.fof",
                "--output",
                "index",
                "--index-mode",
                "minimizers",
                "-k",
                &span.to_string(),
                "-m",
                &size.to_string(),
            ])
            .unwrap()
            .command
            else {
                panic!()
            };
            validate_build_args(&args).unwrap();
            assert_eq!(extraction_span_from_args(&args).unwrap(), span);
            assert_eq!(feature_config_from_args(&args).minimizer_size, size);
        }
    }

    #[test]
    fn minimizer_cli_rejects_ambiguous_or_invalid_windows_before_io() {
        let base = [
            "zorindex",
            "build",
            "--fof",
            "missing.fof",
            "--output",
            "index",
            "--k",
            "31",
        ];
        assert!(
            Cli::try_parse_from(base.into_iter().chain(["--index-mode", "minimizers"])).is_err()
        );
        assert!(Cli::try_parse_from(base.into_iter().chain(["--window", "19"])).is_err());
        for extra in [
            vec!["--minimizer-size", "21"],
            vec!["--index-mode", "minimizers", "--minimizer-size", "0"],
            vec!["--index-mode", "minimizers", "--minimizer-size", "32"],
            vec![
                "--index-mode",
                "minimizers",
                "--minimizer-size",
                "18446744073709551615",
            ],
        ] {
            let Command::Build(args) = Cli::try_parse_from(base.into_iter().chain(extra))
                .unwrap()
                .command
            else {
                panic!()
            };
            assert!(validate_build_args(&args).is_err());
        }
        for (span, size) in [(0, 1), (256, 31), (48, 31), (21, 23)] {
            let Command::Build(args) = Cli::try_parse_from([
                "zorindex",
                "build",
                "-f",
                "missing.fof",
                "-o",
                "index",
                "-i",
                "minimizers",
                "-k",
                &span.to_string(),
                "-m",
                &size.to_string(),
            ])
            .unwrap()
            .command
            else {
                panic!()
            };
            assert!(validate_build_args(&args).is_err());
        }
    }

    #[test]
    fn queries_default_to_one_resident_subindex_with_explicit_all_ram_opt_in() {
        let args = [
            "zorindex", "query", "--index", "index", "--query", "reads.fa",
        ];
        let Command::Query(defaults) = Cli::try_parse_from(args).unwrap().command else {
            panic!()
        };
        assert_eq!(
            defaults.query_memory,
            zorindex::power_merge::QueryMemory::Subindex
        );
        let Command::Query(all) =
            Cli::try_parse_from(args.into_iter().chain(["--query-memory", "all"]))
                .unwrap()
                .command
        else {
            panic!()
        };
        assert_eq!(all.query_memory, zorindex::power_merge::QueryMemory::All);
        assert!(
            Cli::try_parse_from(args.into_iter().chain(["--query-memory", "unknown"])).is_err()
        );
    }

    #[test]
    fn every_option_has_one_short_and_one_long_name() {
        let command = cli_command(true);
        command.clone().debug_assert();
        let readme = include_str!("../../README.md");
        for command in std::iter::once(&command).chain(command.get_subcommands()) {
            for arg in command.get_arguments() {
                if arg.is_positional() {
                    continue;
                }
                // The requested -ah spelling is normalized before Clap parses short clusters.
                let short = if arg.get_id() == "advanced-help" {
                    "-ah".to_owned()
                } else {
                    format!(
                        "-{}",
                        arg.get_short().expect("every option needs a short flag")
                    )
                };
                assert!(
                    arg.get_long().is_some(),
                    "{} -{:?}",
                    command.get_name(),
                    arg.get_short()
                );
                assert!(arg
                    .get_all_aliases()
                    .is_none_or(|aliases| aliases.is_empty()));
                assert!(arg
                    .get_all_short_aliases()
                    .is_none_or(|aliases| aliases.is_empty()));
                let long = format!("--{}", arg.get_long().unwrap());
                if ["--cycle-break", "--tie-scan"].contains(&long.as_str()) {
                    assert!(!readme.contains(&long));
                } else {
                    assert!(
                        readme.contains(&format!("`{short}` | `{long}`")),
                        "undocumented {long}"
                    );
                }
            }
        }
    }

    #[test]
    fn advanced_help_is_explicit_and_preserves_query_threshold_short_flag() {
        let command = cli_command(true);
        for (name, ids) in [
            ("build", ADVANCED_BUILD_OPTIONS),
            ("query", ADVANCED_QUERY_OPTIONS),
        ] {
            let subcommand = command.find_subcommand(name).unwrap();
            for &id in ids {
                let arg = subcommand
                    .get_arguments()
                    .find(|arg| arg.get_id() == id)
                    .unwrap();
                assert!(
                    arg.get_help().is_some() || arg.get_long_help().is_some(),
                    "missing advanced help for {id}"
                );
            }
        }
        for flag in ["-h", "--help"] {
            let error = parse_cli_from(["zorindex", "build", flag]).unwrap_err();
            assert_eq!(error.kind(), clap::error::ErrorKind::DisplayHelp);
            let help = error.to_string();
            assert!(help.contains("--fingerprint-bits"));
            assert!(help.contains("--hashes"));
            assert!(help.contains("--minimizer-size"));
            assert!(help.contains("--max-subindex-size"));
            assert!(help.contains("Mandatory options"));
            assert!(help.contains("Important options"));
            assert!(
                help.find("Mandatory options").unwrap() < help.find("Important options").unwrap()
            );
            assert!(!help.contains("--build-memory-gib"));
            assert!(!help.contains("--cycle-break"));
            assert!(!help.contains("--tie-scan"));
            assert!(!help.contains("--slot-scale"));
            assert!(help.contains("-ah, --advanced-help"));
        }
        for flag in ["-ah", "--advanced-help"] {
            let error = parse_cli_from(["zorindex", "build", flag]).unwrap_err();
            assert_eq!(error.kind(), clap::error::ErrorKind::DisplayHelp);
            let help = error.to_string();
            assert!(
                help.find("Mandatory options").unwrap() < help.find("Important options").unwrap()
            );
            assert!(
                help.find("Important options").unwrap() < help.find("Advanced options").unwrap()
            );
            for expected in [
                "--cycle-break",
                "--tie-scan",
                "--slot-scale",
                "--fingerprint-bits",
                "no-heuristic",
                "lightest",
                "heaviest",
                "most-deg2",
                "min-max-degree",
            ] {
                assert!(help.contains(expected), "missing {expected}: {help}");
            }
            for subcommand in ["query", "info", "fusion"] {
                assert_eq!(
                    parse_cli_from(["zorindex", subcommand, flag])
                        .unwrap_err()
                        .kind(),
                    clap::error::ErrorKind::DisplayHelp
                );
            }
        }
        let Command::Query(args) = parse_cli_from([
            "zorindex", "query", "-i", "index", "-q", "reads.fa", "-a", "0.8",
        ])
        .unwrap()
        .command
        else {
            panic!()
        };
        assert_eq!(args.min_ratio, 0.8);
        let error = parse_cli_from(["zorindex", "build", "--", "-ah"]).unwrap_err();
        assert_ne!(error.kind(), clap::error::ErrorKind::DisplayHelp);
    }

    #[test]
    fn removed_commands_flags_and_long_aliases_are_rejected() {
        for removed in [
            "--build-memory-gib",
            "--max-document-memory-gib",
            "--skip-oversized-documents",
            "--disk-backed",
            "--low-disk",
        ] {
            let error =
                parse_cli_from(["zorindex", "build", "-f", "inputs.fof", "-k", "31", removed])
                    .unwrap_err();
            assert_eq!(
                error.kind(),
                clap::error::ErrorKind::UnknownArgument,
                "{removed}"
            );
        }

        for flag in ["-j", "--other"] {
            let error = parse_cli_from([
                "zorindex", "fusion", "-i", "first", "second", "-o", "output", flag, "third",
            ])
            .unwrap_err();
            assert_eq!(error.kind(), clap::error::ErrorKind::UnknownArgument);
        }
        let query = ["zorindex", "query", "-i", "index", "-q", "reads.fa"];
        for flag in ["--ignore-abandoned", "-R"] {
            let error = parse_cli_from(query.into_iter().chain([flag])).unwrap_err();
            assert_eq!(error.kind(), clap::error::ErrorKind::UnknownArgument);
        }
        for name in ["append", "repack", "repack-subindexes", "merge"] {
            assert!(Cli::try_parse_from(["zorindex", name]).is_err());
        }
        let base = [
            "zorindex",
            "build",
            "-f",
            "inputs.fof",
            "-o",
            "index",
            "-k",
            "31",
        ];
        for name in [
            "--plain",
            "--union-graph",
            "--stack-compression",
            "--technical-layout",
            "--target-technical-features",
            "--max-split-parts",
            "--stack-block-rows",
            "--power-of-two-merge",
            "--build-jobs",
            "--homocompress",
            "--spaced-seed",
            "--randstrobes",
            "--z",
            "--strobe-window-min",
            "--strobe-window-max",
            "--slot-count",
            "--merge-column-order",
            "--merge-max-colors",
            "--merge-min-saving",
            "-C",
            "-O",
            "-A",
            "-E",
        ] {
            assert!(
                Cli::try_parse_from(base.into_iter().chain([name])).is_err(),
                "{name}"
            );
        }
        for planner in ["streaming-unpruned", "all-pairs"] {
            assert!(
                Cli::try_parse_from(base.into_iter().chain(["--merge-planner", planner])).is_err()
            );
        }
        for planner in ["full", "streaming", "naive"] {
            assert!(Cli::try_parse_from(base.into_iter().chain(["-P", planner])).is_ok());
        }
    }

    #[test]
    fn compact_graph_preserves_peeling_decisions() {
        for n in [0usize, 17, 129, 513] {
            let keys: Vec<_> = (0..n).map(|i| splitmix64(i as u64)).collect();
            for hashes in [2usize, 3, 4, 8] {
                for cycle_break in [
                    CycleBreakHeuristic::NoHeuristic,
                    CycleBreakHeuristic::Lightest,
                    CycleBreakHeuristic::Heaviest,
                    CycleBreakHeuristic::MostDeg2,
                    CycleBreakHeuristic::MinMaxDegree,
                ] {
                    for tie_scan in [1, 4, 8] {
                        let slots = n.max(64 * hashes.next_power_of_two()).next_power_of_two();
                        let config = PureZorConfig {
                            layout: calculate_layout_with_segment_length(slots, hashes, 64)
                                .unwrap(),
                            hashes,
                            seed: 69,
                            fingerprint_bits: 16,
                            cycle_break,
                            tie_scan,
                        };
                        let wide = peel_zor_graph_impl::<usize, BinaryHeap<(Reverse<u32>, usize)>>(
                            &keys, config,
                        )
                        .unwrap();
                        let narrow =
                            peel_zor_graph_impl::<u32, IndexedDegreeHeap>(&keys, config).unwrap();
                        assert_eq!(
                            wide.positions,
                            narrow
                                .positions
                                .iter()
                                .map(|&n| n as usize)
                                .collect::<Vec<_>>()
                        );
                        assert_eq!(
                            wide.order,
                            narrow
                                .order
                                .iter()
                                .map(|&(a, b)| (a as usize, b as usize))
                                .collect::<Vec<_>>()
                        );
                        assert_eq!(
                            wide.dropped,
                            narrow
                                .dropped
                                .iter()
                                .map(|&n| n as usize)
                                .collect::<Vec<_>>()
                        );
                    }
                }
            }
        }
        assert!(!compact_graph_eligible(u32::MAX as usize, 16, 4));
        assert!(!compact_graph_eligible(16, u32::MAX as usize + 1, 4));
    }
}
