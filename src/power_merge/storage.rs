use super::{
    abandoned::{self, OriginMap},
    document_guard::{Guard, SkippedDocument},
    planner,
    sketch::Sketch,
    Backend, BuildConfig, ColumnOrder, MergePlanner, PeelGraph, QueryMemory,
};
use crate::{
    create_index_writer, finish_index_writer, open_index_reader, parse_features_from_file,
    read_fof, splitmix64, ColorInput, FeatureConfig, IndexMode, MinimizerScheme, MAX_HASHES,
};
use anyhow::{anyhow, ensure, Context, Result};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::{
    cmp::Reverse,
    collections::{BTreeMap, BinaryHeap, HashMap},
    fs::{self, File},
    hash::{DefaultHasher, Hash, Hasher},
    io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write},
    ops::Range,
    panic::AssertUnwindSafe,
    path::{Path, PathBuf},
    sync::{
        atomic::{AtomicBool, Ordering},
        mpsc,
    },
    time::{Instant, SystemTime, UNIX_EPOCH},
};

const MAGIC: &[u8; 8] = b"ZORPOW1\0";
const PLANE_MAGIC: &[u8; 8] = b"ZORPOW3\0";
const ABANDONED_MAGIC: &[u8; 8] = b"ZORABN1\0";
const ABANDONED_DELTA_MAGIC: &[u8; 8] = b"ZORABN2\0";
const ABANDONED_FILE: &str = "abandoned.zormap.zst";
const FINGERPRINT_SEED: u64 = 0x8E2D_4F6A_9B13_C5D7;

#[path = "bank_io.rs"]
mod bank_io;

#[path = "query_stream.rs"]
mod query_stream;
#[path = "record_output.rs"]
mod record_output;
#[path = "record_query.rs"]
mod record_query;
pub use record_output::RecordFormat;

#[path = "fusion.rs"]
mod fusion;
pub use fusion::{fuse_indexes, fuse_many_indexes};

#[path = "abandoned_io.rs"]
mod abandoned_io;
use abandoned_io::AbandonedReader;

#[derive(Clone, Debug, Serialize, Deserialize)]
struct Color {
    name: String,
    path: String,
    occurrences: u64,
    unique_kmers: usize,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
struct GroupMeta {
    members: Vec<usize>,
    width: usize,
    offset: usize,
    union_cardinality: usize,
    fallback_keys: usize,
    #[serde(default)]
    abandoned_keys: Option<usize>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
struct BankMeta {
    #[serde(default, skip_serializing_if = "is_zero")]
    component: usize,
    file: String,
    slots: usize,
    row_bits: usize,
    groups: Vec<GroupMeta>,
    #[serde(default)]
    fingerprint_first: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    routing: Option<Routing>,
}

fn is_zero(value: &usize) -> bool {
    *value == 0
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
struct Routing {
    source_slots: usize,
    partition: usize,
    partitions: usize,
}

fn partition_id(key: u64, seed: u64, partitions: usize) -> usize {
    (splitmix64(key ^ seed ^ 0x34B6_1D8F_29A7_C053) as usize) & (partitions - 1)
}

impl BankMeta {
    fn source_group(&self) -> (usize, usize) {
        (self.component, self.source_slots())
    }

    fn expected_filename(&self) -> String {
        let name = Self::filename(self.slots, self.routing);
        if self.component == 0 {
            name
        } else {
            format!("component-{}-{name}", self.component)
        }
    }

    fn source_slots(&self) -> usize {
        self.routing.map_or(self.slots, |r| r.source_slots)
    }
    fn filename(slots: usize, routing: Option<Routing>) -> String {
        match routing {
            Some(r) => format!(
                "bank-{}-part-{}-of-{}.zorbank.zst",
                r.source_slots.ilog2(),
                r.partition,
                r.partitions
            ),
            None => format!("bank-{}.zorbank.zst", slots.ilog2()),
        }
    }
}
#[derive(Debug, Serialize, Deserialize)]
struct Manifest {
    version: u32,
    #[serde(default = "recovery_enabled")]
    abandoned_recovery: bool,
    #[serde(default)]
    abandoned_file: Option<String>,
    #[serde(default)]
    abandoned_hashes: usize,
    #[serde(default)]
    abandoned_origins: usize,
    backend: Backend,
    k: usize,
    hashes: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    segment_length: Option<usize>,
    seed: u64,
    fingerprint_bits: usize,
    index_mode: u8,
    minimizer_size: usize,
    #[serde(default)]
    minimizer_scheme: MinimizerScheme,
    #[serde(default)]
    homopolymer_compression: bool,
    #[serde(default)]
    seed_pattern: Option<String>,
    #[serde(default)]
    randstrobemers: Option<crate::RandstrobeConfig>,
    #[serde(default)]
    randstrobe_scheme: Option<String>,
    modimizer_sampling: u64,
    findere_z: usize,
    sketch_size: usize,
    sketch_method: String,
    #[serde(default = "legacy_merge_planner")]
    merge_planner: MergePlanner,
    #[serde(default = "legacy_column_order")]
    column_order: ColumnOrder,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    column_order_method: Option<String>,
    #[serde(default)]
    planning_pair_evaluations: u64,
    #[serde(default)]
    planning_possible_pair_evaluations: u64,
    #[serde(default = "legacy_merge_reverts_enabled")]
    merge_reverts_enabled: bool,
    #[serde(default = "legacy_cell_storage")]
    construction_cell_storage: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    construction_abandoned_storage: Option<String>,
    slot_scale: f64,
    #[serde(default = "legacy_min_subindex_documents")]
    min_subindex_documents: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    max_subindexes: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    max_subindex_size: Option<u64>,
    #[serde(default)]
    promoted_groups: usize,
    singleton_logical_bits: u64,
    predicted_logical_bits: u64,
    planned_groups: usize,
    rejected_groups: usize,
    /// Groups whose exact union needed more slots than planned and were built at the next power of two.
    #[serde(default)]
    overflow_groups: usize,
    exact_group_verifications: u64,
    sketch_seconds: f64,
    plan_seconds: f64,
    build_seconds: f64,
    #[serde(default)]
    input_document_count: usize,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    source_document_ids: Vec<usize>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    skipped_documents: Vec<SkippedDocument>,
    #[serde(default)]
    memory_guard_split_groups: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    merge_memory_budget_bytes: Option<u64>,
    colors: Vec<Color>,
    banks: Vec<BankMeta>,
}
fn recovery_enabled() -> bool {
    true
}
fn legacy_min_subindex_documents() -> usize {
    1
}
fn legacy_merge_reverts_enabled() -> bool {
    true
}
fn legacy_merge_planner() -> MergePlanner {
    MergePlanner::AllPairs
}
fn legacy_column_order() -> ColumnOrder {
    ColumnOrder::Width
}
fn legacy_cell_storage() -> String {
    "in-memory-u64".into()
}
impl Manifest {
    fn collection_status(&self) -> &'static str {
        if self.skipped_documents.is_empty()
            && self.input_document_count.max(self.colors.len()) == self.colors.len()
        {
            "complete"
        } else {
            "incomplete"
        }
    }
    fn features(&self) -> Result<FeatureConfig> {
        ensure!(self.randstrobemers.is_none() || (self.version >= 6 && self.randstrobe_scheme.as_deref() == Some(crate::RANDSTROBE_SCHEME)),
            "unsupported legacy/custom randstrobemer extraction; rebuild this index with the current downstream scheme");
        let config = FeatureConfig::new(
            IndexMode::from_u8(self.index_mode)?,
            self.minimizer_size,
            self.modimizer_sampling,
        )
        .with_findere_z(self.findere_z)
        .with_minimizer_scheme(self.minimizer_scheme)
        .with_homopolymer_compression(self.homopolymer_compression)
        .with_seed_pattern(self.seed_pattern.as_deref().map(str::parse).transpose()?)
        .with_randstrobemers(self.randstrobemers);
        config.validate_anchors(self.k)?;
        Ok(config)
    }
}
struct BuiltGroup {
    ordering_key: Option<std::sync::Arc<Vec<u8>>>,
    meta: GroupMeta,
    slots: usize,
    cells: CompletedCells,
    // Construction-only 32-bit hashes and origin masks. Converted to the single
    // global dataset multimap before publication, never consulted during queries.
    fallback: Vec<(u32, u64)>,
    fallback_path: Option<PathBuf>,
    routing: Option<Routing>,
}

impl BuiltGroup {
    fn discard(self) -> Result<()> {
        self.cells.remove()?;
        if let Some(path) = self.fallback_path {
            fs::remove_file(path)?;
        }
        Ok(())
    }
}

enum CompletedCells {
    Disk(PathBuf),
    Packed {
        data: Vec<u8>,
        width: usize,
        slots: usize,
    },
    // The bytes of `Packed` staged in a file (opt-in disk staging).
    PackedDisk {
        path: PathBuf,
        width: usize,
        slots: usize,
    },
}

// Slots per write when staging packed cells: a multiple of 8, so every chunk is
// byte-aligned and the chunks concatenate into exactly the `Packed` byte image.
const STAGED_CHUNK_SLOTS: usize = 1 << 20;

impl CompletedCells {
    fn packed<T: Copy + Into<u64>>(cells: &[T], width: usize) -> Result<Self> {
        let bits = cells
            .len()
            .checked_mul(width)
            .context("filter size overflow")?;
        let mut data = vec![0; bits.div_ceil(8)];
        for (slot, &value) in cells.iter().enumerate() {
            put_bits(&mut data, slot * width, width, value.into());
        }
        Ok(Self::Packed {
            data,
            width,
            slots: cells.len(),
        })
    }

    // Pack and write a chunk of slots at a time, so staging never needs a
    // second full-size buffer next to the u64 cells.
    fn packed_disk<T: Copy + Into<u64>>(cells: &[T], width: usize, path: PathBuf) -> Result<Self> {
        let mut writer = BufWriter::new(File::create(&path)?);
        for chunk in cells.chunks(STAGED_CHUNK_SLOTS) {
            if let Self::Packed { data, .. } = Self::packed(chunk, width)? {
                writer.write_all(&data)?;
            }
        }
        writer.flush()?;
        Ok(Self::PackedDisk {
            path,
            width,
            slots: cells.len(),
        })
    }

    fn validate(&self, slots: usize, width: usize) -> Result<()> {
        match self {
            Self::Disk(path) => ensure!(
                fs::metadata(path)?.len()
                    == (slots as u64)
                        .checked_mul(8)
                        .context("cell size overflow")?,
                "invalid completed filter cell cache"
            ),
            Self::Packed {
                data,
                width: stored_width,
                slots: stored_slots,
            } => ensure!(
                *stored_width == width
                    && *stored_slots == slots
                    && data.len()
                        == slots
                            .checked_mul(width)
                            .context("filter size overflow")?
                            .div_ceil(8),
                "invalid packed completed filter"
            ),
            Self::PackedDisk {
                path,
                width: stored_width,
                slots: stored_slots,
            } => ensure!(
                *stored_width == width
                    && *stored_slots == slots
                    && fs::metadata(path)?.len()
                        == slots
                            .checked_mul(width)
                            .context("filter size overflow")?
                            .div_ceil(8) as u64,
                "invalid packed completed filter cache"
            ),
        }
        Ok(())
    }

    fn open(&self) -> std::io::Result<Option<File>> {
        match self {
            Self::Disk(path) | Self::PackedDisk { path, .. } => File::open(path).map(Some),
            Self::Packed { .. } => Ok(None),
        }
    }

    // Recreate the original little-endian u64 stream so both storage modes
    // use identical ordering samples and exactly the same final bank packer.
    fn read_block(&self, start: usize, output: &mut [u8], file: Option<&mut File>) -> Result<()> {
        ensure!(output.len() % 8 == 0, "unaligned cell read");
        match self {
            Self::Disk(path) => {
                let mut reopened;
                let reader = if let Some(file) = file {
                    file
                } else {
                    reopened = File::open(path)?;
                    &mut reopened
                };
                reader.seek(SeekFrom::Start(
                    (start as u64)
                        .checked_mul(8)
                        .context("cell offset overflow")?,
                ))?;
                reader.read_exact(output)?;
            }
            Self::Packed { data, width, slots } => {
                ensure!(
                    start
                        .checked_add(output.len() / 8)
                        .is_some_and(|end| end <= *slots),
                    "cell read out of bounds"
                );
                for (i, bytes) in output.chunks_exact_mut(8).enumerate() {
                    bytes.copy_from_slice(
                        &get_bits(data, (start + i) * width, *width).to_le_bytes(),
                    );
                }
            }
            Self::PackedDisk { path, width, slots } => {
                let rows = output.len() / 8;
                ensure!(
                    start.checked_add(rows).is_some_and(|end| end <= *slots),
                    "cell read out of bounds"
                );
                // Read only the bytes covering the requested slots, then unpack
                // them exactly as the in-RAM `Packed` arm does.
                let first_bit = start * width;
                let first_byte = first_bit / 8;
                let mut window = vec![0; ((start + rows) * width).div_ceil(8) - first_byte];
                let mut reopened;
                let reader = if let Some(file) = file {
                    file
                } else {
                    reopened = File::open(path)?;
                    &mut reopened
                };
                reader.seek(SeekFrom::Start(first_byte as u64))?;
                reader.read_exact(&mut window)?;
                for (i, bytes) in output.chunks_exact_mut(8).enumerate() {
                    bytes.copy_from_slice(
                        &get_bits(&window, first_bit % 8 + i * width, *width).to_le_bytes(),
                    );
                }
            }
        }
        Ok(())
    }

    fn remove(self) -> Result<()> {
        if let Self::Disk(path) | Self::PackedDisk { path, .. } = self {
            fs::remove_file(path)?;
        }
        Ok(())
    }
}
struct Bank {
    meta: BankMeta,
    data: Vec<u8>,
    #[cfg(test)]
    fallback: Vec<Vec<(u32, u64)>>,
}
struct Index {
    manifest: Manifest,
    banks: Vec<Bank>,
    abandoned: Option<OriginMap>,
    row_scratch_bytes: usize,
}

// A build writes only to a newly created staging directory. Publishing it at
// the end avoids leaving an apparently usable index after a failed build.
struct Staging(PathBuf);
impl Drop for Staging {
    fn drop(&mut self) {
        if self.0.exists() {
            let _ = fs::remove_dir_all(&self.0);
        }
    }
}
fn cache_path(stage: &Path, id: usize) -> PathBuf {
    stage.join("features").join(format!("{id}.u64"))
}
fn cache_write(path: &Path, keys: &[u64]) -> Result<()> {
    cache_write_values(path, keys)
}
fn cache_write_values<T: Copy + Into<u64>>(path: &Path, keys: &[T]) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    for &key in keys {
        writer.write_all(&key.into().to_le_bytes())?;
    }
    writer.flush()?;
    Ok(())
}
#[cfg(test)]
fn cache_read(path: &Path) -> Result<Vec<u64>> {
    let bytes = fs::read(path)?;
    ensure!(bytes.len() % 8 == 0, "invalid feature cache");
    Ok(bytes
        .chunks_exact(8)
        .map(|b| u64::from_le_bytes(b.try_into().unwrap()))
        .collect())
}

pub fn build(
    fof: &Path,
    output: &Path,
    mut config: BuildConfig,
    graph: impl Fn(&[u64], usize) -> Result<PeelGraph> + Sync,
) -> Result<()> {
    if config.merge.planner == MergePlanner::Naive {
        // Naive means ordered placement only: neither a global size-class
        // promotion pass nor sketch-based physical reordering is performed.
        config.merge.min_subindex_documents = 1;
        config.merge.max_subindexes = None;
    }
    let column_order = if config.merge.planner == MergePlanner::Naive {
        ColumnOrder::Width
    } else {
        ColumnOrder::Compression
    };
    config.validate()?;
    if config.merge.max_ram_gb.is_none() && super::memory::rss_bytes().is_some() {
        config.merge.max_ram_gb = Some(super::memory::budget(None) as f64 / 1_000_000_000.0);
    }
    let rss_guard = super::memory::RssGuard::start(config.merge.max_ram_gb)?;
    config.feature_config.validate_anchors(config.k)?;
    ensure!(
        !output.exists(),
        "output {} already exists; choose a new directory",
        output.display()
    );
    let start = Instant::now();
    let original_inputs = read_fof(fof)?;
    let input_document_count = original_inputs.len();
    let extraction_budget = super::memory::workspace_budget(&config.merge);
    let guard = Guard::new(&config, extraction_budget);
    config.merge.planning_memory_budget_bytes = Some(guard.budget);
    let original_costs = original_inputs
        .par_iter()
        // The file-size heuristic only schedules extraction concurrency. It
        // cannot reject a document: the bounded collector measures its keys.
        .map(|input| {
            super::memory::input_estimate(&input.path)
                .map(|hint| hint.min((guard.budget / 4).max(1)))
        })
        .collect::<Result<Vec<_>>>()?;
    for (id, (input, &cost)) in original_inputs.iter().zip(&original_costs).enumerate() {
        guard
            .check("input-size", cost)
            .with_context(|| format!("dataset {id} {} exceeds --max-ram-gb", input.name))?;
    }
    let inputs = original_inputs;
    let input_costs = original_costs;
    ensure!(!inputs.is_empty(), "{} contains no inputs", fof.display());
    ensure!(
        inputs.len() <= u32::MAX as usize,
        "dataset IDs must fit in 32 bits"
    );
    let parent = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let stamp = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
    let stage = Staging(parent.join(format!(".zor-pow2-building-{}-{stamp}", std::process::id())));
    fs::create_dir(&stage.0)?;
    if config.merge.cache_features() {
        fs::create_dir(stage.0.join("features"))?;
    }
    if config.merge.cache_features() || !config.merge.in_memory() {
        fs::create_dir(stage.0.join("cells"))?;
        if !config.merge.drop_abandoned {
            fs::create_dir(stage.0.join("abandoned"))?;
        }
    }
    eprintln!("build workspace admission: extraction budget={} bytes, workers={}; document guard={}, document budget={} bytes; estimates, not an RSS cap",
        extraction_budget, rayon::current_num_threads(), guard.enabled, guard.budget);
    eprintln!(
        "sketching {} datasets with {} SimdSketch buckets (8-bit fingerprints, 32-bit union minima)",
        inputs.len(),
        config.merge.sketch_size
    );
    let mut colors = Vec::new();
    let mut initial = Vec::new();
    let mut streaming = planner::StreamingPlanner::default();
    let mut singleton_logical_bits = 0u128;
    let mut plan_seconds = 0.0;
    let mut sketch_seconds = 0.0;
    let mut input_signatures = Vec::with_capacity(inputs.len());
    let mut included_inputs = Vec::new();
    let mut included_costs = Vec::new();
    let mut source_document_ids = Vec::new();
    // Every dataset is parsed as a task on the rayon pool while this loop, itself on a pool thread, inserts the
    // results into the planner in FOF order. Idle threads take the oldest queued task first, so they parse while
    // this thread plans, and the planner's own parallel scans split onto the same pool when it is free.
    rayon::scope(|_| {
        pull_jobs_budgeted_with_policy(
            &inputs,
            PLANNER_LOOKAHEAD,
            &input_costs,
            extraction_budget,
            !guard.enabled,
            config.merge.max_ram_gb.is_some(),
            super::memory::total_bytes(config.merge.max_ram_gb),
            |id, input| {
                let (keys, occurrences) = guard.collect(input, &config).with_context(|| {
                    format!("reading dataset {id} {} within --max-ram-gb", input.name)
                })?;
                let sketch_start = Instant::now();
                let sketch = if config.merge.planner == MergePlanner::Naive {
                    Sketch::default()
                } else {
                    Sketch::new(&keys, config.merge.sketch_size, config.seed)
                };
                let sketch_sec = sketch_start.elapsed().as_secs_f64();
                if config.merge.cache_features() {
                    cache_write(&cache_path(&stage.0, id), &keys)?;
                }
                let signature = feature_signature(&keys, occurrences);
                let color = Color {
                    name: input.name.clone(),
                    path: input.path.to_string_lossy().into_owned(),
                    occurrences,
                    unique_kmers: keys.len(),
                };
                Ok((color, sketch, sketch_sec, signature))
            },
            |original_id, result| {
                let (color, sketch, seconds, signature) = result;
                let id = colors.len();
                let mut group = planner::singleton(color.unique_kmers, sketch, id, &config)?;
                group.input_peak_bytes = super::memory::feature_resident_bytes(color.occurrences);
                included_inputs.push(inputs[original_id].clone());
                included_costs.push(super::memory::feature_resident_bytes(color.occurrences));
                source_document_ids.push(original_id);
                input_signatures.push(signature);
                colors.push(color);
                singleton_logical_bits += planner::cost(group.slots, group.width);
                if config.merge.planner.is_streaming() {
                    let plan_start = Instant::now();
                    streaming.insert(group, &config)?;
                    plan_seconds += plan_start.elapsed().as_secs_f64();
                } else {
                    initial.push(group);
                }
                sketch_seconds += seconds;
                Ok(())
            },
        )
    })?;
    let inputs = included_inputs;
    let input_costs = included_costs;
    eprintln!(
        "sketched {} datasets in {:.3}s; the planner was busy for {:.3}s of it",
        colors.len(),
        start.elapsed().as_secs_f64(),
        plan_seconds
    );
    let singleton_logical_bits = u64::try_from(singleton_logical_bits)?;
    let plan_start = Instant::now();
    let plan = if config.merge.planner.is_streaming() {
        streaming.finish()
    } else {
        planner::plan(initial, &config)?
    };
    plan_seconds += plan_start.elapsed().as_secs_f64();
    let planning_pair_evaluations = plan.pair_evaluations;
    let planning_possible_pair_evaluations = plan.possible_pair_evaluations;
    let (planned, memory_guard_split_groups) =
        guard_groups(plan.groups, &colors, &input_costs, &config, guard)?;
    let predicted_logical_bits = u64::try_from(
        planned
            .iter()
            .map(|g| planner::cost(g.slots, g.width))
            .sum::<u128>(),
    )?;
    let planned_groups = planned.len();
    let projected_banks = planned
        .iter()
        .map(|g| g.slots)
        .collect::<std::collections::BTreeSet<_>>()
        .len();
    eprintln!(
        "merge plan ({}; no reverts): {} datasets -> {} groups, predicted bits {} -> {}, planning {:.3}s, pair evaluations {}",
        config.merge.planner.label(),
        colors.len(),
        planned_groups,
        singleton_logical_bits,
        predicted_logical_bits,
        plan_seconds,
        planning_pair_evaluations
    );
    if !config.merge.cache_features() && config.merge.in_memory() {
        eprintln!("two-pass construction: projected {} power-of-two banks, {} logical bits; rereading inputs, no temporary key/cell files; completed filters stay bit-packed in RAM", projected_banks, predicted_logical_bits);
    } else if !config.merge.cache_features() {
        eprintln!("two-pass construction: projected {} power-of-two banks, {} logical bits; rereading inputs, completed filters and abandoned keys staged on disk", projected_banks, predicted_logical_bits);
    }
    let mut built = Vec::new();
    let rejected_groups = 0;
    let mut exact_group_verifications = 0;
    let group_costs = planned
        .iter()
        .map(|group| group_workspace(group, &colors, &input_costs, &config, 0))
        .collect::<Result<Vec<_>>>()?;
    let construction_budget = super::memory::workspace_budget(&config.merge);
    eprintln!("build workspace admission: construction budget={} bytes, workers={}; completed filters and abandoned keys retained in {}",
        construction_budget, rayon::current_num_threads(), if config.merge.in_memory() && !config.merge.cache_features() { "RAM (additional to workspace)" } else { "temporary disk files" });
    pull_jobs_budgeted_with_policy(
        &planned,
        usize::MAX,
        &group_costs,
        construction_budget,
        !guard.enabled,
        config.merge.max_ram_gb.is_some(),
        super::memory::total_bytes(config.merge.max_ram_gb),
        |_, group| {
            if config.merge.cache_features() {
                build_group(group, &stage.0, &config, &graph)
            } else {
                let (keys, values) =
                    union_input_features(group, &inputs, &colors, &input_signatures, &config)?;
                let stage = (!config.merge.in_memory()).then_some(stage.0.as_path());
                construct_group(group, keys, values, stage, &config, &graph, None)
            }
        },
        |_, (group, verified)| {
            built.push(group);
            exact_group_verifications += verified;
            Ok(())
        },
    )?;
    let overflow_groups = planned
        .iter()
        .zip(&built)
        .filter(|(plan, group)| group.slots > plan.slots)
        .count();
    let promoted_groups = promote_small_subindexes(
        &mut built,
        &planned,
        &inputs,
        &colors,
        &input_signatures,
        &input_costs,
        &stage.0,
        &config,
        &graph,
    )?;
    built = partition_large_subindexes(
        built,
        &planned,
        &inputs,
        &colors,
        &input_signatures,
        &input_costs,
        &stage.0,
        &config,
        &graph,
    )?;
    drop(planned);
    let partitioned = built.iter().any(|g| g.routing.is_some());
    if partitioned {
        exact_group_verifications = built
            .iter()
            .map(|g| {
                (g.meta.union_cardinality
                    - if config.merge.drop_abandoned {
                        g.meta.abandoned_keys.unwrap_or(0)
                    } else {
                        0
                    }) as u64
            })
            .sum();
    }
    let (abandoned_hashes, abandoned_origins) = if config.merge.drop_abandoned {
        (0, 0)
    } else {
        save_abandoned_groups(&mut built, &stage.0)?
    };
    let banks = pack_banks_streaming(built, column_order, config.fingerprint_bits, &stage.0)?;
    if config.merge.cache_features() {
        fs::remove_dir(stage.0.join("features"))?;
    }
    if config.merge.cache_features() || !config.merge.in_memory() {
        fs::remove_dir(stage.0.join("cells"))?;
        if !config.merge.drop_abandoned {
            fs::remove_dir(stage.0.join("abandoned"))?;
        }
    }
    if stage.0.join("partitions").exists() {
        fs::remove_dir_all(stage.0.join("partitions"))?;
    }
    let mut manifest = Manifest {
        version: if partitioned { 9 } else if config.segment_length.is_some() { 8 } else { 7 },
        abandoned_recovery: !config.merge.drop_abandoned,
        abandoned_file: (!config.merge.drop_abandoned).then(|| ABANDONED_FILE.into()),
        abandoned_hashes,
        abandoned_origins,
        backend: config.backend,
        k: config.k,
        hashes: config.hashes,
        segment_length: config.segment_length,
        seed: config.seed,
        fingerprint_bits: config.fingerprint_bits,
        index_mode: config.feature_config.mode.as_u8(),
        minimizer_size: config.feature_config.minimizer_size,
        minimizer_scheme: config.feature_config.minimizer_scheme,
        homopolymer_compression: config.feature_config.homopolymer_compression,
        seed_pattern: config.feature_config.seed_pattern.map(|p| p.to_string()),
        randstrobemers: config.feature_config.randstrobemers,
        randstrobe_scheme: config
            .feature_config
            .randstrobemers
            .map(|_| crate::RANDSTROBE_SCHEME.into()),
        modimizer_sampling: config.feature_config.modimizer_sampling,
        findere_z: config.feature_config.findere_z,
        sketch_size: config.merge.sketch_size,
        sketch_method: "simd-sketch-0.3.0 feature-ID/union extensions; 8-bit fingerprints, 32-bit union minima, AVX2 collection/comparison/union; collision-corrected Jaccard".into(),
        merge_planner: config.merge.planner,
        column_order,
        column_order_method: Some(if column_order == ColumnOrder::Compression {
            "union-sketch-fingerprints-lexicographic-v1"
        } else {
            "width-then-dataset-id-v1"
        }.into()),
        planning_pair_evaluations,
        planning_possible_pair_evaluations,
        merge_reverts_enabled: false,
        construction_cell_storage: if config.merge.cache_features() {
            "temporary-disk-u64-tiled-packing"
        } else if config.merge.in_memory() {
            "two-pass-memory-packed-no-scratch"
        } else {
            "two-pass-temporary-disk-packed-tiled-packing"
        }
        .into(),
        construction_abandoned_storage: Some(if config.merge.drop_abandoned {
            "omitted"
        } else if config.merge.in_memory() && !config.merge.cache_features() {
            "memory"
        } else {
            "temporary-disk-sorted-runs"
        }.into()),
        slot_scale: config.slot_scale,
        min_subindex_documents: config.merge.min_subindex_documents,
        max_subindexes: config.merge.max_subindexes,
        max_subindex_size: Some(config.merge.max_subindex_size),
        promoted_groups,
        singleton_logical_bits,
        predicted_logical_bits,
        planned_groups,
        rejected_groups,
        overflow_groups,
        exact_group_verifications,
        sketch_seconds,
        plan_seconds,
        build_seconds: 0.0,
        input_document_count,
        source_document_ids,
        skipped_documents: Vec::new(),
        memory_guard_split_groups,
        merge_memory_budget_bytes: config.merge.planning_memory_budget_bytes,
        colors,
        banks,
    };
    manifest.build_seconds = start.elapsed().as_secs_f64();
    let mut writer = BufWriter::new(File::create(stage.0.join("manifest.json"))?);
    serde_json::to_writer_pretty(&mut writer, &manifest)?;
    writer.flush()?;
    ensure!(
        !output.exists(),
        "output appeared during build: {}",
        output.display()
    );
    rss_guard.check()?;
    fs::rename(&stage.0, output)?;
    eprintln!(
        "collection\tstatus={}\trequested={}\tindexed={}\tskipped={}",
        manifest.collection_status(),
        input_document_count,
        manifest.colors.len(),
        manifest.skipped_documents.len()
    );
    eprintln!("built {} in {:.3}s: banks={}, groups={}, stack_bytes={}, fallback_bytes={}, rejected_groups={}, overflow_groups={}",
        output.display(), start.elapsed().as_secs_f64(), manifest.banks.len(), manifest.banks.iter().map(|b| b.groups.len()).sum::<usize>(),
        manifest.banks.iter().map(|b| b.slots * b.row_bits.div_ceil(8)).sum::<usize>(), abandoned_hashes * 8 + abandoned_origins * 4, rejected_groups, overflow_groups);
    Ok(())
}

/// Maximum datasets parsing may run ahead of the planner. Waiting sketches use
/// 32-bit minima and 8-bit fingerprints per configured cell.
const PLANNER_LOOKAHEAD: usize = 1024;

/// Runs `work` on every item as a task on the rayon pool and hands the results to `consume` in item order on the
/// calling thread, keeping at most `in_flight` items started but not yet consumed. The pool runs the tasks in item
/// order, as many at a time as it has threads, so one slow item never idles the others, and parallel work inside
/// `consume` (the planner) shares those same threads. Nothing more is started after the first error.
#[cfg(test)]
fn pull_jobs<T: Sync, R: Send>(
    items: &[T],
    in_flight: usize,
    work: impl Fn(usize, &T) -> Result<R> + Sync,
    consume: impl FnMut(usize, R) -> Result<()>,
) -> Result<()> {
    pull_jobs_budgeted(
        items,
        in_flight,
        &vec![0; items.len()],
        u64::MAX,
        work,
        consume,
    )
}

#[cfg(test)]
fn pull_jobs_budgeted<T: Sync, R: Send>(
    items: &[T],
    in_flight: usize,
    estimates: &[u64],
    budget: u64,
    work: impl Fn(usize, &T) -> Result<R> + Sync,
    mut consume: impl FnMut(usize, R) -> Result<()>,
) -> Result<()> {
    pull_jobs_budgeted_with_policy(
        items,
        in_flight,
        estimates,
        budget,
        true,
        false,
        None,
        work,
        &mut consume,
    )
}

fn pull_jobs_budgeted_with_policy<T: Sync, R: Send>(
    items: &[T],
    in_flight: usize,
    estimates: &[u64],
    budget: u64,
    allow_oversized: bool,
    explicit_budget: bool,
    total_ram_budget: Option<u64>,
    work: impl Fn(usize, &T) -> Result<R> + Sync,
    mut consume: impl FnMut(usize, R) -> Result<()>,
) -> Result<()> {
    ensure!(
        in_flight > 0 && budget > 0 && estimates.len() == items.len(),
        "invalid admission limits"
    );
    let baseline = super::memory::rss_bytes().unwrap_or(0);
    let failed = AtomicBool::new(false);
    let (tx, rx) = mpsc::channel();
    rayon::in_place_scope(|scope| {
        let mut started = 0;
        let mut reserved = 0u64;
        let mut oversized_warnings = 0;
        // Results arrive in completion order; hold them until their predecessors are consumed.
        let mut pending = BTreeMap::new();
        let outcome = (0..items.len()).try_for_each(|next| {
            while started < items.len().min(next.saturating_add(in_flight)) {
                let i = started;
                let estimate = estimates[i];
                // Admission events also account for retained recovery/filter
                // buffers and other processes consuming RAM since stage start.
                // No timer or monitor thread is involved.
                let limit = if let Some(total) = total_ram_budget {
                    super::memory::total_admission_limit(total, baseline)
                } else if estimate == 0 {
                    budget
                } else {
                    super::memory::admission_limit(budget,explicit_budget)
                };
                // Reserve before spawning, never block a Rayon worker on a
                // semaphore: the same pool may be needed to finish active jobs.
                if started != next && (reserved > limit || estimate > limit.saturating_sub(reserved)) {
                    break;
                }
                ensure!(allow_oversized || estimate <= limit,
                    "build workspace job {i} estimate {estimate} exceeds available admission budget {limit}; memory guard refuses an uncapped job");
                if estimate > limit && oversized_warnings < 8 {
                    eprintln!("build workspace: job {i} estimate={estimate} exceeds admission budget={limit}; running it alone (single-job memory is not capped)");
                    oversized_warnings += 1;
                }
                reserved = reserved.saturating_add(estimate);
                let (tx, failed, work) = (tx.clone(), &failed, &work);
                scope.spawn(move |_| {
                    if failed.load(Ordering::Relaxed) {
                        return;
                    }
                    // A panic must still produce a result, or the loop below would wait for it forever.
                    let result = std::panic::catch_unwind(AssertUnwindSafe(|| work(i, &items[i])))
                        .unwrap_or_else(|_| Err(anyhow!("task {i} panicked")));
                    let _ = tx.send((i, result));
                });
                started += 1;
            }
            while !pending.contains_key(&next) {
                // On a pool thread, run queued tasks rather than wait for the other threads (there may be none).
                let (i, result) = match rx.try_recv() {
                    Ok(message) => message,
                    Err(_) if rayon::yield_now() == Some(rayon::Yield::Executed) => continue,
                    Err(_) => rx.recv()?,
                };
                pending.insert(i, result?);
            }
            consume(next, pending.remove(&next).unwrap())?;
            reserved = reserved.saturating_sub(estimates[next]);
            Ok(())
        });
        if outcome.is_err() {
            failed.store(true, Ordering::Relaxed);
        }
        outcome
    })
}

fn group_workspace(
    group: &planner::Group,
    colors: &[Color],
    input_costs: &[u64],
    config: &BuildConfig,
    slot_floor: usize,
) -> Result<u64> {
    let slot_floor = if config.merge.planner == MergePlanner::Naive {
        slot_floor.max(group.slots)
    } else {
        slot_floor
    };
    if config.merge.planning_memory_budget_bytes.is_some() {
        return Ok(super::memory::group_estimate(
            group.feature_count_upper_bound,
            config
                .slots(group.feature_count_upper_bound as f64)?
                .max(slot_floor),
            config.hashes,
            group.width,
            group.input_peak_bytes,
        ));
    }
    let sum = group.members.iter().fold(0usize, |sum, &id| {
        sum.saturating_add(colors[id].unique_kmers)
    });
    let largest = group
        .members
        .iter()
        .map(|&id| colors[id].unique_kmers)
        .max()
        .unwrap_or(0);
    // The exact union is not available until extraction. Allow a full extra
    // capacity doubling over the sketch plan, bounded by exact input counts.
    let features = sum.min(group.slots.saturating_mul(2)).max(largest);
    let slots = config.slots(features as f64)?.max(slot_floor);
    let input_peak = group
        .members
        .iter()
        .map(|&id| {
            // Include actual first-pass occurrence count (before dedup), not just
            // compressed size, for highly repetitive genomes.
            input_costs[id].max(super::memory::feature_resident_bytes(
                colors[id].occurrences,
            ))
        })
        .max()
        .unwrap_or(0);
    Ok(super::memory::group_estimate(
        features,
        slots,
        config.hashes,
        group.width,
        input_peak,
    ))
}

fn guard_groups(
    groups: Vec<planner::Group>,
    colors: &[Color],
    input_costs: &[u64],
    config: &BuildConfig,
    guard: Guard,
) -> Result<(Vec<planner::Group>, usize)> {
    let mut result = Vec::new();
    let mut split = 0;
    for group in groups {
        // Sum of exact document cardinalities is a safe union bound. Avoid
        // discarding individually admissible documents because a merge is large.
        let sum = group
            .members
            .iter()
            .fold(0usize, |n, &id| n.saturating_add(colors[id].unique_kmers));
        let input_peak = group
            .members
            .iter()
            .map(|&id| {
                input_costs[id].max(super::memory::feature_resident_bytes(
                    colors[id].occurrences,
                ))
            })
            .fold(0u64, u64::saturating_add);
        let cost = super::memory::group_estimate(
            sum,
            config
                .slots(sum as f64)?
                .max(if config.merge.planner == MergePlanner::Naive {
                    group.slots
                } else {
                    0
                }),
            config.hashes,
            group.width,
            input_peak,
        );
        if cost > guard.budget && group.members.len() > 1 {
            split += 1;
            for &id in &group.members {
                result.push(planner::Group {
                    members: vec![id],
                    cardinality: colors[id].unique_kmers as f64,
                    slots: config.slots(colors[id].unique_kmers as f64)?,
                    width: config.width(1),
                    sketch: group.sketch.clone(),
                    ordering_key: group.ordering_key.clone(),
                    feature_count_upper_bound: colors[id].unique_kmers,
                    input_peak_bytes: input_costs[id].max(super::memory::feature_resident_bytes(
                        colors[id].occurrences,
                    )),
                });
            }
        } else {
            result.push(group);
        }
    }
    if split != 0 {
        eprintln!("memory guard: split {split} oversized merged groups into admissible singleton filters; no included documents were dropped");
    }
    Ok((result, split))
}

fn union_features(group: &planner::Group, stage: &Path) -> Result<(Vec<u64>, Vec<u64>)> {
    let mut readers = group
        .members
        .iter()
        .map(|&id| {
            let file = File::open(cache_path(stage, id))?;
            ensure!(file.metadata()?.len() % 8 == 0, "invalid feature cache");
            Ok(BufReader::with_capacity(64 * 1024, file))
        })
        .collect::<Result<Vec<_>>>()?;
    let mut heap = BinaryHeap::new();
    fn next(reader: &mut impl Read) -> Result<Option<u64>> {
        let mut bytes = [0; 8];
        if reader.read(&mut bytes[..1])? == 0 {
            return Ok(None);
        }
        reader.read_exact(&mut bytes[1..])?;
        Ok(Some(u64::from_le_bytes(bytes)))
    }
    for (id, reader) in readers.iter_mut().enumerate() {
        if let Some(key) = next(reader)? {
            heap.push(Reverse((key, id)));
        }
    }
    let mut keys = Vec::new();
    let mut memberships = Vec::new();
    while let Some(Reverse((key, id))) = heap.pop() {
        let mask = 1u64 << id;
        if keys.last() == Some(&key) {
            *memberships.last_mut().unwrap() |= mask;
        } else {
            keys.push(key);
            memberships.push(mask);
        }
        if let Some(next) = next(&mut readers[id])? {
            ensure!(next > key, "unsorted/duplicate feature cache");
            heap.push(Reverse((next, id)));
        }
    }
    Ok((keys, memberships))
}

fn feature_signature(keys: &[u64], occurrences: u64) -> u64 {
    let mut hasher = DefaultHasher::new();
    occurrences.hash(&mut hasher);
    keys.hash(&mut hasher);
    hasher.finish()
}

// Merge backwards in-place, retaining only the current union and one dataset.
// Dataset origin bits and sorted key order match the disk-cache heap union.
#[cfg(test)]
fn extend_union(
    keys: &mut Vec<u64>,
    origins: &mut Vec<u64>,
    incoming: &[u64],
    mask: u64,
) -> Result<()> {
    let mut i = keys.len();
    let mut j = incoming.len();
    let total = i.checked_add(j).context("union size overflow")?;
    keys.resize(total, 0);
    origins.resize(total, 0);
    let mut write = total;
    while i != 0 || j != 0 {
        let (key, value) = if j == 0 || (i != 0 && keys[i - 1] > incoming[j - 1]) {
            i -= 1;
            (keys[i], origins[i])
        } else if i == 0 || incoming[j - 1] > keys[i - 1] {
            j -= 1;
            (incoming[j], mask)
        } else {
            i -= 1;
            j -= 1;
            (keys[i], origins[i] | mask)
        };
        write -= 1;
        keys[write] = key;
        origins[write] = value;
    }
    let length = total - write;
    keys.copy_within(write..total, 0);
    origins.copy_within(write..total, 0);
    keys.truncate(length);
    origins.truncate(length);
    Ok(())
}

fn union_input_features(
    group: &planner::Group,
    inputs: &[ColorInput],
    colors: &[Color],
    signatures: &[u64],
    config: &BuildConfig,
) -> Result<(Vec<u64>, Vec<u64>)> {
    let mut union = super::union::Union::default();
    for (bit, &id) in group.members.iter().enumerate() {
        let mut incoming = Vec::new();
        let mut occurrences = 0;
        parse_features_from_file(
            &inputs[id].path,
            config.k,
            config.seed,
            config.feature_config,
            |key| {
                occurrences += 1;
                incoming.push(key);
            },
        )
        .with_context(|| format!("second-pass extraction of {}", inputs[id].path.display()))?;
        super::sort::sort_unique(&mut incoming);
        ensure!(
            incoming.len() == colors[id].unique_kmers
                && occurrences == colors[id].occurrences
                && feature_signature(&incoming, occurrences) == signatures[id],
            "indexed features changed between passes for {}",
            inputs[id].path.display()
        );
        union.push(incoming, 1u64 << bit);
    }
    Ok(union.finish())
}

fn build_group(
    group: &planner::Group,
    stage: &Path,
    config: &BuildConfig,
    graph: &(impl Fn(&[u64], usize) -> Result<PeelGraph> + Sync),
) -> Result<(BuiltGroup, u64)> {
    let (keys, values) = union_features(group, stage)?;
    let result = construct_group(group, keys, values, Some(stage), config, graph, None)?;
    // The immutable plan partitions datasets: no other group needs these keys.
    for &id in &group.members {
        fs::remove_file(cache_path(stage, id))?;
    }
    Ok(result)
}

fn construct_group(
    group: &planner::Group,
    keys: Vec<u64>,
    values: Vec<u64>,
    cache_stage: Option<&Path>,
    config: &BuildConfig,
    graph: &(impl Fn(&[u64], usize) -> Result<PeelGraph> + Sync),
    slot_floor: Option<usize>,
) -> Result<(BuiltGroup, u64)> {
    let slots = config
        .slots(keys.len() as f64)?
        .max(slot_floor.unwrap_or(0))
        .max(if config.merge.planner == MergePlanner::Naive {
            group.slots
        } else {
            0
        });
    construct_group_at_slots(group, keys, values, cache_stage, config, graph, slots)
}

fn construct_group_at_slots(
    group: &planner::Group,
    keys: Vec<u64>,
    mut values: Vec<u64>,
    cache_stage: Option<&Path>,
    config: &BuildConfig,
    graph: &(impl Fn(&[u64], usize) -> Result<PeelGraph> + Sync),
    slots: usize,
) -> Result<(BuiltGroup, u64)> {
    for (&key, value) in keys.iter().zip(&mut values) {
        *value = fingerprint(key, config.seed, config.fingerprint_bits)
            | if group.members.len() == 1 {
                0
            } else {
                *value << config.fingerprint_bits
            };
    }
    super::memory::trace_memory(
        "graph-start",
        serde_json::json!({
        "members":group.members,"distinct":keys.len(),"hashes":config.hashes,"slots":slots,
        "key_capacity_bytes":keys.capacity().saturating_mul(8),
        "payload_capacity_bytes":values.capacity().saturating_mul(8),
        "key_resident_estimate_bytes":group.input_peak_bytes.max((keys.len() as u64).saturating_mul(8)),
        "prediction_model":super::memory::prediction_kind(keys.len(),slots,config.hashes),
        "predicted_peak_rss_bytes":super::memory::group_prediction(keys.len(),slots,config.hashes,
            group.width,group.input_peak_bytes),
        "admission_workspace_bytes":super::memory::group_estimate(keys.len(),slots,config.hashes,
            group.width,group.input_peak_bytes),
        "conservative_capacity_bound_bytes":super::memory::capacity_bound(keys.len(),slots,config.hashes,
            group.width,group.input_peak_bytes),
        "estimated_workspace_bytes":super::memory::group_estimate(keys.len(),slots,config.hashes,
            group.width,group.input_peak_bytes)}),
    );
    let mut peeled = graph(&keys, slots)?;
    let cells = peeled.assign_verified_narrow(&keys, &values, slots, config.hashes, group.width)?;
    // Recover all dropped origins, including when the unconstrained cells happen
    // to match a fingerprint. Recovery adds hits; it never replaces filter hits.
    let dropped_len = peeled.dropped_len();
    let mut fallback_path = None;
    let mut fallback = Vec::new();
    if !config.merge.drop_abandoned {
        let records = peeled.dropped().map(|edge| {
            (
                abandoned::hash(keys[edge], config.seed),
                if group.members.len() == 1 {
                    1
                } else {
                    values[edge] >> config.fingerprint_bits
                },
            )
        });
        if let Some(stage) = cache_stage {
            if dropped_len != 0 {
                let path = stage
                    .join("abandoned")
                    .join(format!("{}-{slots}.records", group.members[0]));
                let id = u32::try_from(group.members[0]).context("group ID exceeds 32 bits")?;
                super::staged_abandoned::write_sorted(
                    &path,
                    records.map(|(hash, mask)| (hash, id, mask)),
                )?;
                fallback_path = Some(path);
            }
        } else {
            fallback.extend(records);
        }
    }
    let meta = GroupMeta {
        members: group.members.clone(),
        width: group.width,
        offset: 0,
        union_cardinality: keys.len(),
        fallback_keys: if config.merge.drop_abandoned {
            0
        } else {
            dropped_len
        },
        abandoned_keys: Some(dropped_len),
    };
    let verified = (keys.len()
        - if config.merge.drop_abandoned {
            dropped_len
        } else {
            0
        }) as u64;
    // Release graph/union workspaces before retaining completed cells: bit-packed
    // in a staged file (default), in RAM (--build-storage memory), or in the
    // --build-storage cached u64 cell cache.
    drop(peeled);
    drop(values);
    drop(keys);
    fn retain<T: Copy + Into<u64>>(
        cells: &[T],
        cache_stage: Option<&Path>,
        group: &planner::Group,
        config: &BuildConfig,
    ) -> Result<CompletedCells> {
        Ok(match cache_stage {
            Some(stage) if config.merge.cache_features() => {
                let path =
                    stage
                        .join("cells")
                        .join(format!("{}-{}.u64", group.members[0], cells.len()));
                cache_write_values(&path, cells)?;
                CompletedCells::Disk(path)
            }
            Some(stage) => CompletedCells::packed_disk(
                cells,
                group.width,
                stage
                    .join("cells")
                    .join(format!("{}-{}.bits", group.members[0], cells.len())),
            )?,
            None => CompletedCells::packed(cells, group.width)?,
        })
    }
    use super::graph::AssignedCells;
    let cells = match cells {
        AssignedCells::U8(data) if group.width == 8 && cache_stage.is_none() => {
            CompletedCells::Packed {
                data,
                width: 8,
                slots,
            }
        }
        AssignedCells::U8(data) => retain(&data, cache_stage, group, config)?,
        AssignedCells::U16(data) => retain(&data, cache_stage, group, config)?,
        AssignedCells::U32(data) => retain(&data, cache_stage, group, config)?,
        AssignedCells::U64(data) => retain(&data, cache_stage, group, config)?,
    };
    Ok((
        BuiltGroup {
            ordering_key: group.ordering_key.clone(),
            meta,
            slots,
            cells,
            fallback,
            fallback_path,
            routing: None,
        },
        verified,
    ))
}

/// Exact sizes are known only after construction. Rebuild only promoted
/// groups, rereading their own inputs; no key/cell scratch files are introduced.
/// The sketch merge plan, member assignments and fingerprint width stay fixed.
fn promote_small_subindexes(
    built: &mut [BuiltGroup],
    planned: &[planner::Group],
    inputs: &[ColorInput],
    colors: &[Color],
    signatures: &[u64],
    input_costs: &[u64],
    stage: &Path,
    config: &BuildConfig,
    graph: &(impl Fn(&[u64], usize) -> Result<PeelGraph> + Sync),
) -> Result<usize> {
    let mut sizes = BTreeMap::<usize, (usize, usize)>::new();
    for group in built.iter() {
        let bucket = sizes.entry(group.slots).or_default();
        bucket.0 += group.meta.members.len();
        bucket.1 += group.meta.width;
    }
    let buckets: Vec<_> = sizes
        .iter()
        .map(|(&slots, &(docs, width))| (slots, docs, width))
        .collect();
    let targets = super::layout::promotion_targets(
        &buckets,
        config.merge.min_subindex_documents,
        config.merge.max_subindexes,
    )?;
    let mapping: BTreeMap<_, _> = buckets
        .iter()
        .zip(&targets)
        .map(|(b, &target)| (b.0, target))
        .collect();
    let promoted = built
        .iter()
        .filter(|group| mapping[&group.slots] > group.slots)
        .count();
    if promoted == 0 {
        return Ok(0);
    }
    let plans: BTreeMap<_, _> = planned
        .iter()
        .map(|group| (group.members[0], group))
        .collect();
    eprintln!("subindex layout: {} original size buckets -> {} subindexes, minimum={}, maximum={:?}, rebuilding {} promoted groups without temporary files", buckets.len(), targets.iter().collect::<std::collections::BTreeSet<_>>().len(), config.merge.min_subindex_documents, config.merge.max_subindexes, promoted);
    let jobs: Vec<_> = built
        .iter()
        .enumerate()
        .filter_map(|(id, old)| {
            let slots = mapping[&old.slots];
            (slots > old.slots).then(|| (id, slots, plans[&old.meta.members[0]]))
        })
        .collect();
    let costs = jobs
        .iter()
        .map(|(_, slots, group)| group_workspace(group, colors, input_costs, config, *slots))
        .collect::<Result<Vec<_>>>()?;
    if let Some(budget) = config.merge.planning_memory_budget_bytes {
        ensure!(costs.iter().all(|&cost| cost <= budget),
            "requested subindex promotion exceeds the memory guard budget {budget}; reduce layout promotion requirements");
    }
    pull_jobs_budgeted_with_policy(
        &jobs,
        usize::MAX,
        &costs,
        super::memory::workspace_budget(&config.merge),
        false,
        config.merge.max_ram_gb.is_some(),
        super::memory::total_bytes(config.merge.max_ram_gb),
        |_, &(id, slots, group)| {
            let (keys, values) = union_input_features(group, inputs, colors, signatures, config)?;
            let cache =
                (config.merge.cache_features() || !config.merge.in_memory()).then_some(stage);
            let (new, _) = construct_group(group, keys, values, cache, config, graph, Some(slots))?;
            Ok((id, new))
        },
        |_, (id, new)| std::mem::replace(&mut built[id], new).discard(),
    )?;
    Ok(promoted)
}

/// Final geometry includes exact union sizes and promotions. Only oversized
/// classes are reconstructed, using the original merge plan and original keys.
/// Each key belongs to one partition; its H locations use that partition's
/// local power-of-two geometry. Unpartitioned classes are byte-for-byte unchanged.
#[allow(clippy::too_many_arguments)]
fn partition_large_subindexes(
    built: Vec<BuiltGroup>,
    planned: &[planner::Group],
    inputs: &[ColorInput],
    colors: &[Color],
    signatures: &[u64],
    input_costs: &[u64],
    stage: &Path,
    config: &BuildConfig,
    graph: &(impl Fn(&[u64], usize) -> Result<PeelGraph> + Sync),
) -> Result<Vec<BuiltGroup>> {
    let mut widths = BTreeMap::<usize, usize>::new();
    for group in &built {
        *widths.entry(group.slots).or_default() += group.meta.width;
    }
    let mut splits = BTreeMap::new();
    let minimum = config
        .segment_length
        .unwrap_or(1)
        .checked_mul(config.hashes.next_power_of_two())
        .context("partition segment geometry overflow")?;
    for (&slots, &bits) in &widths {
        let bytes = slots
            .checked_mul(bits.div_ceil(8))
            .context("subindex size overflow")?;
        if bytes as u64 <= config.merge.max_subindex_size {
            continue;
        }
        let mut local = slots;
        while (local as u128) * (bits.div_ceil(8) as u128) > config.merge.max_subindex_size as u128
        {
            local /= 2;
        }
        ensure!(local >= minimum,
            "--max-subindex-size={} bytes is too small for {} bytes per slice and minimum {} slots; raise the cap",
            config.merge.max_subindex_size, bits.div_ceil(8), minimum);
        splits.insert(slots, (local, slots / local));
    }
    if splits.is_empty() {
        return Ok(built);
    }
    eprintln!(
        "hash partitioning {} oversized size classes; maximum uncompressed subindex={} bytes",
        splits.len(),
        config.merge.max_subindex_size
    );
    let plans: BTreeMap<_, _> = planned.iter().map(|g| (g.members[0], g)).collect();
    let mut result = Vec::new();
    let mut jobs = Vec::new();
    for old in built {
        if let Some(&(local, partitions)) = splits.get(&old.slots) {
            let source = old.slots;
            let plan = plans[&old.meta.members[0]];
            // Release the superseded cells and abandoned records first.
            old.discard()?;
            jobs.push((plan, source, local, partitions));
        } else {
            result.push(old);
        }
    }
    let costs = jobs
        .iter()
        .map(|(g, source, _, _)| group_workspace(g, colors, input_costs, config, *source))
        .collect::<Result<Vec<_>>>()?;
    pull_jobs_budgeted_with_policy(
        &jobs,
        usize::MAX,
        &costs,
        super::memory::workspace_budget(&config.merge),
        false,
        config.merge.max_ram_gb.is_some(),
        super::memory::total_bytes(config.merge.max_ram_gb),
        |_, &(group, source_slots, local, partitions)| {
            let (keys, masks) = union_input_features(group, inputs, colors, signatures, config)?;
            let mut counts = vec![0usize; partitions];
            for &key in &keys {
                counts[partition_id(key, config.seed, partitions)] += 1;
            }
            let mut keys_by_part: Vec<_> = counts.iter().map(|&n| Vec::with_capacity(n)).collect();
            let mut masks_by_part: Vec<_> = counts.iter().map(|&n| Vec::with_capacity(n)).collect();
            for (key, mask) in keys.into_iter().zip(masks) {
                let id = partition_id(key, config.seed, partitions);
                keys_by_part[id].push(key);
                masks_by_part[id].push(mask);
            }
            let mut built = Vec::with_capacity(partitions);
            for (partition, (keys, masks)) in
                keys_by_part.into_iter().zip(masks_by_part).enumerate()
            {
                let part_stage = stage
                    .join("partitions")
                    .join(format!("{source_slots}-{partition}"));
                let cache = if config.merge.in_memory() && !config.merge.cache_features() {
                    None
                } else {
                    fs::create_dir_all(part_stage.join("cells"))?;
                    if !config.merge.drop_abandoned {
                        fs::create_dir_all(part_stage.join("abandoned"))?;
                    }
                    Some(part_stage.as_path())
                };
                let (mut part, _) =
                    construct_group_at_slots(group, keys, masks, cache, config, graph, local)?;
                part.routing = Some(Routing {
                    source_slots,
                    partition,
                    partitions,
                });
                built.push(part);
            }
            Ok(built)
        },
        |_, parts| {
            result.extend(parts);
            Ok(())
        },
    )?;
    Ok(result)
}

fn order_columns(groups: &mut [BuiltGroup], order: ColumnOrder) {
    // Keep equal-width sections contiguous, then reuse the final union-sketch
    // fingerprints. No filter-cell sampling, feature pass or new sketch is needed.
    groups.sort_by(|a, b| {
        a.meta
            .width
            .cmp(&b.meta.width)
            .then_with(|| {
                if order == ColumnOrder::Compression {
                    a.ordering_key.cmp(&b.ordering_key)
                } else {
                    std::cmp::Ordering::Equal
                }
            })
            .then_with(|| a.meta.members[0].cmp(&b.meta.members[0]))
    });
}

#[cfg(test)]
fn transpose(groups: Vec<BuiltGroup>, order: ColumnOrder) -> Result<Vec<Bank>> {
    let mut by_size = BTreeMap::<usize, Vec<BuiltGroup>>::new();
    for group in groups {
        by_size.entry(group.slots).or_default().push(group);
    }
    let mut banks = Vec::new();
    for (slots, mut groups) in by_size {
        order_columns(&mut groups, order);
        let mut row_bits = 0;
        for group in &mut groups {
            group.meta.offset = row_bits;
            row_bits += group.meta.width;
        }
        let row_bytes = row_bits.div_ceil(8);
        let mut data = vec![0; slots.checked_mul(row_bytes).context("bank size overflow")?];
        for group in &groups {
            let cells = match &group.cells {
                CompletedCells::Disk(path) => cache_read(path)?,
                CompletedCells::Packed { .. } | CompletedCells::PackedDisk { .. } => {
                    let mut bytes = vec![0; slots * 8];
                    group.cells.read_block(0, &mut bytes, None)?;
                    bytes
                        .chunks_exact(8)
                        .map(|b| u64::from_le_bytes(b.try_into().unwrap()))
                        .collect()
                }
            };
            ensure!(cells.len() == slots, "invalid completed filter cell cache");
            for (slot, &value) in cells.iter().enumerate() {
                put_bits(
                    &mut data[slot * row_bytes..(slot + 1) * row_bytes],
                    group.meta.offset,
                    group.meta.width,
                    value,
                );
            }
            if let CompletedCells::Disk(path) = &group.cells {
                fs::remove_file(path)?;
            }
        }
        let meta = BankMeta {
            component: 0,
            file: format!("bank-{}.zorbank.zst", slots.ilog2()),
            slots,
            row_bits,
            groups: groups.iter().map(|g| g.meta.clone()).collect(),
            fingerprint_first: false,
            routing: None,
        };
        banks.push(Bank {
            meta,
            data,
            fallback: groups.into_iter().map(|g| g.fallback).collect(),
        });
    }
    Ok(banks)
}

fn pack_ram_rows(
    groups: &[BuiltGroup],
    fp_bits: usize,
    start: usize,
    output: &mut [u8],
    row_bytes: usize,
) {
    let fp_mask = (1u64 << fp_bits) - 1;
    let rows = output.len() / row_bytes;
    for (index, group) in groups.iter().enumerate() {
        let CompletedCells::Packed { data, width, .. } = &group.cells else {
            unreachable!()
        };
        if *width == fp_bits && matches!(fp_bits, 8 | 16) {
            let bytes = fp_bits / 8;
            let offset = index * bytes;
            for (row, cell) in output
                .chunks_exact_mut(row_bytes)
                .zip(data[start * bytes..(start + rows) * bytes].chunks_exact(bytes))
            {
                row[offset..offset + bytes].copy_from_slice(cell);
            }
            continue;
        }
        let origin_offset = groups.len() * fp_bits + group.meta.offset - index * fp_bits;
        for (n, row) in output.chunks_exact_mut(row_bytes).enumerate() {
            let value = get_bits(data, (start + n) * width, *width);
            put_bits(row, index * fp_bits, fp_bits, value & fp_mask);
            if *width > fp_bits {
                put_bits(row, origin_offset, width - fp_bits, value >> fp_bits);
            }
        }
    }
}

fn pack_banks_streaming(
    groups: Vec<BuiltGroup>,
    order: ColumnOrder,
    fp_bits: usize,
    directory: &Path,
) -> Result<Vec<BankMeta>> {
    let mut by_size = BTreeMap::<(usize, usize), Vec<BuiltGroup>>::new();
    for group in groups {
        let key = group
            .routing
            .map_or((group.slots, 0), |r| (r.source_slots, r.partition));
        by_size.entry(key).or_default().push(group);
    }
    let mut banks = Vec::new();
    for (_, mut groups) in by_size {
        let slots = groups[0].slots;
        let routing = groups[0].routing;
        ensure!(
            groups
                .iter()
                .all(|g| g.slots == slots && g.routing == routing),
            "inconsistent partition geometry"
        );
        order_columns(&mut groups, order);
        let mut row_bits = 0usize;
        for group in &mut groups {
            group.meta.offset = row_bits;
            row_bits = row_bits
                .checked_add(group.meta.width)
                .context("row width overflow")?;
            group.cells.validate(slots, group.meta.width)?;
        }
        let meta = BankMeta {
            component: 0,
            file: BankMeta::filename(slots, routing),
            slots,
            row_bits,
            groups: groups.iter().map(|g| g.meta.clone()).collect(),
            fingerprint_first: true,
            routing,
        };
        let row_bytes = row_bits.div_ceil(8);
        let expected_bytes = slots.checked_mul(row_bytes).context("bank size overflow")?;
        let mut writer = bank_io::Writer::new(&directory.join(&meta.file))?;
        writer.write_all(PLANE_MAGIC)?;
        for n in [slots, row_bits, groups.len(), expected_bytes] {
            writer.write_all(&(n as u64).to_le_bytes())?;
        }
        // Keep all row writes in an 8-MiB tile. Input blocks stay small and each
        // completed bank goes directly to its compressor instead of a Vec<Bank>.
        let tile_rows = (8 * 1024 * 1024 / row_bytes).clamp(1, 4096).min(slots);
        let mut tile = vec![0; tile_rows * row_bytes];
        let mut input = vec![0; tile_rows * 8];
        // Reuse descriptors when allowed. If a host's FD limit is too small,
        // discard partial opens and read one file/block at a time instead.
        let mut files = groups
            .iter()
            .map(|g| g.cells.open())
            .collect::<std::io::Result<Vec<_>>>()
            .ok();
        let fp_mask = (1u64 << fp_bits) - 1;
        let ram = groups
            .iter()
            .all(|g| matches!(g.cells, CompletedCells::Packed { .. }));
        for start in (0..slots).step_by(tile_rows) {
            let rows = tile_rows.min(slots - start);
            let output = &mut tile[..rows * row_bytes];
            output.fill(0);
            if ram {
                // Rows are disjoint ownership units, even when color planes
                // are not byte aligned. Preserve the exact writer tile order.
                if output.len() >= 1024 * 1024 && rows >= 128 && groups.len() >= 32 {
                    let per_task = rows.div_ceil(rayon::current_num_threads());
                    output
                        .par_chunks_mut(per_task * row_bytes)
                        .enumerate()
                        .for_each(|(task, rows)| {
                            pack_ram_rows(
                                &groups,
                                fp_bits,
                                start + task * per_task,
                                rows,
                                row_bytes,
                            );
                        });
                } else {
                    pack_ram_rows(&groups, fp_bits, start, output, row_bytes);
                }
                writer.write_all(output)?;
                continue;
            }
            for (index, group) in groups.iter().enumerate() {
                let block = &mut input[..rows * 8];
                let file = files.as_mut().and_then(|files| files[index].as_mut());
                group.cells.read_block(start, block, file)?;
                let origin_offset = groups.len() * fp_bits + group.meta.offset - index * fp_bits;
                for (row, bytes) in output
                    .chunks_exact_mut(row_bytes)
                    .zip(block.chunks_exact(8))
                {
                    let value = u64::from_le_bytes(bytes.try_into().unwrap());
                    put_bits(row, index * fp_bits, fp_bits, value & fp_mask);
                    if group.meta.width > fp_bits {
                        put_bits(
                            row,
                            origin_offset,
                            group.meta.width - fp_bits,
                            value >> fp_bits,
                        );
                    }
                }
            }
            writer.write_all(output)?;
        }
        writer.finish()?;
        drop(files);
        for group in groups {
            group.cells.remove()?;
        }
        banks.push(meta);
    }
    Ok(banks)
}

type AbandonedRecord = (u32, u32, u64); // hash, construction group, origin mask

fn visit_abandoned(
    records: &mut [AbandonedRecord],
    members: &[Vec<usize>],
    mut emit: impl FnMut(u32, &[u32]) -> Result<()>,
) -> Result<(usize, usize)> {
    // Only hash order matters. Origin lists are sorted/deduplicated below;
    // construction-group order within an equal hash has no observable effect.
    if records.len() >= 16_384 {
        records.par_sort_unstable_by_key(|record| record.0);
    } else {
        records.sort_unstable_by_key(|record| record.0);
    }
    let mut origins = Vec::new();
    let mut hashes = 0usize;
    let mut associations = 0usize;
    for chunk in records.chunk_by(|a, b| a.0 == b.0) {
        origins.clear();
        for &(_, group, mut mask) in chunk {
            let group_members = &members[group as usize];
            while mask != 0 {
                let bit = mask.trailing_zeros() as usize;
                ensure!(
                    bit < group_members.len(),
                    "invalid abandoned membership bit"
                );
                origins
                    .push(u32::try_from(group_members[bit]).context("dataset ID exceeds 32 bits")?);
                mask &= mask - 1;
            }
        }
        origins.sort_unstable();
        origins.dedup();
        ensure!(!origins.is_empty(), "empty abandoned membership");
        emit(chunk[0].0, &origins)?;
        hashes += 1;
        associations = associations
            .checked_add(origins.len())
            .context("abandoned origins overflow")?;
    }
    ensure!(hashes <= u32::MAX as usize, "too many abandoned hashes");
    Ok((hashes, associations))
}

#[cfg(test)]
fn save_abandoned(banks: &mut [Bank], directory: &Path) -> Result<(usize, usize)> {
    let count: usize = banks
        .iter()
        .flat_map(|bank| &bank.fallback)
        .map(Vec::len)
        .sum();
    let mut records = Vec::with_capacity(count);
    let mut members = Vec::new();
    for bank in banks {
        for (group, fallback) in bank.meta.groups.iter().zip(&mut bank.fallback) {
            let id = u32::try_from(members.len()).context("group ID exceeds 32 bits")?;
            members.push(group.members.clone());
            for (hash, mask) in std::mem::take(fallback) {
                records.push((hash, id, mask));
            }
        }
        // Banks have no per-group recovery table in the query representation.
        bank.fallback.clear();
    }
    let mut writer = create_index_writer(&directory.join(ABANDONED_FILE))?;
    writer.write_all(ABANDONED_MAGIC)?;
    let stats = visit_abandoned(&mut records, &members, |hash, origins| {
        writer.write_all(&hash.to_le_bytes())?;
        let count =
            u32::try_from(origins.len()).context("too many origins for an abandoned hash")?;
        writer.write_all(&count.to_le_bytes())?;
        for &origin in origins {
            writer.write_all(&origin.to_le_bytes())?;
        }
        Ok(())
    })?;
    finish_index_writer(writer)?;
    Ok(stats)
}

fn save_abandoned_groups(groups: &mut [BuiltGroup], directory: &Path) -> Result<(usize, usize)> {
    if groups.iter().any(|group| group.fallback_path.is_some()) {
        ensure!(
            groups.iter().all(|group| group.fallback.is_empty()),
            "mixed RAM/disk abandoned-key storage"
        );
        let paths = groups
            .iter_mut()
            .filter_map(|group| group.fallback_path.take())
            .collect();
        let max_id = groups
            .iter()
            .map(|group| group.meta.members[0])
            .max()
            .unwrap_or(0);
        // Only construction metadata is resident: IDs address original groups
        // by their first member. No abandoned-key vector is materialized.
        let mut members = vec![&[][..]; max_id + 1];
        for group in groups.iter() {
            members[group.meta.members[0]] = group.meta.members.as_slice();
        }
        let origin_limit = groups
            .iter()
            .map(|group| group.meta.members.len())
            .sum::<usize>()
            .saturating_mul(2)
            .max(1024);
        return write_abandoned_map(directory, |emit| {
            let mut current = None;
            let mut origins = Vec::new();
            let (mut hashes, mut associations) = (0usize, 0usize);
            let mut flush = |hash, origins: &mut Vec<u32>| -> Result<()> {
                origins.sort_unstable();
                origins.dedup();
                ensure!(!origins.is_empty(), "empty abandoned membership");
                emit(hash, origins)?;
                hashes += 1;
                associations = associations
                    .checked_add(origins.len())
                    .context("abandoned origins overflow")?;
                origins.clear();
                Ok(())
            };
            super::staged_abandoned::visit_sorted(
                paths,
                &directory.join("abandoned"),
                |(hash, group, mut mask)| {
                    if current.is_some_and(|previous| previous != hash) {
                        flush(current.unwrap(), &mut origins)?;
                    }
                    current = Some(hash);
                    let group_members = members
                        .get(group as usize)
                        .context("invalid abandoned group ID")?;
                    while mask != 0 {
                        let bit = mask.trailing_zeros() as usize;
                        let origin = group_members
                            .get(bit)
                            .context("invalid abandoned membership bit")?;
                        origins.push(u32::try_from(*origin).context("dataset ID exceeds 32 bits")?);
                        mask &= mask - 1;
                    }
                    // Even a deliberately pathological single-hash collision cannot
                    // accumulate one origin entry per abandoned key indefinitely.
                    if origins.len() >= origin_limit {
                        origins.sort_unstable();
                        origins.dedup();
                    }
                    Ok(())
                },
            )?;
            if let Some(hash) = current {
                flush(hash, &mut origins)?;
            }
            ensure!(hashes <= u32::MAX as usize, "too many abandoned hashes");
            Ok((hashes, associations))
        });
    }
    let mut records = Vec::with_capacity(groups.iter().map(|g| g.fallback.len()).sum());
    let mut members = Vec::with_capacity(groups.len());
    for group in groups {
        let id = u32::try_from(members.len()).context("group ID exceeds 32 bits")?;
        members.push(group.meta.members.clone());
        for (hash, mask) in std::mem::take(&mut group.fallback) {
            records.push((hash, id, mask));
        }
    }
    write_abandoned_map(directory, |emit| {
        visit_abandoned(&mut records, &members, emit)
    })
}

fn write_abandoned_map(
    directory: &Path,
    visit: impl FnOnce(&mut dyn FnMut(u32, &[u32]) -> Result<()>) -> Result<(usize, usize)>,
) -> Result<(usize, usize)> {
    // Buffer the uncompressed records too: entering the zstd encoder for every
    // tiny integer is expensive on hundreds of millions of abandoned hashes.
    let mut writer = BufWriter::with_capacity(
        1024 * 1024,
        create_index_writer(&directory.join(ABANDONED_FILE))?,
    );
    writer.write_all(ABANDONED_DELTA_MAGIC)?;
    let mut previous_hash = 0u32;
    let stats = visit(&mut |hash, origins| {
        write_var_u32(&mut writer, hash - previous_hash)?;
        previous_hash = hash;
        write_var_u32(
            &mut writer,
            u32::try_from(origins.len()).context("too many abandoned origins")?,
        )?;
        let mut previous_origin = 0u32;
        for &origin in origins {
            write_var_u32(&mut writer, origin - previous_origin)?;
            previous_origin = origin;
        }
        Ok(())
    })?;
    finish_index_writer(writer.into_inner().map_err(|e| e.into_error())?)?;
    Ok(stats)
}

fn write_var_u32(writer: &mut impl Write, mut value: u32) -> Result<()> {
    let mut bytes = [0u8; 5];
    let mut n = 0;
    while value >= 128 {
        bytes[n] = value as u8 | 128;
        value >>= 7;
        n += 1;
    }
    bytes[n] = value as u8;
    writer.write_all(&bytes[..n + 1])?;
    Ok(())
}

#[inline]
fn read_var_u32(reader: &mut impl BufRead) -> Result<u32> {
    // Decode directly in the decompressor's output buffer. Nearly all hash
    // deltas/counts fit one byte and document IDs two; avoid a Read call per
    // byte. The bounded fallback also handles integers crossing buffer edges.
    let bytes = reader.fill_buf()?;
    if bytes.len() >= 5 || bytes.first().is_some_and(|b| b & 128 == 0) {
        let mut value = 0u32;
        for (n, &byte) in bytes.iter().take(5).enumerate() {
            ensure!(n != 4 || byte <= 15, "abandoned variable integer overflow");
            value |= ((byte & 127) as u32) << (n * 7);
            if byte & 128 == 0 {
                reader.consume(n + 1);
                return Ok(value);
            }
        }
        anyhow::bail!("abandoned variable integer overflow");
    }
    let mut value = 0u32;
    for shift in (0..=28).step_by(7) {
        let mut byte = [0];
        reader.read_exact(&mut byte)?;
        ensure!(
            shift != 28 || byte[0] <= 15,
            "abandoned variable integer overflow"
        );
        value |= ((byte[0] & 127) as u32) << shift;
        if byte[0] & 128 == 0 {
            return Ok(value);
        }
    }
    anyhow::bail!("abandoned variable integer overflow")
}

fn read_u32(reader: &mut impl Read) -> Result<u32> {
    let mut bytes = [0; 4];
    reader.read_exact(&mut bytes)?;
    Ok(u32::from_le_bytes(bytes))
}

fn load_abandoned(directory: &Path, manifest: &Manifest) -> Result<OriginMap> {
    load_abandoned_selected(directory, manifest, None)
}

fn load_abandoned_selected(
    directory: &Path,
    manifest: &Manifest,
    remap: Option<&[u32]>,
) -> Result<OriginMap> {
    let mut map = if remap.is_some() {
        OriginMap::default()
    } else {
        OriginMap::with_capacity(manifest.abandoned_hashes, manifest.abandoned_origins)?
    };
    let mut selected = Vec::new();
    visit_abandoned_file(directory, manifest, |hash, origins| {
        if let Some(remap) = remap {
            selected.clear();
            selected.extend(
                origins.iter().filter_map(|&id| {
                    (remap[id as usize] != u32::MAX).then_some(remap[id as usize])
                }),
            );
            if !selected.is_empty() {
                map.push(hash, &selected)?;
            }
        } else {
            map.push(hash, origins)?;
        }
        Ok(())
    })?;
    map.finish()?;
    Ok(map)
}

fn visit_abandoned_file(
    directory: &Path,
    manifest: &Manifest,
    mut visit: impl FnMut(u32, &[u32]) -> Result<()>,
) -> Result<()> {
    let mut reader = AbandonedReader::new(directory, manifest)?;
    while let Some(hash) = reader.advance()? {
        visit(hash, &reader.origins)?;
    }
    Ok(())
}

#[cfg(test)]
fn save_bank(bank: &Bank, directory: &Path) -> Result<()> {
    let mut writer = create_index_writer(&directory.join(&bank.meta.file))?;
    writer.write_all(if bank.meta.fingerprint_first {
        PLANE_MAGIC
    } else {
        MAGIC
    })?;
    for value in [
        bank.meta.slots,
        bank.meta.row_bits,
        bank.meta.groups.len(),
        bank.data.len(),
    ] {
        writer.write_all(&(value as u64).to_le_bytes())?;
    }
    writer.write_all(&bank.data)?;
    finish_index_writer(writer)
}
fn read_u64(reader: &mut impl Read) -> Result<u64> {
    let mut bytes = [0; 8];
    reader.read_exact(&mut bytes)?;
    Ok(u64::from_le_bytes(bytes))
}

fn validate_bank_header(
    reader: &mut impl Read,
    meta: &BankMeta,
    expected_bytes: usize,
) -> Result<()> {
    let mut magic = [0; 8];
    reader.read_exact(&mut magic)?;
    ensure!(
        &magic
            == if meta.fingerprint_first {
                PLANE_MAGIC
            } else {
                MAGIC
            },
        "invalid bank magic"
    );
    for expected in [meta.slots, meta.row_bits, meta.groups.len(), expected_bytes] {
        ensure!(
            read_u64(reader)? == expected as u64,
            "bank header differs from manifest"
        );
    }
    Ok(())
}
fn load_modern_bank(directory: &Path, meta: &BankMeta) -> Result<Bank> {
    let expected_bytes = meta
        .slots
        .checked_mul(meta.row_bits.div_ceil(8))
        .context("bank size overflow")?;
    let path = directory.join(&meta.file);
    let data = if let Some((header, data)) = bank_io::read_parallel(&path, expected_bytes)? {
        validate_bank_header(&mut header.as_slice(), meta, expected_bytes)?;
        data
    } else {
        let mut reader = open_index_reader(&path)?;
        validate_bank_header(&mut reader, meta, expected_bytes)?;
        let mut data = bank_io::zeroed(expected_bytes);
        reader.read_exact(&mut data)?;
        let mut extra = [0];
        ensure!(reader.read(&mut extra)? == 0, "unexpected bank data");
        data
    };
    Ok(Bank {
        meta: meta.clone(),
        data,
        #[cfg(test)]
        fallback: Vec::new(),
    })
}

impl Index {
    fn load(directory: &Path) -> Result<Self> {
        Self::load_with_recovery(directory, true)
    }

    fn load_with_recovery(directory: &Path, recover: bool) -> Result<Self> {
        Self::load_selected(directory, recover, None)
    }

    // Validate all metadata, but materialize only the selected subindex and its
    // abandoned-key associations. Local document IDs keep worker state small.
    fn load_selected(directory: &Path, recover: bool, selected: Option<usize>) -> Result<Self> {
        Self::load_mode(directory, recover, selected, false)
    }

    fn load_mode(
        directory: &Path,
        recover: bool,
        selected: Option<usize>,
        metadata_only: bool,
    ) -> Result<Self> {
        let mut manifest: Manifest =
            serde_json::from_reader(BufReader::new(File::open(directory.join("manifest.json"))?))?;
        let remap = if let Some(id) = selected {
            let meta = manifest
                .banks
                .get(id)
                .context("invalid selected subindex")?;
            let mut ids: Vec<_> = meta
                .groups
                .iter()
                .flat_map(|g| g.members.iter().copied())
                .collect();
            ids.sort_unstable();
            let mut remap = vec![u32::MAX; manifest.colors.len()];
            for (local, global) in ids.into_iter().enumerate() {
                ensure!(global < remap.len(), "invalid dataset mapping");
                remap[global] = u32::try_from(local).context("dataset IDs must fit in 32 bits")?;
            }
            Some(remap)
        } else {
            None
        };
        ensure!(
            (1..=10).contains(&manifest.version) && manifest.backend == Backend::Zor,
            "unsupported power-of-two index format"
        );
        ensure!(
            (2..=MAX_HASHES).contains(&manifest.hashes),
            "invalid hash count"
        );
        ensure!(
            if manifest.version >= 9 {
                manifest
                    .segment_length
                    .is_none_or(|n| n.is_power_of_two() && n <= 4096)
            } else if manifest.version == 8 {
                manifest
                    .segment_length
                    .is_some_and(|n| n.is_power_of_two() && n <= 4096)
            } else {
                manifest.segment_length.is_none()
            },
            "invalid/missing segment length for index version"
        );
        ensure!(
            [1, 2, 4, 8, 16].contains(&manifest.fingerprint_bits),
            "invalid fingerprint width"
        );
        ensure!(
            manifest.colors.len() <= u32::MAX as usize,
            "dataset IDs must fit in 32 bits"
        );
        if manifest.version >= 2 {
            ensure!(
                manifest.abandoned_recovery == manifest.abandoned_file.is_some(),
                "inconsistent abandoned recovery metadata"
            );
            ensure!(
                manifest.abandoned_file.is_none()
                    || manifest.abandoned_file.as_deref() == Some(ABANDONED_FILE),
                "invalid abandoned map filename"
            );
            ensure!(
                manifest.abandoned_recovery
                    || (manifest.abandoned_hashes == 0 && manifest.abandoned_origins == 0),
                "disabled recovery has abandoned records"
            );
        }
        let recover = recover && manifest.abandoned_recovery;
        let mut legacy_records = Vec::new();
        let mut legacy_members = Vec::new();
        let features = manifest.features()?;
        match features.mode {
            IndexMode::Findere => {
                crate::findere::indexed_length(manifest.k, features.findere_z)?;
            }
            IndexMode::Kmers | IndexMode::Modimizers => {
                ensure!(manifest.k > 0, "invalid k")
            }
            IndexMode::Minimizers => ensure!(
                (1..=255).contains(&manifest.k)
                    && manifest.k % 2 == 1
                    && features.minimizer_size <= manifest.k,
                "invalid minimizer window"
            ),
        }
        ensure!(
            (1..=31).contains(&features.minimizer_size) && features.modimizer_sampling > 0,
            "invalid feature parameters"
        );
        let mut seen = vec![None; manifest.colors.len()];
        let mut groups_seen = vec![None; manifest.colors.len()];
        let classes = manifest
            .banks
            .iter()
            .map(BankMeta::source_group)
            .collect::<std::collections::BTreeSet<_>>()
            .len();
        ensure!(
            manifest
                .max_subindexes
                .is_none_or(|limit| limit > 0 && classes <= limit),
            "invalid maximum subindex count"
        );
        let mut banks = Vec::new();
        ensure!(
            manifest.max_subindex_size.is_none_or(u64::is_power_of_two),
            "invalid maximum subindex size"
        );
        let mut previous_group = (0, 0);
        let mut next_partition = 0;
        let mut partition_count = 1;
        let mut first_class_bank = 0;
        for (bank_id, meta) in manifest.banks.iter().enumerate() {
            ensure!(
                meta.component == 0 || manifest.version >= 10,
                "unsupported index component"
            );
            ensure!(
                !meta.fingerprint_first || manifest.version >= 3,
                "unsupported bank layout"
            );
            ensure!(
                meta.slots.is_power_of_two()
                    && meta.slots >= manifest.hashes
                    && meta.source_slots().is_power_of_two(),
                "invalid bank size/order"
            );
            let (partition, partitions) =
                meta.routing.map_or((0, 1), |r| (r.partition, r.partitions));
            ensure!(
                meta.routing.is_none()
                    || (manifest.version >= 9
                        && partitions > 1
                        && partitions.is_power_of_two()
                        && meta.slots.checked_mul(partitions) == Some(meta.source_slots())),
                "invalid hash partition geometry"
            );
            if meta.source_group() != previous_group {
                ensure!(
                    meta.source_group() > previous_group
                        && (previous_group == (0, 0) || next_partition == partition_count),
                    "invalid/incomplete partition order"
                );
                previous_group = meta.source_group();
                next_partition = 0;
                partition_count = partitions;
                first_class_bank = bank_id;
            }
            ensure!(
                partitions == partition_count
                    && partition == next_partition
                    && partition < partitions,
                "invalid/duplicate/missing hash partition"
            );
            next_partition += 1;
            if let Some(length) = manifest.segment_length {
                ensure!(
                    meta.slots % length == 0 && meta.slots / length >= manifest.hashes,
                    "subindex too small for fixed segment geometry"
                );
            }
            ensure!(
                meta.file == meta.expected_filename(),
                "invalid bank file name"
            );
            let mut offset = 0;
            let mut local_seen = std::collections::HashSet::new();
            for (group_id, group) in meta.groups.iter().enumerate() {
                ensure!(
                    !group.members.is_empty()
                        && group.members.len() <= 64 - manifest.fingerprint_bits,
                    "invalid membership width"
                );
                let width = manifest.fingerprint_bits
                    + if group.members.len() == 1 {
                        0
                    } else {
                        group.members.len()
                    };
                ensure!(
                    group.offset == offset && group.width == width,
                    "invalid group bit layout"
                );
                offset = offset.checked_add(width).context("row width overflow")?;
                for &id in &group.members {
                    ensure!(
                        id < seen.len()
                            && local_seen.insert(id)
                            && seen[id].is_none_or(|source| source == meta.source_group()),
                        "invalid/duplicate dataset mapping"
                    );
                    if partition == 0 {
                        seen[id] = Some(meta.source_group());
                        groups_seen[id] = Some((bank_id, group_id));
                    } else {
                        let (first_bank, first_group) =
                            groups_seen[id].context("missing partition dataset")?;
                        ensure!(
                            seen[id] == Some(meta.source_group())
                                && manifest.banks[first_bank].groups[first_group].members
                                    == group.members,
                            "inconsistent partition dataset mapping"
                        );
                    }
                }
                ensure!(
                    group.fallback_keys <= group.union_cardinality,
                    "invalid fallback count"
                );
                if manifest.version >= 2 {
                    let count = group
                        .abandoned_keys
                        .context("missing abandoned equation count")?;
                    ensure!(
                        count <= group.union_cardinality,
                        "invalid abandoned equation count"
                    );
                    ensure!(
                        group.fallback_keys
                            == if manifest.abandoned_recovery {
                                count
                            } else {
                                0
                            },
                        "inconsistent abandoned equation count"
                    );
                }
            }
            ensure!(offset == meta.row_bits && offset > 0, "invalid row width");
            if partition > 0 {
                let first = &manifest.banks[first_class_bank];
                ensure!(
                    first.row_bits == meta.row_bits
                        && first.groups.len() == meta.groups.len()
                        && local_seen.len()
                            == first.groups.iter().map(|g| g.members.len()).sum::<usize>(),
                    "incomplete partition dataset mapping"
                );
            }
            let bytes = meta
                .slots
                .checked_mul(meta.row_bits.div_ceil(8))
                .context("bank size overflow")?;
            ensure!(
                manifest
                    .max_subindex_size
                    .is_none_or(|limit| bytes as u64 <= limit),
                "subindex exceeds maximum uncompressed size"
            );
            // All versions retain the same validation. Version 1 needs its
            // embedded legacy records converted below; newer banks load later.
            if manifest.version >= 2 {
                meta.slots
                    .checked_mul(meta.row_bits.div_ceil(8))
                    .context("bank size overflow")?;
                continue;
            }
            if metadata_only || selected.is_some_and(|id| id != bank_id) {
                continue;
            }
            let expected_bytes = meta
                .slots
                .checked_mul(meta.row_bits.div_ceil(8))
                .context("bank size overflow")?;
            let mut reader = open_index_reader(&directory.join(&meta.file))?;
            let mut magic = [0; 8];
            reader.read_exact(&mut magic)?;
            ensure!(
                &magic
                    == if meta.fingerprint_first {
                        PLANE_MAGIC
                    } else {
                        MAGIC
                    },
                "invalid bank magic"
            );
            for expected in [meta.slots, meta.row_bits, meta.groups.len(), expected_bytes] {
                ensure!(
                    read_u64(&mut reader)? == expected as u64,
                    "bank header differs from manifest"
                );
            }
            let mut data = vec![0; expected_bytes];
            reader.read_exact(&mut data)?;
            if manifest.version == 1 {
                for group in &meta.groups {
                    ensure!(
                        read_u64(&mut reader)? == group.fallback_keys as u64,
                        "invalid fallback length"
                    );
                    let group_id =
                        u32::try_from(legacy_members.len()).context("group ID exceeds 32 bits")?;
                    if recover {
                        legacy_members.push(
                            group
                                .members
                                .iter()
                                .map(|&id| remap.as_ref().map_or(id, |remap| remap[id] as usize))
                                .collect(),
                        );
                    }
                    let mut previous = None;
                    for _ in 0..group.fallback_keys {
                        let key = read_u64(&mut reader)?;
                        let value = read_u64(&mut reader)?;
                        ensure!(
                            group.width == 64 || value < (1u64 << group.width),
                            "invalid fallback payload width"
                        );
                        let fp_mask = (1u64 << manifest.fingerprint_bits) - 1;
                        ensure!(
                            value & fp_mask
                                == fingerprint(key, manifest.seed, manifest.fingerprint_bits),
                            "invalid fallback fingerprint"
                        );
                        ensure!(
                            previous.is_none_or(|last| last < key),
                            "unsorted fallback keys"
                        );
                        previous = Some(key);
                        if recover {
                            let mask = if group.members.len() == 1 {
                                1
                            } else {
                                value >> manifest.fingerprint_bits
                            };
                            ensure!(mask != 0, "empty fallback membership");
                            legacy_records.push((
                                abandoned::hash(key, manifest.seed),
                                group_id,
                                mask,
                            ));
                        }
                    }
                }
            }
            let mut extra = [0];
            ensure!(reader.read(&mut extra)? == 0, "unexpected bank data");
            banks.push(Bank {
                meta: meta.clone(),
                data,
                #[cfg(test)]
                fallback: Vec::new(),
            });
        }
        ensure!(
            !manifest.banks.is_empty()
                && seen.iter().all(Option::is_some)
                && next_partition == partition_count,
            "incomplete dataset mapping"
        );
        if metadata_only {
            return Ok(Self {
                manifest,
                banks: Vec::new(),
                abandoned: None,
                row_scratch_bytes: 0,
            });
        }
        let abandoned = if manifest.version >= 2 {
            let (loaded_banks, loaded_map) = rayon::join(
                || {
                    manifest
                        .banks
                        .par_iter()
                        .enumerate()
                        .filter(|(id, _)| selected.is_none_or(|selected| selected == *id))
                        .map(|(_, meta)| load_modern_bank(directory, meta))
                        .collect::<Result<Vec<_>>>()
                },
                || {
                    if recover {
                        if let Some(remap) = &remap {
                            load_abandoned_selected(directory, &manifest, Some(remap)).map(Some)
                        } else {
                            load_abandoned(directory, &manifest).map(Some)
                        }
                    } else {
                        Ok(None)
                    }
                },
            );
            // Indexed parallel collection preserves the manifest's bank order.
            banks = loaded_banks?;
            loaded_map?
        } else if !recover {
            None
        } else {
            let mut map = OriginMap::default();
            visit_abandoned(&mut legacy_records, &legacy_members, |hash, origins| {
                map.push(hash, origins)
            })?;
            drop(legacy_records);
            map.finish()?;
            Some(map)
        };
        let abandoned = abandoned.filter(|map| map.len() != 0);
        if let Some(remap) = &remap {
            for bank in &mut banks {
                for group in &mut bank.meta.groups {
                    for member in &mut group.members {
                        *member = remap[*member] as usize;
                    }
                }
            }
            manifest.colors = manifest
                .colors
                .into_iter()
                .enumerate()
                .filter_map(|(id, color)| (remap[id] != u32::MAX).then_some(color))
                .collect();
            manifest.banks = banks.iter().map(|bank| bank.meta.clone()).collect();
        }
        let row_scratch_bytes = banks
            .iter()
            // Fingerprint-first banks now compare directly in registers. Only
            // legacy banks require a materialised row before origin decoding.
            .filter(|bank| !bank.meta.fingerprint_first)
            .map(|bank| bank.meta.row_bits.div_ceil(8))
            .max()
            .unwrap_or(0);
        Ok(Self {
            manifest,
            banks,
            abandoned,
            row_scratch_bytes,
        })
    }

    fn lookup(&self, key: u64, scores: &mut [u64], scratch: &mut Vec<u8>) {
        self.lookup_into(key, scores, scratch);
    }

    // Full-RAM indexes contain complete contiguous partition runs. Jump over
    // the run directly: H-position preparation, prefetch and lookup must not
    // grow linearly with the number of physical partitions. A selectively
    // loaded index can contain just one member of such a run.
    fn selected_banks(&self, key: u64) -> impl Iterator<Item = (usize, &Bank)> {
        let mut next = 0;
        std::iter::from_fn(move || {
            while next < self.banks.len() {
                let start = next;
                let bank = &self.banks[start];
                if let Some(routing) = bank.meta.routing {
                    let target = partition_id(key, self.manifest.seed, routing.partitions);
                    if routing.partition == 0 && routing.partitions <= self.banks.len() - start {
                        next += routing.partitions;
                        let id = start + target;
                        debug_assert_eq!(
                            self.banks[id].meta.routing.unwrap().source_slots,
                            routing.source_slots
                        );
                        return Some((id, &self.banks[id]));
                    }
                    next += 1;
                    if target == routing.partition {
                        return Some((start, bank));
                    }
                } else {
                    next += 1;
                    return Some((start, bank));
                }
            }
            None
        })
    }

    fn lookup_into<S: ScoreSink + ?Sized>(&self, key: u64, scores: &mut S, scratch: &mut Vec<u8>) {
        self.lookup_prepared(key, scores, scratch, None);
    }

    fn lookup_prepared<S: ScoreSink + ?Sized>(
        &self,
        key: u64,
        scores: &mut S,
        scratch: &mut Vec<u8>,
        prepared: Option<&[usize]>,
    ) {
        self.lookup_selected(key, scores, scratch, prepared, None);
    }

    fn lookup_selected<S: ScoreSink + ?Sized>(
        &self,
        key: u64,
        scores: &mut S,
        scratch: &mut Vec<u8>,
        prepared: Option<&[usize]>,
        selected: Option<&[Vec<Range<usize>>]>,
    ) {
        // Exactly one map probe for the feature, shared by all groups and banks.
        // Recovery is an additive union, including for 32-bit hash collisions.
        let recovered = self.abandoned.as_ref().map_or(&[][..], |map| {
            map.get(abandoned::hash(key, self.manifest.seed))
        });
        // Save only the recovered datasets' counters beyond the row workspace.
        // After normal filter lookup, add recovery only where the counter did not
        // advance. This avoids a set lookup for every positive filter hit, and the
        // same worker-owned buffer is reused without a per-feature allocation.
        let snapshots_end = self.row_scratch_bytes + recovered.len() * 8;
        if scratch.len() < snapshots_end {
            scratch.resize(snapshots_end, 0);
        }
        for (bytes, &origin) in scratch[self.row_scratch_bytes..snapshots_end]
            .chunks_exact_mut(8)
            .zip(recovered)
        {
            bytes.copy_from_slice(&scores.get(origin as usize).to_le_bytes());
        }
        let fp = fingerprint(key, self.manifest.seed, self.manifest.fingerprint_bits);
        let fp_mask = (1u64 << self.manifest.fingerprint_bits) - 1;
        let position_hash = prepared
            .is_none()
            .then(|| PositionHash::new(key, self.manifest.seed));
        for (bank_id, bank) in self.selected_banks(key) {
            let ranges = selected.map(|banks| banks[bank_id].as_slice());
            if ranges.is_some_and(|ranges| ranges.is_empty()) {
                continue;
            }
            let row_bytes = bank.meta.row_bits.div_ceil(8);
            let lookup_bytes = if bank.meta.fingerprint_first {
                (bank.meta.groups.len() * self.manifest.fingerprint_bits).div_ceil(8)
            } else {
                row_bytes
            };
            let mut computed;
            let slots = if let Some(prepared) = prepared {
                &prepared[bank_id * self.manifest.hashes..(bank_id + 1) * self.manifest.hashes]
            } else {
                computed = [0; MAX_HASHES];
                position_hash.as_ref().unwrap().positions(
                    bank.meta.slots,
                    self.manifest.hashes,
                    self.manifest.segment_length,
                    &mut computed[..self.manifest.hashes],
                );
                &computed[..self.manifest.hashes]
            };
            if bank.meta.fingerprint_first {
                let mut rows = [&[][..]; MAX_HASHES];
                for (row, &slot) in rows.iter_mut().zip(&slots[..self.manifest.hashes]) {
                    *row = &bank.data[slot * row_bytes..slot * row_bytes + lookup_bytes];
                }
                let mut hit = |i: usize| {
                    let group = &bank.meta.groups[i];
                    if group.members.len() == 1 {
                        scores.increment(group.members[0]);
                    } else {
                        let offset = bank.meta.groups.len() * self.manifest.fingerprint_bits
                            + group.offset
                            - i * self.manifest.fingerprint_bits;
                        let mut mask = 0;
                        for &slot in &slots[..self.manifest.hashes] {
                            mask ^= get_bits(
                                &bank.data[slot * row_bytes..(slot + 1) * row_bytes],
                                offset,
                                group.width - self.manifest.fingerprint_bits,
                            );
                        }
                        while mask != 0 {
                            let bit = mask.trailing_zeros() as usize;
                            scores.increment(group.members[bit]);
                            mask &= mask - 1;
                        }
                    }
                };
                if let Some(ranges) = ranges {
                    fingerprint_matches_selected(
                        &rows[..self.manifest.hashes],
                        self.manifest.fingerprint_bits,
                        fp,
                        ranges,
                        &mut hit,
                    );
                } else {
                    fingerprint_matches_rows(
                        &rows[..self.manifest.hashes],
                        self.manifest.fingerprint_bits,
                        fp,
                        bank.meta.groups.len(),
                        &mut hit,
                    );
                }
                continue;
            }
            if let Some(ranges) = ranges {
                // Legacy interleaved layouts can also skip retired groups:
                // decode only selected payloads, without materialising a row.
                for range in ranges {
                    for group in &bank.meta.groups[range.clone()] {
                        let mut value = 0;
                        for &slot in slots {
                            value ^= get_bits(
                                &bank.data[slot * row_bytes..(slot + 1) * row_bytes],
                                group.offset,
                                group.width,
                            );
                        }
                        if value & fp_mask != fp {
                            continue;
                        }
                        if group.members.len() == 1 {
                            scores.increment(group.members[0]);
                        } else {
                            let mut mask = value >> self.manifest.fingerprint_bits;
                            while mask != 0 {
                                scores.increment(group.members[mask.trailing_zeros() as usize]);
                                mask &= mask - 1;
                            }
                        }
                    }
                }
                continue;
            }
            // Legacy interleaved fingerprints/origins still need a complete row.
            let row = &mut scratch[..lookup_bytes];
            row.copy_from_slice(&bank.data[slots[0] * row_bytes..(slots[0] + 1) * row_bytes]);
            for &slot in &slots[1..self.manifest.hashes] {
                xor_row(row, &bank.data[slot * row_bytes..(slot + 1) * row_bytes]);
            }
            for group in &bank.meta.groups {
                let value = get_bits(row, group.offset, group.width);
                if value & fp_mask != fp {
                    continue;
                }
                if group.members.len() == 1 {
                    scores.increment(group.members[0]);
                } else {
                    let mut mask = value >> self.manifest.fingerprint_bits;
                    while mask != 0 {
                        let bit = mask.trailing_zeros() as usize;
                        scores.increment(group.members[bit]);
                        mask &= mask - 1;
                    }
                }
            }
        }
        for (bytes, &origin) in scratch[self.row_scratch_bytes..snapshots_end]
            .chunks_exact(8)
            .zip(recovered)
        {
            let before = u64::from_le_bytes(bytes.try_into().unwrap());
            if scores.get(origin as usize) == before {
                scores.increment(origin as usize);
            }
        }
    }

    /// Small per-worker blocks expose independent row accesses before lookup.
    /// Geometry is prepared once per key, not recomputed for prefetch and XOR.
    fn lookup_batch<S: ScoreSink + ?Sized>(
        &self,
        keys: &[u64],
        scores: &mut S,
        scratch: &mut Vec<u8>,
        positions: &mut Vec<usize>,
    ) {
        const BLOCK: usize = 32;
        const AHEAD: usize = 4;
        if keys.len() < AHEAD {
            for &key in keys {
                self.lookup_into(key, scores, scratch);
            }
            return;
        }
        let hashes = self.manifest.hashes;
        let stride = self.banks.len() * hashes;
        for block in keys.chunks(BLOCK) {
            positions.resize(block.len() * stride, 0);
            for (&key, slots) in block.iter().zip(positions.chunks_exact_mut(stride)) {
                let hash = PositionHash::new(key, self.manifest.seed);
                for (id, bank) in self.selected_banks(key) {
                    hash.positions(
                        bank.meta.slots,
                        hashes,
                        self.manifest.segment_length,
                        &mut slots[id * hashes..(id + 1) * hashes],
                    );
                }
            }
            for i in 0..block.len().min(AHEAD) {
                self.prefetch_lookup(block[i], &positions[i * stride..(i + 1) * stride]);
            }
            for (i, &key) in block.iter().enumerate() {
                if i + AHEAD < block.len() {
                    let next = i + AHEAD;
                    self.prefetch_lookup(
                        block[next],
                        &positions[next * stride..(next + 1) * stride],
                    );
                }
                self.lookup_prepared(
                    key,
                    scores,
                    scratch,
                    Some(&positions[i * stride..(i + 1) * stride]),
                );
            }
        }
    }

    #[inline]
    fn prefetch_record(
        &self,
        key: u64,
        positions: &[usize],
        selected: Option<&[Vec<Range<usize>>]>,
    ) {
        #[cfg(target_arch = "x86_64")]
        {
            if let Some(map) = &self.abandoned {
                map.prefetch(abandoned::hash(key, self.manifest.seed));
            }
            for (id, bank) in self.selected_banks(key) {
                let slots = &positions[id * self.manifest.hashes..(id + 1) * self.manifest.hashes];
                let first = match selected {
                    Some(ranges) => match ranges[id].first() {
                        Some(range) => range.start,
                        None => continue,
                    },
                    None => 0,
                };
                let row_bytes = bank.meta.row_bits.div_ceil(8);
                let offset = if bank.meta.fingerprint_first {
                    first * self.manifest.fingerprint_bits / 8
                } else {
                    0
                };
                // Hint only the start of each independent row stream; hardware
                // handles its sequential continuation. Never prefetch every line.
                for &slot in slots {
                    unsafe {
                        std::arch::x86_64::_mm_prefetch(
                            bank.data.as_ptr().add(slot * row_bytes + offset).cast(),
                            std::arch::x86_64::_MM_HINT_T0,
                        );
                    }
                }
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        let _ = (key, positions, selected);
    }

    #[inline]
    fn prefetch_lookup(&self, key: u64, positions: &[usize]) {
        #[cfg(target_arch = "x86_64")]
        {
            use std::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
            if let Some(map) = &self.abandoned {
                map.prefetch(abandoned::hash(key, self.manifest.seed));
            }
            for (id, bank) in self.selected_banks(key) {
                let slots = &positions[id * self.manifest.hashes..(id + 1) * self.manifest.hashes];
                let row_bytes = bank.meta.row_bits.div_ceil(8);
                let bytes = if bank.meta.fingerprint_first {
                    (bank.meta.groups.len() * self.manifest.fingerprint_bits).div_ceil(8)
                } else {
                    row_bytes
                };
                for &slot in slots {
                    for byte in (0..bytes).step_by(64) {
                        // Every prefetched address is inside a validated row.
                        unsafe {
                            _mm_prefetch(
                                bank.data.as_ptr().add(slot * row_bytes + byte).cast(),
                                _MM_HINT_T0,
                            );
                        }
                    }
                }
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        let _ = (key, positions);
    }
}

/// The dense aggregate path and sparse record path share the exact same lookup,
/// including additive abandoned-key handling. Monomorphisation avoids dynamic
/// dispatch or touched-list bookkeeping in the aggregate query hot path.
trait ScoreSink {
    fn get(&self, origin: usize) -> u64;
    fn increment(&mut self, origin: usize);
    fn add(&mut self, origin: usize, count: u64);
}

impl ScoreSink for [u64] {
    #[inline]
    fn get(&self, origin: usize) -> u64 {
        self[origin]
    }

    #[inline]
    fn increment(&mut self, origin: usize) {
        self[origin] += 1;
    }

    #[inline]
    fn add(&mut self, origin: usize, count: u64) {
        self[origin] += count;
    }
}

struct WeightedScores<'a, S: ScoreSink + ?Sized> {
    scores: &'a mut S,
    weight: u64,
}

impl<S: ScoreSink + ?Sized> ScoreSink for WeightedScores<'_, S> {
    #[inline]
    fn get(&self, origin: usize) -> u64 {
        self.scores.get(origin)
    }
    #[inline]
    fn increment(&mut self, origin: usize) {
        self.scores.add(origin, self.weight);
    }
    #[inline]
    fn add(&mut self, origin: usize, count: u64) {
        self.scores.add(origin, count * self.weight);
    }
}

struct RecordScores {
    values: Vec<u64>,
    touched: Vec<usize>,
}

impl RecordScores {
    fn write_record(
        &mut self,
        name: &[u8],
        total: u64,
        min_ratio: f64,
        format: RecordFormat,
        colors: &[Color],
        out: &mut impl Write,
    ) -> Result<()> {
        let all = if min_ratio == 0.0 {
            self.values.len()
        } else {
            0
        };
        let sparse = if min_ratio > 0.0 && total != 0 {
            self.touched.sort_unstable();
            self.touched.len()
        } else {
            0
        };
        let hits = self.values[..all]
            .iter()
            .copied()
            .enumerate()
            .chain(
                self.touched[..sparse]
                    .iter()
                    .map(|&id| (id, self.values[id])),
            )
            .filter(|&(_, count)| {
                let ratio = if total == 0 {
                    0.0
                } else {
                    count as f64 / total as f64
                };
                ratio >= min_ratio
            });
        format.write_record(out, name, total, hits, colors)
    }
    fn new(colors: usize) -> Self {
        Self {
            values: vec![0; colors],
            touched: Vec::new(),
        }
    }

    fn clear(&mut self) {
        for origin in self.touched.drain(..) {
            self.values[origin] = 0;
        }
    }

    #[cfg(test)]
    fn write(&mut self, total: u64, min_ratio: f64, out: &mut impl Write) -> Result<()> {
        write!(out, "{total}\t")?;
        let mut separator = "";
        let mut emit = |origin: usize, count: u64| -> Result<()> {
            // Retain the original floating-point threshold semantics exactly.
            let ratio = if total == 0 {
                0.0
            } else {
                count as f64 / total as f64
            };
            if ratio >= min_ratio {
                write!(out, "{separator}{origin}:{count}")?;
                separator = ",";
            }
            Ok(())
        };
        if min_ratio == 0.0 {
            // At threshold zero, zero-count documents must still be reported.
            for (origin, &count) in self.values.iter().enumerate() {
                emit(origin, count)?;
            }
        } else if total != 0 {
            self.touched.sort_unstable();
            for &origin in &self.touched {
                emit(origin, self.values[origin])?;
            }
        }
        writeln!(out)?;
        Ok(())
    }
}

impl ScoreSink for RecordScores {
    #[inline]
    fn get(&self, origin: usize) -> u64 {
        self.values[origin]
    }

    #[inline]
    fn increment(&mut self, origin: usize) {
        self.add(origin, 1);
    }

    #[inline]
    fn add(&mut self, origin: usize, count: u64) {
        if self.values[origin] == 0 {
            self.touched.push(origin);
        }
        self.values[origin] += count;
    }
}

fn fingerprint(key: u64, seed: u64, bits: usize) -> u64 {
    splitmix64(key ^ seed ^ FINGERPRINT_SEED) & ((1u64 << bits) - 1)
}

/// XOR and compare directly in registers, without materialising an intermediate
/// fingerprint row or repeatedly reading/writing it for every selected slot.
fn fingerprint_matches_selected(
    rows: &[&[u8]],
    bits: usize,
    fp: u64,
    ranges: &[Range<usize>],
    hit: &mut impl FnMut(usize),
) {
    let mut selected = [&[][..]; MAX_HASHES];
    for range in ranges {
        debug_assert!(range.start * bits % 8 == 0);
        let start = range.start * bits / 8;
        let end = (range.end * bits).div_ceil(8);
        for (out, row) in selected.iter_mut().zip(rows) {
            *out = &row[start..end];
        }
        fingerprint_matches_rows(&selected[..rows.len()], bits, fp, range.len(), |group| {
            hit(range.start + group)
        });
    }
}

fn fingerprint_matches_rows(
    rows: &[&[u8]],
    bits: usize,
    fp: u64,
    count: usize,
    mut hit: impl FnMut(usize),
) {
    debug_assert!(!rows.is_empty());
    debug_assert!(rows
        .iter()
        .all(|row| row.len() == (count * bits).div_ceil(8)));
    if matches!(bits, 2 | 4) {
        crate::packed_fingerprints::matches_rows(rows, bits, fp, count, hit);
        return;
    }
    if bits == 1 {
        for base in (0..count).step_by(64) {
            let byte = base / 8;
            let len = (count - base).min(64).div_ceil(8);
            let mut value = 0u64;
            for row in rows {
                let mut word = [0; 8];
                word[..len].copy_from_slice(&row[byte..byte + len]);
                value ^= u64::from_le_bytes(word);
            }
            let remaining = (count - base).min(64);
            let valid = u64::MAX >> (64 - remaining);
            let mut mask = (if fp == 1 { value } else { !value }) & valid;
            while mask != 0 {
                hit(base + mask.trailing_zeros() as usize);
                mask &= mask - 1;
            }
        }
        return;
    }
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx2") {
        // Check once per bank, rather than dispatching for every selected row.
        unsafe { fingerprint_matches_rows_avx2(rows, bits, fp, &mut hit) };
        return;
    }
    fingerprint_matches_rows_scalar(rows, bits, fp, 0, &mut hit);
}

fn fingerprint_matches_rows_scalar(
    rows: &[&[u8]],
    bits: usize,
    fp: u64,
    byte_start: usize,
    hit: &mut impl FnMut(usize),
) {
    let width = bits / 8;
    for byte in (byte_start..rows[0].len()).step_by(width) {
        let mut value = 0u64;
        for row in rows {
            value ^= if bits == 8 {
                row[byte] as u64
            } else {
                u16::from_le_bytes([row[byte], row[byte + 1]]) as u64
            };
        }
        if value == fp {
            hit(byte / width);
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn fingerprint_matches_rows_avx2(
    rows: &[&[u8]],
    bits: usize,
    fp: u64,
    hit: &mut impl FnMut(usize),
) {
    use std::arch::x86_64::*;
    let target = if bits == 8 {
        _mm256_set1_epi8(fp as i8)
    } else {
        _mm256_set1_epi16(fp as i16)
    };
    let end = rows[0].len() / 32 * 32;
    for byte in (0..end).step_by(32) {
        let mut value = unsafe { _mm256_loadu_si256(rows[0].as_ptr().add(byte).cast()) };
        for row in &rows[1..] {
            let next = unsafe { _mm256_loadu_si256(row.as_ptr().add(byte).cast()) };
            value = _mm256_xor_si256(value, next);
        }
        let equal = if bits == 8 {
            _mm256_cmpeq_epi8(value, target)
        } else {
            _mm256_cmpeq_epi16(value, target)
        };
        let mut mask = _mm256_movemask_epi8(equal) as u32;
        if bits == 16 {
            mask &= 0x5555_5555;
        }
        while mask != 0 {
            hit((byte + mask.trailing_zeros() as usize) / (bits / 8));
            mask &= mask - 1;
        }
    }
    fingerprint_matches_rows_scalar(rows, bits, fp, end, hit);
}

#[cfg(test)]
fn fingerprint_matches(row: &[u8], bits: usize, fp: u64, count: usize, mut hit: impl FnMut(usize)) {
    if matches!(bits, 2 | 4) {
        crate::packed_fingerprints::matches_rows(&[row], bits, fp, count, hit);
        return;
    }
    if bits == 1 {
        // Each XORed bit is already a comparison result. Screen 64 groups at
        // once, and never treat padding or adjacent origin bits as groups.
        for (chunk, bytes) in row.chunks(8).enumerate() {
            let mut word = [0u8; 8];
            word[..bytes.len()].copy_from_slice(bytes);
            let value = u64::from_le_bytes(word);
            let remaining = count.saturating_sub(chunk * 64).min(64);
            let valid = if remaining == 64 {
                u64::MAX
            } else {
                (1u64 << remaining) - 1
            };
            let mut mask = (if fp == 1 { value } else { !value }) & valid;
            while mask != 0 {
                hit(chunk * 64 + mask.trailing_zeros() as usize);
                mask &= mask - 1;
            }
        }
        return;
    }
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx2") {
        // Runtime dispatch keeps indexes portable to machines without AVX2.
        unsafe {
            fingerprint_matches_avx2(row, bits, fp, &mut hit);
        }
        return;
    }
    fingerprint_matches_scalar(row, bits, fp, 0, &mut hit);
}

#[cfg(test)]
fn fingerprint_matches_scalar(
    row: &[u8],
    bits: usize,
    fp: u64,
    base: usize,
    hit: &mut impl FnMut(usize),
) {
    for (i, value) in row.chunks_exact(bits / 8).enumerate() {
        let value = if bits == 8 {
            value[0] as u64
        } else {
            u16::from_le_bytes(value.try_into().unwrap()) as u64
        };
        if value == fp {
            hit(base + i);
        }
    }
}

#[cfg(all(test, target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn fingerprint_matches_avx2(row: &[u8], bits: usize, fp: u64, hit: &mut impl FnMut(usize)) {
    use std::arch::x86_64::*;
    let target = if bits == 8 {
        _mm256_set1_epi8(fp as i8)
    } else {
        _mm256_set1_epi16(fp as i16)
    };
    let mut chunks = row.chunks_exact(32);
    for (chunk, bytes) in chunks.by_ref().enumerate() {
        let value = unsafe { _mm256_loadu_si256(bytes.as_ptr().cast()) };
        let equal = if bits == 8 {
            _mm256_cmpeq_epi8(value, target)
        } else {
            _mm256_cmpeq_epi16(value, target)
        };
        let mut mask = _mm256_movemask_epi8(equal) as u32;
        if bits == 16 {
            mask &= 0x5555_5555;
        }
        while mask != 0 {
            hit((chunk * 32 + mask.trailing_zeros() as usize) / (bits / 8));
            mask &= mask - 1;
        }
    }
    let tail = chunks.remainder();
    fingerprint_matches_scalar(tail, bits, fp, (row.len() - tail.len()) / (bits / 8), hit);
}
/// Same segmented mapping as production ZOR. Total cell count is a power of
/// two; the number of overlapping segment starts is not necessarily one.
#[cfg(test)]
fn positions(key: u64, cells: usize, hashes: usize, seed: u64, out: &mut [usize]) {
    PositionHash::new(key, seed).positions(cells, hashes, None, out);
}

struct PositionHash {
    first: u64,
    step: u64,
}

impl PositionHash {
    #[inline]
    fn new(key: u64, seed: u64) -> Self {
        let first = crate::mixsplit(key, seed);
        Self {
            first,
            step: splitmix64(first),
        }
    }

    #[inline]
    fn positions(&self, cells: usize, hashes: usize, fixed: Option<usize>, out: &mut [usize]) {
        let length = fixed.unwrap_or_else(|| {
            let n = (cells / hashes).clamp(1, 4096);
            1usize << n.ilog2()
        });
        let count = cells / length - hashes + 1;
        let base = (((self.first as u128 * count as u128) >> 64) as usize) * length;
        let mut hash = self.first;
        for (i, slot) in out.iter_mut().enumerate() {
            *slot = base + i * length + ((hash ^ (hash >> 33)) as usize & (length - 1));
            hash = hash.wrapping_add(self.step);
        }
    }
}

fn xor_row(acc: &mut [u8], row: &[u8]) {
    let mut a = acc.chunks_exact_mut(8);
    let mut b = row.chunks_exact(8);
    for (x, y) in a.by_ref().zip(b.by_ref()) {
        let value =
            u64::from_le_bytes(x.try_into().unwrap()) ^ u64::from_le_bytes(y.try_into().unwrap());
        x.copy_from_slice(&value.to_le_bytes());
    }
    for (x, &y) in a.into_remainder().iter_mut().zip(b.remainder()) {
        *x ^= y;
    }
}
fn put_bits(row: &mut [u8], offset: usize, width: usize, value: u64) {
    let byte = offset / 8;
    let shift = offset % 8;
    let value = if width == 64 {
        value
    } else {
        value & ((1u64 << width) - 1)
    };
    if row.len() - byte >= 8 {
        let previous = u64::from_le_bytes(row[byte..byte + 8].try_into().unwrap());
        row[byte..byte + 8].copy_from_slice(&(previous | (value << shift)).to_le_bytes());
        if shift + width > 64 {
            row[byte + 8] |= (value >> (64 - shift)) as u8;
        }
        return;
    }
    let mut left = width;
    let mut bit = offset;
    let mut consumed = 0;
    while left > 0 {
        let chunk = left.min(8 - bit % 8);
        let mask = (1u64 << chunk) - 1;
        row[bit / 8] |= (((value >> consumed) & mask) as u8) << (bit % 8);
        left -= chunk;
        bit += chunk;
        consumed += chunk;
    }
}
fn get_bits(row: &[u8], offset: usize, width: usize) -> u64 {
    let byte = offset / 8;
    let shift = offset % 8;
    let count = (shift + width).div_ceil(8);
    let mut buffer = [0; 8];
    let first = count.min(8);
    buffer[..first].copy_from_slice(&row[byte..byte + first]);
    let mut value = u64::from_le_bytes(buffer) >> shift;
    if count == 9 {
        value |= (row[byte + 8] as u64) << (64 - shift);
    }
    if width < 64 {
        value &= (1u64 << width) - 1;
    }
    value
}

struct QueryWorker {
    scores: Vec<u64>,
    scratch: Vec<u8>,
    positions: Vec<usize>,
}
fn flush(index: &Index, batch: &mut Vec<u64>, workers: &mut Vec<QueryWorker>) {
    if batch.is_empty() {
        return;
    }
    // Global batch deduplication happens before partitioning, so a repeated key
    // is looked up once even when its occurrences would cross worker boundaries.
    batch.sort_unstable();
    let mut runs: Vec<(u64, u64)> = Vec::new();
    for &key in batch.iter() {
        if let Some(last) = runs.last_mut().filter(|last| last.0 == key) {
            last.1 += 1;
        } else {
            runs.push((key, 1));
        }
    }
    let threads = rayon::current_num_threads().min(runs.len());
    while workers.len() < threads {
        workers.push(QueryWorker {
            scores: vec![0; index.manifest.colors.len()],
            scratch: Vec::new(),
            positions: Vec::new(),
        });
    }
    if runs.len() != batch.len() {
        workers[..threads]
            .par_iter_mut()
            .zip(runs.par_chunks(runs.len().div_ceil(threads)))
            .for_each(|(worker, runs)| {
                for &(key, weight) in runs {
                    index.lookup_into(
                        key,
                        &mut WeightedScores {
                            scores: worker.scores.as_mut_slice(),
                            weight,
                        },
                        &mut worker.scratch,
                    );
                }
            });
    } else {
        workers[..threads]
            .par_iter_mut()
            .zip(batch.par_chunks(batch.len().div_ceil(threads)))
            .for_each(|(worker, keys)| {
                index.lookup_batch(
                    keys,
                    worker.scores.as_mut_slice(),
                    &mut worker.scratch,
                    &mut worker.positions,
                );
            });
    }
    batch.clear();
}

pub fn query_command(directory: &Path, query: &Path, min_ratio: f64) -> Result<()> {
    query_command_with_recovery(directory, query, min_ratio, true)
}

/// Record-level reporting for ordinary feature modes. Uses the native lookup
/// with named or compact TSV output, retaining counters/scratch per worker
/// and formats only touched origins (unless threshold zero requires all).
pub fn query_records_command(
    directory: &Path,
    query: &Path,
    min_ratio: f64,
    recover: bool,
    batch_size: usize,
    format: RecordFormat,
) -> Result<()> {
    query_records_command_with_memory(
        directory,
        query,
        min_ratio,
        recover,
        batch_size,
        QueryMemory::default(),
        format,
    )
}

pub fn query_records_command_with_memory(
    directory: &Path,
    query: &Path,
    min_ratio: f64,
    recover: bool,
    batch_size: usize,
    memory: QueryMemory,
    format: RecordFormat,
) -> Result<()> {
    let load_start = Instant::now();
    ensure!(
        (0.0..=1.0).contains(&min_ratio),
        "--min-ratio must be in [0, 1]"
    );
    ensure!(batch_size > 0, "record batch size must be positive");
    if memory == QueryMemory::Subindex {
        return query_stream::query_records(
            directory, query, min_ratio, recover, batch_size, format,
        );
    }
    let index = Index::load_with_recovery(directory, recover)
        .with_context(|| format!("loading {}", directory.display()))?;
    let load_seconds = load_start.elapsed().as_secs_f64();
    let features = index.manifest.features()?;
    ensure!(
        features.mode != IndexMode::Findere,
        "--records does not support findere indexes"
    );
    let query_start = Instant::now();
    super::records::query_file(
        query,
        batch_size,
        format.header(),
        || record_query::Worker::new(&index),
        |worker, records, out| {
            worker.query_records(&index, features, min_ratio, format, records, out)
        },
    )?;
    eprintln!(
        "query_timings\tload_s={:.6}\tquery_s={:.6}",
        load_seconds,
        query_start.elapsed().as_secs_f64()
    );
    Ok(())
}

pub fn query_command_with_recovery(
    directory: &Path,
    query: &Path,
    min_ratio: f64,
    recover: bool,
) -> Result<()> {
    query_command_with_memory(directory, query, min_ratio, recover, QueryMemory::default())
}

pub fn query_command_with_memory(
    directory: &Path,
    query: &Path,
    min_ratio: f64,
    recover: bool,
    memory: QueryMemory,
) -> Result<()> {
    ensure!(
        (0.0..=1.0).contains(&min_ratio),
        "--min-ratio must be in [0, 1]"
    );
    if memory == QueryMemory::Subindex {
        return query_stream::query_aggregate(directory, query, min_ratio, recover);
    }
    let index = Index::load_with_recovery(directory, recover)
        .with_context(|| format!("loading {}", directory.display()))?;
    let start = Instant::now();
    let (scores, total, probes) = query_index(&index, query)?;
    print_aggregate(
        &index.manifest,
        query,
        &scores,
        total,
        probes,
        probes
            * index
                .manifest
                .banks
                .iter()
                .map(BankMeta::source_group)
                .collect::<std::collections::BTreeSet<_>>()
                .len() as u64,
        min_ratio,
        start.elapsed().as_secs_f64(),
    );
    Ok(())
}

fn query_index(index: &Index, query: &Path) -> Result<(Vec<u64>, u64, u64)> {
    let features = index.manifest.features()?;
    let (scores, total, probes) = if features.mode == IndexMode::Findere {
        let result = crate::findere::query_file_with_homopolymer_compression(
            query,
            index.manifest.k,
            features.findere_z,
            index.manifest.colors.len(),
            features.homopolymer_compression,
            Vec::<u8>::new,
            |scratch, key, hits| {
                index.lookup(key, hits, scratch);
                true
            },
        )?;
        (result.scores, result.total_kmers, result.smer_queries)
    } else if features.randstrobemers.is_some() {
        query_strobe_file(&index, query, features)?
    } else {
        let mut workers = Vec::new();
        let mut batch = Vec::with_capacity(262_144);
        let mut total = 0;
        parse_features_from_file(
            query,
            index.manifest.k,
            index.manifest.seed,
            features,
            |key| {
                total += 1;
                batch.push(key);
                if batch.len() == 262_144 {
                    flush(&index, &mut batch, &mut workers);
                }
            },
        )?;
        flush(&index, &mut batch, &mut workers);
        let mut scores = vec![0; index.manifest.colors.len()];
        for worker in workers {
            for (score, v) in scores.iter_mut().zip(worker.scores) {
                *score += v;
            }
        }
        (scores, total, total)
    };
    Ok((scores, total, probes))
}

#[allow(clippy::too_many_arguments)]
fn print_aggregate(
    manifest: &Manifest,
    query: &Path,
    scores: &[u64],
    total: u64,
    probes: u64,
    bank_probes: u64,
    min_ratio: f64,
    seconds: f64,
) {
    let mode = IndexMode::from_u8(manifest.index_mode).expect("validated index mode");
    eprintln!("#query\t{}\tk={}\tindex_mode={}\ttotal_features={}\tstack_probes={}\tbank_probes={}\telapsed_s={:.6}",
        query.display(), manifest.k, mode.label(), total, probes, bank_probes, seconds);
    println!("color\tmatches\tratio\tunique_kmers");
    for (color, &matches) in manifest.colors.iter().zip(scores) {
        let ratio = if total == 0 {
            0.0
        } else {
            matches as f64 / total as f64
        };
        if ratio >= min_ratio {
            println!(
                "{}\t{}\t{:.8}\t{}",
                color.name, matches, ratio, color.unique_kmers
            );
        }
    }
}

fn query_strobe_file(
    index: &Index,
    path: &Path,
    features: FeatureConfig,
) -> Result<(Vec<u64>, u64, u64)> {
    let n = rayon::current_num_threads();
    let mut workers: Vec<_> = (0..n)
        .map(|_| {
            (
                record_query::Worker::new(index),
                vec![0u64; index.manifest.colors.len()],
                0u64,
            )
        })
        .collect();
    let mut batch = Vec::new();
    let mut bases = 0;
    let mut flush = |batch: &mut Vec<Vec<u8>>| {
        if batch.is_empty() {
            return;
        }
        let chunk = batch.len().div_ceil(n);
        workers
            .par_iter_mut()
            .zip(batch.par_chunks(chunk))
            .for_each(|((worker, sums, total), reads)| {
                for seq in reads {
                    *total += worker.query_sequence(index, features, 0.0, seq);
                    for &origin in &worker.scores.touched {
                        sums[origin] += worker.scores.values[origin];
                    }
                }
            });
        batch.clear();
    };
    crate::parse_sequences_from_file(path, |seq| {
        bases += seq.len();
        batch.push(seq.to_vec());
        if bases >= 1024 * 1024 || batch.len() == 4096 {
            flush(&mut batch);
            bases = 0;
        }
    })?;
    flush(&mut batch);
    let mut scores = vec![0; index.manifest.colors.len()];
    let mut total = 0;
    for (_, sums, count) in workers {
        total += count;
        for (score, value) in scores.iter_mut().zip(sums) {
            *score += value;
        }
    }
    Ok((scores, total, total * 2))
}

pub fn info_command(directory: &Path) -> Result<()> {
    let mut index = Index::load_mode(directory, false, None, true)?;
    if index.manifest.version == 1 {
        index = Index::load(directory)?;
    } else {
        // Validate payloads one at a time; info must not defeat the query cap.
        for meta in &index.manifest.banks {
            drop(load_modern_bank(directory, meta)?);
        }
        if index.manifest.abandoned_recovery {
            index.abandoned = Some(load_abandoned(directory, &index.manifest)?);
        }
    }
    let m = &index.manifest;
    println!(
        "collection_status\t{}\nindexed_documents\t{}",
        m.collection_status(),
        m.colors.len()
    );
    let copied = m.construction_cell_storage == "copied-compressed-banks";
    if copied {
        println!(
            "components\t{}",
            m.banks
                .iter()
                .map(|b| b.component)
                .collect::<std::collections::BTreeSet<_>>()
                .len()
        );
    }
    let features = m.features()?;
    let groups: usize = m.banks.iter().map(|b| b.groups.len()).sum();
    let logical_bits: usize = m.banks.iter().map(|b| b.slots * b.row_bits).sum();
    let stack_bytes: usize = m
        .banks
        .iter()
        .map(|b| b.slots * b.row_bits.div_ceil(8))
        .sum();
    let fallback_bytes: usize = if m.version == 1 {
        m.banks
            .iter()
            .flat_map(|b| &b.groups)
            .map(|g| g.fallback_keys * 16)
            .sum()
    } else {
        m.abandoned_hashes * 8 + m.abandoned_origins * 4
    };
    let metadata_bytes = fs::metadata(directory.join("manifest.json"))?.len() as usize
        + m.banks.len() * 40
        + if m.version == 1 {
            groups * 8
        } else if m.abandoned_file.is_some() {
            8
        } else {
            0
        };
    let disk_bytes: u64 = fs::read_dir(directory)?
        .map(|entry| entry.and_then(|entry| entry.metadata()).map(|m| m.len()))
        .collect::<std::io::Result<Vec<_>>>()?
        .into_iter()
        .sum();
    println!("index\t{}\nk\t{}\nhashes\t{}\nseed\t{}\nfingerprint_bits\t{}\nindex_mode\t{}\ntechnical_layout\tmerge-pow2\ncolors\t{}\ngroups\t{}\nbanks\t{}",
        directory.display(), m.k, m.hashes, m.seed, m.fingerprint_bits, features.mode.label(), m.colors.len(), groups, m.banks.len());
    println!(
        "size_classes\t{}\nmax_subindex_size\t{}\nlargest_subindex_bytes\t{}\nhash_partitions\t{}",
        m.banks
            .iter()
            .map(BankMeta::source_group)
            .collect::<std::collections::BTreeSet<_>>()
            .len(),
        m.max_subindex_size
            .map_or_else(|| "unlimited".into(), |v| v.to_string()),
        m.banks
            .iter()
            .map(|b| b.slots * b.row_bits.div_ceil(8))
            .max()
            .unwrap_or(0),
        m.banks.iter().filter(|b| b.routing.is_some()).count()
    );
    println!(
        "input_documents\t{}\nskipped_documents\t{}\nmemory_guard_split_groups\t{}",
        m.input_document_count.max(m.colors.len()),
        m.skipped_documents.len(),
        m.memory_guard_split_groups
    );
    if features.mode == IndexMode::Findere {
        println!(
            "findere_z\t{}\nindexed_k\t{}",
            features.findere_z,
            m.k - features.findere_z
        );
    }
    println!(
        "minimizer_size\t{}\nmodimizer_sampling\t{}",
        features.minimizer_size, features.modimizer_sampling
    );
    println!("minimizer_scheme\t{}", features.minimizer_scheme.label());
    println!(
        "homopolymer_compression\t{}",
        features.homopolymer_compression
    );
    crate::anchors::print_settings(features);
    println!(
        "merge_planner\t{}\nplanning_pair_evaluations\t{}\nplanning_possible_pair_evaluations\t{}\nmerge_reverts_enabled\t{}",
        if copied { "preserved" } else { m.merge_planner.label() },
        m.planning_pair_evaluations,
        m.planning_possible_pair_evaluations,
        m.merge_reverts_enabled
    );
    println!("construction_cell_storage\t{}", m.construction_cell_storage);
    if let Some(budget) = m.merge_memory_budget_bytes {
        println!("merge_memory_budget_bytes\t{budget}");
    }
    println!(
        "construction_abandoned_storage\t{}",
        m.construction_abandoned_storage
            .as_deref()
            .unwrap_or("legacy-memory")
    );
    if let Some(length) = m.segment_length {
        println!("segment_length\t{length}");
    }
    println!("overflow_groups\t{}", m.overflow_groups);
    println!(
        "merge_column_order\t{}",
        if copied {
            "preserved"
        } else {
            m.column_order.label()
        }
    );
    if let Some(method) = &m.column_order_method {
        println!("column_order_method\t{method}");
    }
    println!("abandoned_recovery\t{}\nabandoned_hash_bits\t32\nabandoned_keys\t{}\nabandoned_unique_hashes\t{}\nabandoned_dataset_associations\t{}\nabandoned_map_ram_bytes\t{}",
        m.abandoned_recovery,
        m.banks.iter().flat_map(|b| &b.groups).map(|g| g.abandoned_keys.unwrap_or(g.fallback_keys)).sum::<usize>(),
        index.abandoned.as_ref().map_or(0, OriginMap::len),
        index.abandoned.as_ref().map_or(0, OriginMap::origins_len),
        index.abandoned.as_ref().map_or(0, OriginMap::memory_bytes));
    println!("logical_bits\t{}\npadding_bits\t{}\nstack_bytes\t{}\nstack_stored_bytes\t{}\nfallback_bytes\t{}\nstructure_bytes\t{}\nsingleton_logical_bits\t{}\npredicted_logical_bits\t{}\nplanned_groups\t{}\nrejected_groups\t{}\nexact_group_verifications\t{}\nsketch_seconds\t{:.6}\nplan_seconds\t{:.6}\nbuild_seconds\t{:.6}",
        logical_bits, stack_bytes * 8 - logical_bits, stack_bytes, stack_bytes, fallback_bytes, stack_bytes + fallback_bytes,
        m.singleton_logical_bits, m.predicted_logical_bits, m.planned_groups, m.rejected_groups, m.exact_group_verifications,
        m.sketch_seconds, m.plan_seconds, m.build_seconds);
    println!(
        "metadata_bytes\t{}\nuncompressed_index_bytes\t{}\ndisk_bytes\t{}",
        metadata_bytes,
        stack_bytes + fallback_bytes + metadata_bytes,
        disk_bytes
    );
    println!(
        "min_subindex_documents\t{}\npromoted_groups\t{}",
        m.min_subindex_documents, m.promoted_groups
    );
    if let Some(limit) = m.max_subindexes {
        println!("max_subindexes\t{limit}");
    }
    println!("bank\tslots\tgroups\trow_bits\tpadding_bits_per_row\tfallback_keys\tdocuments");
    for b in &m.banks {
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            b.file,
            b.slots,
            b.groups.len(),
            b.row_bits,
            b.row_bits.div_ceil(8) * 8 - b.row_bits,
            b.groups.iter().map(|g| g.fallback_keys).sum::<usize>(),
            b.groups.iter().map(|g| g.members.len()).sum::<usize>()
        );
    }
    println!("group\tslots\twidth\tunion_cardinality\tfallback_keys\tcolors");
    let mut id = 0;
    for b in &m.banks {
        for g in &b.groups {
            println!(
                "{}\t{}\t{}\t{}\t{}\t{}",
                id,
                b.slots,
                g.width,
                g.union_cardinality,
                g.fallback_keys,
                g.members
                    .iter()
                    .map(|&c| m.colors[c].name.as_str())
                    .collect::<Vec<_>>()
                    .join(",")
            );
            id += 1;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn guarded_admission_never_runs_an_oversized_job_alone() {
        let called = std::sync::atomic::AtomicBool::new(false);
        let error = pull_jobs_budgeted_with_policy(
            &[1],
            1,
            &[101],
            100,
            false,
            false,
            None,
            |_, _| {
                called.store(true, Ordering::Relaxed);
                Ok(())
            },
            |_, _| Ok(()),
        )
        .unwrap_err();
        assert!(!called.load(Ordering::Relaxed));
        assert!(error.to_string().contains("refuses an uncapped job"));
    }

    #[test]
    fn oversized_merge_is_split_without_skipping_admissible_documents() {
        let mut config = scratch_free_config(false);
        config.hashes = 8;
        let colors: Vec<_> = (0..2)
            .map(|id| Color {
                name: id.to_string(),
                path: id.to_string(),
                occurrences: 100_000,
                unique_kmers: 100_000,
            })
            .collect();
        let group = planner::Group {
            members: vec![0, 1],
            cardinality: 200_000.0,
            slots: config.slots(200_000.0).unwrap(),
            width: 10,
            sketch: Sketch::new(&[], 64, 69),
            ordering_key: None,
            feature_count_upper_bound: 200_000,
            input_peak_bytes: 16 << 20,
        };
        let (groups, count) = guard_groups(
            vec![group],
            &colors,
            &[8 << 20; 2],
            &config,
            Guard {
                enabled: true,
                budget: 32 << 20,
            },
        )
        .unwrap();
        assert_eq!(count, 1);
        assert_eq!(
            groups.iter().map(|g| g.members.clone()).collect::<Vec<_>>(),
            vec![vec![0], vec![1]]
        );
        assert!(groups
            .iter()
            .all(|g| group_workspace(g, &colors, &[8 << 20; 2], &config, 0).unwrap() <= 32 << 20));
    }

    #[test]
    fn workspace_admission_limits_live_and_pending_results_without_deadlock() {
        use std::sync::atomic::AtomicU64;
        let items = [4u64, 6, 20, 3, 7, 1, 9, 2];
        for threads in [1, 4] {
            let live = AtomicU64::new(0);
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            struct Reservation<'a> {
                bytes: u64,
                live: &'a AtomicU64,
            }
            impl Drop for Reservation<'_> {
                fn drop(&mut self) {
                    self.live.fetch_sub(self.bytes, Ordering::SeqCst);
                }
            }
            let mut order = Vec::new();
            pool.install(|| {
                pull_jobs_budgeted(
                    &items,
                    usize::MAX,
                    &items,
                    10,
                    |id, &bytes| {
                        let previous = live.fetch_add(bytes, Ordering::SeqCst);
                        assert!(if bytes > 10 {
                            previous == 0
                        } else {
                            previous + bytes <= 10
                        });
                        if id % 2 == 0 {
                            std::thread::yield_now();
                        }
                        Ok(Reservation { bytes, live: &live })
                    },
                    |id, value| {
                        order.push(id);
                        drop(value);
                        Ok(())
                    },
                )
            })
            .unwrap();
            assert_eq!(order, (0..items.len()).collect::<Vec<_>>());
            assert_eq!(live.load(Ordering::SeqCst), 0);
        }
    }

    #[test]
    fn workspace_admission_releases_results_when_worker_or_consumer_fails() {
        for fail_worker in [true, false] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(4)
                .build()
                .unwrap();
            let result = pool.install(|| {
                pull_jobs_budgeted(
                    &[0, 1, 2, 3],
                    4,
                    &[4, 4, 11, 4],
                    8,
                    |id, _| {
                        ensure!(id != 2 || !fail_worker, "worker failed");
                        Ok(id)
                    },
                    |id, _| {
                        ensure!(id != 2 || fail_worker, "consumer failed");
                        Ok(())
                    },
                )
            });
            assert!(result.is_err());
        }
    }

    #[test]
    fn task_pipeline_preserves_input_order_with_one_or_multiple_workers() {
        let items: Vec<_> = (0..65usize).collect();
        for threads in [1, 4] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            for ahead in [1, 16, PLANNER_LOOKAHEAD, usize::MAX] {
                let mut consumed = Vec::new();
                pool.install(|| {
                    pull_jobs(
                        &items,
                        ahead,
                        |id, value| {
                            if id % 7 == 0 {
                                std::thread::yield_now();
                            }
                            Ok(*value)
                        },
                        |id, value| {
                            assert_eq!(id, consumed.len());
                            consumed.push(value);
                            Ok(())
                        },
                    )
                    .unwrap();
                });
                assert_eq!(consumed, items);
            }
        }
    }

    #[test]
    fn task_pipeline_propagates_worker_consumer_and_panic_errors() {
        let items: Vec<_> = (0..16usize).collect();
        for threads in [1, 4] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            for failure in ["worker", "consumer", "panic"] {
                let error = pool.install(|| {
                    pull_jobs(
                        &items,
                        4,
                        |id, value| {
                            if id == 3 {
                                match failure {
                                    "worker" => anyhow::bail!("worker failure"),
                                    "panic" => panic!("test worker panic"),
                                    _ => {}
                                }
                            }
                            Ok(*value)
                        },
                        |id, _| {
                            ensure!(id != 3 || failure != "consumer", "consumer failure");
                            Ok(())
                        },
                    )
                    .unwrap_err()
                });
                assert!(error.to_string().contains(if failure == "panic" {
                    "task 3 panicked"
                } else if failure == "worker" {
                    "worker failure"
                } else {
                    "consumer failure"
                }));
            }
        }
    }

    #[test]
    fn in_place_union_matches_exact_origin_map() {
        let mut keys = Vec::new();
        let mut masks = Vec::new();
        let mut expected = BTreeMap::<u64, u64>::new();
        for bit in 0..63 {
            let mut incoming: Vec<_> = (0..bit * 3).map(|i| splitmix64(i) % 500).collect();
            incoming.extend([0, u64::MAX]);
            incoming.sort_unstable();
            incoming.dedup();
            let mask = 1u64 << bit;
            for &key in &incoming {
                *expected.entry(key).or_default() |= mask;
            }
            extend_union(&mut keys, &mut masks, &incoming, mask).unwrap();
            assert_eq!(keys, expected.keys().copied().collect::<Vec<_>>());
            assert_eq!(masks, expected.values().copied().collect::<Vec<_>>());
        }
        extend_union(&mut keys, &mut masks, &[], 1).unwrap();
        assert_eq!(keys, expected.keys().copied().collect::<Vec<_>>());
    }

    #[test]
    fn packed_cells_preserve_reads_for_all_payload_widths() {
        let stage = temporary_stage();
        for width in 1..=64 {
            for slots in [4, 17, 65] {
                let mask = if width == 64 {
                    u64::MAX
                } else {
                    (1u64 << width) - 1
                };
                let cells: Vec<_> = (0..slots).map(|i| splitmix64(i as u64) & mask).collect();
                let mut group = test_group(&stage, width * 100 + slots, width, &cells);
                group.cells.remove().unwrap();
                group.cells = CompletedCells::packed(&cells, width).unwrap();
                if let CompletedCells::Packed { data, .. } = &group.cells {
                    assert_eq!(data.len(), (width * slots).div_ceil(8));
                }
                let staged = stage.0.join(format!("{width}-{slots}.bits"));
                for cells_in_file in [false, true] {
                    if cells_in_file {
                        group.cells =
                            CompletedCells::packed_disk(&cells, width, staged.clone()).unwrap();
                        assert_eq!(
                            fs::metadata(&staged).unwrap().len(),
                            (width * slots).div_ceil(8) as u64
                        );
                    }
                    let mut block = vec![0; (slots - 1) * 8];
                    group.cells.read_block(1, &mut block, None).unwrap();
                    let actual: Vec<_> = block
                        .chunks_exact(8)
                        .map(|b| u64::from_le_bytes(b.try_into().unwrap()))
                        .collect();
                    assert_eq!(actual, cells[1..]);
                    assert!(group.cells.read_block(slots, &mut [0; 8], None).is_err());
                }
                group.cells.remove().unwrap();
                assert!(!staged.exists());
            }
        }
    }

    fn scratch_free_config(low_disk: bool) -> BuildConfig {
        BuildConfig {
            backend: Backend::Zor,
            k: 3,
            hashes: 2,
            seed: 69,
            fingerprint_bits: 8,
            feature_config: FeatureConfig::legacy_kmers(),
            slot_scale: 1.0,
            minimum_slots: None,
            segment_length: None,
            merge: super::super::MergeArgs {
                disabled: true,
                sketch_size: 64,
                build_storage: if low_disk {
                    super::super::BuildStorage::Memory
                } else {
                    super::super::BuildStorage::Disk
                },
                ..Default::default()
            },
        }
    }

    fn dropped_test_graph(keys: &[u64], slots: usize) -> PeelGraph {
        let mut mapped = Vec::new();
        for &key in keys {
            let mut p = [0; 2];
            positions(key, slots, 2, 69, &mut p);
            mapped.extend(p.map(|v| v as u32));
        }
        PeelGraph::Compact(super::super::TypedPeelGraph {
            positions: mapped,
            order: Vec::new(),
            dropped: (0..keys.len() as u32).collect(),
        })
    }

    #[test]
    fn default_two_pass_stages_packed_cells_and_abandoned_keys_during_construction() {
        assert!(!super::super::MergeArgs::default().in_memory());
        two_pass_stages_only_packed_cells(false);
    }

    #[test]
    fn low_disk_two_pass_never_creates_scratch_files_during_construction() {
        two_pass_stages_only_packed_cells(true);
    }

    #[test]
    fn disk_group_retains_no_abandoned_records_and_promotion_uses_distinct_files() {
        let root = temporary_stage();
        fs::create_dir(root.0.join("cells")).unwrap();
        fs::create_dir(root.0.join("abandoned")).unwrap();
        let config = scratch_free_config(false);
        let keys = vec![1, 2, 3, 4, 5];
        let group = planner::Group {
            members: vec![7],
            cardinality: 5.0,
            slots: 8,
            width: 8,
            sketch: Sketch::new(&keys, 64, 69),
            ordering_key: None,
            feature_count_upper_bound: keys.len(),
            input_peak_bytes: 0,
        };
        let (first, _) = construct_group(
            &group,
            keys.clone(),
            vec![1; keys.len()],
            Some(&root.0),
            &config,
            &|keys, slots| Ok(dropped_test_graph(keys, slots)),
            None,
        )
        .unwrap();
        let (promoted, _) = construct_group(
            &group,
            keys.clone(),
            vec![1; keys.len()],
            Some(&root.0),
            &config,
            &|keys, slots| Ok(dropped_test_graph(keys, slots)),
            Some(64),
        )
        .unwrap();
        assert!(first.fallback.is_empty() && promoted.fallback.is_empty());
        assert_eq!(first.meta.fallback_keys, keys.len());
        assert_ne!(first.fallback_path, promoted.fallback_path);
        assert_eq!(
            fs::metadata(first.fallback_path.as_ref().unwrap())
                .unwrap()
                .len(),
            16 * keys.len() as u64
        );
        let promoted_path = promoted.fallback_path.clone().unwrap();
        first.discard().unwrap();
        assert!(promoted_path.exists());
        promoted.discard().unwrap();
        assert!(!promoted_path.exists());
        assert_eq!(fs::read_dir(root.0.join("cells")).unwrap().count(), 0);
        assert_eq!(fs::read_dir(root.0.join("abandoned")).unwrap().count(), 0);
    }

    #[test]
    fn global_disk_abandoned_merge_is_byte_identical_to_ram_with_hash_collisions() {
        let root = temporary_stage();
        let ram = root.0.join("ram");
        let disk = root.0.join("disk");
        fs::create_dir(&ram).unwrap();
        fs::create_dir(&disk).unwrap();
        fs::create_dir(disk.join("abandoned")).unwrap();
        let mut ram_groups = Vec::new();
        let mut disk_groups = Vec::new();
        // >64 group files forces an external merge pass. Sparse group IDs,
        // shared hashes, repeated masks and 0/MAX hashes exercise all unions.
        for id in 0..70usize {
            let mut records = vec![(u32::MAX, 1), (0, 2), (43, 3), (id as u32, 1)];
            records.extend(std::iter::repeat_n((19, 3), 1100));
            let group = || BuiltGroup {
                ordering_key: None,
                meta: GroupMeta {
                    members: vec![id * 3, id * 3 + 1],
                    width: 10,
                    offset: 0,
                    union_cardinality: records.len(),
                    fallback_keys: records.len(),
                    abandoned_keys: Some(records.len()),
                },
                slots: 8,
                cells: CompletedCells::Packed {
                    data: vec![0; 10],
                    width: 10,
                    slots: 8,
                },
                fallback: Vec::new(),
                fallback_path: None,
                routing: None,
            };
            let mut memory = group();
            memory.fallback = records.clone();
            ram_groups.push(memory);
            let mut staged = group();
            let path = disk.join("abandoned").join(format!("group-{id}.records"));
            super::super::staged_abandoned::write_sorted(
                &path,
                records
                    .into_iter()
                    .map(|(hash, mask)| (hash, (id * 3) as u32, mask)),
            )
            .unwrap();
            staged.fallback_path = Some(path);
            disk_groups.push(staged);
        }
        let expected = save_abandoned_groups(&mut ram_groups, &ram).unwrap();
        let actual = save_abandoned_groups(&mut disk_groups, &disk).unwrap();
        assert_eq!(actual, expected);
        assert_eq!(
            fs::read(ram.join(ABANDONED_FILE)).unwrap(),
            fs::read(disk.join(ABANDONED_FILE)).unwrap()
        );
        assert_eq!(fs::read_dir(disk.join("abandoned")).unwrap().count(), 0);
        assert!(disk_groups
            .iter()
            .all(|group| group.fallback.is_empty() && group.fallback_path.is_none()));
    }

    fn two_pass_stages_only_packed_cells(low_disk: bool) {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let root = temporary_stage();
        let first = root.0.join("a.fa");
        let second = root.0.join("b.fa");
        fs::write(&first, ">a\nACGTACGT\n").unwrap();
        fs::write(&second, ">b\nAAAAAAAA\n").unwrap();
        let fof = root.0.join("inputs.fof");
        fs::write(
            &fof,
            format!("a\t{}\nb\t{}\n", first.display(), second.display()),
        )
        .unwrap();
        let output = root.0.join("index");
        let calls = AtomicUsize::new(0);
        build(
            &fof,
            &output,
            scratch_free_config(low_disk),
            |keys, slots| {
                calls.fetch_add(1, Ordering::Relaxed);
                let stages: Vec<_> = fs::read_dir(&root.0)
                    .unwrap()
                    .map(|e| e.unwrap().path())
                    .filter(|p| {
                        p.file_name()
                            .unwrap()
                            .to_string_lossy()
                            .starts_with(".zor-pow2-building-")
                    })
                    .collect();
                assert_eq!(stages.len(), 1);
                // Before final packing, the publication directory holds no
                // features/, cache files or other construction artifacts: nothing
                // at all with --build-storage memory, packed cells AND abandoned-key runs
                // of completed groups otherwise.
                let entries: Vec<_> = fs::read_dir(&stages[0])
                    .unwrap()
                    .map(|e| e.unwrap().path())
                    .collect();
                if low_disk {
                    assert!(entries.is_empty());
                } else {
                    let entries: std::collections::BTreeSet<_> = entries.into_iter().collect();
                    assert_eq!(
                        entries,
                        [stages[0].join("cells"), stages[0].join("abandoned")].into()
                    );
                    for cell in fs::read_dir(stages[0].join("cells")).unwrap() {
                        let name = cell.unwrap().file_name();
                        assert!(name.to_string_lossy().ends_with(".bits"), "{name:?}");
                    }
                }
                Ok(dropped_test_graph(keys, slots))
            },
        )
        .unwrap();
        assert_eq!(calls.load(Ordering::Relaxed), 2);
        let index = Index::load(&output).unwrap();
        assert_eq!(
            index.manifest.construction_abandoned_storage.as_deref(),
            Some(if low_disk {
                "memory"
            } else {
                "temporary-disk-sorted-runs"
            })
        );
        assert_eq!(
            index.manifest.construction_cell_storage,
            if low_disk {
                "two-pass-memory-packed-no-scratch"
            } else {
                "two-pass-temporary-disk-packed-tiled-packing"
            }
        );
        for entry in fs::read_dir(&output).unwrap() {
            let entry = entry.unwrap();
            assert!(entry.file_type().unwrap().is_file());
            let name = entry.file_name().to_str().unwrap().to_owned();
            assert!(
                name == "manifest.json" || name == ABANDONED_FILE || name.ends_with(".zorbank.zst")
            );
        }
    }

    #[test]
    fn two_pass_rejects_changed_features_and_removes_unpublished_output() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let root = temporary_stage();
        let first = root.0.join("a.fa");
        let second = root.0.join("b.fa");
        fs::write(&first, ">a\nACGTACGT\n").unwrap();
        fs::write(&second, ">b\nAAAAAAAA\n").unwrap();
        let fof = root.0.join("inputs.fof");
        fs::write(
            &fof,
            format!("a\t{}\nb\t{}\n", first.display(), second.display()),
        )
        .unwrap();
        let calls = AtomicUsize::new(0);
        let output = root.0.join("index");
        // One thread builds the groups one after another, so the second group rereads its dataset after the
        // first group's callback has changed it (retaining cardinality and occurrence count, so only the feature
        // signature can detect the change).
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap();
        let result = pool.install(|| {
            build(&fof, &output, scratch_free_config(false), |keys, slots| {
                calls.fetch_add(1, Ordering::Relaxed);
                fs::write(&first, ">a\nTTTTTTTT\n").unwrap();
                fs::write(&second, ">b\nCCCCCCCC\n").unwrap();
                Ok(dropped_test_graph(keys, slots))
            })
        });
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("changed between passes"));
        assert_eq!(calls.load(Ordering::Relaxed), 1);
        assert!(!output.exists());
        assert_eq!(fs::read_dir(&root.0).unwrap().count(), 3);
    }

    #[test]
    fn one_bit_screening_ignores_padding_and_adjacent_origin_bits() {
        for count in [1usize, 7, 8, 9, 63, 64, 65, 257] {
            for padding in [0u8, 255] {
                let mut row = vec![padding; count.div_ceil(8)];
                for i in 0..count {
                    let mask = 1 << (i % 8);
                    row[i / 8] =
                        (row[i / 8] & !mask) | (((splitmix64(i as u64) & 1) as u8) << (i % 8));
                }
                for fp in [0, 1] {
                    let mut actual = Vec::new();
                    fingerprint_matches(&row, 1, fp, count, |i| actual.push(i));
                    let expected: Vec<_> = (0..count)
                        .filter(|&i| splitmix64(i as u64) & 1 == fp)
                        .collect();
                    assert_eq!(actual, expected);
                }
            }
        }
    }

    #[test]
    fn selected_subbyte_ranges_match_independent_field_reference() {
        for bits in [2usize, 4] {
            let unit = 8 / bits;
            let count = 513;
            let ranges = [unit..unit * 5, unit * 9..unit * 17, unit * 33..count];
            for hashes in [2, 5, 6, 7, 32] {
                let data: Vec<Vec<u8>> = (0..hashes)
                    .map(|h| {
                        (0..(count * bits).div_ceil(8))
                            .map(|i| splitmix64((h * 3333 + i) as u64) as u8)
                            .collect()
                    })
                    .collect();
                let rows: Vec<_> = data.iter().map(Vec::as_slice).collect();
                for fp in 0..1u64 << bits {
                    let expected: Vec<_> = ranges
                        .iter()
                        .flat_map(|range| range.clone())
                        .filter(|&i| {
                            rows.iter().fold(0u64, |v, r| {
                                v ^ ((r[i * bits / 8] >> (i * bits % 8)) as u64 & ((1 << bits) - 1))
                            }) == fp
                        })
                        .collect();
                    let mut actual = Vec::new();
                    fingerprint_matches_selected(&rows, bits, fp, &ranges, &mut |i| actual.push(i));
                    assert_eq!(actual, expected, "bits={bits} hashes={hashes} fp={fp}");
                }
            }
        }
    }

    #[test]
    fn staging_guard_cleans_partial_artifacts_on_error_and_unwind() {
        let parent = temporary_stage();
        let sentinel = parent.0.join("sentinel");
        fs::write(&sentinel, b"keep").unwrap();
        for panic in [false, true] {
            let path = parent.0.join(format!("partial-{panic}"));
            let result = std::panic::catch_unwind(|| -> Result<()> {
                fs::create_dir(&path)?;
                let _stage = Staging(path.clone());
                for directory in ["features", "cells", "abandoned"] {
                    fs::create_dir(path.join(directory))?;
                    fs::write(path.join(directory).join("partial"), b"unfinished")?;
                }
                if panic {
                    panic!("injected worker unwind");
                }
                anyhow::bail!("injected build error")
            });
            if panic {
                assert!(result.is_err());
            } else {
                assert!(result.unwrap().is_err());
            }
            assert!(!path.exists());
            assert_eq!(fs::read(&sentinel).unwrap(), b"keep");
        }
    }

    #[test]
    fn fused_screening_matches_materialized_rows_for_all_widths_and_hash_counts() {
        for bits in [1usize, 2, 4, 8, 16] {
            for count in [1usize, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 257] {
                for hashes in [2, 3, 4, 8, MAX_HASHES] {
                    let storage: Vec<Vec<u8>> = (0..hashes)
                        .map(|slot| {
                            (0..(count * bits).div_ceil(8))
                                .map(|byte| splitmix64((slot * 1024 + byte) as u64) as u8)
                                .collect()
                        })
                        .collect();
                    let rows: Vec<_> = storage.iter().map(Vec::as_slice).collect();
                    let mut materialized = vec![0; rows[0].len()];
                    for row in &rows {
                        xor_row(&mut materialized, row);
                    }
                    for fp in [0, 1, (1u64 << bits) - 1] {
                        let mut expected = Vec::new();
                        fingerprint_matches(&materialized, bits, fp, count, |i| expected.push(i));
                        let mut actual = Vec::new();
                        fingerprint_matches_rows(&rows, bits, fp, count, |i| actual.push(i));
                        assert_eq!(
                            actual, expected,
                            "bits={bits} count={count} hashes={hashes}"
                        );
                        if bits >= 8 {
                            actual.clear();
                            fingerprint_matches_rows_scalar(&rows, bits, fp, 0, &mut |i| {
                                actual.push(i)
                            });
                            assert_eq!(actual, expected);
                        }
                    }
                    // Guarantee a full set of hits, including zero fingerprints.
                    let zero = vec![0; rows[0].len()];
                    let mut actual = Vec::new();
                    fingerprint_matches_rows(&[&zero, &zero], bits, 0, count, |i| actual.push(i));
                    assert_eq!(actual, (0..count).collect::<Vec<_>>());
                }
            }
        }
    }

    #[test]
    fn shared_position_hash_matches_original_segment_mapping() {
        for hashes in [2, 3, 4, 8, MAX_HASHES] {
            for exponent in 5..=24 {
                let cells = 1usize << exponent;
                for key in [0, 1, 123, u64::MAX] {
                    let seed = 69;
                    let length = 1usize << (cells / hashes).clamp(1, 4096).ilog2();
                    let count = cells / length - hashes + 1;
                    let hash1 = crate::mixsplit(key, seed);
                    let hash2 = splitmix64(hash1);
                    let base = (((hash1 as u128 * count as u128) >> 64) as usize) * length;
                    let mut hash = hash1;
                    let expected: Vec<_> = (0..hashes)
                        .map(|i| {
                            let slot =
                                base + i * length + ((hash ^ (hash >> 33)) as usize & (length - 1));
                            hash = hash.wrapping_add(hash2);
                            slot
                        })
                        .collect();
                    let mut actual = vec![0; hashes];
                    PositionHash::new(key, seed).positions(cells, hashes, None, &mut actual);
                    assert_eq!(actual, expected);
                }
            }
        }
    }

    #[test]
    fn sparse_record_scores_preserve_thresholds_order_and_zero_count_documents() {
        let mut sparse = RecordScores::new(8);
        for origin in [7, 3, 7, 1, 3, 7] {
            sparse.increment(origin);
        }
        for total in [0, 3, 6] {
            for threshold in [0.0, 0.1, 1.0 / 3.0, 0.5, 1.0] {
                let mut expected = format!("{total}\t");
                let mut separator = "";
                for (origin, &count) in sparse.values.iter().enumerate() {
                    let ratio = if total == 0 {
                        0.0
                    } else {
                        count as f64 / total as f64
                    };
                    if ratio >= threshold {
                        expected.push_str(&format!("{separator}{origin}:{count}"));
                        separator = ",";
                    }
                }
                expected.push('\n');
                let mut actual = Vec::new();
                sparse.write(total, threshold, &mut actual).unwrap();
                assert_eq!(actual, expected.as_bytes());
            }
        }
        sparse.clear();
        assert!(sparse.values.iter().all(|&value| value == 0));
        assert!(sparse.touched.is_empty());
        sparse.increment(3);
        assert_eq!(sparse.touched, [3]);
        assert_eq!(sparse.values[3], 1);
    }

    #[test]
    fn sparse_lookup_preserves_abandoned_union_and_repeated_feature_counts() {
        for recovered in [&[][..], &[0][..], &[1][..], &[0, 1][..]] {
            let index = union_index(123, recovered);
            let mut dense = vec![0u64; index.manifest.colors.len()];
            let mut sparse = RecordScores::new(dense.len());
            let mut scratch = Vec::new();
            let mut sparse_scratch = Vec::new();
            for _record in 0..3 {
                dense.fill(0);
                sparse.clear();
                for key in [123, 123, 7, 19, u64::MAX, 123] {
                    index.lookup(key, &mut dense, &mut scratch);
                    index.lookup_into(key, &mut sparse, &mut sparse_scratch);
                    assert_eq!(dense, sparse.values);
                }
                assert_eq!(
                    sparse.touched.len(),
                    dense.iter().filter(|&&value| value != 0).count()
                );
            }
        }
    }

    #[test]
    fn batched_prefetch_preserves_dense_and_sparse_scores_and_repeated_keys() {
        for recovered in [&[][..], &[0][..], &[1][..], &[0, 1][..]] {
            let index = union_index(123, recovered);
            let mut positions = Vec::new();
            for length in [0, 1, 3, 4, 31, 32, 33, 63, 65, 257] {
                let keys: Vec<_> = (0..length)
                    .map(|i| {
                        if i % 3 == 0 {
                            123
                        } else {
                            splitmix64(i as u64)
                        }
                    })
                    .collect();
                let mut scalar = vec![0u64; index.manifest.colors.len()];
                let mut dense = scalar.clone();
                let mut sparse = RecordScores::new(scalar.len());
                for &key in &keys {
                    index.lookup(key, &mut scalar, &mut Vec::new());
                }
                index.lookup_batch(&keys, dense.as_mut_slice(), &mut Vec::new(), &mut positions);
                index.lookup_batch(&keys, &mut sparse, &mut Vec::new(), &mut positions);
                assert_eq!(scalar, dense);
                assert_eq!(scalar, sparse.values);
            }
        }
    }

    #[test]
    fn aggregate_deduplication_preserves_weighted_counts_across_workers_and_flushes() {
        for recovered in [&[][..], &[0][..], &[1][..], &[0, 1][..]] {
            let index = union_index(123, recovered);
            for threads in [1, 2, 4] {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap()
                    .install(|| {
                        let mut workers = Vec::new();
                        let mut expected = vec![0; index.manifest.colors.len()];
                        for keys in [
                            vec![123; 100],
                            (0..91).collect(),
                            vec![83_843, 95_900, 123, 83_843, 123, 95_900, 123],
                            Vec::new(),
                        ] {
                            for &key in &keys {
                                index.lookup(key, &mut expected, &mut Vec::new());
                            }
                            let mut batch = keys;
                            flush(&index, &mut batch, &mut workers);
                            assert!(batch.is_empty());
                            let actual: Vec<u64> = (0..expected.len())
                                .map(|origin| workers.iter().map(|w| w.scores[origin]).sum())
                                .collect();
                            assert_eq!(actual, expected);
                        }
                    });
            }
        }
    }

    #[test]
    fn simd_screening_matches_scalar_including_zero_and_tails() {
        for bits in [8, 16] {
            for count in [1, 15, 16, 17, 31, 32, 33, 257] {
                let row: Vec<_> = (0..count * bits / 8)
                    .map(|i| splitmix64(i as u64) as u8)
                    .collect();
                for fp in [0, 1, 255, if bits == 8 { 128 } else { 65535 }] {
                    let mut scalar = Vec::new();
                    fingerprint_matches_scalar(&row, bits, fp, 0, &mut |i| scalar.push(i));
                    let mut simd = Vec::new();
                    fingerprint_matches(&row, bits, fp, count, |i| simd.push(i));
                    assert_eq!(scalar, simd);
                }
                let mut hits = Vec::new();
                fingerprint_matches(&vec![0; row.len()], bits, 0, count, |i| hits.push(i));
                assert_eq!(hits, (0..count).collect::<Vec<_>>());
            }
        }
    }

    #[test]
    fn wide_bit_packing_matches_bit_reference() {
        for offset in 0usize..64 {
            for width in 1..=64 {
                let value = splitmix64((offset * 64 + width) as u64);
                for bytes in [(offset + width).div_ceil(8), 24] {
                    let mut reference = vec![0; bytes];
                    for bit in 0..width {
                        reference[(offset + bit) / 8] |=
                            (((value >> bit) & 1) as u8) << ((offset + bit) % 8);
                    }
                    let mut actual = vec![0; bytes];
                    put_bits(&mut actual, offset, width, value);
                    assert_eq!(actual, reference);
                }
            }
        }
    }

    fn temporary_stage() -> Staging {
        let stamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path =
            std::env::temp_dir().join(format!("zor-column-test-{}-{stamp}", std::process::id()));
        fs::create_dir(&path).unwrap();
        Staging(path)
    }

    fn test_group(stage: &Staging, id: usize, width: usize, cells: &[u64]) -> BuiltGroup {
        let path = stage.0.join(format!("{id}.u64"));
        cache_write(&path, cells).unwrap();
        BuiltGroup {
            ordering_key: None,
            meta: GroupMeta {
                members: vec![id],
                width,
                offset: 0,
                union_cardinality: 1,
                fallback_keys: 1,
                abandoned_keys: Some(1),
            },
            slots: cells.len(),
            cells: CompletedCells::Disk(path),
            fallback: vec![(id as u32, cells[0])],
            fallback_path: None,
            routing: None,
        }
    }

    #[test]
    fn bounded_union_combines_origins_and_rejects_invalid_caches() {
        let stage = temporary_stage();
        fs::create_dir(stage.0.join("features")).unwrap();
        let group = planner::Group {
            members: vec![0, 1, 2],
            cardinality: 5.0,
            slots: 8,
            width: 11,
            sketch: Sketch::new(&[], 4096, 69),
            ordering_key: None,
            feature_count_upper_bound: 5,
            input_peak_bytes: 0,
        };
        cache_write(&cache_path(&stage.0, 0), &[1, 3, 5]).unwrap();
        cache_write(&cache_path(&stage.0, 1), &[2, 3, 5]).unwrap();
        cache_write(&cache_path(&stage.0, 2), &[]).unwrap();
        assert_eq!(
            union_features(&group, &stage.0).unwrap(),
            (vec![1, 2, 3, 5], vec![1, 2, 3, 3])
        );
        cache_write(&cache_path(&stage.0, 2), &[3, 5]).unwrap();
        assert_eq!(
            union_features(&group, &stage.0).unwrap().1,
            vec![1, 2, 7, 7]
        );
        cache_write(&cache_path(&stage.0, 2), &[5, 3]).unwrap();
        assert!(union_features(&group, &stage.0).is_err());
        cache_write(&cache_path(&stage.0, 2), &[3, 3]).unwrap();
        assert!(union_features(&group, &stage.0).is_err());
        fs::write(cache_path(&stage.0, 2), [1, 2, 3]).unwrap();
        assert!(union_features(&group, &stage.0).is_err());
    }

    #[test]
    fn sketch_fingerprints_order_keeps_widths_contiguous_and_stable_ties() {
        let stage = temporary_stage();
        let mut groups = vec![
            test_group(&stage, 0, 18, &[9; 32]),
            test_group(&stage, 1, 16, &[250; 32]),
            test_group(&stage, 2, 18, &[3; 32]),
            test_group(&stage, 3, 18, &[9; 32]),
            test_group(&stage, 4, 16, &[2; 32]),
        ];
        for (group, key) in groups.iter_mut().zip([9, 250, 3, 9, 2]) {
            group.ordering_key = Some(std::sync::Arc::new(vec![key]));
            if let CompletedCells::Disk(path) = &group.cells {
                fs::remove_file(path).unwrap();
            }
        }
        order_columns(&mut groups, ColumnOrder::Compression);
        assert_eq!(
            groups.iter().map(|g| g.meta.members[0]).collect::<Vec<_>>(),
            [4, 1, 2, 0, 3]
        );
        order_columns(&mut groups, ColumnOrder::Width);
        assert_eq!(
            groups.iter().map(|g| g.meta.members[0]).collect::<Vec<_>>(),
            [1, 4, 0, 2, 3]
        );
    }

    #[test]
    fn tiled_planes_preserve_every_payload_across_banks() {
        for fp_bits in [1, 2, 4, 8, 16] {
            let stage = temporary_stage();
            let mut groups = Vec::new();
            for (id, count) in [1, 2, 3, 64 - fp_bits].into_iter().enumerate() {
                let width = fp_bits + if count == 1 { 0 } else { count };
                let cells: Vec<_> = (0..8192)
                    .map(|slot| {
                        let value = splitmix64((id * 8192 + slot) as u64);
                        if width == 64 {
                            value
                        } else {
                            value & ((1 << width) - 1)
                        }
                    })
                    .collect();
                let mut group = test_group(&stage, id, width, &cells);
                group.meta.members = (0..count).map(|m| id * 64 + m).collect();
                groups.push(group);
            }
            let banks =
                pack_banks_streaming(groups, ColumnOrder::Width, fp_bits, &stage.0).unwrap();
            for meta in banks {
                let mut reader = open_index_reader(&stage.0.join(&meta.file)).unwrap();
                let mut header = [0; 40];
                reader.read_exact(&mut header).unwrap();
                assert_eq!(&header[..8], PLANE_MAGIC);
                let mut data = Vec::new();
                reader.read_to_end(&mut data).unwrap();
                let stride = meta.row_bits.div_ceil(8);
                for (slot, row) in data.chunks_exact(stride).enumerate() {
                    for (i, group) in meta.groups.iter().enumerate() {
                        let mut actual = get_bits(row, i * fp_bits, fp_bits);
                        if group.width > fp_bits {
                            actual |= get_bits(
                                row,
                                meta.groups.len() * fp_bits + group.offset - i * fp_bits,
                                group.width - fp_bits,
                            ) << fp_bits;
                        }
                        let id = group.members[0] / 64;
                        let expected = splitmix64((id * 8192 + slot) as u64);
                        let expected = if group.width == 64 {
                            expected
                        } else {
                            expected & ((1 << group.width) - 1)
                        };
                        assert_eq!(actual, expected);
                    }
                }
            }
        }
    }

    #[test]
    fn compression_order_clusters_repeated_columns_and_preserves_payloads() {
        let stage = temporary_stage();
        let slots = 1024;
        let mut patterns: Vec<_> = (0..64).map(|id| (id, splitmix64(id as u64))).collect();
        patterns.sort_by_key(|&(_, hash)| hash);
        let patterns: Vec<_> = patterns.into_iter().map(|(id, _)| id / 4).collect();
        let make = |prefix: usize| {
            patterns
                .iter()
                .enumerate()
                .map(|(id, &pattern)| {
                    let cells: Vec<_> = (0..slots)
                        .map(|slot| splitmix64(slot as u64 ^ (pattern as u64 * 10007)) & 0xffff)
                        .collect();
                    let mut group = test_group(&stage, prefix + id, 16, &cells);
                    group.meta.members = vec![id];
                    group.ordering_key = Some(std::sync::Arc::new(vec![pattern as u8]));
                    group.fallback = vec![(id as u32, cells[0])];
                    group
                })
                .collect()
        };
        let width = transpose(make(0), ColumnOrder::Width).unwrap().remove(0);
        let compression = transpose(make(1000), ColumnOrder::Compression)
            .unwrap()
            .remove(0);
        assert_eq!(width.data.len(), compression.data.len());
        for (index, group) in compression.meta.groups.iter().enumerate() {
            let id = group.members[0];
            assert_eq!(compression.fallback[index], width.fallback[id]);
            for slot in 0..slots {
                let row = &compression.data[slot * 128..(slot + 1) * 128];
                assert_eq!(
                    get_bits(row, group.offset, group.width),
                    get_bits(&width.data[slot * 128..(slot + 1) * 128], id * 16, 16)
                );
            }
        }
        let old_size = zstd::bulk::compress(&width.data, 3).unwrap().len();
        let new_size = zstd::bulk::compress(&compression.data, 3).unwrap().len();
        assert!(
            new_size < old_size,
            "width={old_size}, compression={new_size}"
        );
        eprintln!("repeated-column fixture zstd bytes: {old_size} -> {new_size}");
    }
    #[test]
    fn bit_slices_cross_bytes_and_words_without_wasting_payload_bits() {
        let mut row = vec![0; 48];
        let mut offset = 0;
        for width in [8, 18, 19, 64, 9, 63, 16] {
            let value = if width == 64 {
                u64::MAX
            } else {
                splitmix64(width as u64) & ((1u64 << width) - 1)
            };
            put_bits(&mut row, offset, width, value);
            assert_eq!(get_bits(&row, offset, width), value);
            offset += width;
        }
    }

    #[test]
    fn fixed_segments_have_the_same_length_at_every_filter_size() {
        let length = super::super::DEFAULT_SEGMENT_LENGTH;
        for hashes in [4, 6, 8] {
            for cells in [512usize, 1024, 65536, 524288, 4194304] {
                for key in [0, 1, 69, 123456, u64::MAX] {
                    let hash = PositionHash::new(key, 69);
                    let mut positions = [0; MAX_HASHES];
                    hash.positions(cells, hashes, Some(length), &mut positions[..hashes]);
                    let base = positions[0] / length;
                    for (i, &slot) in positions[..hashes].iter().enumerate() {
                        assert!(slot < cells);
                        assert_eq!(slot / length, base + i);
                    }
                }
            }
        }
        let mut config = scratch_free_config(true);
        config.segment_length = Some(length);
        for hashes in [4usize, 6, 8] {
            config.hashes = hashes;
            config.validate().unwrap();
            assert_eq!(
                config.slots(0.0).unwrap(),
                length * hashes.next_power_of_two()
            );
            assert_eq!(config.slots(60000.0).unwrap(), 65536);
        }
    }

    pub(super) fn union_index(key: u64, recovered: &[u32]) -> Index {
        let seed = 69;
        let group = GroupMeta {
            members: vec![0, 1],
            width: 10,
            offset: 0,
            union_cardinality: 1,
            fallback_keys: usize::from(!recovered.is_empty()),
            abandoned_keys: Some(usize::from(!recovered.is_empty())),
        };
        let meta = BankMeta {
            component: 0,
            file: "bank-2.zorbank.zst".into(),
            fingerprint_first: false,
            routing: None,
            slots: 4,
            row_bits: 10,
            groups: vec![group],
        };
        let mut data = vec![0; 8];
        let mut slots = [0; 2];
        positions(key, 4, 2, seed, &mut slots);
        // Arrange for the XOR filter to return only dataset 0.
        put_bits(
            &mut data[slots[0] * 2..(slots[0] + 1) * 2],
            0,
            10,
            fingerprint(key, seed, 8) | (1 << 8),
        );
        let abandoned = if recovered.is_empty() {
            None
        } else {
            let mut map = OriginMap::default();
            map.push(abandoned::hash(key, seed), recovered).unwrap();
            map.finish().unwrap();
            Some(map)
        };
        let manifest = Manifest {
            version: 2,
            abandoned_recovery: !recovered.is_empty(),
            abandoned_file: (!recovered.is_empty()).then(|| ABANDONED_FILE.into()),
            abandoned_hashes: usize::from(!recovered.is_empty()),
            abandoned_origins: recovered.len(),
            backend: Backend::Zor,
            k: 31,
            hashes: 2,
            segment_length: None,
            seed,
            fingerprint_bits: 8,
            index_mode: IndexMode::Kmers.as_u8(),
            minimizer_size: 31,
            minimizer_scheme: MinimizerScheme::Simd,
            homopolymer_compression: false,
            seed_pattern: None,
            randstrobemers: None,
            randstrobe_scheme: None,
            modimizer_sampling: 10,
            findere_z: 10,
            sketch_size: 64,
            sketch_method: "test".into(),
            merge_planner: MergePlanner::Streaming,
            column_order: ColumnOrder::Width,
            column_order_method: None,
            planning_pair_evaluations: 0,
            planning_possible_pair_evaluations: 0,
            merge_reverts_enabled: false,
            construction_cell_storage: "temporary-disk-u64".into(),
            construction_abandoned_storage: None,
            slot_scale: 1.0,
            singleton_logical_bits: 64,
            min_subindex_documents: 1,
            max_subindexes: None,
            max_subindex_size: None,
            promoted_groups: 0,
            predicted_logical_bits: 40,
            planned_groups: 1,
            rejected_groups: 0,
            overflow_groups: 0,
            exact_group_verifications: 1,
            sketch_seconds: 0.0,
            plan_seconds: 0.0,
            build_seconds: 0.0,
            input_document_count: 2,
            source_document_ids: Vec::new(),
            skipped_documents: Vec::new(),
            memory_guard_split_groups: 0,
            merge_memory_budget_bytes: None,
            colors: ["a", "b"]
                .into_iter()
                .map(|name| Color {
                    name: name.into(),
                    path: name.into(),
                    occurrences: 1,
                    unique_kmers: 1,
                })
                .collect(),
            banks: vec![meta.clone()],
        };
        Index {
            manifest,
            banks: vec![Bank {
                meta,
                data,
                fallback: Vec::new(),
            }],
            abandoned,
            row_scratch_bytes: 2,
        }
    }

    #[test]
    fn fingerprint_planes_preserve_legacy_lookup_and_abandoned_union() {
        for bits in [1usize, 2, 4, 8, 16] {
            check_fingerprint_planes(bits);
        }
    }

    fn check_fingerprint_planes(bits: usize) {
        let mut legacy = union_index(123, &[1]);
        legacy.manifest.fingerprint_bits = bits;
        let mut next = 0;
        let mut offset = 0;
        let groups: Vec<_> = (0..40)
            .map(|i| {
                let count = i % 3 + 1;
                let members = (next..next + count).collect();
                next += count;
                let width = bits + if count == 1 { 0 } else { count };
                let group = GroupMeta {
                    members,
                    width,
                    offset,
                    union_cardinality: 0,
                    fallback_keys: 0,
                    abandoned_keys: Some(0),
                };
                offset += width;
                group
            })
            .collect();
        let stride = offset.div_ceil(8);
        legacy.banks[0].meta.groups = groups;
        legacy.banks[0].meta.row_bits = offset;
        legacy.banks[0].data = vec![0; 4 * stride];
        legacy.row_scratch_bytes = stride;
        let mut planes = union_index(123, &[1]);
        planes.manifest.fingerprint_bits = bits;
        planes.banks[0].meta = legacy.banks[0].meta.clone();
        planes.banks[0].meta.fingerprint_first = true;
        planes.banks[0].data = vec![0; 4 * stride];
        planes.row_scratch_bytes = (40 * bits).div_ceil(8);
        let groups = legacy.banks[0].meta.groups.clone();
        for slot in 0..4 {
            for (i, g) in groups.iter().enumerate() {
                let value = splitmix64((slot * 40 + i) as u64) & ((1 << g.width) - 1);
                let old_row = &mut legacy.banks[0].data[slot * stride..(slot + 1) * stride];
                let new_row = &mut planes.banks[0].data[slot * stride..(slot + 1) * stride];
                put_bits(old_row, g.offset, g.width, value);
                put_bits(new_row, i * bits, bits, value & ((1 << bits) - 1));
                if g.width > bits {
                    put_bits(
                        new_row,
                        40 * bits + g.offset - i * bits,
                        g.width - bits,
                        value >> bits,
                    );
                }
            }
        }
        let mut old_scratch = Vec::new();
        let mut new_scratch = Vec::new();
        let mut sparse_scratch = Vec::new();
        let mut sparse = RecordScores::new(next);
        for key in 0..20000 {
            let mut old = vec![0; next];
            let mut new = vec![0; next];
            legacy.lookup(key, &mut old, &mut old_scratch);
            planes.lookup(key, &mut new, &mut new_scratch);
            assert_eq!(old, new, "key {key}");
            sparse.clear();
            planes.lookup_into(key, &mut sparse, &mut sparse_scratch);
            assert_eq!(old, sparse.values, "sparse key {key}");
        }
        let keys: Vec<_> = (0..2000).collect();
        let mut old = vec![0; next];
        for &key in &keys {
            legacy.lookup(key, &mut old, &mut old_scratch);
        }
        sparse.clear();
        planes.lookup_batch(&keys, &mut sparse, &mut sparse_scratch, &mut Vec::new());
        assert_eq!(old, sparse.values, "batched fingerprint width {bits}");
    }

    #[test]
    fn recovery_adds_origins_never_overwrites_and_counts_each_origin_once() {
        for (recovered, expected) in [
            (&[][..], [1, 0]),
            (&[0][..], [1, 0]),
            (&[1][..], [1, 1]),
            (&[0, 1][..], [1, 1]),
        ] {
            let index = union_index(123, recovered);
            let mut scores = [0; 2];
            let mut scratch = Vec::new();
            index.lookup(123, &mut scores, &mut scratch);
            assert_eq!(scores, expected);
            index.lookup(123, &mut scores, &mut scratch);
            assert_eq!(scores, expected.map(|n| n * 2));
        }
    }

    #[test]
    fn hash_collisions_union_origins_across_groups_and_deduplicate() {
        // Two distinct input keys with the same *actual* 32-bit hash at seed 69.
        let a = 83_843;
        let b = 95_900;
        let hash = abandoned::hash(a, 69);
        assert_eq!(hash, abandoned::hash(b, 69));
        let mut records = vec![(hash, 0, 1), (hash, 1, 1), (hash, 0, 1), (hash, 0, 2)];
        let mut map = OriginMap::default();
        let stats = visit_abandoned(&mut records, &[vec![2, 0], vec![1]], |hash, origins| {
            map.push(hash, origins)
        })
        .unwrap();
        map.finish().unwrap();
        assert_eq!(stats, (1, 3));
        assert_eq!(map.get(abandoned::hash(a, 69)), [0, 1, 2]);
        assert_eq!(map.get(abandoned::hash(b, 69)), [0, 1, 2]);
        // There is no original 64-bit key check gating recovery.
        let mut index = union_index(a, &[]);
        let mut colliding_map = OriginMap::default();
        colliding_map.push(abandoned::hash(b, 69), &[1]).unwrap();
        colliding_map.finish().unwrap();
        index.abandoned = Some(colliding_map);
        let mut scores = [0; 2];
        index.lookup(a, &mut scores, &mut Vec::new());
        assert_eq!(scores, [1, 1]);
    }

    #[test]
    fn parallel_packed_ram_transpose_matches_disk_tiles_byte_for_byte() {
        let stage = temporary_stage();
        for fp in [1usize, 2, 4, 8, 16] {
            let make = |ram, prefix| {
                (0..256)
                    .map(|id| {
                        let count = [1, 2, 3, 7, 64 - fp][id % 5];
                        let width = fp + if count == 1 { 0 } else { count };
                        let mask = u64::MAX >> (64 - width);
                        let cells: Vec<_> = (0..4096)
                            .map(|slot| splitmix64((id * 4096 + slot) as u64) & mask)
                            .collect();
                        let mut group = test_group(&stage, prefix + id, width, &cells);
                        group.meta.members = (0..count).map(|i| id * 64 + i).collect();
                        if ram {
                            group.cells.remove().unwrap();
                            group.cells = CompletedCells::packed(&cells, width).unwrap();
                        }
                        group
                    })
                    .collect::<Vec<_>>()
            };
            let reference =
                pack_banks_streaming(make(false, 0), ColumnOrder::Compression, fp, &stage.0)
                    .unwrap();
            let bytes = fs::read(stage.0.join(&reference[0].file)).unwrap();
            for threads in [1, 4] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                let groups = make(true, threads * 1000);
                let actual = pool
                    .install(|| {
                        pack_banks_streaming(groups, ColumnOrder::Compression, fp, &stage.0)
                    })
                    .unwrap();
                assert_eq!(bytes, fs::read(stage.0.join(&actual[0].file)).unwrap());
                assert_eq!(reference[0].row_bits, actual[0].row_bits);
            }
        }
    }

    #[test]
    fn variable_integers_and_delta_recovery_preserve_full_u32_values() {
        for value in [0, 1, 127, 128, 16383, 16384, 1 << 28, u32::MAX] {
            let mut bytes = Vec::new();
            write_var_u32(&mut bytes, value).unwrap();
            assert_eq!(read_var_u32(&mut bytes.as_slice()).unwrap(), value);
        }
        assert!(read_var_u32(&mut &[255u8; 5][..]).is_err());
        assert!(read_var_u32(&mut &[128u8][..]).is_err());
        let values = [0, 1, 127, 128, 16383, 16384, 1 << 28, u32::MAX];
        let mut encoded = Vec::new();
        for &value in &values {
            write_var_u32(&mut encoded, value).unwrap();
        }
        for capacity in 1..=8 {
            let mut reader = BufReader::with_capacity(capacity, encoded.as_slice());
            for value in values {
                assert_eq!(read_var_u32(&mut reader).unwrap(), value);
            }
            assert!(read_var_u32(&mut reader).is_err());
        }
        let stage = temporary_stage();
        let mut index = union_index(123, &[]);
        index
            .manifest
            .colors
            .resize(65, index.manifest.colors[0].clone());
        index.manifest.abandoned_file = Some(ABANDONED_FILE.into());
        let expected = [
            (0u32, vec![0u32]),
            (7, vec![0, 32, 64]),
            (u32::MAX, vec![64]),
        ];
        index.manifest.abandoned_hashes = expected.len();
        index.manifest.abandoned_origins = 5;
        let mut data = ABANDONED_DELTA_MAGIC.to_vec();
        let mut previous = 0;
        for (hash, origins) in &expected {
            write_var_u32(&mut data, hash - previous).unwrap();
            previous = *hash;
            write_var_u32(&mut data, origins.len() as u32).unwrap();
            let mut prior = 0;
            for &id in origins {
                write_var_u32(&mut data, id - prior).unwrap();
                prior = id;
            }
        }
        for mode in 0..3 {
            let mut writer = create_index_writer(&stage.0.join(ABANDONED_FILE)).unwrap();
            writer
                .write_all(if mode == 1 {
                    &data[..data.len() - 1]
                } else {
                    &data
                })
                .unwrap();
            if mode == 2 {
                writer.write_all(&[0]).unwrap();
            }
            finish_index_writer(writer).unwrap();
            let result = load_abandoned(&stage.0, &index.manifest);
            if mode == 0 {
                let map = result.unwrap();
                for (hash, origins) in &expected {
                    assert_eq!(map.get(*hash), origins);
                }
            } else {
                assert!(result.is_err());
            }
        }
    }

    #[test]
    fn global_map_roundtrip_and_skip_loading() {
        let stage = temporary_stage();
        let key = 123;
        let mut index = union_index(key, &[1]);
        index.banks[0].fallback = vec![vec![(abandoned::hash(key, 69), 2)]];
        let stats = save_abandoned(&mut index.banks, &stage.0).unwrap();
        assert_eq!(stats, (1, 1));
        assert!(index.banks[0].fallback.is_empty());
        save_bank(&index.banks[0], &stage.0).unwrap();
        serde_json::to_writer(
            File::create(stage.0.join("manifest.json")).unwrap(),
            &index.manifest,
        )
        .unwrap();
        let bytes = zstd::decode_all(File::open(stage.0.join(ABANDONED_FILE)).unwrap()).unwrap();
        assert_eq!(bytes.len(), 8 + 4 + 4 + 4); // magic, hash32, count32, dataset32
        assert_eq!(&bytes[8..12], &abandoned::hash(key, 69).to_le_bytes());
        let loaded = Index::load(&stage.0).unwrap();
        assert!(loaded.banks.iter().all(|bank| bank.fallback.is_empty()));
        let mut scores = [0; 2];
        loaded.lookup(key, &mut scores, &mut Vec::new());
        assert_eq!(scores, [1, 1]);
        fs::remove_file(stage.0.join(ABANDONED_FILE)).unwrap();
        let skipped = Index::load_with_recovery(&stage.0, false).unwrap();
        assert!(skipped.abandoned.is_none());
        let mut scores = [0; 2];
        skipped.lookup(key, &mut scores, &mut Vec::new());
        assert_eq!(scores, [1, 0]);
        assert!(Index::load(&stage.0).is_err());
    }

    #[test]
    fn bulk_abandoned_postings_cross_buffers_and_reject_truncation_and_trailing_data() {
        let stage = temporary_stage();
        let mut index = union_index(123, &[]);
        // More than one 64KiB postings block, with a partial final word block.
        let count = 16_387;
        index
            .manifest
            .colors
            .resize(count, index.manifest.colors[0].clone());
        index.manifest.abandoned_file = Some(ABANDONED_FILE.into());
        index.manifest.abandoned_hashes = 1;
        index.manifest.abandoned_origins = count;
        let hash = abandoned::hash(123, 69);
        let mut data = Vec::new();
        data.extend_from_slice(ABANDONED_MAGIC);
        data.extend_from_slice(&hash.to_le_bytes());
        data.extend_from_slice(&(count as u32).to_le_bytes());
        for origin in 0..count as u32 {
            data.extend_from_slice(&origin.to_le_bytes());
        }
        for length in [data.len(), data.len() - 1, data.len() - 4] {
            let mut writer = create_index_writer(&stage.0.join(ABANDONED_FILE)).unwrap();
            writer.write_all(&data[..length]).unwrap();
            finish_index_writer(writer).unwrap();
            let loaded = load_abandoned(&stage.0, &index.manifest);
            if length == data.len() {
                assert_eq!(
                    loaded.unwrap().get(hash),
                    (0..count as u32).collect::<Vec<_>>()
                );
            } else {
                assert!(loaded.is_err());
            }
        }
        data.push(0);
        let mut writer = create_index_writer(&stage.0.join(ABANDONED_FILE)).unwrap();
        writer.write_all(&data).unwrap();
        finish_index_writer(writer).unwrap();
        assert!(load_abandoned(&stage.0, &index.manifest).is_err());
    }

    #[test]
    fn map_postings_are_checked_for_invalid_or_duplicate_dataset_ids() {
        let stage = temporary_stage();
        let mut index = union_index(123, &[1]);
        save_bank(&index.banks[0], &stage.0).unwrap();
        let corrupt_origins: [&[u32]; 4] = [&[2], &[1, 1], &[1, 0], &[]];
        for origins in corrupt_origins {
            index.manifest.abandoned_origins = origins.len();
            serde_json::to_writer(
                File::create(stage.0.join("manifest.json")).unwrap(),
                &index.manifest,
            )
            .unwrap();
            let mut writer = create_index_writer(&stage.0.join(ABANDONED_FILE)).unwrap();
            writer.write_all(ABANDONED_MAGIC).unwrap();
            writer
                .write_all(&abandoned::hash(123, 69).to_le_bytes())
                .unwrap();
            writer
                .write_all(&(origins.len() as u32).to_le_bytes())
                .unwrap();
            for &origin in origins {
                writer.write_all(&origin.to_le_bytes()).unwrap();
            }
            finish_index_writer(writer).unwrap();
            assert!(Index::load(&stage.0).is_err(), "origins={origins:?}");
        }
    }

    #[test]
    fn one_map_recovers_across_banks_with_different_row_widths() {
        let stage = temporary_stage();
        let key = 123;
        let mut index = union_index(key, &[1]);
        index.banks[0].fallback = vec![vec![(abandoned::hash(key, 69), 2)]];
        let meta = BankMeta {
            component: 0,
            file: "bank-3.zorbank.zst".into(),
            fingerprint_first: false,
            routing: None,
            slots: 8,
            row_bits: 8,
            groups: vec![GroupMeta {
                members: vec![2],
                width: 8,
                offset: 0,
                union_cardinality: 1,
                fallback_keys: 1,
                abandoned_keys: Some(1),
            }],
        };
        let mut slots = [0; 2];
        positions(key, 8, 2, 69, &mut slots);
        let mut data = vec![0; 8];
        data[slots[0]] = fingerprint(key, 69, 8) as u8;
        index.banks.push(Bank {
            meta: meta.clone(),
            data,
            fallback: vec![vec![(abandoned::hash(key, 69), 1)]],
        });
        index.manifest.banks.push(meta);
        index.manifest.colors.push(Color {
            name: "c".into(),
            path: "c".into(),
            occurrences: 1,
            unique_kmers: 1,
        });
        let (hashes, origins) = save_abandoned(&mut index.banks, &stage.0).unwrap();
        assert_eq!((hashes, origins), (1, 2));
        index.manifest.abandoned_hashes = hashes;
        index.manifest.abandoned_origins = origins;
        for bank in &index.banks {
            save_bank(bank, &stage.0).unwrap();
        }
        serde_json::to_writer(
            File::create(stage.0.join("manifest.json")).unwrap(),
            &index.manifest,
        )
        .unwrap();
        let loaded = Index::load(&stage.0).unwrap();
        assert_eq!(
            loaded
                .abandoned
                .as_ref()
                .unwrap()
                .get(abandoned::hash(key, 69)),
            [1, 2]
        );
        let mut scores = [0; 3];
        let mut scratch = Vec::new();
        loaded.lookup(key, &mut scores, &mut scratch);
        assert_eq!(scores, [1, 1, 1]);
        loaded.lookup(key, &mut scores, &mut scratch);
        assert_eq!(scores, [2, 2, 2]);
        let skipped = Index::load_with_recovery(&stage.0, false).unwrap();
        let mut raw = [0; 3];
        skipped.lookup(key, &mut raw, &mut scratch);
        assert_eq!(raw, [1, 0, 1]);

        // A shared abandoned hash spans both subindexes. Loading either one
        // must retain only its own origins, remap IDs, and never double count.
        for (selected, expected, expected_raw) in
            [(0, vec![1, 1], vec![1, 0]), (1, vec![1], vec![1])]
        {
            let other_path = stage.0.join(&index.manifest.banks[1 - selected].file);
            let other_data = fs::read(&other_path).unwrap();
            fs::remove_file(&other_path).unwrap();
            // The other subindex is absent: this also proves it is not read.
            for (recover, expected) in [(true, &expected), (false, &expected_raw)] {
                let partial = Index::load_selected(&stage.0, recover, Some(selected)).unwrap();
                assert_eq!(partial.banks.len(), 1);
                assert_eq!(partial.manifest.colors.len(), expected.len());
                if recover {
                    assert_eq!(partial.abandoned.as_ref().unwrap().origins_len(), 1);
                } else {
                    assert!(partial.abandoned.is_none());
                }
                let mut actual = vec![0; expected.len()];
                partial.lookup(key, actual.as_mut_slice(), &mut Vec::new());
                assert_eq!(&actual, expected);
            }
            fs::write(other_path, other_data).unwrap();
        }

        // Version 1 embeds abandoned keys after each subindex rather than in
        // the global file. In particular, subindex 1's original ID 2 must map
        // to local ID 0, including the abandoned-key association.
        index.manifest.version = 1;
        index.manifest.abandoned_file = None;
        index.manifest.abandoned_hashes = 0;
        index.manifest.abandoned_origins = 0;
        for bank in &index.banks {
            let mut writer = create_index_writer(&stage.0.join(&bank.meta.file)).unwrap();
            writer.write_all(MAGIC).unwrap();
            for n in [
                bank.meta.slots,
                bank.meta.row_bits,
                bank.meta.groups.len(),
                bank.data.len(),
            ] {
                writer.write_all(&(n as u64).to_le_bytes()).unwrap();
            }
            writer.write_all(&bank.data).unwrap();
            writer.write_all(&1u64.to_le_bytes()).unwrap();
            writer.write_all(&key.to_le_bytes()).unwrap();
            let value = fingerprint(key, 69, 8) | if bank.meta.row_bits == 10 { 2 << 8 } else { 0 };
            writer.write_all(&value.to_le_bytes()).unwrap();
            finish_index_writer(writer).unwrap();
        }
        serde_json::to_writer(
            File::create(stage.0.join("manifest.json")).unwrap(),
            &index.manifest,
        )
        .unwrap();
        for selected in 0..2 {
            let partial = Index::load_selected(&stage.0, true, Some(selected)).unwrap();
            let mut actual = vec![0; partial.manifest.colors.len()];
            partial.lookup(key, actual.as_mut_slice(), &mut Vec::new());
            assert_eq!(actual, vec![1; if selected == 0 { 2 } else { 1 }]);
        }
    }

    #[test]
    fn legacy_per_group_fallback_is_converted_to_one_32_bit_map() {
        let stage = temporary_stage();
        let key = 123;
        let mut index = union_index(key, &[1]);
        index.manifest.version = 1;
        index.manifest.abandoned_file = None;
        index.manifest.abandoned_hashes = 0;
        index.manifest.abandoned_origins = 0;
        index.manifest.banks[0].groups[0].abandoned_keys = None;
        let bank = &index.banks[0];
        let mut writer = create_index_writer(&stage.0.join(&bank.meta.file)).unwrap();
        writer.write_all(MAGIC).unwrap();
        for n in [
            bank.meta.slots,
            bank.meta.row_bits,
            bank.meta.groups.len(),
            bank.data.len(),
        ] {
            writer.write_all(&(n as u64).to_le_bytes()).unwrap();
        }
        writer.write_all(&bank.data).unwrap();
        writer.write_all(&1u64.to_le_bytes()).unwrap();
        writer.write_all(&key.to_le_bytes()).unwrap();
        writer
            .write_all(&(fingerprint(key, 69, 8) | (2 << 8)).to_le_bytes())
            .unwrap();
        finish_index_writer(writer).unwrap();
        serde_json::to_writer(
            File::create(stage.0.join("manifest.json")).unwrap(),
            &index.manifest,
        )
        .unwrap();
        let loaded = Index::load(&stage.0).unwrap();
        assert_eq!(
            loaded
                .abandoned
                .as_ref()
                .unwrap()
                .get(abandoned::hash(key, 69)),
            [1]
        );
        assert!(loaded.banks[0].fallback.is_empty());
        let mut scores = [0; 2];
        loaded.lookup(key, &mut scores, &mut Vec::new());
        assert_eq!(scores, [1, 1]);
        let skipped = Index::load_with_recovery(&stage.0, false).unwrap();
        assert!(skipped.abandoned.is_none());
        let partial = Index::load_selected(&stage.0, true, Some(0)).unwrap();
        let mut scores = [0; 2];
        partial.lookup(key, &mut scores, &mut Vec::new());
        assert_eq!(scores, [1, 1]);
    }
}
