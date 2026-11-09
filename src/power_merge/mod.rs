//! Merge-only colored indexes with one slot-major bank per power of two.
mod abandoned;
mod document_guard;
mod graph;
mod layout;
mod planner;
mod records;
mod sketch;
mod sort;
mod staged_abandoned;
mod storage;
mod union;

use crate::FeatureConfig;
use anyhow::{ensure, Context, Result};
use clap::{Args, ValueEnum};
pub use graph::{DegreeQueue, GraphIndex, IndexedDegreeHeap, PeelGraph, TypedPeelGraph};
pub use memory::trace_memory;
use serde::{Deserialize, Serialize};
/// Fixed segment size for newly built merged indexes, independent of filter size.
pub const DEFAULT_SEGMENT_LENGTH: usize = 64;
pub const DEFAULT_MAX_SUBINDEX_SIZE: u64 = 8 * 1024 * 1024 * 1024;

fn parse_subindex_size(value: &str) -> std::result::Result<u64, String> {
    let value = value.trim();
    let (number, scale) = [
        ("GiB", 1u64 << 30),
        ("MiB", 1u64 << 20),
        ("KiB", 1u64 << 10),
        ("B", 1),
    ]
    .into_iter()
    .find_map(|(suffix, scale)| value.strip_suffix(suffix).map(|n| (n, scale)))
    .unwrap_or((value, 1));
    let bytes = number
        .parse::<u64>()
        .ok()
        .and_then(|n| n.checked_mul(scale))
        .ok_or_else(|| "expected bytes or an integer followed by KiB, MiB or GiB".to_owned())?;
    if !bytes.is_power_of_two() {
        return Err("maximum subindex size must be a non-zero power of two".into());
    }
    Ok(bytes)
}
pub use storage::{
    build, fuse_indexes, fuse_many_indexes, info_command, query_command, query_command_with_memory,
    query_command_with_recovery, query_records_command, query_records_command_with_memory,
    RecordFormat,
};
mod memory;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub enum QueryMemory {
    /// Prepare query seeds once, scan abandoned keys once, and load one subindex at a time.
    #[default]
    Subindex,
    /// Load all subindexes and the global abandoned-key map together for maximum throughput.
    All,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum MergePlanner {
    /// Full maximum-saving agglomeration.
    #[value(name = "full")]
    #[serde(alias = "full")]
    AllPairs,
    #[default]
    Streaming,
    /// Exhaustive reference implementation used only by planner tests.
    #[cfg(test)]
    #[value(skip)]
    StreamingUnpruned,
    /// No dataset unions: append each dataset to the smallest sufficient existing subindex in FOF order.
    Naive,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ColumnOrder {
    /// Keep width sections contiguous, then use input dataset IDs.
    Width,
    /// Within width sections, order columns by retained union-sketch fingerprints.
    #[default]
    Compression,
}
impl ColumnOrder {
    pub fn label(self) -> &'static str {
        match self {
            Self::Width => "width",
            Self::Compression => "compression",
        }
    }
}
impl MergePlanner {
    pub fn label(self) -> &'static str {
        match self {
            Self::AllPairs => "full",
            Self::Streaming => "streaming",
            #[cfg(test)]
            Self::StreamingUnpruned => "streaming-unpruned",
            Self::Naive => "naive",
        }
    }

    pub fn is_streaming(self) -> bool {
        match self {
            Self::Streaming | Self::Naive => true,
            #[cfg(test)]
            Self::StreamingUnpruned => true,
            _ => false,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub enum BuildStorage {
    /// Reread inputs and stage packed filters on disk.
    #[default]
    Disk,
    /// Cache extracted seeds on disk to avoid rereading inputs.
    Cached,
    /// Reread inputs and keep completed filters in RAM, without scratch files.
    Memory,
}

#[derive(Args, Clone, Debug)]
pub struct MergeArgs {
    /// streaming (default), full, or naive (ordered placement without unions).
    #[arg(short = 'P', long = "merge-planner", value_enum, default_value_t = MergePlanner::Streaming)]
    pub planner: MergePlanner,
    /// Number of SimdSketch buckets (8-bit comparison fingerprints + 32-bit union minima).
    ///
    /// A power of two in 64..=1048576, default 128. Larger sketches improve
    /// estimation precision at the cost of sketch memory and comparison work.
    /// Use this to tune planning accuracy against construction resources.
    #[arg(short = 'K', long = "merge-sketch-size", default_value_t = 128)]
    pub sketch_size: usize,
    /// Plan each merged group as if its estimated union were this fraction larger, so sketch noise does not push
    /// the exact union past the planned power of two at build time (which doubles the group's slots).
    ///
    /// A finite fraction in [0,0.5), default 0.03. This provides headroom for
    /// sketch estimation error when choosing a capacity. A larger margin can
    /// reject merges whose estimated savings depend on a tight capacity fit.
    #[arg(short = 'G', long = "merge-plan-margin", default_value_t = 0.03)]
    pub plan_margin: f64,
    /// Minimum original documents per subindex. Promote smaller filters, never fold.
    /// Collections smaller than this use one subindex. Set 1 to retain small size classes.
    ///
    /// Positive count, default 128. Promote undersized classes to larger filter
    /// capacities to reduce the number of size classes probed during queries.
    /// This trades filter space for fewer class probes. The naive planner ignores it.
    #[arg(short = 'D', long, default_value_t = 128)]
    pub min_subindex_documents: usize,
    /// Maximum size-class count before hash partitioning. Promote smaller filters. Omit for unlimited.
    ///
    /// Positive count. Use this to put an explicit bound on query size classes,
    /// potentially increasing filter capacity through promotion. Physical hash
    /// partitions are outside this count. The naive planner ignores it.
    #[arg(short = 'N', long)]
    pub max_subindexes: Option<usize>,
    /// Maximum uncompressed merged-subindex payload. Hash-partition larger size classes.
    /// A power of two, in bytes or with KiB/MiB/GiB suffix (default: 8 GiB).
    #[arg(short = 'X', long, default_value = "8GiB", value_parser = parse_subindex_size)]
    pub max_subindex_size: u64,
    /// Keep singleton power-of-two banks, for measuring sizing independently of merging.
    ///
    /// Disable dataset unions while retaining the regular construction and
    /// size-class-promotion pipeline. This is an experimental comparison control,
    /// distinct from the naive planner's input-order placement.
    #[arg(short = 'U', long = "merge-disabled")]
    pub disabled: bool,
    /// Omit abandoned-key recovery storage and lookup (merged indexes only).
    /// Abandoned features may then be missed.
    ///
    /// Save recovery-file space and its query work by omitting those records
    /// during construction. This deliberately allows false negatives. Restoring
    /// recovery later requires rebuilding the index.
    #[arg(short = 'R', long)]
    pub drop_abandoned: bool,
    /// Construction storage: disk staging, cached seeds, or in-memory filters.
    #[arg(short = 'e', long, value_enum, default_value = "disk")]
    pub build_storage: BuildStorage,

    /// Total process memory limit in decimal GB. Default: 80% of available RAM.
    ///
    /// The same limit controls workspace admission and sampled process RSS.
    /// Every dataset must fit. Construction fails rather than excluding inputs.
    #[arg(short = 'g', long)]
    pub max_ram_gb: Option<f64>,

    /// Runtime-only frozen workspace limit shared with merge candidate checks.
    #[arg(skip)]
    pub planning_memory_budget_bytes: Option<u64>,
}
impl MergeArgs {
    fn cache_features(&self) -> bool {
        self.build_storage == BuildStorage::Cached
    }
    fn in_memory(&self) -> bool {
        self.build_storage == BuildStorage::Memory
    }
}
impl Default for MergeArgs {
    fn default() -> Self {
        Self {
            planner: MergePlanner::Streaming,
            sketch_size: 128,
            plan_margin: 0.03,
            min_subindex_documents: 128,
            max_subindexes: None,
            max_subindex_size: DEFAULT_MAX_SUBINDEX_SIZE,
            disabled: false,
            drop_abandoned: false,
            build_storage: BuildStorage::Disk,
            max_ram_gb: None,
            planning_memory_budget_bytes: None,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Backend {
    Zor,
}

#[derive(Clone, Debug)]
pub struct BuildConfig {
    pub backend: Backend,
    pub k: usize,
    pub hashes: usize,
    pub seed: u64,
    pub fingerprint_bits: usize,
    pub feature_config: FeatureConfig,
    pub slot_scale: f64,
    /// Cell-count floor for internal layout fixtures. Production sizing is automatic.
    #[cfg(test)]
    pub minimum_slots: Option<usize>,
    /// Fixed segment size. None is retained for legacy library/layout fixtures.
    pub segment_length: Option<usize>,
    pub merge: MergeArgs,
}
impl BuildConfig {
    pub fn validate(&self) -> Result<()> {
        ensure!(
            self.merge.max_subindex_size.is_power_of_two(),
            "--max-subindex-size must be a non-zero power of two"
        );
        ensure!(
            self.merge.max_ram_gb.is_none_or(|v| v.is_finite()
                && v > 0.0
                && v * 1_000_000_000.0 >= 1.0
                && v * 1_000_000_000.0 < u64::MAX as f64),
            "--max-ram-gb must be a positive finite decimal GB value of at least one byte"
        );
        ensure!(
            self.merge.max_subindexes != Some(0),
            "--max-subindexes must be positive"
        );
        ensure!(
            self.merge.min_subindex_documents > 0,
            "--min-subindex-documents must be positive"
        );
        ensure!(
            [1, 2, 4, 8, 16].contains(&self.fingerprint_bits),
            "fingerprint bits must be 1, 2, 4, 8 or 16"
        );
        match self.feature_config.mode {
            crate::IndexMode::Kmers | crate::IndexMode::Modimizers => {
                ensure!(self.k > 0, "k must be positive")
            }
            crate::IndexMode::Minimizers => ensure!(
                (1..=255).contains(&self.k)
                    && self.k % 2 == 1
                    && (1..=31).contains(&self.feature_config.minimizer_size)
                    && self.feature_config.minimizer_size <= self.k,
                "invalid minimizer window/length"
            ),
            crate::IndexMode::Findere => {
                crate::findere::indexed_length(self.k, self.feature_config.findere_z)?;
            }
        }
        ensure!(
            self.merge.sketch_size.is_power_of_two()
                && (64..=1_048_576).contains(&self.merge.sketch_size),
            "--merge-sketch-size must be a power of two in 64..=1048576"
        );
        ensure!(
            self.merge.plan_margin.is_finite() && (0.0..0.5).contains(&self.merge.plan_margin),
            "--merge-plan-margin must be in [0, 0.5)"
        );
        ensure!(
            (2..=crate::MAX_HASHES).contains(&self.hashes),
            "invalid hash count"
        );
        ensure!(
            self.slot_scale.is_finite() && self.slot_scale > 0.0,
            "invalid slot scale"
        );
        ensure!(
            self.segment_length
                .is_none_or(|n| n.is_power_of_two() && n <= 4096),
            "segment length must be a power of two in 1..=4096"
        );
        Ok(())
    }
    /// A multi-dataset payload stores one membership bit per dataset beside its fingerprint.
    pub fn max_group_members(&self) -> usize {
        64usize.saturating_sub(self.fingerprint_bits)
    }
    pub fn width(&self, colors: usize) -> usize {
        self.fingerprint_bits + if colors == 1 { 0 } else { colors }
    }
    pub fn slots(&self, cardinality: f64) -> Result<usize> {
        let needed = (cardinality.max(1.0) * self.slot_scale).ceil();
        ensure!(
            needed.is_finite() && needed < usize::MAX as f64,
            "filter capacity overflow"
        );
        let needed = needed as usize;
        #[cfg(test)]
        let needed = needed.max(self.minimum_slots.unwrap_or(0));
        let needed = needed.max(self.hashes.next_power_of_two()).max(
            self.segment_length
                .unwrap_or(1)
                .checked_mul(self.hashes.next_power_of_two())
                .context("segment geometry overflow")?,
        );
        needed
            .checked_next_power_of_two()
            .context("power-of-two capacity overflow")
    }
}
