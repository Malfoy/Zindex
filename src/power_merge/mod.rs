//! Merge-only colored indexes with one slot-major bank per power of two.
mod abandoned;
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
use serde::{Deserialize, Serialize};
/// Fixed segment size for newly built merged indexes, independent of filter size.
pub const DEFAULT_SEGMENT_LENGTH: usize = 64;
pub use storage::{
    build, info_command, query_command, query_command_with_recovery, query_records_command,
};
mod memory;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum MergePlanner {
    AllPairs,
    #[default]
    Streaming,
    /// Original streaming planner, retained as an unpruned benchmark baseline.
    StreamingUnpruned,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ColumnOrder {
    /// Keep width sections contiguous, then use input dataset IDs.
    Width,
    /// Within width sections, cluster columns by sampled filter cell contents.
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
            Self::AllPairs => "all-pairs",
            Self::Streaming => "streaming",
            Self::StreamingUnpruned => "streaming-unpruned",
        }
    }

    pub fn is_streaming(self) -> bool {
        matches!(self, Self::Streaming | Self::StreamingUnpruned)
    }
}

#[derive(Args, Clone, Debug)]
pub struct MergeArgs {
    /// Compatibility flag: merged power-of-two banks are already the default.
    #[arg(long, action = clap::ArgAction::SetTrue, default_value_t = true)]
    pub power_of_two_merge: bool,
    /// Streaming uses lossless cost-bound pruning; streaming-unpruned retains the exhaustive scan.
    #[arg(long = "merge-planner", value_enum, default_value_t = MergePlanner::Streaming)]
    pub planner: MergePlanner,
    /// Physical group order within each bank; does not change planning or query output order.
    #[arg(long = "merge-column-order", value_enum, default_value_t = ColumnOrder::Compression)]
    pub column_order: ColumnOrder,
    /// Number of SimdSketch buckets (8-bit comparison fingerprints + 32-bit union minima).
    /// Must be a power of two.
    #[arg(long = "merge-sketch-size", default_value_t = 128)]
    pub sketch_size: usize,
    /// Maximum datasets sharing a ZOR payload (fingerprint + membership must fit 64 bits).
    #[arg(long = "merge-max-colors", default_value_t = 48)]
    pub max_colors: usize,
    /// Require this fraction of predicted bit savings before accepting a merge.
    #[arg(long = "merge-min-saving", default_value_t = 0.0)]
    pub min_saving: f64,
    /// Plan each merged group as if its estimated union were this fraction larger, so sketch noise does not push
    /// the exact union past the planned power of two at build time (which doubles the group's slots).
    #[arg(long = "merge-plan-margin", default_value_t = 0.0)]
    pub plan_margin: f64,
    /// Minimum original documents per subindex; promote smaller filters, never fold.
    /// Collections smaller than this use one subindex. Set 1 for the legacy layout.
    #[arg(long, default_value_t = 128)]
    pub min_subindex_documents: usize,
    /// Maximum subindex count; promote smaller filters, never fold. Omit for unlimited.
    #[arg(long)]
    pub max_subindexes: Option<usize>,
    /// Keep singleton power-of-two banks, for measuring sizing independently of merging.
    #[arg(long = "merge-disabled")]
    pub disabled: bool,
    /// Omit abandoned-key recovery storage and lookup (merged indexes only).
    /// Abandoned features may then be missed.
    #[arg(long, conflicts_with = "plain")]
    pub drop_abandoned: bool,
    /// Use temporary feature/cell/abandoned-key files instead of the scratch-free two-pass build.
    /// The default keeps completed filters bit-packed in RAM until compression.
    #[arg(long, conflicts_with = "plain")]
    pub disk_backed: bool,
    /// Keep completed filters and abandoned-key records in RAM (default).
    /// Set --low-disk=false to stage both on disk, with a bounded abandoned-key merge.
    #[arg(long, action = clap::ArgAction::Set, num_args = 0..=1,
        require_equals = true, default_missing_value = "true", default_value_t = true,
        conflicts_with_all = ["plain", "disk_backed"])]
    pub low_disk: bool,

    /// Estimated concurrent build-workspace budget in GiB. Default: 80% of
    /// available RAM (including visible cgroup limits). Oversized jobs run alone.
    /// This is admission control, not an RSS cap; retained filters/recovery and
    /// a single oversized genome still need enough memory. Feature sorting stays in RAM.
    #[arg(long)]
    pub build_memory_gib: Option<f64>,
}
impl Default for MergeArgs {
    fn default() -> Self {
        Self {
            power_of_two_merge: true,
            planner: MergePlanner::Streaming,
            column_order: ColumnOrder::Compression,
            sketch_size: 128,
            max_colors: 48,
            min_saving: 0.0,
            plan_margin: 0.0,
            min_subindex_documents: 128,
            max_subindexes: None,
            disabled: false,
            drop_abandoned: false,
            disk_backed: false,
            low_disk: true,
            build_memory_gib: None,
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
    pub minimum_slots: Option<usize>,
    /// Fixed segment size. None is retained for legacy library/layout fixtures.
    pub segment_length: Option<usize>,
    pub merge: MergeArgs,
}
impl BuildConfig {
    pub fn validate(&self) -> Result<()> {
        ensure!(
            self.merge
                .build_memory_gib
                .is_none_or(|v| v.is_finite() && v > 0.0 && v * 1073741824.0 < u64::MAX as f64),
            "--build-memory-gib must be a positive finite GiB value"
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
            crate::IndexMode::Kmers | crate::IndexMode::Modimizers => ensure!(
                (1..=31).contains(&self.k),
                "k must be in 1..=31 for k-mers/modimizers"
            ),
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
            (1..=64).contains(&self.merge.max_colors),
            "--merge-max-colors must be in 1..=64"
        );
        ensure!(
            self.merge.min_saving.is_finite() && (0.0..1.0).contains(&self.merge.min_saving),
            "--merge-min-saving must be in [0, 1)"
        );
        ensure!(
            self.merge.plan_margin.is_finite() && (0.0..0.5).contains(&self.merge.plan_margin),
            "--merge-plan-margin must be in [0, 0.5)"
        );
        ensure!(
            self.merge.max_colors + self.fingerprint_bits <= 64,
            "fingerprint bits + --merge-max-colors must be <= 64"
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
    pub fn width(&self, colors: usize) -> usize {
        self.fingerprint_bits + if colors == 1 { 0 } else { colors }
    }
    pub fn slots(&self, cardinality: f64) -> Result<usize> {
        let needed = (cardinality.max(1.0) * self.slot_scale).ceil();
        ensure!(
            needed.is_finite() && needed < usize::MAX as f64,
            "filter capacity overflow"
        );
        let needed = (needed as usize)
            .max(self.minimum_slots.unwrap_or(0))
            .max(self.hashes.next_power_of_two())
            .max(
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
