//! Merge-only colored indexes with one slot-major bank per power of two.
mod abandoned;
mod graph;
mod planner;
mod sketch;
mod storage;

use crate::FeatureConfig;
use anyhow::{ensure, Context, Result};
use clap::{Args, ValueEnum};
pub use graph::{DegreeQueue, GraphIndex, IndexedDegreeHeap, PeelGraph, TypedPeelGraph};
use serde::{Deserialize, Serialize};
pub use storage::{build, info_command, query_command, query_command_with_recovery};

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
    /// Streaming uses lossless size/origin pruning; streaming-unpruned retains the original scan.
    #[arg(long = "merge-planner", value_enum, default_value_t = MergePlanner::Streaming)]
    pub planner: MergePlanner,
    /// Physical group order within each bank; does not change planning or query output order.
    #[arg(long = "merge-column-order", value_enum, default_value_t = ColumnOrder::Compression)]
    pub column_order: ColumnOrder,
    /// Number of full-minimum buckets in the mergeable sketch. Must be a power of two.
    #[arg(long = "merge-sketch-size", default_value_t = 4096)]
    pub sketch_size: usize,
    /// Maximum datasets sharing a ZOR payload (fingerprint + membership must fit 64 bits).
    #[arg(long = "merge-max-colors", default_value_t = 48)]
    pub max_colors: usize,
    /// Require this fraction of predicted bit savings before accepting a merge.
    #[arg(long = "merge-min-saving", default_value_t = 0.0)]
    pub min_saving: f64,
    /// Keep singleton power-of-two banks, for measuring sizing independently of merging.
    #[arg(long = "merge-disabled")]
    pub disabled: bool,
    /// Omit abandoned-key recovery storage and lookup (merged indexes only).
    /// Abandoned features may then be missed.
    #[arg(long, conflicts_with = "plain")]
    pub drop_abandoned: bool,
    /// Use temporary key/cell files instead of the default scratch-free two-pass build.
    /// The default keeps completed filters bit-packed in RAM until compression.
    #[arg(long, conflicts_with = "plain")]
    pub disk_backed: bool,
    /// Keep completed filters bit-packed in RAM until bank compression instead of
    /// staging them as bit-packed cell files next to the output. Needs no scratch
    /// disk, but RAM must then hold every completed filter.
    #[arg(long, conflicts_with_all = ["plain", "disk_backed"])]
    pub low_disk: bool,
}
impl Default for MergeArgs {
    fn default() -> Self {
        Self {
            power_of_two_merge: true,
            planner: MergePlanner::Streaming,
            column_order: ColumnOrder::Compression,
            sketch_size: 4096,
            max_colors: 48,
            min_saving: 0.0,
            disabled: false,
            drop_abandoned: false,
            disk_backed: false,
            low_disk: false,
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
    pub merge: MergeArgs,
}
impl BuildConfig {
    pub fn validate(&self) -> Result<()> {
        ensure!(
            [1, 8, 16].contains(&self.fingerprint_bits),
            "fingerprint bits must be 1, 8 or 16"
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
            .max(self.hashes.next_power_of_two());
        needed
            .checked_next_power_of_two()
            .context("power-of-two capacity overflow")
    }
}
