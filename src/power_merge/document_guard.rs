//! Bound extraction and graph workspace without excluding datasets.
use super::{memory, sort, BuildConfig};
use crate::{try_parse_features_from_file, ColorInput};
use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::fmt;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub(super) struct SkippedDocument {
    pub source_document_id: usize,
    pub name: String,
    pub path: String,
    pub stage: String,
    pub estimated_workspace_bytes: u64,
    pub memory_budget_bytes: u64,
    pub reason: String,
}

#[derive(Debug)]
pub(super) struct DocumentTooLarge {
    pub stage: &'static str,
    pub estimated: u64,
    pub budget: u64,
}
impl fmt::Display for DocumentTooLarge {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} workspace estimate {} bytes exceeds document memory budget {} bytes",
            self.stage, self.estimated, self.budget
        )
    }
}
impl std::error::Error for DocumentTooLarge {}

#[derive(Clone, Copy)]
pub(super) struct Guard {
    pub enabled: bool,
    pub budget: u64,
}
impl Guard {
    pub fn new(_config: &BuildConfig, build_budget: u64) -> Self {
        Self {
            enabled: true,
            budget: build_budget,
        }
    }
    pub fn check(&self, stage: &'static str, estimated: u64) -> Result<()> {
        if self.enabled && estimated > self.budget {
            return Err(DocumentTooLarge {
                stage,
                estimated,
                budget: self.budget,
            }
            .into());
        }
        Ok(())
    }
    pub fn collect(&self, input: &ColorInput, config: &BuildConfig) -> Result<(Vec<u64>, u64)> {
        let mut keys = Vec::new();
        let mut occurrences = 0u64;
        // Reserve at most a quarter of the per-document budget for extraction
        // keys, leaving room for the parser and later union/payload work. Exact
        // reservation avoids a last Vec doubling beyond this allocation bound.
        let max_keys = (self.budget / 32).min(usize::MAX as u64) as usize;
        try_parse_features_from_file(
            &input.path,
            config.k,
            config.seed,
            config.feature_config,
            self.enabled.then_some(self.budget / 4),
            |key| {
                occurrences = occurrences.saturating_add(1);
                if self.enabled {
                    if keys.len() >= max_keys {
                        return Err(DocumentTooLarge {
                            stage: "feature-vector",
                            estimated: occurrences.saturating_mul(32),
                            budget: self.budget,
                        }
                        .into());
                    }
                    if keys.len() == keys.capacity() {
                        let next = keys.capacity().max(4096).saturating_mul(2).min(max_keys);
                        keys.try_reserve_exact(next - keys.len())
                            .context("reserving document feature vector")?;
                    }
                }
                keys.push(key);
                Ok(())
            },
        )?;
        sort::sort_unique(&mut keys);
        let slots = config.slots(keys.len() as f64)?;
        let input_resident = memory::feature_resident_bytes(occurrences);
        let workspace = memory::group_estimate(
            keys.len(),
            slots,
            config.hashes,
            config.width(1),
            input_resident,
        );
        memory::trace_memory(
            "document-extracted",
            serde_json::json!({
            "name":input.name,"occurrences":occurrences,"distinct":keys.len(),
            "key_capacity_bytes":(keys.capacity() as u64).saturating_mul(8),
            "key_resident_bytes":input_resident,
            "hashes":config.hashes,"slots":slots,
            "prediction_model":memory::prediction_kind(keys.len(),slots,config.hashes),
            "predicted_peak_rss_bytes":memory::group_prediction(keys.len(),slots,config.hashes,config.width(1),input_resident),
            "admission_workspace_bytes":workspace,
            "conservative_capacity_bound_bytes":memory::capacity_bound(keys.len(),slots,config.hashes,config.width(1),input_resident),
            "estimated_graph_workspace_bytes":workspace,"document_budget_bytes":self.budget}),
        );
        self.check("singleton-graph", workspace)?;
        Ok((keys, occurrences))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn oversized_documents_fail_with_a_typed_error() {
        let guard = Guard {
            enabled: true,
            budget: 100,
        };
        assert!(guard.check("test", 100).is_ok());
        let error = guard
            .check("test", 101)
            .unwrap_err()
            .context("outer context");
        let limit = error.downcast_ref::<DocumentTooLarge>().unwrap();
        assert_eq!(limit.budget, 100);
        assert_eq!(limit.estimated, 101);
    }
}
