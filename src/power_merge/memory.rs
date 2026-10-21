//! Admission estimates, not a virtual-memory limit or an external sorter.
use anyhow::Result;
use std::{fs, path::Path};

const GIB: u64 = 1 << 30;

pub(super) fn available_bytes() -> Option<u64> {
    let mut available = fs::read_to_string("/proc/meminfo").ok().and_then(|text| {
        text.lines().find_map(|line| {
            line.strip_prefix("MemAvailable:")?
                .split_whitespace()
                .next()?
                .parse::<u64>()
                .ok()
                .map(|n| n.saturating_mul(1024))
        })
    });
    // Respect visible cgroup-v2 ancestor limits as well as host RAM. Some
    // systems hide the cgroup filesystem; in that case retain MemAvailable.
    if let Ok(cgroup) = fs::read_to_string("/proc/self/cgroup") {
        if let Some(relative) = cgroup.lines().find_map(|s| s.strip_prefix("0::/")) {
            let root = Path::new("/sys/fs/cgroup");
            let mut path = root.join(relative);
            while path.starts_with(root) {
                let read = |name| {
                    fs::read_to_string(path.join(name))
                        .ok()?
                        .trim()
                        .parse::<u64>()
                        .ok()
                };
                if let (Some(limit), Some(current)) = (read("memory.max"), read("memory.current")) {
                    let free = limit.saturating_sub(current);
                    available = Some(available.map_or(free, |n| n.min(free)));
                }
                if path == root || !path.pop() {
                    break;
                }
            }
        }
    }
    available
}

pub(super) fn budget(gib: Option<f64>) -> u64 {
    gib.map(|n| (n * GIB as f64) as u64)
        .unwrap_or_else(|| available_bytes().unwrap_or(8 * GIB).saturating_mul(4) / 5)
        .max(1)
}

pub(super) fn input_estimate(path: &Path) -> Result<u64> {
    // Extraction stores occurrences before deduplication. Compressed size is
    // only a heuristic (especially for very repetitive data), never a bound.
    let bytes = fs::metadata(path)?.len();
    Ok(bytes
        .saturating_mul(if crate::is_compressed_path(path) {
            64
        } else {
            24
        })
        .max(8 << 20))
}

pub(super) fn group_estimate(
    features: usize,
    slots: usize,
    hashes: usize,
    width: usize,
    input_peak: u64,
) -> u64 {
    let narrow = features <= u32::MAX as usize
        && slots <= u32::MAX as usize
        && features
            .checked_mul(hashes)
            .is_some_and(|n| n <= u32::MAX as usize);
    let id = if narrow { 4u64 } else { 8 };
    // Union keys/payload capacity, both incidence arrays, peel order, recovery;
    // then degrees, offsets/cursor, queue/heap and assigned/packed cells.
    let per_key = 32 + 2 * hashes as u64 * id + 2 * id + 16;
    let per_slot = 4
        + 3 * id
        + if narrow { 12 } else { 64 }
        + (width.div_ceil(8).next_power_of_two() + width.div_ceil(8)) as u64;
    (features as u64)
        .saturating_mul(per_key)
        .saturating_add((slots as u64).saturating_mul(per_slot))
        .saturating_add(input_peak)
        .max(8 << 20)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn estimates_account_for_wide_graph_and_do_not_overflow() {
        let narrow = group_estimate(1000, 2048, 8, 8, 0);
        let wide = group_estimate(536_870_912, 1 << 30, 8, 8, 0);
        assert!(wide > narrow);
        assert!(wide > 536_870_912u64 * 128);
        assert_eq!(
            group_estimate(usize::MAX, usize::MAX, 32, 64, u64::MAX),
            u64::MAX
        );
        assert_eq!(budget(Some(0.5)), GIB / 2);
    }
}
