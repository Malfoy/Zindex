//! Minimum-cost contiguous partition of exact power-of-two size buckets.
//! Each interval is promoted to its largest size; no documents/groups are merged.
use anyhow::{ensure, Context, Result};

/// Buckets are (slots, original document count, sum of payload widths).
/// Costs include the byte padding of the final stacked rows.
pub(super) fn promotion_targets(
    buckets: &[(usize, usize, usize)],
    minimum: usize,
) -> Result<Vec<usize>> {
    ensure!(
        minimum > 0,
        "minimum subindex document count must be positive"
    );
    ensure!(
        buckets.windows(2).all(|b| b[0].0 < b[1].0),
        "unsorted slot buckets"
    );
    ensure!(
        buckets
            .iter()
            .all(|&(slots, docs, width)| slots.is_power_of_two() && docs > 0 && width > 0),
        "invalid slot bucket"
    );
    if minimum == 1 || buckets.iter().all(|b| b.1 >= minimum) {
        return Ok(buckets.iter().map(|b| b.0).collect());
    }
    let total: usize = buckets.iter().map(|b| b.1).sum();
    let minimum = minimum.min(total);
    let n = buckets.len();
    let mut best = vec![None::<(u128, usize, usize)>; n + 1];
    best[0] = Some((0, 0, 0));
    for end in 1..=n {
        let mut docs = 0usize;
        let mut width = 0usize;
        for start in (0..end).rev() {
            docs = docs
                .checked_add(buckets[start].1)
                .context("document count overflow")?;
            width = width
                .checked_add(buckets[start].2)
                .context("row width overflow")?;
            if docs < minimum {
                continue;
            }
            if let Some((cost, banks, _)) = best[start] {
                let candidate = (
                    cost + buckets[end - 1].0 as u128 * width.div_ceil(8) as u128 * 8,
                    banks + 1,
                    start,
                );
                if best[end].is_none_or(|old| candidate < old) {
                    best[end] = Some(candidate);
                }
            }
        }
    }
    let mut targets = vec![0; n];
    let mut end = n;
    while end > 0 {
        let (_, _, start) = best[end].context("no feasible minimum-document layout")?;
        targets[start..end].fill(buckets[end - 1].0);
        end = start;
    }
    Ok(targets)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn meets_document_minimum_without_folding_or_modifying_full_layouts() {
        let buckets = [(8, 2, 16), (16, 130, 1040), (32, 10, 80), (64, 140, 1120)];
        assert_eq!(promotion_targets(&buckets, 1).unwrap(), [8, 16, 32, 64]);
        let targets = promotion_targets(&buckets, 128).unwrap();
        let mut counts = std::collections::BTreeMap::<usize, usize>::new();
        for (&(old, docs, _), &target) in buckets.iter().zip(&targets) {
            assert!(target >= old && target.is_power_of_two());
            *counts.entry(target).or_default() += docs;
        }
        assert!(counts.values().all(|&count| count >= 128));
        assert_eq!(
            promotion_targets(&[(8, 2, 16), (64, 3, 24)], 128).unwrap(),
            [64, 64]
        );
        assert_eq!(
            promotion_targets(&[(8, 128, 1024), (64, 128, 1024)], 128).unwrap(),
            [8, 64]
        );
        assert!(promotion_targets(&buckets, 0).is_err());
    }

    #[test]
    fn dynamic_program_matches_exhaustive_contiguous_partitions() {
        for seed in 0..64 {
            let buckets: Vec<_> = (0..7)
                .map(|i| {
                    let docs = (crate::splitmix64(seed + i) % 150 + 1) as usize;
                    (8 << i, docs, docs * 9)
                })
                .collect();
            let targets = promotion_targets(&buckets, 128).unwrap();
            let cost = |targets: &[usize]| {
                let mut widths = std::collections::BTreeMap::<usize, usize>::new();
                for (&(_, _, width), &slots) in buckets.iter().zip(targets) {
                    *widths.entry(slots).or_default() += width;
                }
                widths
                    .into_iter()
                    .map(|(slots, width)| slots as u128 * width.div_ceil(8) as u128 * 8)
                    .sum::<u128>()
            };
            let mut expected = u128::MAX;
            for cuts in 0..1u32 << 6 {
                let mut trial = Vec::new();
                let mut start = 0;
                let mut valid = true;
                for end in 1..=7 {
                    if end == 7 || cuts & (1 << (end - 1)) != 0 {
                        valid &= buckets[start..end].iter().map(|b| b.1).sum::<usize>() >= 128;
                        trial.resize(end, buckets[end - 1].0);
                        start = end;
                    }
                }
                if valid {
                    expected = expected.min(cost(&trial));
                }
            }
            assert_eq!(cost(&targets), expected);
        }
    }
}
