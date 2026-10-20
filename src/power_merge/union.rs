//! Bounded, balanced unions of sorted feature sets and their original color bits.
#[derive(Default)]
struct Run {
    keys: Vec<u64>,
    origins: Vec<u64>,
}

impl Run {
    fn merge(mut self, mut other: Self) -> Self {
        if self.keys.len() < other.keys.len() {
            std::mem::swap(&mut self, &mut other);
        }
        let (mut i, mut j) = (self.keys.len(), other.keys.len());
        let total = i.checked_add(j).expect("union size overflow");
        self.keys.resize(total, 0);
        self.origins.resize(total, 0);
        let mut write = total;
        while i != 0 || j != 0 {
            let (key, mask) = if j == 0 || (i != 0 && self.keys[i - 1] > other.keys[j - 1]) {
                i -= 1;
                (self.keys[i], self.origins[i])
            } else if i == 0 || other.keys[j - 1] > self.keys[i - 1] {
                j -= 1;
                (other.keys[j], other.origins[j])
            } else {
                i -= 1;
                j -= 1;
                (self.keys[i], self.origins[i] | other.origins[j])
            };
            write -= 1;
            self.keys[write] = key;
            self.origins[write] = mask;
        }
        self.keys.copy_within(write..total, 0);
        self.origins.copy_within(write..total, 0);
        self.keys.truncate(total - write);
        self.origins.truncate(total - write);
        self
    }
}

#[derive(Default)]
pub(super) struct Union {
    levels: Vec<Option<Run>>,
}

impl Union {
    pub(super) fn push(&mut self, keys: Vec<u64>, mask: u64) {
        let mut run = Run {
            origins: vec![mask; keys.len()],
            keys,
        };
        let mut level = 0;
        loop {
            if level == self.levels.len() {
                self.levels.push(None);
            }
            if let Some(previous) = self.levels[level].take() {
                run = previous.merge(run);
                level += 1;
            } else {
                self.levels[level] = Some(run);
                break;
            }
        }
        // Similar genomes have nearly constant union size: do not retain one
        // full copy per level. Bound retained entries to twice the largest run.
        let total: usize = self.levels.iter().flatten().map(|r| r.keys.len()).sum();
        let largest = self
            .levels
            .iter()
            .flatten()
            .map(|r| r.keys.len())
            .max()
            .unwrap_or(0);
        if total > largest.saturating_mul(2) {
            let last = self.levels.iter().rposition(Option::is_some).unwrap();
            let run = self.collapse();
            self.levels[last] = Some(run);
        }
    }

    fn collapse(&mut self) -> Run {
        let mut runs: Vec<_> = self.levels.iter_mut().filter_map(Option::take).collect();
        while runs.len() > 1 {
            runs.sort_unstable_by_key(|r| std::cmp::Reverse(r.keys.len()));
            let a = runs.pop().unwrap();
            let b = runs.pop().unwrap();
            runs.push(a.merge(b));
        }
        runs.pop().unwrap_or_default()
    }

    pub(super) fn finish(mut self) -> (Vec<u64>, Vec<u64>) {
        let run = self.collapse();
        (run.keys, run.origins)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn balanced_union_keeps_every_origin_for_disjoint_and_overlapping_sets() {
        for overlap in [0, 7, 999] {
            let mut union = Union::default();
            let mut expected = std::collections::BTreeMap::<u64, u64>::new();
            for bit in 0..64 {
                let mut keys: Vec<_> = (0..1000)
                    .map(|i| crate::splitmix64(i + (1000 - overlap) * bit))
                    .collect();
                keys.sort_unstable();
                keys.dedup();
                for &key in &keys {
                    *expected.entry(key).or_default() |= 1u64 << bit;
                }
                union.push(keys, 1u64 << bit);
                let sizes: Vec<_> = union
                    .levels
                    .iter()
                    .flatten()
                    .map(|r| r.keys.len())
                    .collect();
                assert!(
                    sizes.iter().sum::<usize>() <= 2 * sizes.iter().max().copied().unwrap_or(0)
                );
            }
            let (keys, origins) = union.finish();
            assert_eq!(keys, expected.keys().copied().collect::<Vec<_>>());
            assert_eq!(origins, expected.values().copied().collect::<Vec<_>>());
        }
        assert_eq!(Union::default().finish(), (vec![], vec![]));
    }
}
