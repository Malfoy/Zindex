use super::{sketch::Sketch, BuildConfig, MergePlanner};
use anyhow::{Context, Result};
use rayon::prelude::*;
use std::{
    cmp::Ordering,
    collections::{BTreeMap, BinaryHeap},
};

#[derive(Clone, Debug)]
pub struct Group {
    pub members: Vec<usize>,
    pub cardinality: f64,
    pub slots: usize,
    pub width: usize,
    pub sketch: Sketch,
}

pub struct Plan {
    pub groups: Vec<Group>,
    /// Candidate pair evaluations, including pairs rejected by the origin limit.
    pub pair_evaluations: u64,
    /// Comparisons an unpruned planner would have considered (including origin rejects).
    pub possible_pair_evaluations: u64,
}

#[derive(Clone, Debug)]
struct Candidate {
    saving: u128,
    a: usize,
    b: usize,
    cardinality: f64,
    slots: usize,
}
impl PartialEq for Candidate {
    fn eq(&self, b: &Self) -> bool {
        (self.saving, self.a, self.b) == (b.saving, b.a, b.b)
    }
}
impl Eq for Candidate {}
impl PartialOrd for Candidate {
    fn partial_cmp(&self, b: &Self) -> Option<Ordering> {
        Some(self.cmp(b))
    }
}
impl Ord for Candidate {
    fn cmp(&self, b: &Self) -> Ordering {
        self.saving
            .cmp(&b.saving)
            .then_with(|| b.a.cmp(&self.a))
            .then_with(|| b.b.cmp(&self.b))
    }
}

pub fn cost(slots: usize, width: usize) -> u128 {
    slots as u128 * width as u128
}

pub fn singleton(count: usize, sketch: Sketch, id: usize, config: &BuildConfig) -> Result<Group> {
    Ok(Group {
        members: vec![id],
        cardinality: count as f64,
        slots: config.slots(count as f64)?,
        width: config.width(1),
        sketch,
    })
}

fn candidate(
    a: usize,
    b: usize,
    ga: &Group,
    gb: &Group,
    config: &BuildConfig,
) -> Result<Option<Candidate>> {
    let count = ga.members.len() + gb.members.len();
    if count > config.merge.max_colors {
        return Ok(None);
    }
    let j = ga.sketch.jaccard(&gb.sketch);
    let cardinality = ((ga.cardinality + gb.cardinality) / (1.0 + j)).clamp(
        ga.cardinality.max(gb.cardinality),
        ga.cardinality + gb.cardinality,
    );
    // Size the merged group with a safety margin; the chained estimate itself stays uninflated.
    let slots = config.slots(cardinality * (1.0 + config.merge.plan_margin))?;
    let old = cost(ga.slots, ga.width) + cost(gb.slots, gb.width);
    let new = cost(slots, config.width(count));
    if new >= old || (old - new) as f64 <= old as f64 * config.merge.min_saving {
        return Ok(None);
    }
    Ok(Some(Candidate {
        saving: old - new,
        a,
        b,
        cardinality,
        slots,
    }))
}

fn merge(ga: Group, gb: Group, c: &Candidate, config: &BuildConfig) -> Group {
    let mut members = ga.members;
    members.extend(gb.members);
    members.sort_unstable();
    Group {
        width: config.width(members.len()),
        members,
        cardinality: c.cardinality,
        slots: c.slots,
        sketch: ga.sketch.union(&gb.sketch),
    }
}

/// Optimistic cost: union capacity cannot be smaller than either input's
/// capacity. If even this lower bound cannot save bits, Jaccard cannot help.
/// This is valid for singleton arrivals AND arrivals containing multiple origins.
fn bucket_can_save(slots: usize, origins: usize, incoming: &Group, config: &BuildConfig) -> bool {
    let count = origins + incoming.members.len();
    if count > config.merge.max_colors {
        return false;
    }
    let old = cost(slots, config.width(origins)) + cost(incoming.slots, incoming.width);
    let lower_new = cost(slots.max(incoming.slots), config.width(count));
    lower_new < old && (old - lower_new) as f64 > old as f64 * config.merge.min_saving
}

/// Incremental planning state for one batch. Only sketches/metadata are merged;
/// no union keys or filters are constructed here. Previously indexed batches
/// must not be inserted into this state when planning an update.
#[derive(Default)]
pub struct StreamingPlanner {
    groups: Vec<Group>,
    pair_evaluations: u64,
    possible_pair_evaluations: u64,
    // Sub-bucket by origin count, so both width and slot cost are constant.
    // IDs never move; ties therefore match the original FOF-ordered scan.
    buckets: BTreeMap<(u32, usize), Vec<usize>>,
    bucket_positions: Vec<usize>,
}
/// Smallest slice of candidate evaluations handed to a pool thread: a few thousand sketch comparisons are
/// about a millisecond, so smaller scans run on the calling thread instead of paying for thread hand-offs.
const CANDIDATES_PER_TASK: usize = 1024;

impl StreamingPlanner {
    fn add_to_bucket(&mut self, id: usize) {
        let group = &self.groups[id];
        let key = (group.slots.trailing_zeros(), group.members.len());
        let bucket = self.buckets.entry(key).or_default();
        self.bucket_positions[id] = bucket.len();
        bucket.push(id);
    }

    fn remove_from_bucket(&mut self, id: usize) {
        let group = &self.groups[id];
        let key = (group.slots.trailing_zeros(), group.members.len());
        let bucket = self
            .buckets
            .get_mut(&key)
            .expect("missing streaming bucket");
        let position = self.bucket_positions[id];
        assert_eq!(bucket.swap_remove(position), id);
        if position < bucket.len() {
            self.bucket_positions[bucket[position]] = position;
        }
        if bucket.is_empty() {
            self.buckets.remove(&key);
        }
    }

    pub fn insert(&mut self, incoming: Group, config: &BuildConfig) -> Result<()> {
        let unpruned = config.merge.planner == MergePlanner::StreamingUnpruned;
        if config.merge.disabled {
            self.groups.push(incoming);
            self.bucket_positions.push(0);
            if !unpruned {
                self.add_to_bucket(self.groups.len() - 1);
            }
            return Ok(());
        }
        self.possible_pair_evaluations += self.groups.len() as u64;
        // Incoming groups are processed in FOF order. Independent comparisons
        // can run in parallel; Candidate's total order makes ties deterministic.
        let best = if unpruned {
            // Keep the previous algorithm as a selectable, unpruned baseline.
            self.pair_evaluations += self.groups.len() as u64;
            self.groups
                .par_iter()
                .with_min_len(CANDIDATES_PER_TASK)
                .enumerate()
                .map(|(id, group)| candidate(id, self.groups.len(), group, &incoming, config))
                .try_reduce(|| None, |a, b| Ok(a.max(b)))?
        } else {
            let eligible: Vec<usize> = self
                .buckets
                .iter()
                .filter(|(&(exponent, origins), _)| {
                    bucket_can_save(1usize << exponent, origins, &incoming, config)
                })
                .flat_map(|(_, ids)| ids.iter().copied())
                .collect();
            self.pair_evaluations += eligible.len() as u64;
            eligible
                .par_iter()
                .with_min_len(CANDIDATES_PER_TASK)
                .map(|&id| candidate(id, self.groups.len(), &self.groups[id], &incoming, config))
                .try_reduce(|| None, |a, b| Ok(a.max(b)))?
        };
        if let Some(c) = best {
            // Do not search for additional old-group merges after insertion.
            if unpruned {
                let previous = self.groups.remove(c.a);
                self.groups
                    .insert(c.a, merge(previous, incoming, &c, config));
            } else {
                self.remove_from_bucket(c.a);
                let group = &mut self.groups[c.a];
                group.members.extend(incoming.members);
                group.members.sort_unstable();
                group.width = config.width(group.members.len());
                group.cardinality = c.cardinality;
                group.slots = c.slots;
                group.sketch = group.sketch.union(&incoming.sketch);
                self.add_to_bucket(c.a);
            }
        } else {
            self.groups.push(incoming);
            self.bucket_positions.push(0);
            if !unpruned {
                self.add_to_bucket(self.groups.len() - 1);
            }
        }
        Ok(())
    }

    pub fn finish(self) -> Plan {
        Plan {
            groups: self.groups,
            pair_evaluations: self.pair_evaluations,
            possible_pair_evaluations: self.possible_pair_evaluations,
        }
    }
}

/// Maximum-saving agglomeration. Only candidates involving the new group need
/// recalculation; heap entries involving retired groups are ignored. Pairwise
/// sketch comparisons remain O(N²), and the heap can also occupy O(N²) space.
pub fn plan(initial: Vec<Group>, config: &BuildConfig) -> Result<Plan> {
    if config.merge.planner.is_streaming() {
        let mut state = StreamingPlanner::default();
        for group in initial {
            state.insert(group, config)?;
        }
        return Ok(state.finish());
    }
    if config.merge.disabled {
        return Ok(Plan {
            groups: initial,
            pair_evaluations: 0,
            possible_pair_evaluations: 0,
        });
    }
    let mut pair_evaluations = initial.len() as u64 * initial.len().saturating_sub(1) as u64 / 2;
    let rows = (0..initial.len())
        .into_par_iter()
        .map(|a| {
            let mut row = Vec::new();
            for b in a + 1..initial.len() {
                if let Some(c) = candidate(a, b, &initial[a], &initial[b], config)? {
                    row.push(c);
                }
            }
            Ok(row)
        })
        .collect::<Result<Vec<_>>>()?;
    let mut heap = BinaryHeap::new();
    for row in rows {
        heap.extend(row);
    }
    let mut groups: Vec<_> = initial.into_iter().map(Some).collect();
    while let Some(c) = heap.pop() {
        if groups[c.a].is_none() || groups[c.b].is_none() {
            continue;
        }
        let ga = groups[c.a].take().context("missing merge group")?;
        let gb = groups[c.b].take().context("missing merge group")?;
        let group = merge(ga, gb, &c, config);
        let id = groups.len();
        for (other_id, other) in groups.iter().enumerate() {
            if let Some(other) = other {
                pair_evaluations += 1;
                if let Some(next) = candidate(other_id, id, other, &group, config)? {
                    heap.push(next);
                }
            }
        }
        groups.push(Some(group));
    }
    let mut result: Vec<_> = groups.into_iter().flatten().collect();
    result.sort_by_key(|g| g.members[0]);
    Ok(Plan {
        groups: result,
        pair_evaluations,
        possible_pair_evaluations: pair_evaluations,
    })
}

#[cfg(test)]
mod tests {
    use super::super::{Backend, MergeArgs};
    use super::*;
    use crate::FeatureConfig;
    fn config(backend: Backend) -> BuildConfig {
        BuildConfig {
            backend,
            k: 31,
            hashes: 4,
            seed: 69,
            fingerprint_bits: 16,
            feature_config: FeatureConfig::legacy_kmers(),
            slot_scale: 1.0,
            minimum_slots: None,
            // Exact sizes: these tests check the cost rule itself, not the sketch-noise margin.
            merge: MergeArgs {
                plan_margin: 0.0,
                ..MergeArgs::default()
            },
        }
    }
    #[test]
    fn rounded_cost_controls_merges() {
        let keys: Vec<_> = (0..1000).collect();
        {
            let backend = Backend::Zor;
            let c = config(backend);
            let initial: Vec<_> = (0..8)
                .map(|i| singleton(keys.len(), Sketch::new(&keys, 4096, 69), i, &c).unwrap())
                .collect();
            let old: u128 = initial.iter().map(|g| cost(g.slots, g.width)).sum();
            let groups = plan(initial, &c).unwrap().groups;
            let new: u128 = groups.iter().map(|g| cost(g.slots, g.width)).sum();
            assert!(groups.iter().all(|g| g.slots.is_power_of_two()));
            assert_eq!(groups.len(), 1);
            assert!(new < old);
        }
    }

    fn initial(sets: &[Vec<u64>], c: &BuildConfig) -> Vec<Group> {
        sets.iter()
            .enumerate()
            .map(|(id, keys)| {
                singleton(
                    keys.len(),
                    Sketch::new(keys, c.merge.sketch_size, c.seed),
                    id,
                    c,
                )
                .unwrap()
            })
            .collect()
    }

    #[test]
    fn streaming_batch_and_incremental_insert_are_identical() {
        let mut c = config(Backend::Zor);
        c.merge.planner = MergePlanner::Streaming;
        let sets: Vec<_> = (0..12)
            .map(|i| (i * 200..i * 200 + 1000).collect())
            .collect();
        let groups = initial(&sets, &c);
        let batch = plan(groups.clone(), &c).unwrap();
        let mut state = StreamingPlanner::default();
        // A sketch can arrive later: state is retained, not a whole pair table.
        for group in groups {
            state.insert(group, &c).unwrap();
        }
        let incremental = state.finish();
        assert_eq!(batch.pair_evaluations, incremental.pair_evaluations);
        assert_eq!(batch.groups.len(), incremental.groups.len());
        for (a, b) in batch.groups.iter().zip(incremental.groups) {
            assert_eq!(a.members, b.members);
            assert_eq!(a.slots, b.slots);
            assert_eq!(a.cardinality, b.cardinality);
            assert_eq!(a.sketch.jaccard(&b.sketch), 1.0);
        }
    }

    #[test]
    fn streaming_comparison_count_tracks_current_groups() {
        let mut c = config(Backend::Zor);
        c.merge.planner = MergePlanner::StreamingUnpruned;
        let identical = vec![(0..1000).collect(); 12];
        let merged = plan(initial(&identical, &c), &c).unwrap();
        assert_eq!(merged.groups.len(), 1);
        assert_eq!(merged.pair_evaluations, 11);
        let disjoint: Vec<_> = (0..12)
            .map(|i| (i * 1000..(i + 1) * 1000).collect())
            .collect();
        let separate = plan(initial(&disjoint, &c), &c).unwrap();
        assert_eq!(separate.groups.len(), 12);
        assert_eq!(separate.pair_evaluations, 66);
        c.merge.disabled = true;
        let disabled = plan(initial(&identical, &c), &c).unwrap();
        assert_eq!(disabled.groups.len(), 12);
        assert_eq!(disabled.pair_evaluations, 0);
    }

    #[test]
    fn streaming_respects_origin_limit_and_deterministic_ties() {
        let mut c = config(Backend::Zor);
        c.merge.planner = MergePlanner::Streaming;
        c.merge.max_colors = 2;
        let identical = vec![(0..1000).collect(); 5];
        let result = plan(initial(&identical, &c), &c).unwrap();
        assert_eq!(
            result
                .groups
                .iter()
                .map(|g| g.members.clone())
                .collect::<Vec<_>>(),
            [vec![0, 1], vec![2, 3], vec![4]]
        );
        assert!(result.groups.iter().all(|g| g.slots.is_power_of_two()));
    }

    #[test]
    fn streaming_chooses_largest_saving_not_first_eligible_group() {
        let mut c = config(Backend::Zor);
        c.merge.planner = MergePlanner::Streaming;
        let small = (0..200).collect::<Vec<_>>();
        let big = (0..1000).collect::<Vec<_>>();
        // Establish two groups independently, then enable merging for the third
        // arrival: it must choose the most profitable group, not the first one.
        c.merge.disabled = true;
        let mut state = StreamingPlanner::default();
        for g in initial(&[small.clone(), big.clone()], &c) {
            state.insert(g, &c).unwrap();
        }
        c.merge.disabled = false;
        let incoming = singleton(
            big.len(),
            Sketch::new(&big, c.merge.sketch_size, c.seed),
            2,
            &c,
        )
        .unwrap();
        state.insert(incoming, &c).unwrap();
        let result = state.finish();
        assert_eq!(
            result
                .groups
                .iter()
                .map(|g| g.members.clone())
                .collect::<Vec<_>>(),
            [vec![0], vec![1, 2]]
        );
    }

    #[test]
    fn streaming_never_revisits_pairs_of_existing_groups() {
        let mut c = config(Backend::Zor);
        c.merge.planner = MergePlanner::Streaming;
        let old = (0..1000).collect::<Vec<_>>();
        let incoming = (1000..2000).collect::<Vec<_>>();
        c.merge.disabled = true;
        let mut state = StreamingPlanner::default();
        for g in initial(&[old.clone(), old], &c) {
            state.insert(g, &c).unwrap();
        }
        c.merge.disabled = false;
        let g = singleton(
            incoming.len(),
            Sketch::new(&incoming, c.merge.sketch_size, c.seed),
            2,
            &c,
        )
        .unwrap();
        state.insert(g, &c).unwrap();
        let result = state.finish();
        // The old pair is profitable, but only arrival-vs-current pairs are considered.
        assert_eq!(result.groups.len(), 3);
        assert_eq!(result.pair_evaluations, 2);
    }

    #[test]
    fn streaming_parallel_comparisons_have_thread_independent_ties() {
        let mut c = config(Backend::Zor);
        c.merge.planner = MergePlanner::Streaming;
        c.merge.max_colors = 3;
        let sets = vec![(0..1000).collect(); 16];
        let run = |threads| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| {
                    plan(initial(&sets, &c), &c)
                        .unwrap()
                        .groups
                        .into_iter()
                        .map(|g| g.members)
                        .collect::<Vec<_>>()
                })
        };
        assert_eq!(run(1), run(4));
    }

    fn assert_same_groups(a: &[Group], b: &[Group]) {
        assert_eq!(a.len(), b.len());
        for (a, b) in a.iter().zip(b) {
            assert_eq!(a.members, b.members);
            assert_eq!(a.slots, b.slots);
            assert_eq!(a.width, b.width);
            assert_eq!(a.cardinality.to_bits(), b.cardinality.to_bits());
            assert_eq!(a.sketch, b.sketch);
        }
    }

    #[test]
    fn pruned_streaming_matches_original_across_widths_rounding_limits_and_ties() {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap();
        pool.install(|| {
            for fingerprint_bits in [1, 8, 16] {
                for max_colors in [2, 5, 16] {
                    for min_saving in [0.0, 0.05, 0.5] {
                        for minimum_slots in [None, Some(256)] {
                            let mut c = config(Backend::Zor);
                            c.fingerprint_bits = fingerprint_bits;
                            c.merge.max_colors = max_colors;
                            c.merge.min_saving = min_saving;
                            c.minimum_slots = minimum_slots;
                            c.merge.sketch_size = 64;
                            let sets: Vec<Vec<u64>> = (0..60)
                                .map(|i| {
                                    let h = crate::splitmix64(i);
                                    let count = 1usize << ((h % 11) as u32);
                                    // Overlap, disjoint families, repeated identical sets, and
                                    // cardinalities exactly at power-of-two boundaries.
                                    let start = ((h >> 8) % 4) * 4096;
                                    (start..start + count as u64).collect()
                                })
                                .collect();
                            let original = initial(&sets, &c);
                            c.merge.planner = MergePlanner::StreamingUnpruned;
                            let old = plan(original.clone(), &c).unwrap();
                            c.merge.planner = MergePlanner::Streaming;
                            let new = plan(original, &c).unwrap();
                            assert_same_groups(&old.groups, &new.groups);
                            assert_eq!(new.possible_pair_evaluations, old.pair_evaluations);
                            assert!(new.pair_evaluations <= old.pair_evaluations);
                        }
                    }
                }
            }
        });
    }

    #[test]
    fn bucket_bound_handles_singletons_exact_equality_and_origin_limits() {
        let c = config(Backend::Zor);
        let keys: Vec<_> = (0..128).collect();
        let incoming = singleton(keys.len(), Sketch::new(&keys, 64, c.seed), 0, &c).unwrap();
        assert_eq!(incoming.slots, 128);
        assert!(!bucket_can_save(2048, 2, &incoming, &c)); // 16*S == L
        assert!(bucket_can_save(1024, 2, &incoming, &c));
        assert!(!bucket_can_save(1024, 1, &incoming, &c)); // first merge adds TWO bits
        assert!(!bucket_can_save(128, c.merge.max_colors, &incoming, &c));
        assert!(!bucket_can_save(
            4,
            2,
            &singleton(
                4096,
                Sketch::new(&(0..4096).collect::<Vec<_>>(), 64, c.seed),
                1,
                &c
            )
            .unwrap(),
            &c
        ));
    }

    #[test]
    fn full_origin_groups_are_not_compared_but_larger_future_arrivals_remain_eligible() {
        let mut c = config(Backend::Zor);
        c.merge.max_colors = 2;
        let sets = vec![(0..1000).collect(); 12];
        let new = plan(initial(&sets, &c), &c).unwrap();
        c.merge.planner = MergePlanner::StreamingUnpruned;
        let old = plan(initial(&sets, &c), &c).unwrap();
        assert_same_groups(&old.groups, &new.groups);
        assert_eq!(new.pair_evaluations, 6);
        assert!(new.pair_evaluations < old.pair_evaluations);

        // A big group skipped for a small arrival must not be permanently retired.
        c.merge.planner = MergePlanner::Streaming;
        c.merge.max_colors = 48;
        let big: Vec<_> = (0..4096).collect();
        let small: Vec<_> = (0..64).collect();
        let sets = [big.clone(), small, big];
        let new = plan(initial(&sets, &c), &c).unwrap();
        c.merge.planner = MergePlanner::StreamingUnpruned;
        let old = plan(initial(&sets, &c), &c).unwrap();
        assert_same_groups(&old.groups, &new.groups);
        assert!(new.groups.iter().any(|g| g.members == [0, 2]));
    }

    #[test]
    fn pruning_is_valid_for_multi_origin_arrivals_too() {
        let mut c = config(Backend::Zor);
        c.merge.sketch_size = 64;
        let sets: Vec<Vec<u64>> = (0..12).map(|i| (i * 50..i * 50 + 1000).collect()).collect();
        let mut groups = initial(&sets, &c);
        for (id, group) in groups.iter_mut().enumerate() {
            group.members = vec![id * 2, id * 2 + 1];
            group.width = c.width(2);
        }
        c.merge.planner = MergePlanner::StreamingUnpruned;
        let old = plan(groups.clone(), &c).unwrap();
        c.merge.planner = MergePlanner::Streaming;
        let new = plan(groups, &c).unwrap();
        assert_same_groups(&old.groups, &new.groups);
    }
}
