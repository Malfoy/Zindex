use super::{sketch::Sketch, BuildConfig, MergePlanner};
use anyhow::{Context, Result};
use rayon::prelude::*;
use std::{
    cmp::Ordering,
    collections::{BTreeMap, BTreeSet, BinaryHeap},
    sync::{
        atomic::{AtomicU64, AtomicUsize, Ordering as AtomicOrdering},
        Arc, Mutex,
    },
};
mod packed;
use packed::{PackedSketches, LANES};

#[derive(Clone, Debug)]
pub struct Group {
    pub members: Vec<usize>,
    pub cardinality: f64,
    pub slots: usize,
    pub width: usize,
    pub sketch: Sketch,
    /// Final union-sketch bytes, shared with rebuilt and hash-partitioned filters.
    pub ordering_key: Option<Arc<Vec<u8>>>,
    /// Sum of exact input cardinalities: a safe union bound, unlike the sketch.
    pub feature_count_upper_bound: usize,
    /// Sum of source workspaces, covering retained vectors during a union.
    pub input_peak_bytes: u64,
}

impl Group {
    fn retire_sketch(&mut self) {
        if !self.sketch.is_released() {
            self.ordering_key = Some(Arc::new(self.sketch.take_fingerprints()));
        }
    }
}

pub struct Plan {
    pub groups: Vec<Group>,
    /// Candidate pairs visited, including partial SIMD scans rejected early.
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
        ordering_key: None,
        feature_count_upper_bound: count,
        input_peak_bytes: 0,
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
    if count > config.max_group_members() || !merge_fits_memory(ga, gb, config)? {
        return Ok(None);
    }
    let j = ga.sketch.jaccard(&gb.sketch);
    candidate_with_similarity(a, b, ga, gb, config, j)
}

fn candidate_with_similarity(
    a: usize,
    b: usize,
    ga: &Group,
    gb: &Group,
    config: &BuildConfig,
    j: f64,
) -> Result<Option<Candidate>> {
    let count = ga.members.len() + gb.members.len();
    if count > config.max_group_members() || !merge_fits_memory(ga, gb, config)? {
        return Ok(None);
    }
    let cardinality = ((ga.cardinality + gb.cardinality) / (1.0 + j)).clamp(
        ga.cardinality.max(gb.cardinality),
        ga.cardinality + gb.cardinality,
    );
    // Size the merged group with a safety margin; the chained estimate itself stays uninflated.
    let slots = config.slots(cardinality * (1.0 + config.merge.plan_margin))?;
    let old = cost(ga.slots, ga.width) + cost(gb.slots, gb.width);
    let new = cost(slots, config.width(count));
    if new >= old {
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

fn merge_fits_memory(ga: &Group, gb: &Group, config: &BuildConfig) -> Result<bool> {
    let Some(budget) = config.merge.planning_memory_budget_bytes else {
        return Ok(true);
    };
    let features = ga
        .feature_count_upper_bound
        .saturating_add(gb.feature_count_upper_bound);
    // Use a cardinality upper bound and include the plan's safety margin. A
    // noisy/high-similarity sketch must not approve an unbounded graph job.
    let Ok(slots) = config.slots(features as f64 * (1.0 + config.merge.plan_margin)) else {
        return Ok(false);
    };
    Ok(super::memory::group_estimate(
        features,
        slots,
        config.hashes,
        config.width(ga.members.len() + gb.members.len()),
        ga.input_peak_bytes.saturating_add(gb.input_peak_bytes),
    ) <= budget)
}

fn merge(mut ga: Group, gb: Group, c: &Candidate, config: &BuildConfig) -> Group {
    ga.sketch.merge_in_place(&gb.sketch);
    let mut members = ga.members;
    members.extend(gb.members);
    members.sort_unstable();
    let mut group = Group {
        width: config.width(members.len()),
        members,
        cardinality: c.cardinality,
        slots: c.slots,
        sketch: ga.sketch,
        ordering_key: None,
        feature_count_upper_bound: ga
            .feature_count_upper_bound
            .saturating_add(gb.feature_count_upper_bound),
        input_peak_bytes: ga.input_peak_bytes.saturating_add(gb.input_peak_bytes),
    };
    if group.members.len() >= config.max_group_members() {
        group.retire_sketch();
    }
    group
}

/// Optimistic cost: union capacity cannot be smaller than either input's
/// capacity. If even this lower bound cannot save bits, Jaccard cannot help.
/// This is valid for singleton arrivals AND arrivals containing multiple origins.
#[cfg(test)]
fn bucket_can_save(slots: usize, origins: usize, incoming: &Group, config: &BuildConfig) -> bool {
    bucket_max_saving(slots, origins, incoming, config).is_some()
}

fn bucket_max_saving(
    slots: usize,
    origins: usize,
    incoming: &Group,
    config: &BuildConfig,
) -> Option<u128> {
    let count = origins + incoming.members.len();
    if count > config.max_group_members() {
        return None;
    }
    let old = cost(slots, config.width(origins)) + cost(incoming.slots, incoming.width);
    let lower_new = cost(slots.max(incoming.slots), config.width(count));
    if lower_new < old {
        Some(old - lower_new)
    } else {
        None
    }
}

/// Equality is prunable only when no earlier ID could win the existing tie.
fn can_beat(upper: u128, first_id: usize, best: &Option<Candidate>) -> bool {
    best.as_ref()
        .is_none_or(|best| upper > best.saving || (upper == best.saving && first_id < best.a))
}

#[derive(Default)]
struct Bucket {
    blocks: Vec<usize>,
    // Conservative earliest logical ID in each physical block. Stale-low
    // values after removals only reduce pruning, never change tie decisions.
    block_first: Vec<usize>,
    len: usize,
    first: usize,
    // A stale lower bound after removal is conservative, never unsafe.
    min_cardinality: f64,
}

impl Bucket {
    fn slot(&self, position: usize) -> usize {
        self.blocks[position / LANES] * LANES + position % LANES
    }

    /// Minimum matching bytes that could improve the incumbent. All groups in
    /// a cost bucket have identical old cost and merged width. Its smallest
    /// cardinality gives the most optimistic union for every possible match
    /// count. Binary search uses the *same* rounded estimator and slot sizing
    /// as the exact candidate calculation, including strict saving thresholds.
    fn required_matches(
        &self,
        slots: usize,
        origins: usize,
        incoming: &Group,
        config: &BuildConfig,
        best: &Option<Candidate>,
    ) -> Result<usize> {
        let size = incoming.sketch.fingerprints().len();
        let old = cost(slots, config.width(origins)) + cost(incoming.slots, incoming.width);
        let width = config.width(origins + incoming.members.len());
        let can_win = |matches| -> Result<bool> {
            let j = Sketch::similarity_from_matches(matches, size);
            let card = ((self.min_cardinality + incoming.cardinality) / (1.0 + j)).clamp(
                self.min_cardinality.max(incoming.cardinality),
                self.min_cardinality + incoming.cardinality,
            );
            let new = cost(
                config.slots(card * (1.0 + config.merge.plan_margin))?,
                width,
            );
            Ok(new < old && can_beat(old - new, self.first, best))
        };
        let (mut low, mut high) = (0, size + 1);
        while low < high {
            let mid = (low + high) / 2;
            if can_win(mid)? {
                high = mid;
            } else {
                low = mid + 1;
            }
        }
        Ok(low)
    }
}

/// Incremental planning state for one batch. Only sketches/metadata are merged;
/// no union keys or filters are constructed here. Previously indexed batches
/// must not be inserted into this state when planning an update.
#[derive(Default)]
pub struct StreamingPlanner {
    groups: Vec<Group>,
    // Naive placement never expands or revisits an existing size class.
    naive_sizes: BTreeSet<usize>,
    pair_evaluations: u64,
    possible_pair_evaluations: u64,
    // Sub-bucket by origin count, so both width and slot cost are constant.
    // Logical IDs never move; physical columns swap-remove and recycle storage.
    buckets: BTreeMap<(u32, usize), Bucket>,
    packed: PackedSketches,
    positions: Vec<usize>,
    bucket_work: Vec<(u128, usize, (u32, usize))>,
    // Short scans cannot amortize waking a large global pool. Keep their
    // workers isolated so unused global workers sleep rather than steal/spin.
    scan_pools: [Option<rayon::ThreadPool>; 2],
}
/// Smallest slice of candidate evaluations handed to a pool thread: a few thousand sketch comparisons are
/// about a millisecond, so smaller scans run on the calling thread instead of paying for thread hand-offs.
const CANDIDATES_PER_TASK: usize = 1024;

impl StreamingPlanner {
    fn add_to_bucket(&mut self, id: usize) {
        let group = &self.groups[id];
        if group.sketch.is_released() {
            return;
        }
        let key = (group.slots.trailing_zeros(), group.members.len());
        let bucket = self.buckets.entry(key).or_default();
        if bucket.len % LANES == 0 {
            bucket
                .blocks
                .push(self.packed.allocate(group.sketch.fingerprints().len()));
            bucket.block_first.push(id);
        }
        bucket.block_first[bucket.len / LANES] = bucket.block_first[bucket.len / LANES].min(id);
        self.positions.resize(self.groups.len(), 0);
        self.positions[id] = bucket.len;
        self.packed
            .write(bucket.slot(bucket.len), id, group.sketch.fingerprints());
        if bucket.len == 0 {
            bucket.first = id;
            bucket.min_cardinality = group.cardinality;
        } else {
            bucket.first = bucket.first.min(id);
            bucket.min_cardinality = bucket.min_cardinality.min(group.cardinality);
        }
        bucket.len += 1;
    }

    fn remove_from_bucket(&mut self, id: usize) {
        let group = &self.groups[id];
        let key = (group.slots.trailing_zeros(), group.members.len());
        let bucket = self
            .buckets
            .get_mut(&key)
            .expect("missing streaming bucket");
        let position = self.positions[id];
        let moved = self
            .packed
            .move_slot(bucket.slot(bucket.len - 1), bucket.slot(position));
        self.positions[moved] = position;
        bucket.block_first[position / LANES] = bucket.block_first[position / LANES].min(moved);
        bucket.len -= 1;
        if bucket.len % LANES == 0 {
            self.packed.free(bucket.blocks.pop().unwrap());
            bucket.block_first.pop();
        }
        if bucket.len == 0 {
            self.buckets.remove(&key);
        } else if bucket.first == id {
            bucket.first = (0..bucket.len)
                .map(|i| self.packed.ids[bucket.slot(i)])
                .min()
                .unwrap();
        }
    }

    fn best_bounded(
        &mut self,
        incoming: &Group,
        config: &BuildConfig,
    ) -> Result<Option<Candidate>> {
        self.bucket_work.clear();
        for (&key @ (exponent, origins), bucket) in &self.buckets {
            if let Some(upper) = bucket_max_saving(1usize << exponent, origins, incoming, config) {
                self.bucket_work.push((upper, bucket.first, key));
            }
        }
        self.bucket_work
            .sort_unstable_by(|a, b| b.0.cmp(&a.0).then_with(|| a.1.cmp(&b.1)));
        let mut best = None;
        for index in 0..self.bucket_work.len() {
            let (upper, first, key) = self.bucket_work[index];
            if !can_beat(upper, first, &best) {
                // Sorted by decreasing bound, then increasing earliest ID:
                // every remaining bucket is also unable to win.
                break;
            }
            let bytes = self.buckets[&key].len.saturating_mul(self.packed.size);
            let pool_index = if bytes < 16 * 1024 * 1024 {
                Some(0)
            } else if bytes < 32 * 1024 * 1024 {
                Some(1)
            } else {
                None
            };
            let pool_index = pool_index.filter(|&i| rayon::current_num_threads() > (8 << i));
            // Delay pool creation until after the cheap earliest-ID probe.
            let bucket = &self.buckets[&key];
            // Preserve the O(1) hot-cluster case and the earliest-ID tie.
            self.pair_evaluations += 1;
            best = best.max(candidate(
                first,
                self.groups.len(),
                &self.groups[first],
                incoming,
                config,
            )?);
            if !can_beat(upper, first, &best) {
                continue;
            }
            let required =
                bucket.required_matches(1usize << key.0, key.1, incoming, config, &best)?;
            if required > self.packed.size {
                continue;
            }
            if bucket.blocks.len() >= 1024 {
                if let Some(i) = pool_index {
                    if self.scan_pools[i].is_none() {
                        self.scan_pools[i] = Some(
                            rayon::ThreadPoolBuilder::new()
                                .num_threads(8 << i)
                                .thread_name(move |n| format!("layout-{}-{n}", 8 << i))
                                .build()
                                .context("creating streaming scan worker pool")?,
                        );
                    }
                }
            }
            let incumbent = Mutex::new(best.clone());
            let revision = AtomicUsize::new(0);
            let visited = AtomicU64::new(0);
            let scan = |start: usize, blocks: &[usize]| -> Result<Option<Candidate>> {
                let mut found = best.clone();
                let mut required = required;
                let mut observed_revision = 0;
                let mut comparisons = 0u64;
                let mut counts = [0u32; 4 * LANES];
                for (offset, chunk) in blocks.chunks(4).enumerate() {
                    let current_revision = revision.load(AtomicOrdering::Acquire);
                    if current_revision != observed_revision {
                        found = found.max(incumbent.lock().unwrap().clone());
                        observed_revision = current_revision;
                        required = bucket.required_matches(
                            1usize << key.0,
                            key.1,
                            incoming,
                            config,
                            &found,
                        )?;
                    }
                    if required > self.packed.size {
                        break;
                    }
                    let chunk_start = start + offset * 4;
                    let first_id = bucket.block_first[chunk_start..chunk_start + chunk.len()]
                        .iter()
                        .copied()
                        .min()
                        .unwrap();
                    if !can_beat(upper, first_id, &found) {
                        continue;
                    }
                    let valid =
                        (bucket.len - (start + offset * 4) * LANES).min(chunk.len() * LANES);
                    comparisons += valid as u64
                        - u64::from(
                            (chunk_start * LANES..chunk_start * LANES + valid)
                                .contains(&self.positions[first]),
                        );
                    let previous = found.clone();
                    let query = incoming.sketch.fingerprints();
                    let mut live = if valid == 4 * LANES {
                        packed::matches_quad(
                            [
                                self.packed.block(chunk[0]),
                                self.packed.block(chunk[1]),
                                self.packed.block(chunk[2]),
                                self.packed.block(chunk[3]),
                            ],
                            query,
                            required,
                            &mut counts,
                        )
                    } else {
                        let mut live = 0u128;
                        for (i, &block) in chunk.iter().enumerate() {
                            let out = (&mut counts[i * LANES..(i + 1) * LANES])
                                .try_into()
                                .unwrap();
                            let mask = packed::matches(
                                self.packed.block(block),
                                query,
                                required,
                                (valid - i * LANES).min(LANES),
                                out,
                            );
                            live |= (mask as u128) << (i * LANES);
                        }
                        live
                    };
                    while live != 0 {
                        let lane = live.trailing_zeros() as usize;
                        live &= live - 1;
                        let id = self.packed.ids[chunk[lane / LANES] * LANES + lane % LANES];
                        if id != first && can_beat(upper, id, &found) {
                            let j = Sketch::similarity_from_matches(
                                counts[lane] as usize,
                                self.packed.size,
                            );
                            found = found.max(candidate_with_similarity(
                                id,
                                self.groups.len(),
                                &self.groups[id],
                                incoming,
                                config,
                                j,
                            )?);
                        }
                    }
                    if found != previous {
                        required = bucket.required_matches(
                            1usize << key.0,
                            key.1,
                            incoming,
                            config,
                            &found,
                        )?;
                        let mut shared = incumbent.lock().unwrap();
                        if found > *shared {
                            *shared = found.clone();
                            revision.fetch_add(1, AtomicOrdering::Release);
                        }
                    }
                }
                visited.fetch_add(comparisons, AtomicOrdering::Relaxed);
                Ok(found)
            };
            // One coarse dispatch, not hundreds of speculative batches. Small
            // scans stay local to avoid making scheduling the bottleneck again.
            let found = if bucket.blocks.len() >= 1024 {
                // One pinned job per worker, with a small shared block queue.
                // Fast and slow cores claim different numbers of chunks, but
                // avoid recursively splitting/stealing millions of tiny jobs.
                let next = AtomicUsize::new(0);
                let worker = |_: rayon::BroadcastContext<'_>| -> Result<Option<Candidate>> {
                    let mut found = None;
                    loop {
                        let start = next.fetch_add(256, AtomicOrdering::Relaxed);
                        if start >= bucket.blocks.len() {
                            break;
                        }
                        let end = (start + 256).min(bucket.blocks.len());
                        found = found.max(scan(start, &bucket.blocks[start..end])?);
                    }
                    Ok(found)
                };
                let results = if let Some(i) = pool_index {
                    self.scan_pools[i].as_ref().unwrap().broadcast(worker)
                } else {
                    rayon::broadcast(worker)
                };
                results
                    .into_iter()
                    .try_fold(None, |a, b| Ok::<_, anyhow::Error>(a.max(b?)))?
            } else {
                scan(0, &bucket.blocks)?
            };
            best = best.max(found);
            self.pair_evaluations += visited.load(AtomicOrdering::Relaxed);
        }
        Ok(best)
    }

    pub fn insert(&mut self, mut incoming: Group, config: &BuildConfig) -> Result<()> {
        if config.merge.planner == MergePlanner::Naive {
            let required = incoming.slots;
            incoming.slots = match self.naive_sizes.range(required..).next().copied() {
                Some(slots) => slots,
                None => {
                    self.naive_sizes.insert(required);
                    required
                }
            };
            incoming.sketch.release();
            self.groups.push(incoming);
            return Ok(());
        }
        #[cfg(test)]
        let unpruned = config.merge.planner == MergePlanner::StreamingUnpruned;
        #[cfg(not(test))]
        let unpruned = false;
        if config.merge.disabled {
            if incoming.members.len() >= config.max_group_members() {
                incoming.retire_sketch();
            }
            self.groups.push(incoming);
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
            self.best_bounded(&incoming, config)?
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
                group.sketch.merge_in_place(&incoming.sketch);
                group.feature_count_upper_bound = group
                    .feature_count_upper_bound
                    .saturating_add(incoming.feature_count_upper_bound);
                group.input_peak_bytes = group
                    .input_peak_bytes
                    .saturating_add(incoming.input_peak_bytes);
                if group.members.len() >= config.max_group_members() {
                    group.retire_sketch();
                }
                self.add_to_bucket(c.a);
            }
        } else {
            if incoming.members.len() >= config.max_group_members() {
                incoming.retire_sketch();
            }
            self.groups.push(incoming);
            if !unpruned {
                self.add_to_bucket(self.groups.len() - 1);
            }
        }
        Ok(())
    }

    pub fn finish(mut self) -> Plan {
        for group in &mut self.groups {
            group.retire_sketch();
        }
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
pub fn plan(mut initial: Vec<Group>, config: &BuildConfig) -> Result<Plan> {
    if config.merge.planner.is_streaming() {
        let mut state = StreamingPlanner::default();
        for group in initial {
            state.insert(group, config)?;
        }
        return Ok(state.finish());
    }
    if config.merge.disabled {
        for group in &mut initial {
            group.retire_sketch();
        }
        return Ok(Plan {
            groups: initial,
            pair_evaluations: 0,
            possible_pair_evaluations: 0,
        });
    }
    for group in &mut initial {
        if group.members.len() >= config.max_group_members() {
            group.retire_sketch();
        }
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
    for group in &mut result {
        group.retire_sketch();
    }
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
            segment_length: None,
            // Exact sizes: these tests check the cost rule itself, not the sketch-noise margin.
            merge: MergeArgs {
                plan_margin: 0.0,
                ..MergeArgs::default()
            },
        }
    }

    #[test]
    fn naive_places_in_smallest_existing_capacity_without_revisiting_or_merging() {
        let mut c = config(Backend::Zor);
        c.merge.planner = MergePlanner::Naive;
        let mut state = StreamingPlanner::default();
        for (id, slots) in [16, 8, 64, 32, 128, 16].into_iter().enumerate() {
            let mut incoming = singleton(3, Sketch::default(), id, &c).unwrap();
            incoming.slots = slots;
            state.insert(incoming, &c).unwrap();
        }
        let plan = state.finish();
        assert_eq!(plan.pair_evaluations, 0);
        assert_eq!(plan.possible_pair_evaluations, 0);
        assert_eq!(
            plan.groups.iter().map(|g| g.slots).collect::<Vec<_>>(),
            [16, 16, 64, 64, 128, 16]
        );
        for (id, group) in plan.groups.iter().enumerate() {
            assert_eq!(group.members, [id]);
            assert!(group.sketch.is_released());
            assert_eq!(group.width, c.fingerprint_bits);
        }
        // Changing arrival order may change padding, never the set of datasets.
        let mut reversed = StreamingPlanner::default();
        for (id, slots) in [128, 64, 16].into_iter().enumerate() {
            let mut incoming = singleton(3, Sketch::default(), id, &c).unwrap();
            incoming.slots = slots;
            reversed.insert(incoming, &c).unwrap();
        }
        assert!(reversed.finish().groups.iter().all(|g| g.slots == 128));
    }

    #[test]
    fn every_planner_keeps_safe_pairs_instead_of_forming_an_oversized_group() {
        for mode in [
            MergePlanner::Streaming,
            MergePlanner::StreamingUnpruned,
            MergePlanner::AllPairs,
        ] {
            let mut c = config(Backend::Zor);
            c.hashes = 8;
            c.fingerprint_bits = 8;
            c.merge.planner = mode;
            c.merge.planning_memory_budget_bytes = Some(40 << 20);
            let initial: Vec<_> = (0..4)
                .map(|id| {
                    let mut group =
                        singleton(50_000, Sketch::new(&[1, 2, 3], 64, 69), id, &c).unwrap();
                    group.input_peak_bytes = 8 << 20;
                    group
                })
                .collect();
            let groups = if mode.is_streaming() {
                let mut planner = StreamingPlanner::default();
                for group in initial {
                    planner.insert(group, &c).unwrap();
                }
                planner.finish().groups
            } else {
                plan(initial, &c).unwrap().groups
            };
            assert_eq!(groups.len(), 2, "{mode:?}");
            assert!(groups.iter().all(|g| g.members.len() == 2
                && g.feature_count_upper_bound == 100_000
                && g.input_peak_bytes == 16 << 20));
            let mut members: Vec<_> = groups
                .iter()
                .flat_map(|g| g.members.iter().copied())
                .collect();
            members.sort_unstable();
            assert_eq!(members, vec![0, 1, 2, 3]);
            assert!(groups.iter().all(|g| super::super::memory::group_estimate(
                g.feature_count_upper_bound,
                c.slots(g.feature_count_upper_bound as f64).unwrap(),
                c.hashes,
                g.width,
                g.input_peak_bytes
            ) <= 40 << 20));
        }
    }

    #[test]
    fn memory_check_uses_exact_count_bound_even_when_sketches_are_identical() {
        let mut c = config(Backend::Zor);
        c.hashes = 8;
        c.fingerprint_bits = 8;
        c.merge.planning_memory_budget_bytes = Some(200 << 30);
        let a = singleton(300_000_000, Sketch::new(&[1, 2, 3], 64, 69), 0, &c).unwrap();
        let b = singleton(300_000_000, Sketch::new(&[1, 2, 3], 64, 69), 1, &c).unwrap();
        assert!(candidate_with_similarity(0, 1, &a, &b, &c, 1.0)
            .unwrap()
            .is_none());
        c.merge.planning_memory_budget_bytes = Some(1 << 40);
        let candidate = candidate_with_similarity(0, 1, &a, &b, &c, 1.0)
            .unwrap()
            .unwrap();
        let merged = merge(a, b, &candidate, &c);
        assert_eq!(merged.feature_count_upper_bound, 600_000_000);
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

    #[test]
    fn one_and_two_bit_singletons_cannot_save_bits_by_pairwise_merging() {
        for bits in [1, 2] {
            let mut c = config(Backend::Zor);
            c.fingerprint_bits = bits;
            for a in [1, 63, 64, 65, 127, 128, 129, 4096] {
                for b in [1, 63, 64, 65, 127, 128, 129, 4096] {
                    let ga = singleton(a, Sketch::new(&[1, 2, 3], 64, 69), 0, &c).unwrap();
                    let gb = singleton(b, Sketch::new(&[1, 2, 3], 64, 69), 1, &c).unwrap();
                    for similarity in [0.0, 0.5, 1.0] {
                        assert!(candidate_with_similarity(0, 1, &ga, &gb, &c, similarity)
                            .unwrap()
                            .is_none());
                    }
                }
            }
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
            assert_eq!(a.sketch, b.sketch);
            assert_eq!(a.ordering_key, b.ordering_key);
            assert!(a.sketch.is_released());
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
        let cap = c.max_group_members();
        let identical = vec![(0..1000).collect(); 2 * cap + 1];
        let result = plan(initial(&identical, &c), &c).unwrap();
        assert_eq!(
            result
                .groups
                .iter()
                .map(|g| g.members.clone())
                .collect::<Vec<_>>(),
            [
                (0..cap).collect::<Vec<_>>(),
                (cap..2 * cap).collect::<Vec<_>>(),
                vec![2 * cap]
            ]
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
            assert_eq!(a.ordering_key, b.ordering_key);
        }
    }

    #[test]
    fn packed_columns_stay_exact_through_merges_retirement_and_reuse() {
        for size in [64, 256, 4096] {
            let mut c = config(Backend::Zor);
            c.merge.sketch_size = size;
            c.merge.disabled = true;
            let sets: Vec<Vec<u64>> = (0..100)
                .map(|i| (i * 8192..i * 8192 + 1000).collect())
                .collect();
            let mut packed = StreamingPlanner::default();
            let mut old = StreamingPlanner::default();
            for group in initial(&sets, &c) {
                packed.insert(group.clone(), &c).unwrap();
                let mut legacy = c.clone();
                legacy.merge.planner = MergePlanner::StreamingUnpruned;
                old.insert(group, &legacy).unwrap();
            }
            c.merge.disabled = false;
            for i in 0..210 {
                let target = (i * 37) % sets.len();
                let group = singleton(
                    1000,
                    Sketch::new(&sets[target], size, c.seed),
                    sets.len() + i,
                    &c,
                )
                .unwrap();
                let mut legacy = c.clone();
                legacy.merge.planner = MergePlanner::StreamingUnpruned;
                old.insert(group.clone(), &legacy).unwrap();
                packed.insert(group, &c).unwrap();
                assert_same_groups(&packed.groups, &old.groups);
                let mut seen = vec![false; packed.groups.len()];
                for bucket in packed.buckets.values() {
                    for position in 0..bucket.len {
                        let slot = bucket.slot(position);
                        let id = packed.packed.ids[slot];
                        assert!(!seen[id]);
                        seen[id] = true;
                        assert_eq!(packed.positions[id], position);
                        for (row, &byte) in
                            packed.groups[id].sketch.fingerprints().iter().enumerate()
                        {
                            assert_eq!(
                                packed.packed.block(slot / LANES)[row * LANES + slot % LANES],
                                byte
                            );
                        }
                    }
                }
                for (group, active) in packed.groups.iter().zip(seen) {
                    assert_eq!(active, !group.sketch.is_released());
                }
            }
        }
    }

    #[test]
    fn pruned_streaming_matches_original_across_widths_rounding_limits_and_ties() {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap();
        pool.install(|| {
            for fingerprint_bits in [1, 2, 4, 8, 16, 32] {
                for minimum_slots in [None, Some(256)] {
                    let mut c = config(Backend::Zor);
                    c.fingerprint_bits = fingerprint_bits;
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
        });
    }

    #[test]
    fn default_sketch_size_and_optional_margin_preserve_exhaustive_streaming_decisions() {
        let sets: Vec<Vec<u64>> = (0..48)
            .map(|i| {
                let start = (i % 4) * 40_000;
                let count = [8000, 8192, 8500, 16000, 16384, 17000][i as usize % 6];
                (start..start + count).collect()
            })
            .collect();
        for plan_margin in [0.0, 0.03, 0.1, 0.49] {
            let mut c = config(Backend::Zor);
            assert_eq!(c.merge.sketch_size, 128);
            c.merge.plan_margin = plan_margin;
            let groups = initial(&sets, &c);
            c.merge.planner = MergePlanner::StreamingUnpruned;
            let old = plan(groups.clone(), &c).unwrap();
            c.merge.planner = MergePlanner::Streaming;
            let new = plan(groups, &c).unwrap();
            assert_same_groups(&old.groups, &new.groups);
            assert_eq!(new.possible_pair_evaluations, old.pair_evaluations);
            assert!(new.pair_evaluations <= old.pair_evaluations);
        }
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
        assert!(!bucket_can_save(128, c.max_group_members(), &incoming, &c));
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
        let sets = vec![(0..1000).collect(); c.max_group_members() + 12];
        let new = plan(initial(&sets, &c), &c).unwrap();
        c.merge.planner = MergePlanner::StreamingUnpruned;
        let old = plan(initial(&sets, &c), &c).unwrap();
        assert_same_groups(&old.groups, &new.groups);
        assert_eq!(new.groups[0].members.len(), c.max_group_members());
        assert!(new.pair_evaluations < old.pair_evaluations);

        // A big group skipped for a small arrival must not be permanently retired.
        c.merge.planner = MergePlanner::Streaming;
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

    #[test]
    fn saturating_match_skips_the_rest_of_a_hot_bucket_without_changing_ties() {
        let mut c = config(Backend::Zor);
        c.merge.sketch_size = 64;
        c.merge.disabled = true;
        let keys: Vec<_> = (0..1000).collect();
        let sets = vec![keys.clone(); 600];
        let mut state = StreamingPlanner::default();
        for g in initial(&sets, &c) {
            state.insert(g, &c).unwrap();
        }
        let mut exhaustive = StreamingPlanner::default();
        c.merge.planner = MergePlanner::StreamingUnpruned;
        for g in initial(&sets, &c) {
            exhaustive.insert(g, &c).unwrap();
        }
        c.merge.disabled = false;
        for id in 600..603 {
            let g = singleton(keys.len(), Sketch::new(&keys, 64, c.seed), id, &c).unwrap();
            exhaustive.insert(g.clone(), &c).unwrap();
            c.merge.planner = MergePlanner::Streaming;
            state.insert(g, &c).unwrap();
            c.merge.planner = MergePlanner::StreamingUnpruned;
            assert_same_groups(&state.groups, &exhaustive.groups);
        }
        assert_eq!(state.pair_evaluations, 3);
        assert_eq!(exhaustive.pair_evaluations, 1800);
        assert_eq!(state.groups[0].members, [0, 600, 601, 602]);
    }

    #[test]
    fn late_optimal_candidate_prunes_later_blocks_without_losing_earlier_ties() {
        let mut c = config(Backend::Zor);
        c.merge.disabled = true;
        let sets: Vec<Vec<u64>> = (0..1024)
            .map(|id| (id * 2000..id * 2000 + 1000).collect())
            .collect();
        let mut fast = StreamingPlanner::default();
        let mut old = StreamingPlanner::default();
        for group in initial(&sets, &c) {
            fast.insert(group.clone(), &c).unwrap();
            c.merge.planner = MergePlanner::StreamingUnpruned;
            old.insert(group, &c).unwrap();
            c.merge.planner = MergePlanner::Streaming;
        }
        c.merge.disabled = false;
        let incoming = singleton(
            1000,
            Sketch::new(&sets[128], c.merge.sketch_size, c.seed),
            1024,
            &c,
        )
        .unwrap();
        fast.insert(incoming.clone(), &c).unwrap();
        c.merge.planner = MergePlanner::StreamingUnpruned;
        old.insert(incoming, &c).unwrap();
        assert_same_groups(&fast.groups, &old.groups);
        assert_eq!(fast.groups[128].members, [128, 1024]);
        assert!(fast.pair_evaluations < old.pair_evaluations / 2);
    }

    #[test]
    fn planning_allocations_are_freed_while_ordering_fingerprints_are_retained() {
        let mut c = config(Backend::Zor);
        c.merge.sketch_size = 64;
        let sets = vec![(0..1000).collect(); c.max_group_members() + 2];
        let mut state = StreamingPlanner::default();
        for g in initial(&sets, &c) {
            state.insert(g, &c).unwrap();
        }
        assert!(state.groups[0].sketch.is_released());
        assert_eq!(state.groups[0].sketch.allocated_bytes(), 0);
        let key = state.groups[0].ordering_key.clone().unwrap();
        assert_eq!(key.len(), c.merge.sketch_size);
        assert!(state
            .buckets
            .values()
            .all(|bucket| (0..bucket.len).all(|i| state.packed.ids[bucket.slot(i)] != 0)));
        assert!(!state.groups[1].sketch.is_released());
        let plan = state.finish();
        assert!(plan.groups.iter().all(|g| g.sketch.allocated_bytes() == 0));
        assert!(Arc::ptr_eq(
            &key,
            plan.groups[0].ordering_key.as_ref().unwrap()
        ));
        assert!(plan
            .groups
            .iter()
            .all(|g| g.ordering_key.as_ref().unwrap().len() == c.merge.sketch_size));
        c.merge.planner = MergePlanner::AllPairs;
        assert!(super::plan(initial(&sets, &c), &c)
            .unwrap()
            .groups
            .iter()
            .all(|g| g.sketch.allocated_bytes() == 0));
    }

    #[test]
    fn ordering_fingerprints_equal_the_final_group_union_sketch() {
        let sets: Vec<Vec<u64>> = vec![
            (0..256).collect(),
            (128..384).collect(),
            (4000..4256).collect(),
            (0..256).collect(),
            Vec::new(),
        ];
        for planner in [MergePlanner::Streaming, MergePlanner::AllPairs] {
            for disabled in [false, true] {
                let mut c = config(Backend::Zor);
                c.merge.planner = planner;
                c.merge.disabled = disabled;
                let result = plan(initial(&sets, &c), &c).unwrap();
                for group in result.groups {
                    let mut union: Vec<_> = group
                        .members
                        .iter()
                        .flat_map(|&id| sets[id].iter().copied())
                        .collect();
                    union.sort_unstable();
                    union.dedup();
                    let expected = Sketch::new(&union, c.merge.sketch_size, c.seed);
                    assert_eq!(
                        group.ordering_key.as_ref().unwrap().as_slice(),
                        expected.fingerprints()
                    );
                    assert_eq!(group.sketch.allocated_bytes(), 0);
                }
            }
        }
    }

    #[test]
    fn ordered_parallel_batches_match_exhaustive_before_and_after_every_arrival() {
        let mut c = config(Backend::Zor);
        c.merge.sketch_size = 64;
        c.merge.disabled = true;
        let sets: Vec<Vec<u64>> = (0..520)
            .map(|i| (i * 4096..i * 4096 + 1000).collect())
            .collect();
        let mut bounded = StreamingPlanner::default();
        let mut exhaustive = StreamingPlanner::default();
        for g in initial(&sets, &c) {
            bounded.insert(g.clone(), &c).unwrap();
            c.merge.planner = MergePlanner::StreamingUnpruned;
            exhaustive.insert(g, &c).unwrap();
            c.merge.planner = MergePlanner::Streaming;
        }
        c.merge.disabled = false;
        for (offset, target) in [400, 1, 519, 20, 400, 300, 1].into_iter().enumerate() {
            let g = singleton(
                sets[target].len(),
                Sketch::new(&sets[target], 64, c.seed),
                520 + offset,
                &c,
            )
            .unwrap();
            bounded.insert(g.clone(), &c).unwrap();
            c.merge.planner = MergePlanner::StreamingUnpruned;
            exhaustive.insert(g, &c).unwrap();
            c.merge.planner = MergePlanner::Streaming;
            assert_same_groups(&bounded.groups, &exhaustive.groups);
        }
    }

    #[test]
    fn parallel_cancellation_keeps_the_earliest_equal_saving_id() {
        let mut c = config(Backend::Zor);
        c.merge.sketch_size = 64;
        c.merge.disabled = true;
        let mut sets: Vec<Vec<u64>> = (0..600)
            .map(|i| (i * 4096..i * 4096 + 1000).collect())
            .collect();
        sets[450] = sets[300].clone();
        let groups = initial(&sets, &c);
        for threads in [1, 4, 32] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            pool.install(|| {
                let mut state = StreamingPlanner::default();
                for group in groups.clone() {
                    state.insert(group, &c).unwrap();
                }
                let mut enabled = c.clone();
                enabled.merge.disabled = false;
                let incoming =
                    singleton(1000, Sketch::new(&sets[300], 64, c.seed), 600, &c).unwrap();
                state.insert(incoming, &enabled).unwrap();
                assert_eq!(state.groups[300].members, [300, 600]);
                assert_eq!(state.groups[450].members, [450]);
            });
        }
    }

    #[test]
    fn packed_parallel_chunks_choose_the_earliest_profitable_winner() {
        let mut c = config(Backend::Zor);
        c.merge.sketch_size = 64;
        c.merge.disabled = true;
        let ordinary = singleton(
            1000,
            Sketch::new(&(0..1000).collect::<Vec<_>>(), 64, c.seed),
            0,
            &c,
        )
        .unwrap();
        let target = Sketch::new(&(100000..101000).collect::<Vec<_>>(), 64, c.seed);
        // Exceed the coarse dispatch threshold; put equal-saving winners in
        // different tasks, well away from the earliest-ID serial probe.
        for threads in [1, 4, 32] {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| {
                    let mut packed = StreamingPlanner::default();
                    let mut legacy = StreamingPlanner::default();
                    let mut old_config = c.clone();
                    old_config.merge.planner = MergePlanner::StreamingUnpruned;
                    for id in 0..33000 {
                        let mut group = ordinary.clone();
                        group.members[0] = id;
                        if id == 300 || id == 32900 {
                            group.sketch = target.clone();
                        }
                        packed.insert(group.clone(), &c).unwrap();
                        legacy.insert(group, &old_config).unwrap();
                    }
                    let mut enabled = c.clone();
                    enabled.merge.disabled = false;
                    old_config.merge.disabled = false;
                    let arrival = singleton(1000, target.clone(), 33000, &enabled).unwrap();
                    packed.insert(arrival.clone(), &enabled).unwrap();
                    legacy.insert(arrival, &old_config).unwrap();
                    assert_same_groups(&packed.groups, &legacy.groups);
                    assert_eq!(packed.groups[300].members, [300, 33000]);
                    assert_eq!(packed.groups[32900].members, [32900]);
                });
        }
    }
}
