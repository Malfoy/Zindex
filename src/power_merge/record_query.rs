//! Exact threshold pruning and bounded, chunk-local feature-result reuse.
//! Nothing in this module is serialized. A cached result always covers ALL
//! documents, independently of the candidate set of the read that requested it.
use super::*;
use crate::power_merge::records::Record;

const FEATURE_BUDGET: usize = 32 * 1024;
const CACHED_ORIGIN_BUDGET: usize = 256 * 1024;

// Match the reporter's division/comparison, including floating-point boundaries,
// rather than assuming ceil(total * ratio) has exactly the same rounding.
fn minimum_matches(total: u64, ratio: f64) -> u64 {
    if total == 0 || ratio == 0.0 {
        return 0;
    }
    let (mut lo, mut hi) = (0, total);
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if mid as f64 / total as f64 >= ratio {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    lo
}

#[derive(Default)]
struct Candidates {
    ranges: Vec<Vec<Range<usize>>>,
    spare: Vec<Vec<Range<usize>>>,
    required: u64,
    remaining: u64,
    until_check: u64,
    enabled: bool,
    alive: bool,
}

impl Candidates {
    fn reset(&mut self, index: &Index, total: u64, ratio: f64) {
        self.ranges.resize_with(index.banks.len(), Vec::new);
        self.spare.resize_with(index.banks.len(), Vec::new);
        for (ranges, bank) in self.ranges.iter_mut().zip(&index.banks) {
            ranges.clear();
            ranges.push(0..bank.meta.groups.len());
        }
        self.required = minimum_matches(total, ratio);
        self.remaining = total;
        self.until_check = (total - self.required).saturating_add(1);
        self.enabled = ratio > 0.0 && total != 0;
        self.alive = total != 0;
    }

    fn advance(&mut self, index: &Index, scores: &RecordScores, count: u64) {
        self.remaining -= count;
        if !self.enabled || self.remaining == 0 {
            return;
        }
        if count < self.until_check {
            self.until_check -= count;
            return;
        }
        self.alive = false;
        self.until_check = u64::MAX;
        for (bank_id, bank) in index.banks.iter().enumerate() {
            let mut next = std::mem::take(&mut self.spare[bank_id]);
            next.clear();
            // Keep packed range starts byte aligned. Retain a whole packing
            // unit if any member can pass (one-bit uses its existing word unit).
            let bits = index.manifest.fingerprint_bits;
            let unit = if bank.meta.fingerprint_first && bits < 8 {
                if bits == 1 {
                    64
                } else {
                    8 / bits
                }
            } else {
                1
            };
            for range in &self.ranges[bank_id] {
                for start in (range.start..range.end).step_by(unit) {
                    let end = (start + unit).min(range.end);
                    let best = bank.meta.groups[start..end]
                        .iter()
                        .flat_map(|group| &group.members)
                        .map(|&id| scores.values[id])
                        .max()
                        .unwrap_or(0);
                    let upper = best + self.remaining;
                    if upper < self.required {
                        continue;
                    }
                    self.alive = true;
                    // No block can become impossible sooner than this many
                    // additional misses. Avoid rescanning candidates per key.
                    self.until_check = self.until_check.min(upper - self.required + 1);
                    if let Some(last) = next.last_mut().filter(|last| last.end == start) {
                        last.end = end;
                    } else {
                        next.push(start..end);
                    }
                }
            }
            self.spare[bank_id] = std::mem::replace(&mut self.ranges[bank_id], next);
        }
    }
}

#[derive(Clone, Copy)]
enum Cached {
    Pending,
    Ready { start: usize, end: usize },
    TooLarge,
}

// One membership lookup can contribute independently to both orientations.
// get() observes a strand with nonzero weight, preserving the abandoned-key
// union's before/after counter test without a temporary origin list.
struct StrandPairScores<'a> {
    forward: &'a mut RecordScores,
    reverse: &'a mut RecordScores,
    forward_weight: u64,
    reverse_weight: u64,
}

impl ScoreSink for StrandPairScores<'_> {
    #[inline]
    fn get(&self, origin: usize) -> u64 {
        if self.forward_weight != 0 {
            self.forward.values[origin]
        } else {
            self.reverse.values[origin]
        }
    }
    #[inline]
    fn increment(&mut self, origin: usize) {
        self.add(origin, 1);
    }
    #[inline]
    fn add(&mut self, origin: usize, count: u64) {
        if self.forward_weight != 0 {
            self.forward.add(origin, count * self.forward_weight);
        }
        if self.reverse_weight != 0 {
            self.reverse.add(origin, count * self.reverse_weight);
        }
    }
}

pub(super) struct Worker {
    pub(super) scores: RecordScores,
    strand_scores: Option<RecordScores>,
    full_hits: RecordScores,
    scratch: Vec<u8>,
    keys: Vec<u64>,
    sorted: Vec<u64>,
    ends: Vec<usize>,
    cache: HashMap<u64, Cached>,
    origins: Vec<u32>,
    candidates: Candidates,
    positions: Vec<usize>,
    #[cfg(test)]
    full_lookups: usize,
    #[cfg(test)]
    cache_hits: usize,
    #[cfg(test)]
    selected_lookups: usize,
}

impl Worker {
    pub(super) fn new(index: &Index) -> Self {
        Self {
            scores: RecordScores::new(index.manifest.colors.len()),
            strand_scores: None,
            full_hits: RecordScores::new(index.manifest.colors.len()),
            scratch: Vec::new(),
            keys: Vec::new(),
            sorted: Vec::new(),
            ends: Vec::new(),
            cache: HashMap::new(),
            origins: Vec::new(),
            candidates: Candidates::default(),
            positions: Vec::new(),
            #[cfg(test)]
            full_lookups: 0,
            #[cfg(test)]
            cache_hits: 0,
            #[cfg(test)]
            selected_lookups: 0,
        }
    }

    fn prepare_cache(&mut self) {
        self.cache.clear();
        self.origins.clear();
        self.sorted.clear();
        self.sorted.extend_from_slice(&self.keys);
        self.sorted.sort_unstable();
        let mut i = 0;
        while i < self.sorted.len() {
            let key = self.sorted[i];
            let start = i;
            while i < self.sorted.len() && self.sorted[i] == key {
                i += 1;
            }
            if i - start > 1 {
                self.cache.insert(key, Cached::Pending);
            }
        }
    }

    fn query_key(&mut self, index: &Index, key: u64, weight: u64) {
        self.query_key_prepared(index, key, weight, None);
    }

    fn query_key_prepared(
        &mut self,
        index: &Index,
        key: u64,
        weight: u64,
        prepared: Option<&[usize]>,
    ) {
        match self.cache.get(&key).copied() {
            Some(Cached::Ready { start, end }) => {
                for &origin in &self.origins[start..end] {
                    self.scores.add(origin as usize, weight);
                }
                #[cfg(test)]
                {
                    self.cache_hits += 1;
                }
            }
            Some(Cached::Pending) => {
                // Cache only complete results. In particular, never memoize a
                // lookup that skipped groups retired for a different record.
                self.full_hits.clear();
                index.lookup_prepared(key, &mut self.full_hits, &mut self.scratch, prepared);
                #[cfg(test)]
                {
                    self.full_lookups += 1;
                }
                let start = self.origins.len();
                let fits = self.full_hits.touched.len() <= CACHED_ORIGIN_BUDGET - start;
                for &origin in &self.full_hits.touched {
                    debug_assert_eq!(self.full_hits.values[origin], 1);
                    self.scores.add(origin, weight);
                    if fits {
                        self.origins.push(origin as u32);
                    }
                }
                self.cache.insert(
                    key,
                    if fits {
                        Cached::Ready {
                            start,
                            end: self.origins.len(),
                        }
                    } else {
                        Cached::TooLarge
                    },
                );
            }
            _ => {
                let selected = self
                    .candidates
                    .enabled
                    .then_some(self.candidates.ranges.as_slice());
                index.lookup_selected(
                    key,
                    &mut WeightedScores {
                        scores: &mut self.scores,
                        weight,
                    },
                    &mut self.scratch,
                    prepared,
                    selected,
                );
                #[cfg(test)]
                {
                    self.selected_lookups += 1;
                }
            }
        }
        self.candidates.advance(index, &self.scores, weight);
    }

    fn query_features(&mut self, index: &Index, range: Range<usize>, ratio: f64) -> u64 {
        let total = range.len() as u64;
        self.scores.clear();
        self.candidates.reset(index, total, ratio);
        self.keys[range.clone()].sort_unstable();
        let mut i = range.start;
        // Prepare independent addresses ahead, but commit scores/cache/pruning
        // strictly in the original key order. Future prefetches are only hints.
        const BLOCK: usize = 16;
        const AHEAD: usize = 4;
        let stride = index.banks.len() * index.manifest.hashes;
        let mut positions = std::mem::take(&mut self.positions);
        let mut keys = [0u64; BLOCK];
        let mut weights = [0u64; BLOCK];
        while i < range.end && self.candidates.alive {
            let mut count = 0;
            while count < BLOCK && i < range.end {
                let key = self.keys[i];
                let start = i;
                while i < range.end && self.keys[i] == key {
                    i += 1;
                }
                keys[count] = key;
                weights[count] = (i - start) as u64;
                count += 1;
            }
            positions.resize(count * stride, 0);
            for (&key, out) in keys[..count].iter().zip(positions.chunks_exact_mut(stride)) {
                let hash = PositionHash::new(key, index.manifest.seed);
                for (bank, out) in index
                    .banks
                    .iter()
                    .zip(out.chunks_exact_mut(index.manifest.hashes))
                {
                    hash.positions(
                        bank.meta.slots,
                        index.manifest.hashes,
                        index.manifest.segment_length,
                        out,
                    );
                }
            }
            for n in 0..count.min(AHEAD) {
                index.prefetch_record(
                    keys[n],
                    &positions[n * stride..(n + 1) * stride],
                    self.candidates
                        .enabled
                        .then_some(self.candidates.ranges.as_slice()),
                );
            }
            for n in 0..count {
                if !self.candidates.alive {
                    break;
                }
                if n + AHEAD < count {
                    let next = n + AHEAD;
                    index.prefetch_record(
                        keys[next],
                        &positions[next * stride..(next + 1) * stride],
                        self.candidates
                            .enabled
                            .then_some(self.candidates.ranges.as_slice()),
                    );
                }
                self.query_key_prepared(
                    index,
                    keys[n],
                    weights[n],
                    Some(&positions[n * stride..(n + 1) * stride]),
                );
            }
        }
        self.positions = positions;
        total
    }

    // Very long records do not allocate an unbounded feature list or result
    // cache. The first scan already counted features; rescan in source order.
    fn query_long(
        &mut self,
        index: &Index,
        features: FeatureConfig,
        seq: &[u8],
        total: u64,
        ratio: f64,
    ) {
        self.cache.clear();
        self.origins.clear();
        self.scores.clear();
        self.candidates.reset(index, total, ratio);
        let mut previous = None;
        let mut weight = 0;
        crate::scan_index_features(
            seq,
            index.manifest.k,
            index.manifest.seed,
            features,
            &mut |key| {
                if !self.candidates.alive {
                    return;
                }
                if previous == Some(key) {
                    weight += 1;
                } else {
                    if let Some(previous) = previous {
                        self.query_key(index, previous, weight);
                    }
                    previous = Some(key);
                    weight = 1;
                }
            },
        );
        if self.candidates.alive {
            if let Some(key) = previous {
                self.query_key(index, key, weight);
            }
        }
    }

    pub(super) fn query_records(
        &mut self,
        index: &Index,
        features: FeatureConfig,
        ratio: f64,
        records: &[Record],
        out: &mut Vec<u8>,
    ) -> Result<()> {
        if features.randstrobemers.is_some() {
            return self.query_strobe_records(index, features, ratio, records, out);
        }
        self.query_records_forward(index, features, ratio, records, out)
    }

    /// Merge sorted strand keys for full counts; thresholded queries retain
    /// separate pruning. Neither path needs a second cache-preparation sort.
    fn query_strobe_records(
        &mut self,
        index: &Index,
        features: FeatureConfig,
        ratio: f64,
        records: &[Record],
        out: &mut Vec<u8>,
    ) -> Result<()> {
        let mut first = 0;
        while first < records.len() {
            // Do not perform a speculative expensive extraction before the
            // existing long-record fallback. Small paired lists are bounded
            // by twice the input length even before counting their anchors.
            if records[first].1.len() > FEATURE_BUDGET / 2 {
                let total = self.query_sequence(index, features, ratio, &records[first].1);
                out.extend_from_slice(&records[first].0);
                out.push(b'\t');
                self.scores.write(total, ratio, out)?;
                first += 1;
                continue;
            }
            self.keys.clear();
            self.ends.clear();
            let mut end = first;
            while end < records.len() {
                let seq = &records[end].1;
                if seq.len() > FEATURE_BUDGET / 2 {
                    break;
                }
                let start = self.keys.len();
                let available = FEATURE_BUDGET - start;
                let (middle, stop) = if features.mode == IndexMode::Kmers
                    && features.seed_pattern.is_none()
                    && seq.len() <= available / 2
                {
                    let total = crate::anchors::append_dense_randstrobe_orientations(
                        seq,
                        index.manifest.k,
                        index.manifest.seed,
                        features,
                        &mut self.keys,
                    );
                    (start + total, self.keys.len())
                } else {
                    let mut count = 0;
                    crate::scan_index_features(
                        seq,
                        index.manifest.k,
                        index.manifest.seed,
                        features,
                        &mut |key| {
                            count += 1;
                            if self.keys.len() < FEATURE_BUDGET {
                                self.keys.push(key);
                            }
                        },
                    );
                    if count > available / 2 {
                        self.keys.truncate(start);
                        break;
                    }
                    let middle = self.keys.len();
                    let mut reverse_count = 0;
                    crate::scan_index_features(
                        &crate::reverse_complement(seq),
                        index.manifest.k,
                        index.manifest.seed,
                        features,
                        &mut |key| {
                            reverse_count += 1;
                            if self.keys.len() < FEATURE_BUDGET {
                                self.keys.push(key);
                            }
                        },
                    );
                    assert_eq!(
                        count, reverse_count,
                        "canonical anchor counts must agree across strands"
                    );
                    (middle, self.keys.len())
                };
                self.ends.extend_from_slice(&[middle, stop]);
                end += 1;
                if self.keys.len() == FEATURE_BUDGET {
                    break;
                }
            }
            if end == first {
                // Oversized reads retain the bounded two-pass strand fallback.
                let total = self.query_sequence(index, features, ratio, &records[first].1);
                out.extend_from_slice(&records[first].0);
                out.push(b'\t');
                self.scores.write(total, ratio, out)?;
                first += 1;
                continue;
            }
            self.cache.clear();
            self.origins.clear();
            let mut start = 0;
            for (offset, (id, _)) in records[first..end].iter().enumerate() {
                let middle = self.ends[offset * 2];
                let stop = self.ends[offset * 2 + 1];
                let total = if ratio == 0.0 {
                    self.query_strobe_pair_full(index, start, middle, stop)
                } else {
                    let total = self.query_features(index, start..middle, ratio);
                    let strand = self
                        .strand_scores
                        .get_or_insert_with(|| RecordScores::new(index.manifest.colors.len()));
                    std::mem::swap(&mut self.scores, strand);
                    let reverse_total = self.query_features(index, middle..stop, ratio);
                    assert_eq!(total, reverse_total);
                    self.merge_forward_scores();
                    total
                };
                out.extend_from_slice(id);
                out.push(b'\t');
                self.scores.write(total, ratio, out)?;
                start = stop;
            }
            first = end;
        }
        Ok(())
    }

    fn query_strobe_pair_full(
        &mut self,
        index: &Index,
        start: usize,
        middle: usize,
        stop: usize,
    ) -> u64 {
        assert_eq!(middle - start, stop - middle);
        self.keys[start..middle].sort_unstable();
        self.keys[middle..stop].sort_unstable();
        self.scores.clear();
        let reverse = self
            .strand_scores
            .get_or_insert_with(|| RecordScores::new(index.manifest.colors.len()));
        reverse.clear();
        let (mut a, mut b) = (start, middle);
        while a < middle || b < stop {
            let key = if a == middle {
                self.keys[b]
            } else if b == stop {
                self.keys[a]
            } else {
                self.keys[a].min(self.keys[b])
            };
            let (first_a, first_b) = (a, b);
            while a < middle && self.keys[a] == key {
                a += 1;
            }
            while b < stop && self.keys[b] == key {
                b += 1;
            }
            index.lookup_selected(
                key,
                &mut StrandPairScores {
                    forward: &mut self.scores,
                    reverse: &mut *reverse,
                    forward_weight: (a - first_a) as u64,
                    reverse_weight: (b - first_b) as u64,
                },
                &mut self.scratch,
                None,
                None,
            );
        }
        self.merge_forward_scores();
        (middle - start) as u64
    }

    fn merge_forward_scores(&mut self) {
        let forward = self.strand_scores.as_ref().unwrap();
        for &origin in &forward.touched {
            let extra = forward.values[origin].saturating_sub(self.scores.values[origin]);
            if extra != 0 {
                self.scores.add(origin, extra);
            }
        }
    }

    /// Query both orientations independently, then take the per-document
    /// maximum. One seed per selected anchor in either orientation gives the
    /// same denominator; never pool both strands into one count/denominator.
    pub(super) fn query_sequence(
        &mut self,
        index: &Index,
        features: FeatureConfig,
        ratio: f64,
        seq: &[u8],
    ) -> u64 {
        let total = self.query_strand(index, features, ratio, seq);
        if features.randstrobemers.is_none() {
            return total;
        }
        let strand = self
            .strand_scores
            .get_or_insert_with(|| RecordScores::new(index.manifest.colors.len()));
        std::mem::swap(&mut self.scores, strand);
        let reverse = crate::reverse_complement(seq);
        let reverse_total = self.query_strand(index, features, ratio, &reverse);
        assert_eq!(
            total, reverse_total,
            "canonical anchor counts must agree across strands"
        );
        self.merge_forward_scores();
        total
    }

    fn query_strand(
        &mut self,
        index: &Index,
        features: FeatureConfig,
        ratio: f64,
        seq: &[u8],
    ) -> u64 {
        self.keys.clear();
        let mut total = 0;
        crate::scan_index_features(
            seq,
            index.manifest.k,
            index.manifest.seed,
            features,
            &mut |key| {
                total += 1;
                if self.keys.len() < FEATURE_BUDGET {
                    self.keys.push(key);
                }
            },
        );
        if total > FEATURE_BUDGET as u64 {
            self.query_long(index, features, seq, total, ratio);
        } else {
            // query_features already folds all repeated keys into weights;
            // a cache local to this ONE strand cannot produce any reuse.
            self.cache.clear();
            self.origins.clear();
            self.query_features(index, 0..self.keys.len(), ratio);
        }
        total
    }

    fn query_records_forward(
        &mut self,
        index: &Index,
        features: FeatureConfig,
        ratio: f64,
        records: &[Record],
        out: &mut Vec<u8>,
    ) -> Result<()> {
        let mut first = 0;
        while first < records.len() {
            self.keys.clear();
            self.ends.clear();
            let mut end = first;
            let mut long = None;
            while end < records.len() && self.keys.len() < FEATURE_BUDGET {
                let start = self.keys.len();
                let mut count = 0u64;
                crate::scan_index_features(
                    &records[end].1,
                    index.manifest.k,
                    index.manifest.seed,
                    features,
                    &mut |key| {
                        count += 1;
                        if self.keys.len() < FEATURE_BUDGET {
                            self.keys.push(key);
                        }
                    },
                );
                if count > (FEATURE_BUDGET - start) as u64 {
                    self.keys.truncate(start);
                    if end == first {
                        long = Some(count);
                    }
                    break;
                }
                self.ends.push(self.keys.len());
                end += 1;
            }
            if let Some(total) = long {
                self.query_long(index, features, &records[first].1, total, ratio);
                out.extend_from_slice(&records[first].0);
                out.push(b'\t');
                self.scores.write(total, ratio, out)?;
                first += 1;
                continue;
            }
            self.prepare_cache();
            let mut start = 0;
            for (offset, (id, _)) in records[first..end].iter().enumerate() {
                let next = self.ends[offset];
                let total = self.query_features(index, start..next, ratio);
                out.extend_from_slice(id);
                out.push(b'\t');
                self.scores.write(total, ratio, out)?;
                start = next;
            }
            first = end;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random_index(bits: usize, planes: bool) -> Index {
        let mut index = super::super::tests::union_index(123, &[]);
        index.manifest.fingerprint_bits = bits;
        index.banks.clear();
        let mut id = 0;
        for (bank_id, count) in [73, 35].into_iter().enumerate() {
            let mut groups = Vec::new();
            let mut offset = 0;
            for i in 0..count {
                let n = 1 + i % 3;
                let members = (id..id + n).collect();
                id += n;
                let width = bits + if n == 1 { 0 } else { n };
                groups.push(GroupMeta {
                    members,
                    width,
                    offset,
                    union_cardinality: 100,
                    fallback_keys: 0,
                    abandoned_keys: Some(0),
                });
                offset += width;
            }
            let slots = 32 << bank_id;
            let data = (0..slots * offset.div_ceil(8))
                .map(|i| splitmix64((bank_id * 100_000 + i) as u64) as u8)
                .collect();
            index.banks.push(Bank {
                meta: BankMeta {
                    file: format!("bank-{}.zorbank.zst", slots.ilog2()),
                    slots,
                    row_bits: offset,
                    groups,
                    fingerprint_first: planes,
                },
                data,
                fallback: Vec::new(),
            });
        }
        index
            .manifest
            .colors
            .resize(id, index.manifest.colors[0].clone());
        index.manifest.banks = index.banks.iter().map(|bank| bank.meta.clone()).collect();
        index.row_scratch_bytes = if planes {
            0
        } else {
            index
                .banks
                .iter()
                .map(|bank| bank.meta.row_bits.div_ceil(8))
                .max()
                .unwrap()
        };
        let mut map = OriginMap::default();
        let mut entries: Vec<_> = [123u64, 83_843, 95_900, 456, 999]
            .into_iter()
            .map(|key| abandoned::hash(key, index.manifest.seed))
            .collect();
        entries.sort_unstable();
        entries.dedup();
        for hash in entries {
            map.push(hash, &[0, 1, 101]).unwrap();
        }
        map.finish().unwrap();
        index.abandoned = Some(map);
        index
    }

    fn reference(index: &Index, keys: &[u64], ratio: f64) -> Vec<u8> {
        let mut scores = RecordScores::new(index.manifest.colors.len());
        for &key in keys {
            index.lookup_into(key, &mut scores, &mut Vec::new());
        }
        let mut out = Vec::new();
        scores.write(keys.len() as u64, ratio, &mut out).unwrap();
        out
    }

    #[test]
    fn integer_threshold_matches_original_floating_point_comparison() {
        for total in 1..=151 {
            for ratio in [0.0, f64::MIN_POSITIVE, 0.1, 1.0 / 3.0, 0.5, 0.8, 1.0] {
                let required = minimum_matches(total, ratio);
                for count in 0..=total {
                    assert_eq!(count >= required, count as f64 / total as f64 >= ratio);
                }
            }
        }
        for total in [(1u64 << 53) - 1, (1u64 << 53) + 1, u64::MAX] {
            for ratio in [0.1, 1.0 / 3.0, 0.8, 1.0] {
                let required = minimum_matches(total, ratio);
                assert!(required as f64 / total as f64 >= ratio);
                assert!(required == 0 || ((required - 1) as f64 / total as f64) < ratio);
            }
        }
        assert_eq!(minimum_matches(0, 0.8), 0);
    }

    #[test]
    fn pruning_and_cached_weighted_hits_match_full_lookup_for_all_layouts() {
        for bits in [1, 2, 4, 8, 16] {
            for planes in [false, true] {
                let mut index = random_index(bits, planes);
                for recovery in [true, false] {
                    if !recovery {
                        index.abandoned = None;
                    }
                    let mut worker = Worker::new(&index);
                    for ratio in [0.0, 0.1, 1.0 / 3.0, 0.5, 0.8, 1.0] {
                        for record in 0..12 {
                            let keys: Vec<_> = (0..113)
                                .map(|i| match i % 7 {
                                    0 | 1 => 123,
                                    2 => 83_843,
                                    3 => 95_900,
                                    _ => splitmix64((record * 1000 + i) as u64),
                                })
                                .collect();
                            let expected = reference(&index, &keys, ratio);
                            worker.keys = keys;
                            worker.prepare_cache();
                            let total = worker.query_features(&index, 0..worker.keys.len(), ratio);
                            let mut actual = Vec::new();
                            worker.scores.write(total, ratio, &mut actual).unwrap();
                            assert_eq!(
                                actual, expected,
                                "bits={bits} planes={planes} recovery={recovery} ratio={ratio}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn cache_is_complete_after_first_read_retires_all_its_groups() {
        let index = super::super::tests::union_index(123, &[1]);
        // Find a true nonmatch, not merely an assumed absent key in a filter.
        let bad = (0..123)
            .find(|&key| {
                let mut scores = [0; 2];
                index.lookup(key, &mut scores, &mut Vec::new());
                scores == [0, 0]
            })
            .unwrap();
        let mut worker = Worker::new(&index);
        worker.keys = [vec![bad; 40], vec![123; 60], vec![123; 100], vec![123; 100]].concat();
        worker.prepare_cache();
        let mut start = 0;
        for end in [100, 200, 300] {
            let expected = reference(&index, &worker.keys[start..end], 0.8);
            let total = worker.query_features(&index, start..end, 0.8);
            let mut actual = Vec::new();
            worker.scores.write(total, 0.8, &mut actual).unwrap();
            assert_eq!(actual, expected);
            start = end;
        }
        assert_eq!(worker.full_lookups, 2); // bad once, 123 once
        assert_eq!(worker.cache_hits, 1); // third read reuses the second's result
        assert_eq!(worker.selected_lookups, 0);
    }

    #[test]
    fn unrelated_read_retires_banks_after_first_impossible_threshold() {
        let index = super::super::tests::union_index(123, &[]);
        let keys: Vec<_> = (0..10_000)
            .filter(|&key| {
                let mut scores = [0; 2];
                index.lookup(key, &mut scores, &mut Vec::new());
                scores == [0, 0]
            })
            .take(100)
            .collect();
        let mut worker = Worker::new(&index);
        worker.keys = keys.clone();
        worker.prepare_cache();
        worker.query_features(&index, 0..100, 0.8);
        assert!(!worker.candidates.alive);
        assert!(worker.candidates.ranges.iter().all(Vec::is_empty));
        assert_eq!(worker.selected_lookups, 21);
        let mut actual = Vec::new();
        worker.scores.write(100, 0.8, &mut actual).unwrap();
        assert_eq!(actual, reference(&index, &keys, 0.8));
    }

    #[test]
    fn oversized_cached_postings_fall_back_without_losing_hits() {
        let index = super::super::tests::union_index(123, &[1]);
        let mut worker = Worker::new(&index);
        worker.origins.resize(CACHED_ORIGIN_BUDGET - 1, 0);
        worker.cache.insert(123, Cached::Pending);
        worker.candidates.reset(&index, 2, 0.0);
        worker.query_key(&index, 123, 1);
        assert!(matches!(worker.cache[&123], Cached::TooLarge));
        assert_eq!(worker.origins.len(), CACHED_ORIGIN_BUDGET - 1);
        worker.query_key(&index, 123, 1);
        assert_eq!(worker.scores.values, [2, 2]);
        assert_eq!(worker.full_lookups, 1);
        assert_eq!(worker.selected_lookups, 1);
    }

    #[test]
    fn paired_sink_preserves_weighted_counts_and_abandoned_union() {
        for bits in [1, 2, 4, 8, 16, 32] {
            for planes in [false, true] {
                let index = random_index(bits, planes);
                for (forward_weight, reverse_weight) in
                    [(0, 1), (1, 0), (1, 1), (7, 11), (0, 9), (6, 0)]
                {
                    let mut forward = RecordScores::new(index.manifest.colors.len());
                    let mut reverse = RecordScores::new(index.manifest.colors.len());
                    let mut expected = RecordScores::new(index.manifest.colors.len());
                    for key in [123, 83_843, 95_900, 0, u64::MAX, 123, 83_843] {
                        index.lookup_into(key, &mut expected, &mut Vec::new());
                        index.lookup_selected(
                            key,
                            &mut StrandPairScores {
                                forward: &mut forward,
                                reverse: &mut reverse,
                                forward_weight,
                                reverse_weight,
                            },
                            &mut Vec::new(),
                            None,
                            None,
                        );
                        for origin in 0..expected.values.len() {
                            assert_eq!(
                                forward.values[origin],
                                expected.values[origin] * forward_weight
                            );
                            assert_eq!(
                                reverse.values[origin],
                                expected.values[origin] * reverse_weight
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn record_chunks_and_long_record_fallback_preserve_all_feature_modes() {
        let mut index = random_index(8, true);
        let pattern = ("#_".repeat(15) + "#").parse().unwrap();
        let base = FeatureConfig::new(IndexMode::Kmers, 31, 10);
        let records: Vec<Record> = [
            Vec::new(),
            b"ACGTN".to_vec(),
            (0..1900)
                .map(|i| b"ACGT"[(splitmix64(i) & 3) as usize])
                .collect(),
            b"A".repeat(FEATURE_BUDGET + 35),
            b"ACGT".repeat(FEATURE_BUDGET / 2),
            b"ACGTACGTNNNNACGTACGT".repeat(30),
        ]
        .into_iter()
        .enumerate()
        .map(|(id, seq)| (format!("r{id}").into_bytes(), seq))
        .collect();
        let strobe = crate::RandstrobeConfig {
            min_distance: 35,
            max_distance: 75,
        };
        for features in [
            base,
            base.with_homopolymer_compression(true),
            base.with_seed_pattern(Some(pattern)),
            base.with_randstrobemers(Some(strobe)),
            base.with_randstrobemers(Some(strobe))
                .with_homopolymer_compression(true),
            base.with_randstrobemers(Some(strobe))
                .with_seed_pattern(Some(pattern)),
            FeatureConfig::new(IndexMode::Modimizers, 31, 10),
            FeatureConfig::new(IndexMode::Minimizers, 31, 10)
                .with_minimizer_scheme(MinimizerScheme::OpenClosedMod),
            FeatureConfig::new(IndexMode::Minimizers, 31, 10)
                .with_minimizer_scheme(MinimizerScheme::Simd),
            FeatureConfig::new(IndexMode::Minimizers, 31, 10)
                .with_minimizer_scheme(MinimizerScheme::OpenClosedMod)
                .with_randstrobemers(Some(strobe)),
            FeatureConfig::new(IndexMode::Minimizers, 31, 10)
                .with_minimizer_scheme(MinimizerScheme::Simd)
                .with_randstrobemers(Some(strobe)),
        ] {
            index.manifest.k = if features.mode == IndexMode::Minimizers {
                49
            } else {
                31
            };
            for ratio in [0.0, 0.8, 1.0] {
                let mut expected = Vec::new();
                for (id, seq) in &records {
                    let mut keys = Vec::new();
                    crate::scan_index_features(
                        seq,
                        index.manifest.k,
                        index.manifest.seed,
                        features,
                        &mut |key| keys.push(key),
                    );
                    let mut result = reference(&index, &keys, ratio);
                    if features.randstrobemers.is_some() {
                        let mut scores = RecordScores::new(index.manifest.colors.len());
                        let mut reverse = RecordScores::new(index.manifest.colors.len());
                        for &key in &keys {
                            index.lookup_into(key, &mut scores, &mut Vec::new());
                        }
                        let mut reverse_total = 0;
                        crate::scan_index_features(
                            &crate::reverse_complement(seq),
                            index.manifest.k,
                            index.manifest.seed,
                            features,
                            &mut |key| {
                                reverse_total += 1;
                                index.lookup_into(key, &mut reverse, &mut Vec::new());
                            },
                        );
                        assert_eq!(reverse_total, keys.len());
                        for &origin in &reverse.touched {
                            let extra =
                                reverse.values[origin].saturating_sub(scores.values[origin]);
                            if extra > 0 {
                                scores.add(origin, extra);
                            }
                        }
                        result.clear();
                        scores.write(keys.len() as u64, ratio, &mut result).unwrap();
                    }
                    expected.extend_from_slice(id);
                    expected.push(b'\t');
                    expected.extend(result);
                }
                let mut worker = Worker::new(&index);
                let mut actual = Vec::new();
                worker
                    .query_records(&index, features, ratio, &records, &mut actual)
                    .unwrap();
                assert_eq!(actual, expected, "features={features:?} ratio={ratio}");
                assert!(worker.keys.len() <= FEATURE_BUDGET);
                assert!(worker.origins.len() <= CACHED_ORIGIN_BUDGET);
            }
        }
    }
}
