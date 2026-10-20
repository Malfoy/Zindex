//! One reusable arena of transposed 8-bit fingerprints. Each block contains
//! 32 groups: fingerprint zero of all groups, fingerprint one of all groups, ...
//! Full minima remain cold union state; comparisons never follow their pointers.
pub(super) const LANES: usize = 32;

#[derive(Default)]
pub(super) struct PackedSketches {
    pub bytes: Vec<u8>,
    pub ids: Vec<usize>,
    free: Vec<usize>,
    pub size: usize,
}

impl PackedSketches {
    pub fn allocate(&mut self, size: usize) -> usize {
        assert!(
            self.size == 0 || self.size == size,
            "inconsistent sketch sizes"
        );
        self.size = size;
        self.free.pop().unwrap_or_else(|| {
            let block = self.ids.len() / LANES;
            self.ids.resize(self.ids.len() + LANES, usize::MAX);
            self.bytes.resize(self.bytes.len() + size * LANES, 0);
            block
        })
    }

    pub fn write(&mut self, slot: usize, id: usize, fp: &[u8]) {
        assert_eq!(fp.len(), self.size);
        self.ids[slot] = id;
        let base = slot / LANES * self.size * LANES + slot % LANES;
        for (i, &byte) in fp.iter().enumerate() {
            self.bytes[base + i * LANES] = byte;
        }
    }

    pub fn move_slot(&mut self, from: usize, to: usize) -> usize {
        let id = self.ids[from];
        if from != to {
            self.ids[to] = id;
            let a = from / LANES * self.size * LANES + from % LANES;
            let b = to / LANES * self.size * LANES + to % LANES;
            for row in 0..self.size {
                self.bytes[b + row * LANES] = self.bytes[a + row * LANES];
            }
        }
        id
    }

    pub fn free(&mut self, block: usize) {
        self.free.push(block);
    }

    pub fn block(&self, block: usize) -> &[u8] {
        &self.bytes[block * self.size * LANES..(block + 1) * self.size * LANES]
    }
}

/// Return full match counts only for lanes that can reach `required`. No
/// fingerprint value is reserved; the tail mask describes groups, not cells.
pub(super) fn matches(
    block: &[u8],
    query: &[u8],
    required: usize,
    valid: usize,
    counts: &mut [u32; LANES],
) -> u32 {
    debug_assert_eq!(block.len(), query.len() * LANES);
    debug_assert!(valid <= LANES);
    if required > query.len() || valid == 0 {
        return 0;
    }
    let mask = u32::MAX >> (LANES - valid);
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx2") {
        // SAFETY: AVX2 checked; every vector is one complete allocated row.
        return unsafe { matches_avx2(block, query, required, mask, counts) };
    }
    counts.fill(0);
    let mut live = mask;
    for (row, &q) in query.iter().enumerate() {
        for (lane, count) in counts.iter_mut().enumerate().take(valid) {
            *count += u32::from(block[row * LANES + lane] == q);
            if *count as usize + query.len() - row - 1 < required {
                live &= !(1 << lane);
            }
        }
        if live == 0 {
            return 0;
        }
    }
    live
}

/// Two independent vectors share each query broadcast and loop bookkeeping.
/// This keeps the compact 32-column arena blocks while scanning 64 groups.
pub(super) fn matches_pair(
    a: &[u8],
    b: &[u8],
    query: &[u8],
    required: usize,
    counts: &mut [u32; 2 * LANES],
) -> u64 {
    debug_assert_eq!(a.len(), query.len() * LANES);
    debug_assert_eq!(b.len(), a.len());
    if required > query.len() {
        return 0;
    }
    #[cfg(target_arch = "x86_64")]
    if query.len() <= 128 && std::is_x86_feature_detected!("avx2") {
        // SAFETY: complete rows, equal sketch lengths, checked CPU support.
        return unsafe { matches_pair_avx2(a, b, query, required, counts) };
    }
    let (low, high) = counts.split_at_mut(LANES);
    let lo = matches(a, query, required, LANES, low.try_into().unwrap());
    let hi = matches(b, query, required, LANES, high.try_into().unwrap());
    lo as u64 | ((hi as u64) << LANES)
}

/// Four independent vectors amortize broadcasts and loop bookkeeping across
/// 128 groups. The arena still consists of reusable 32-column blocks.
pub(super) fn matches_quad(
    blocks: [&[u8]; 4],
    query: &[u8],
    required: usize,
    counts: &mut [u32; 4 * LANES],
) -> u128 {
    debug_assert!(blocks.iter().all(|b| b.len() == query.len() * LANES));
    if required > query.len() {
        return 0;
    }
    #[cfg(target_arch = "x86_64")]
    if query.len() <= 128 && std::is_x86_feature_detected!("avx2") {
        // SAFETY: full equal-length rows, checked CPU support and byte counters
        // cannot overflow for the dispatched sketch lengths.
        return unsafe { matches_quad_avx2(blocks, query, required, counts) };
    }
    let (low, high) = counts.split_at_mut(2 * LANES);
    let lo = matches_pair(
        blocks[0],
        blocks[1],
        query,
        required,
        low.try_into().unwrap(),
    );
    let hi = matches_pair(
        blocks[2],
        blocks[3],
        query,
        required,
        high.try_into().unwrap(),
    );
    lo as u128 | ((hi as u128) << (2 * LANES))
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn matches_quad_avx2(
    blocks: [&[u8]; 4],
    query: &[u8],
    required: usize,
    counts: &mut [u32; 4 * LANES],
) -> u128 {
    use std::arch::x86_64::*;
    unsafe {
        let mut sums = [_mm256_setzero_si256(); 4];
        let mut row = 0;
        // Check when three prefix matches would be needed. Checks needing
        // only one or two usually leave a live lane in a wide random block.
        // This changes when we check, never the bound or required match count.
        let mut next_check = (query.len() - required + 3).min(query.len());
        loop {
            while row < next_check {
                let q = _mm256_set1_epi8(*query.get_unchecked(row) as i8);
                for lane in 0..4 {
                    sums[lane] = _mm256_sub_epi8(
                        sums[lane],
                        _mm256_cmpeq_epi8(
                            _mm256_loadu_si256(blocks[lane].as_ptr().add(row * LANES).cast()),
                            q,
                        ),
                    );
                }
                row += 1;
            }
            let bound = _mm256_set1_epi8(required.saturating_sub(query.len() - row) as i8);
            let mut live = 0u128;
            for (lane, &sum) in sums.iter().enumerate() {
                let mask = _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_max_epu8(sum, bound), sum))
                    as u32;
                live |= (mask as u128) << (lane * LANES);
            }
            if live == 0 {
                return 0;
            }
            if row == query.len() {
                let mut bytes = [0u8; 4 * LANES];
                for (lane, &sum) in sums.iter().enumerate() {
                    _mm256_storeu_si256(bytes.as_mut_ptr().add(lane * LANES).cast(), sum);
                }
                for i in 0..4 * LANES {
                    counts[i] = bytes[i] as u32;
                }
                return live;
            }
            next_check = (next_check + 2).min(query.len());
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn matches_pair_avx2(
    a: &[u8],
    b: &[u8],
    query: &[u8],
    required: usize,
    counts: &mut [u32; 2 * LANES],
) -> u64 {
    use std::arch::x86_64::*;
    unsafe {
        let mut low = _mm256_setzero_si256();
        let mut high = _mm256_setzero_si256();
        let mut row = 0;
        let mut next_check = (query.len() - required + 3).min(query.len());
        loop {
            while row < next_check {
                let q = _mm256_set1_epi8(*query.get_unchecked(row) as i8);
                low = _mm256_sub_epi8(
                    low,
                    _mm256_cmpeq_epi8(_mm256_loadu_si256(a.as_ptr().add(row * LANES).cast()), q),
                );
                high = _mm256_sub_epi8(
                    high,
                    _mm256_cmpeq_epi8(_mm256_loadu_si256(b.as_ptr().add(row * LANES).cast()), q),
                );
                row += 1;
            }
            let bound = _mm256_set1_epi8(required.saturating_sub(query.len() - row) as i8);
            let lo =
                _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_max_epu8(low, bound), low)) as u32;
            let hi =
                _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_max_epu8(high, bound), high)) as u32;
            let live = lo as u64 | ((hi as u64) << LANES);
            if live == 0 {
                return 0;
            }
            if row == query.len() {
                let mut bytes = [0u8; 2 * LANES];
                _mm256_storeu_si256(bytes.as_mut_ptr().cast(), low);
                _mm256_storeu_si256(bytes.as_mut_ptr().add(LANES).cast(), high);
                for i in 0..2 * LANES {
                    counts[i] = bytes[i] as u32;
                }
                return live;
            }
            next_check = (next_check + 2).min(query.len());
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn matches_avx2(
    block: &[u8],
    query: &[u8],
    required: usize,
    mask: u32,
    counts: &mut [u32; LANES],
) -> u32 {
    use std::arch::x86_64::*;
    unsafe {
        if query.len() > u16::MAX as usize {
            // The CLI permits up to a million cells. Keep those comparisons
            // vectorized too, without overflowing the compact word counters.
            counts.fill(0);
            let mut live = mask;
            for start in (0..query.len()).step_by(16384) {
                let end = (start + 16384).min(query.len());
                let mut partial = [0u32; LANES];
                matches_avx2(
                    &block[start * LANES..end * LANES],
                    &query[start..end],
                    0,
                    mask,
                    &mut partial,
                );
                for lane in (0..LANES).step_by(8) {
                    let a = _mm256_loadu_si256(counts.as_ptr().add(lane).cast());
                    let b = _mm256_loadu_si256(partial.as_ptr().add(lane).cast());
                    _mm256_storeu_si256(
                        counts.as_mut_ptr().add(lane).cast(),
                        _mm256_add_epi32(a, b),
                    );
                }
                live = mask;
                for (lane, &count) in counts.iter().enumerate() {
                    if count as usize + query.len() - end < required {
                        live &= !(1 << lane);
                    }
                }
                if live == 0 {
                    break;
                }
            }
            return live;
        }
        if query.len() <= 128 {
            let mut sum = _mm256_setzero_si256();
            let mut row = 0;
            // Delay the first check until three matches would be needed, then
            // check every two rows. Every rejection still uses the exact bound.
            let mut next_check = (query.len() - required + 3).min(query.len());
            loop {
                while row < next_check {
                    let equal = _mm256_cmpeq_epi8(
                        _mm256_loadu_si256(block.as_ptr().add(row * LANES).cast()),
                        _mm256_set1_epi8(*query.get_unchecked(row) as i8),
                    );
                    sum = _mm256_sub_epi8(sum, equal);
                    row += 1;
                }
                let need = required.saturating_sub(query.len() - row);
                let bound = _mm256_set1_epi8(need as i8);
                let live = _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_max_epu8(sum, bound), sum))
                    as u32
                    & mask;
                if live == 0 {
                    return 0;
                }
                if row == query.len() {
                    let mut bytes = [0u8; LANES];
                    _mm256_storeu_si256(bytes.as_mut_ptr().cast(), sum);
                    for i in 0..LANES {
                        counts[i] = bytes[i] as u32;
                    }
                    return live;
                }
                next_check = (next_check + 2).min(query.len());
            }
        }
        // Long sketches use 16-bit counters, accumulating bytes in chunks of
        // at most 64 rows to avoid overflow without widening each comparison.
        let mut low = _mm256_setzero_si256();
        let mut high = _mm256_setzero_si256();
        for start in (0..query.len()).step_by(64) {
            let end = (start + 64).min(query.len());
            let mut sum = _mm256_setzero_si256();
            for row in start..end {
                let equal = _mm256_cmpeq_epi8(
                    _mm256_loadu_si256(block.as_ptr().add(row * LANES).cast()),
                    _mm256_set1_epi8(*query.get_unchecked(row) as i8),
                );
                sum = _mm256_sub_epi8(sum, equal);
            }
            low = _mm256_add_epi16(low, _mm256_cvtepu8_epi16(_mm256_castsi256_si128(sum)));
            high = _mm256_add_epi16(
                high,
                _mm256_cvtepu8_epi16(_mm256_extracti128_si256::<1>(sum)),
            );
            let need = required.saturating_sub(query.len() - end);
            let bound = _mm256_set1_epi16(need as i16);
            let lo = _mm256_cmpeq_epi16(_mm256_max_epu16(low, bound), low);
            let hi = _mm256_cmpeq_epi16(_mm256_max_epu16(high, bound), high);
            if ((_mm256_movemask_epi8(lo) as u32 & expand_mask(mask as u16))
                | (_mm256_movemask_epi8(hi) as u32 & expand_mask((mask >> 16) as u16)))
                == 0
            {
                return 0;
            }
        }
        let mut words = [0u16; LANES];
        _mm256_storeu_si256(words.as_mut_ptr().cast(), low);
        _mm256_storeu_si256(words.as_mut_ptr().add(16).cast(), high);
        for (count, word) in counts.iter_mut().zip(words) {
            *count = word as u32;
        }
        let mut live = 0;
        for (i, &count) in counts.iter().enumerate() {
            if count as usize >= required {
                live |= 1 << i;
            }
        }
        live & mask
    }
}

#[cfg(target_arch = "x86_64")]
fn expand_mask(mask: u16) -> u32 {
    let mut x = mask as u32;
    x = (x | (x << 8)) & 0x00ff00ff;
    x = (x | (x << 4)) & 0x0f0f0f0f;
    x = (x | (x << 2)) & 0x33333333;
    x = (x | (x << 1)) & 0x55555555;
    x | (x << 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transposed_simd_matches_scalar_with_exact_thresholds_and_tails() {
        for size in [64, 128, 256, 4096, 65536] {
            let query: Vec<_> = (0..size).map(|i| i as u8).collect();
            let mut block = vec![0u8; size * LANES];
            let mut expected = [0u32; LANES];
            for row in 0..size {
                for lane in 0..LANES {
                    let value = if row % LANES <= lane {
                        query[row]
                    } else {
                        query[row].wrapping_add(1)
                    };
                    block[row * LANES + lane] = value;
                    expected[lane] += u32::from(value == query[row]);
                }
            }
            // All 256 fingerprint values occur, including zero and 255.
            for valid in [0, 1, 7, 16, 31, 32] {
                for required in [0, 1, size / 2, size - 1, size, size + 1] {
                    let mut counts = [u32::MAX; LANES];
                    let live = matches(&block, &query, required, valid, &mut counts);
                    for lane in 0..LANES {
                        let pass = lane < valid && expected[lane] as usize >= required;
                        assert_eq!(
                            live >> lane & 1 != 0,
                            pass,
                            "size={size} valid={valid} required={required} lane={lane}"
                        );
                        if pass {
                            assert_eq!(counts[lane], expected[lane]);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn packed_arena_overwrites_reuses_and_does_not_reserve_fingerprints() {
        let mut arena = PackedSketches::default();
        let a = arena.allocate(256);
        let b = arena.allocate(256);
        let forward: Vec<_> = (0..=255u8).collect();
        let backward: Vec<_> = forward.iter().rev().copied().collect();
        arena.write(a * LANES + 3, 73, &forward);
        arena.write(b * LANES + 31, 12, &backward);
        assert_eq!(arena.move_slot(b * LANES + 31, a * LANES + 3), 12);
        for row in 0..256 {
            assert_eq!(arena.block(a)[row * LANES + 3], backward[row]);
        }
        let bytes = arena.bytes.len();
        arena.free(b);
        assert_eq!(arena.allocate(256), b);
        assert_eq!(arena.bytes.len(), bytes);
        arena.write(b * LANES, 99, &forward);
        assert_eq!(arena.ids[b * LANES], 99);
    }

    #[test]
    fn paired_kernel_has_exact_counts_and_thresholds() {
        for size in [64, 128, 256, 4096] {
            let query: Vec<_> = (0..size).map(|i| i as u8).collect();
            let mut a = vec![0u8; size * LANES];
            let mut b = a.clone();
            let mut expected = [0u32; 2 * LANES];
            for row in 0..size {
                for lane in 0..LANES {
                    a[row * LANES + lane] = if row % 32 <= lane {
                        query[row]
                    } else {
                        query[row].wrapping_add(1)
                    };
                    b[row * LANES + lane] = if row % 32 >= lane {
                        query[row]
                    } else {
                        query[row].wrapping_sub(1)
                    };
                    expected[lane] += u32::from(a[row * LANES + lane] == query[row]);
                    expected[LANES + lane] += u32::from(b[row * LANES + lane] == query[row]);
                }
            }
            for required in 0..=size + 1 {
                let mut counts = [u32::MAX; 2 * LANES];
                let live = matches_pair(&a, &b, &query, required, &mut counts);
                for lane in 0..2 * LANES {
                    let pass = expected[lane] as usize >= required;
                    assert_eq!(live >> lane & 1 != 0, pass);
                    if pass {
                        assert_eq!(counts[lane], expected[lane]);
                    }
                }
            }
        }
    }

    #[test]
    fn quadruple_kernel_matches_exact_counts_at_every_threshold() {
        for size in [64, 128, 256, 4096] {
            let query: Vec<_> = (0..size).map(|i| i as u8).collect();
            let mut blocks: [Vec<u8>; 4] = std::array::from_fn(|_| vec![0; size * LANES]);
            let mut expected = [0u32; 4 * LANES];
            for row in 0..size {
                for lane in 0..4 * LANES {
                    let wanted = lane * size / (4 * LANES - 1);
                    let value = if (row * 17 + lane * 5) % size < wanted {
                        query[row]
                    } else {
                        query[row].wrapping_add(1)
                    };
                    blocks[lane / LANES][row * LANES + lane % LANES] = value;
                    expected[lane] += u32::from(value == query[row]);
                }
            }
            let thresholds: Vec<_> = if size <= 128 {
                (0..=size + 1).collect()
            } else {
                vec![0, 1, size / 2, size - 1, size, size + 1]
            };
            for required in thresholds {
                let mut counts = [u32::MAX; 4 * LANES];
                let live = matches_quad(
                    std::array::from_fn(|i| blocks[i].as_slice()),
                    &query,
                    required,
                    &mut counts,
                );
                for lane in 0..4 * LANES {
                    let pass = expected[lane] as usize >= required;
                    assert_eq!(
                        live >> lane & 1 != 0,
                        pass,
                        "size={size} threshold={required} lane={lane}"
                    );
                    if pass {
                        assert_eq!(counts[lane], expected[lane]);
                    }
                }
            }
            // Exercise the early-all-rejected path too, not only surviving lanes.
            for block in &mut blocks {
                for (row, &q) in query.iter().enumerate() {
                    block[row * LANES..(row + 1) * LANES].fill(q.wrapping_add(1));
                }
            }
            let mut counts = [u32::MAX; 4 * LANES];
            assert_eq!(
                matches_quad(
                    std::array::from_fn(|i| blocks[i].as_slice()),
                    &query,
                    size / 2,
                    &mut counts
                ),
                0
            );
        }
    }
}
