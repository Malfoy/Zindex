//! Feature-ID and merge extensions to Ragnar Groot Koerkamp's SimdSketch.
//! Compare bytes, but retain full minima for exact unions (never min the bytes).
use super::BucketSketch;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MergeableSketch {
    minima: Vec<u32>,
    fingerprints: Vec<u8>,
}

impl MergeableSketch {
    pub fn new(keys: &[u64], size: usize, seed: u64) -> Self {
        assert!(size.is_power_of_two() && size >= 64);
        let mut minima = vec![u32::MAX; size];
        if !keys.is_empty() {
            let fraction = (2.0 * (size as f64).ln() * size as f64 / keys.len() as f64).min(1.0);
            let mut bound = (fraction * u32::MAX as f64) as u32;
            let mut selected = Vec::with_capacity((size * 32).min(keys.len()).max(8));
            loop {
                collect(keys, seed, bound, &mut selected);
                for &hash in &selected {
                    let bucket = hash as usize & (size - 1);
                    minima[bucket] = minima[bucket].min(hash);
                }
                if bound == u32::MAX || !minima.contains(&u32::MAX) {
                    break;
                }
                bound = bound.saturating_mul(2).max(1);
            }
        }
        let fingerprints = minima
            .iter()
            .map(|x| (x >> size.trailing_zeros()) as u8)
            .collect();
        Self {
            minima,
            fingerprints,
        }
    }

    pub fn jaccard(&self, other: &Self) -> f64 {
        debug_assert!(!self.is_released() && !other.is_released());
        BucketSketch::dense_b8_similarity(&self.fingerprints, &other.fingerprints) as f64
    }

    /// Dense bytes for a planner's packed/transposed comparison arena.
    pub fn fingerprints(&self) -> &[u8] {
        &self.fingerprints
    }

    /// Keep the estimator's f32 rounding identical for row- and column-wise scans.
    pub fn similarity_from_matches(matches: usize, size: usize) -> f64 {
        ((256.0 * matches as f64 / size as f64 - 1.0).max(0.0) / 255.0) as f32 as f64
    }

    pub fn merge_in_place(&mut self, other: &Self) {
        debug_assert!(!self.is_released() && self.minima.len() == other.minima.len());
        merge_minima(&mut self.minima, &other.minima, &mut self.fingerprints);
    }

    /// Actually free allocations, rather than clear() and retain capacity.
    pub fn release(&mut self) {
        self.minima = Vec::new();
        self.fingerprints = Vec::new();
    }
    pub fn is_released(&self) -> bool {
        self.minima.is_empty()
    }
    pub fn allocated_bytes(&self) -> usize {
        self.minima.capacity() * 4 + self.fingerprints.capacity()
    }
}

#[inline]
fn feature_hash(key: u64, seed: u64) -> u32 {
    let mut x = (key ^ seed).wrapping_add(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
    x ^= x >> 31;
    (x >> 32) as u32
}
fn collect(keys: &[u64], seed: u64, bound: u32, out: &mut Vec<u32>) {
    out.clear();
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx2") {
        // SAFETY: checked CPU support; kernel uses only allocated ranges.
        unsafe {
            collect_avx2(keys, seed, bound, out);
        }
        return;
    }
    collect_scalar(keys, seed, bound, out);
}
fn collect_scalar(keys: &[u64], seed: u64, bound: u32, out: &mut Vec<u32>) {
    out.extend(
        keys.iter()
            .map(|&key| feature_hash(key, seed))
            .filter(|&x| x <= bound),
    );
}

pub(crate) fn byte_matches(a: &[u8], b: &[u8]) -> usize {
    debug_assert_eq!(a.len(), b.len());
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx2") {
        // SAFETY: checked CPU support; tails handled explicitly.
        return unsafe { byte_matches_avx2(a, b) };
    }
    byte_matches_scalar(a, b)
}
fn byte_matches_scalar(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b).filter(|(a, b)| a == b).count()
}
fn merge_minima(a: &mut [u32], b: &[u32], fp: &mut [u8]) {
    let shift = a.len().trailing_zeros();
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx2") {
        // SAFETY: equal >=64 power-of-two lengths are enforced by the caller.
        unsafe {
            merge_avx2(a, b, fp, shift);
        }
        return;
    }
    for ((a, &b), fp) in a.iter_mut().zip(b).zip(fp) {
        *a = (*a).min(b);
        *fp = (*a >> shift) as u8;
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn collect_avx2(keys: &[u64], seed: u64, bound: u32, out: &mut Vec<u32>) {
    use std::arch::x86_64::*;
    unsafe {
        #[inline]
        #[target_feature(enable = "avx2")]
        unsafe fn mul64(x: __m256i, c: u64) -> __m256i {
            let cv = _mm256_set1_epi64x(c as i64);
            let lo = _mm256_mul_epu32(x, cv);
            let cross = _mm256_add_epi64(
                _mm256_mul_epu32(_mm256_srli_epi64::<32>(x), cv),
                _mm256_mul_epu32(x, _mm256_srli_epi64::<32>(cv)),
            );
            _mm256_add_epi64(lo, _mm256_slli_epi64::<32>(cross))
        }
        #[inline]
        #[target_feature(enable = "avx2")]
        unsafe fn hashes(ptr: *const u64, seed: u64) -> __m128i {
            unsafe {
                let x = _mm256_xor_si256(
                    _mm256_loadu_si256(ptr.cast()),
                    _mm256_set1_epi64x(seed as i64),
                );
                let x = _mm256_add_epi64(x, _mm256_set1_epi64x(0x9e3779b97f4a7c15u64 as i64));
                let x = mul64(
                    _mm256_xor_si256(x, _mm256_srli_epi64::<30>(x)),
                    0xbf58476d1ce4e5b9,
                );
                let x = mul64(
                    _mm256_xor_si256(x, _mm256_srli_epi64::<27>(x)),
                    0x94d049bb133111eb,
                );
                let x = _mm256_xor_si256(x, _mm256_srli_epi64::<31>(x));
                let x = _mm256_permutevar8x32_epi32(x, _mm256_setr_epi32(1, 3, 5, 7, 0, 0, 0, 0));
                _mm256_castsi256_si128(x)
            }
        }
        let mut write = 0;
        let sign = _mm256_set1_epi32(i32::MIN);
        let limit = _mm256_xor_si256(_mm256_set1_epi32(bound as i32), sign);
        let end = keys.len() / 8 * 8;
        for i in (0..end).step_by(8) {
            let values = _mm256_inserti128_si256::<1>(
                _mm256_castsi128_si256(hashes(keys.as_ptr().add(i), seed)),
                hashes(keys.as_ptr().add(i + 4), seed),
            );
            let mask = _mm256_xor_si256(
                _mm256_cmpgt_epi32(_mm256_xor_si256(values, sign), limit),
                _mm256_set1_epi32(-1),
            );
            if write + 8 > out.capacity() {
                out.reserve(out.capacity() + 8);
            }
            // Upstream Ragnar lane-compaction implementation.
            super::intrinsics::append_from_mask(
                std::mem::transmute::<__m256i, packed_seq::u32x8>(values),
                std::mem::transmute::<__m256i, packed_seq::u32x8>(mask),
                out,
                &mut write,
            );
        }
        out.set_len(write);
        collect_scalar(&keys[end..], seed, bound, out);
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn byte_matches_avx2(a: &[u8], b: &[u8]) -> usize {
    use std::arch::x86_64::*;
    unsafe {
        // Dense-only comparison: no mask storage, reads, or empty checks.
        // Four independent counters hide the popcount/add dependency chain.
        let mut counts = [0usize; 4];
        let end = a.len() / 128 * 128;
        for i in (0..end).step_by(128) {
            for (lane, count) in counts.iter_mut().enumerate() {
                let offset = i + lane * 32;
                let equal = _mm256_cmpeq_epi8(
                    _mm256_loadu_si256(a.as_ptr().add(offset).cast()),
                    _mm256_loadu_si256(b.as_ptr().add(offset).cast()),
                );
                *count += (_mm256_movemask_epi8(equal) as u32).count_ones() as usize;
            }
        }
        let mut matched = counts.into_iter().sum::<usize>();
        let vectors_end = a.len() / 32 * 32;
        for i in (end..vectors_end).step_by(32) {
            let equal = _mm256_cmpeq_epi8(
                _mm256_loadu_si256(a.as_ptr().add(i).cast()),
                _mm256_loadu_si256(b.as_ptr().add(i).cast()),
            );
            matched += (_mm256_movemask_epi8(equal) as u32).count_ones() as usize;
        }
        matched += a[vectors_end..]
            .iter()
            .zip(&b[vectors_end..])
            .filter(|(a, b)| a == b)
            .count();

        matched
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn merge_avx2(a: &mut [u32], b: &[u32], fp: &mut [u8], shift: u32) {
    use std::arch::x86_64::*;
    unsafe {
        let shuffle = _mm256_setr_epi8(
            0, 4, 8, 12, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 4, 8, 12, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -1,
        );
        for i in (0..a.len()).step_by(8) {
            let m = _mm256_min_epu32(
                _mm256_loadu_si256(a.as_ptr().add(i).cast()),
                _mm256_loadu_si256(b.as_ptr().add(i).cast()),
            );
            _mm256_storeu_si256(a.as_mut_ptr().add(i).cast(), m);
            let low = _mm256_shuffle_epi8(
                _mm256_srl_epi32(m, _mm_cvtsi32_si128(shift as i32)),
                shuffle,
            );
            let bytes = _mm_unpacklo_epi32(
                _mm256_castsi256_si128(low),
                _mm256_extracti128_si256::<1>(low),
            );
            _mm_storel_epi64(fp.as_mut_ptr().add(i).cast(), bytes);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exact_unions_and_release() {
        for size in [64, 4096] {
            let mut all = Vec::new();
            let mut union = MergeableSketch::new(&[], size, 69);
            for (start, len) in [
                (0, 0),
                (0, 1),
                (10, 7),
                (100, 300),
                (0, 10000),
                (5000, 10000),
            ] {
                let keys: Vec<_> = (start..start + len).collect();
                let sketch = MergeableSketch::new(&keys, size, 69);
                let ptr = union.minima.as_ptr();
                let fp_ptr = union.fingerprints.as_ptr();
                union.merge_in_place(&sketch);
                assert_eq!(ptr, union.minima.as_ptr());
                assert_eq!(fp_ptr, union.fingerprints.as_ptr());
                all.extend(keys);
                all.sort_unstable();
                all.dedup();
                assert_eq!(union, MergeableSketch::new(&all, size, 69));
                assert_eq!(union.jaccard(&union), 1.0);
            }
            union.release();
            assert!(union.is_released());
            assert_eq!(union.allocated_bytes(), 0);
        }
    }
    #[test]
    fn collectors_match_scalar() {
        let keys: Vec<_> = (0..2057)
            .map(|i: u64| i.wrapping_mul(0xa4093822299f31d0))
            .collect();
        for len in [0, 1, 7, 8, 9, 31, 32, 33, 2057] {
            for seed in [0, 69, u64::MAX] {
                for bound in [0, 1, 256, 1 << 24, 1 << 31, u32::MAX - 1, u32::MAX] {
                    let mut a = Vec::new();
                    let mut b = Vec::new();
                    collect_scalar(&keys[..len], seed, bound, &mut a);
                    collect(&keys[..len], seed, bound, &mut b);
                    assert_eq!(a, b, "len={len} seed={seed} bound={bound}");
                }
            }
        }
    }
    #[test]
    fn byte_matches_scalar_and_all_256_fingerprints_are_valid() {
        for len in [0usize, 1, 31, 32, 33, 63, 64, 65, 128, 4096] {
            let a: Vec<_> = (0..len).map(|i| (i * 37) as u8).collect();
            let b: Vec<_> = (0..len)
                .map(|i| if i % 3 == 0 { a[i] } else { 255 })
                .collect();
            assert_eq!(byte_matches(&a, &b), byte_matches_scalar(&a, &b));
        }
        for value in 0..=255u8 {
            assert_eq!(byte_matches(&[value; 64], &[value; 64]), 64);
            assert_eq!(
                BucketSketch::dense_b8_similarity(&[value; 64], &[value; 64]),
                1.0
            );
        }
    }
    #[test]
    fn accuracy_and_memory() {
        let a: Vec<_> = (0..100_000).collect();
        let b: Vec<_> = (50_000..150_000).collect();
        let sa = MergeableSketch::new(&a, 4096, 69);
        let sb = MergeableSketch::new(&b, 4096, 69);
        assert!((sa.jaccard(&sb) - 1.0 / 3.0).abs() < 0.04);
        assert_eq!(sa.allocated_bytes(), 4096 * 5);
    }
}
