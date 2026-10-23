//! Four independent overlapping sequence regions per AVX2 vector. Chunk
//! boundaries own window starts, not selected positions; ordered emission and
//! cross-chunk deduplication therefore exactly match the portable iterator.
//! Prefix/suffix block minima avoid data-dependent deque operations. All ranks
//! retain the complete 64-bit hash plus a separate syncmer-priority class.
use super::{inner_hash_table, parameters, scan_run_portable};
use std::arch::x86_64::*;

const TILE_WINDOWS: usize = 1024;

#[target_feature(enable = "avx2")]
pub(super) unsafe fn first_invalid(seq: &[u8]) -> usize {
    let uppercase = _mm256_set1_epi8(!32);
    let mut offset = 0;
    while seq.len() - offset >= 32 {
        let bytes = _mm256_and_si256(
            _mm256_loadu_si256(seq.as_ptr().add(offset).cast()),
            uppercase,
        );
        let ac = _mm256_or_si256(
            _mm256_cmpeq_epi8(bytes, _mm256_set1_epi8(b'A' as i8)),
            _mm256_cmpeq_epi8(bytes, _mm256_set1_epi8(b'C' as i8)),
        );
        let gt = _mm256_or_si256(
            _mm256_cmpeq_epi8(bytes, _mm256_set1_epi8(b'G' as i8)),
            _mm256_cmpeq_epi8(bytes, _mm256_set1_epi8(b'T' as i8)),
        );
        let valid = _mm256_movemask_epi8(_mm256_or_si256(ac, gt)) as u32;
        if valid != u32::MAX {
            return offset + (!valid).trailing_zeros() as usize;
        }
        offset += 32;
    }
    offset
        + seq[offset..]
            .iter()
            .position(|&base| crate::encode_base(base).is_none())
            .unwrap_or(seq.len() - offset)
}

#[derive(Clone, Copy)]
struct Rank {
    class: __m256i,
    hash: __m256i,
    left: __m256i,
    right: __m256i,
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn unsigned_gt(a: __m256i, b: __m256i) -> __m256i {
    let sign = _mm256_set1_epi64x(i64::MIN);
    _mm256_cmpgt_epi64(_mm256_xor_si256(a, sign), _mm256_xor_si256(b, sign))
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn minimum(a: __m256i, b: __m256i) -> __m256i {
    _mm256_blendv_epi8(a, b, unsigned_gt(a, b))
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn combine(a: Rank, b: Rank) -> Rank {
    let class_equal = _mm256_cmpeq_epi64(a.class, b.class);
    let a_less = _mm256_or_si256(
        _mm256_cmpgt_epi64(b.class, a.class),
        _mm256_and_si256(class_equal, unsigned_gt(b.hash, a.hash)),
    );
    let equal = _mm256_and_si256(class_equal, _mm256_cmpeq_epi64(a.hash, b.hash));
    // a is the earlier interval. Leftmost ties come from a; rightmost from b.
    Rank {
        class: _mm256_blendv_epi8(b.class, a.class, a_less),
        hash: _mm256_blendv_epi8(b.hash, a.hash, a_less),
        left: _mm256_blendv_epi8(b.left, a.left, _mm256_or_si256(a_less, equal)),
        right: _mm256_blendv_epi8(b.right, a.right, a_less),
    }
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn multiply(x: __m256i, c: u64) -> __m256i {
    let low = _mm256_set1_epi64x((c as u32) as i64);
    let high = _mm256_set1_epi64x((c >> 32) as i64);
    let cross = _mm256_add_epi64(
        _mm256_mul_epu32(_mm256_srli_epi64::<32>(x), low),
        _mm256_mul_epu32(x, high),
    );
    _mm256_add_epi64(_mm256_mul_epu32(x, low), _mm256_slli_epi64::<32>(cross))
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn hash(mut x: __m256i, seed: u64) -> __m256i {
    x = _mm256_add_epi64(
        x,
        _mm256_set1_epi64x(seed.wrapping_add(0x9E37_79B9_7F4A_7C15) as i64),
    );
    x = multiply(
        _mm256_xor_si256(x, _mm256_srli_epi64::<30>(x)),
        0xBF58_476D_1CE4_E5B9,
    );
    x = multiply(
        _mm256_xor_si256(x, _mm256_srli_epi64::<27>(x)),
        0x94D0_49BB_1331_11EB,
    );
    _mm256_xor_si256(x, _mm256_srli_epi64::<31>(x))
}

struct Window {
    ring: Vec<__m256i>,
    index: usize,
    prefix: __m256i,
}
impl Window {
    #[target_feature(enable = "avx2")]
    unsafe fn new(length: usize) -> Self {
        let maximum = _mm256_set1_epi64x(-1);
        Self {
            ring: vec![maximum; length],
            index: 0,
            prefix: maximum,
        }
    }
    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn push(&mut self, value: __m256i) -> __m256i {
        self.ring[self.index] = value;
        self.prefix = minimum(self.prefix, value);
        self.index += 1;
        if self.index == self.ring.len() {
            self.index = 0;
            for i in (0..self.ring.len() - 1).rev() {
                self.ring[i] = minimum(self.ring[i], self.ring[i + 1]);
            }
            self.prefix = value;
        }
        minimum(self.prefix, self.ring[self.index])
    }
    #[target_feature(enable = "avx2")]
    unsafe fn reset(&mut self) {
        let maximum = _mm256_set1_epi64x(-1);
        self.ring.fill(maximum);
        self.index = 0;
        self.prefix = maximum;
    }
}

struct Anchors {
    ring: Vec<Rank>,
    index: usize,
    prefix: Rank,
    empty: Rank,
}
impl Anchors {
    #[target_feature(enable = "avx2")]
    unsafe fn new(length: usize) -> Self {
        let empty = Rank {
            class: _mm256_set1_epi64x(3),
            hash: _mm256_set1_epi64x(-1),
            left: _mm256_setzero_si256(),
            right: _mm256_setzero_si256(),
        };
        Self {
            ring: vec![empty; length],
            index: 0,
            prefix: empty,
            empty,
        }
    }
    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn push(&mut self, value: Rank) -> Rank {
        self.ring[self.index] = value;
        self.prefix = combine(self.prefix, value);
        self.index += 1;
        if self.index == self.ring.len() {
            self.index = 0;
            for i in (0..self.ring.len() - 1).rev() {
                self.ring[i] = combine(self.ring[i], self.ring[i + 1]);
            }
            self.prefix = value;
        }
        combine(self.ring[self.index], self.prefix)
    }
    fn reset(&mut self) {
        self.ring.fill(self.empty);
        self.index = 0;
        self.prefix = self.empty;
    }
}

/// Caller must check AVX2 support and provide a DNA-only run with len>=k.
#[target_feature(enable = "avx2")]
pub(super) unsafe fn scan_run<F: FnMut(usize, u64)>(
    seq: &[u8],
    k: usize,
    m: usize,
    seed: u64,
    emit: &mut F,
) {
    scan_run_at_threshold(seq, k, m, seed, k.max(48), emit);
}

#[target_feature(enable = "avx2")]
pub(super) unsafe fn scan_run_at_threshold<F: FnMut(usize, u64)>(
    seq: &[u8],
    k: usize,
    m: usize,
    seed: u64,
    minimum_windows: usize,
    emit: &mut F,
) {
    debug_assert!(minimum_windows >= 4);
    debug_assert!(seq.len() >= k && m > 0 && m <= 31 && m <= k && k & 1 != 0);
    let (w, t, s) = parameters(k, m);
    let hashes = inner_hash_table(s, seed);
    let remainders: Vec<i64> = (0..k - t + 1).map(|i| (i % w) as i64).collect();
    let ring_length = k.next_power_of_two();
    let ring_mask = ring_length - 1;
    let zero = _mm256_setzero_si256();
    let mut keys = vec![zero; ring_length];
    let mut inner_hashes = vec![zero; ring_length];
    let mut inner_min = Window::new(t - s + 1);
    let mut anchors = Anchors::new(k - t + 1);
    let mut outputs: [Vec<(usize, u64)>; 4] = std::array::from_fn(|_| Vec::new());
    let mask_m = _mm256_set1_epi64x(((1u64 << (2 * m)) - 1) as i64);
    let mask_t = _mm256_set1_epi64x(((1u64 << (2 * t)) - 1) as i64);
    let mask_s = _mm256_set1_epi64x(((1u64 << (2 * s)) - 1) as i64);
    let rc_shift = _mm256_set1_epi64x((2 * (m - 1)) as i64);
    let t_shift = _mm256_set1_epi64x((2 * (m - t)) as i64);
    let lanes = _mm256_set_epi64x(3, 2, 1, 0);
    let ring_mask_v = _mm256_set1_epi64x(ring_mask as i64);
    let two = _mm256_set1_epi64x(2);
    let three = _mm256_set1_epi64x(3);
    let gt_half = _mm256_set1_epi64x((k / 2) as i64);
    let all_windows = seq.len() - k + 1;
    let mut consumed = 0;
    let mut previous = None;
    while all_windows - consumed >= minimum_windows {
        let owned = ((all_windows - consumed) / 4).min(TILE_WINDOWS);
        let mut forward = zero;
        let mut reverse = zero;
        let mut gt_count = zero;
        let mut previous_vector = _mm256_set1_epi64x(-1);
        inner_min.reset();
        anchors.reset();
        for output in &mut outputs {
            output.clear();
        }
        for end in 0..owned + k - 1 {
            // All four windows are in bounds by the ownership calculation.
            // ASCII DNA maps to A=0,C=1,G=2,T=3 even for lowercase.
            let bytes = _mm256_set_epi64x(
                seq[consumed + 3 * owned + end] as i64,
                seq[consumed + 2 * owned + end] as i64,
                seq[consumed + owned + end] as i64,
                seq[consumed + end] as i64,
            );
            let code = _mm256_and_si256(
                _mm256_xor_si256(_mm256_srli_epi64::<1>(bytes), _mm256_srli_epi64::<2>(bytes)),
                three,
            );
            gt_count = _mm256_add_epi64(gt_count, _mm256_srli_epi64::<1>(code));
            if end >= k {
                let old = _mm256_set_epi64x(
                    seq[consumed + 3 * owned + end - k] as i64,
                    seq[consumed + 2 * owned + end - k] as i64,
                    seq[consumed + owned + end - k] as i64,
                    seq[consumed + end - k] as i64,
                );
                let gt = _mm256_and_si256(
                    _mm256_xor_si256(_mm256_srli_epi64::<2>(old), _mm256_srli_epi64::<3>(old)),
                    _mm256_set1_epi64x(1),
                );
                gt_count = _mm256_sub_epi64(gt_count, gt);
            }
            forward = _mm256_and_si256(
                _mm256_or_si256(_mm256_slli_epi64::<2>(forward), code),
                mask_m,
            );
            reverse = _mm256_or_si256(
                _mm256_srli_epi64::<2>(reverse),
                _mm256_sllv_epi64(_mm256_xor_si256(code, three), rc_shift),
            );
            if end + 1 >= m {
                keys[(end + 1 - m) & ring_mask] = minimum(forward, reverse);
            }
            if end + 1 < s {
                continue;
            }
            let h = _mm256_i64gather_epi64::<8>(
                hashes.as_ptr().cast(),
                _mm256_and_si256(forward, mask_s),
            );
            inner_hashes[(end + 1 - s) & ring_mask] = h;
            let min = inner_min.push(h);
            if end + 1 < t {
                continue;
            }
            let anchor_pos = end + 1 - t;
            let open =
                _mm256_cmpeq_epi64(inner_hashes[(anchor_pos + (t - s) / 2) & ring_mask], min);
            let closed = _mm256_or_si256(
                _mm256_cmpeq_epi64(inner_hashes[anchor_pos & ring_mask], min),
                _mm256_cmpeq_epi64(h, min),
            );
            let class =
                _mm256_andnot_si256(open, _mm256_blendv_epi8(two, _mm256_set1_epi64x(1), closed));
            let anchor_key = minimum(
                _mm256_and_si256(forward, mask_t),
                _mm256_srlv_epi64(reverse, t_shift),
            );
            let pos_v = _mm256_set1_epi64x(anchor_pos as i64);
            let selected = anchors.push(Rank {
                class,
                hash: hash(anchor_key, seed),
                left: pos_v,
                right: pos_v,
            });
            if end + 1 < k {
                continue;
            }
            let start_v = _mm256_set1_epi64x((end + 1 - k) as i64);
            let selected_pos = _mm256_blendv_epi8(
                selected.left,
                selected.right,
                _mm256_cmpgt_epi64(gt_count, gt_half),
            );
            let offset = _mm256_sub_epi64(selected_pos, start_v);
            let remainder = if w == 1 {
                zero
            } else if k - t < w {
                offset
            } else if (k - t) / w == 1 {
                let wv = _mm256_set1_epi64x(w as i64);
                let ge = _mm256_cmpgt_epi64(offset, _mm256_set1_epi64x((w - 1) as i64));
                _mm256_sub_epi64(offset, _mm256_and_si256(ge, wv))
            } else {
                _mm256_i64gather_epi64::<8>(remainders.as_ptr(), offset)
            };
            let pos = _mm256_add_epi64(start_v, remainder);
            let mut changed = (!(_mm256_movemask_pd(_mm256_castsi256_pd(_mm256_cmpeq_epi64(
                pos,
                previous_vector,
            ))) as u32))
                & 15;
            previous_vector = pos;
            // Most windows repeat all four selections. Avoid the representative
            // gather, vector spills and scalar emission bookkeeping in that case.
            if changed == 0 {
                continue;
            }
            let indices = _mm256_add_epi64(
                _mm256_slli_epi64::<2>(_mm256_and_si256(pos, ring_mask_v)),
                lanes,
            );
            let key = _mm256_i64gather_epi64::<8>(keys.as_ptr().cast(), indices);
            let mut positions = [0u64; 4];
            let mut values = [0u64; 4];
            _mm256_storeu_si256(positions.as_mut_ptr().cast(), pos);
            _mm256_storeu_si256(values.as_mut_ptr().cast(), key);
            while changed != 0 {
                let lane = changed.trailing_zeros() as usize;
                changed &= changed - 1;
                let p = positions[lane] as usize;
                outputs[lane].push((consumed + lane * owned + p, values[lane]));
            }
        }
        for output in &outputs {
            for &(pos, key) in output {
                if previous != Some(pos) {
                    emit(pos, key);
                    previous = Some(pos);
                }
            }
        }
        consumed += 4 * owned;
    }
    if consumed < all_windows {
        scan_run_portable(&seq[consumed..], k, m, seed, &mut |pos, key| {
            let pos = consumed + pos;
            if previous != Some(pos) {
                emit(pos, key);
                previous = Some(pos);
            }
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn full_width_vector_hashes_equal_scalar_including_unsigned_high_bits() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }
        // SAFETY: feature check above; fixed-size buffers match vector stores.
        unsafe {
            for seed in [0, 69, u64::MAX] {
                for i in 0..1000 {
                    let keys = [
                        crate::splitmix64(i),
                        crate::splitmix64(i + 1),
                        u64::MAX - i,
                        i,
                    ];
                    let input = _mm256_loadu_si256(keys.as_ptr().cast());
                    let mut actual = [0u64; 4];
                    _mm256_storeu_si256(actual.as_mut_ptr().cast(), hash(input, seed));
                    assert_eq!(actual, keys.map(|key| crate::mixsplit(key, seed)));
                }
            }
        }
    }
}
