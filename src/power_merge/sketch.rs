//! Mergeable mod-bucket sketches of the *indexed feature set*. Following
//! SimdSketch, threshold filtering avoids most random bucket writes. Keeping
//! full minima permits exact sketch unions, including sparse/empty sketches.
use crate::splitmix64;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Sketch {
    minima: Vec<u64>,
}

impl Sketch {
    pub fn new(keys: &[u64], size: usize, seed: u64) -> Self {
        assert!(size.is_power_of_two());
        let mut minima = vec![u64::MAX; size];
        // A full pass for sparse sets. For dense sets this threshold usually
        // fills every bucket; retry with a larger threshold if necessary.
        let fraction =
            (2.0 * (size as f64).ln().max(1.0) * size as f64 / keys.len().max(1) as f64).min(1.0);
        let mut threshold = (fraction * u64::MAX as f64) as u64;
        loop {
            for &key in keys {
                let hash = splitmix64(key ^ seed);
                if hash <= threshold {
                    let bucket = hash as usize & (size - 1);
                    minima[bucket] = minima[bucket].min(hash >> size.trailing_zeros());
                }
            }
            if threshold == u64::MAX || !minima.contains(&u64::MAX) {
                break;
            }
            threshold = threshold.saturating_mul(2).max(1);
        }
        Self { minima }
    }

    pub fn union(&self, other: &Self) -> Self {
        Self {
            minima: self
                .minima
                .iter()
                .zip(&other.minima)
                .map(|(&a, &b)| a.min(b))
                .collect(),
        }
    }

    pub fn jaccard(&self, other: &Self) -> f64 {
        let (matches, empty) = counts(&self.minima, &other.minima);
        let occupied = self.minima.len() - empty;
        if occupied == 0 {
            1.0
        } else {
            (matches - empty) as f64 / occupied as f64
        }
    }
}

fn counts(a: &[u64], b: &[u64]) -> (usize, usize) {
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx2") {
        // SAFETY: runtime AVX2 check; the routine only loads complete chunks.
        return unsafe { counts_avx2(a, b) };
    }
    counts_scalar(a, b)
}

fn counts_scalar(a: &[u64], b: &[u64]) -> (usize, usize) {
    a.iter().zip(b).fold((0, 0), |(matches, empty), (&x, &y)| {
        (
            matches + usize::from(x == y),
            empty + usize::from(x == u64::MAX && y == u64::MAX),
        )
    })
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn counts_avx2(a: &[u64], b: &[u64]) -> (usize, usize) {
    use std::arch::x86_64::*;
    let mut matches = 0;
    let mut empty = 0;
    let sentinel = _mm256_set1_epi64x(-1);
    let end = a.len() / 4 * 4;
    for i in (0..end).step_by(4) {
        let x = _mm256_loadu_si256(a.as_ptr().add(i).cast());
        let y = _mm256_loadu_si256(b.as_ptr().add(i).cast());
        let equal = _mm256_cmpeq_epi64(x, y);
        matches += _mm256_movemask_pd(_mm256_castsi256_pd(equal)).count_ones() as usize;
        let both_empty = _mm256_and_si256(equal, _mm256_cmpeq_epi64(x, sentinel));
        empty += _mm256_movemask_pd(_mm256_castsi256_pd(both_empty)).count_ones() as usize;
    }
    let tail = counts_scalar(&a[end..], &b[end..]);
    (matches + tail.0, empty + tail.1)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn unions_and_sparse_empty_buckets_are_correct() {
        for size in [64, 4096] {
            let a: Vec<_> = (0..10_000).collect();
            let b: Vec<_> = (5000..15_000).collect();
            let all: Vec<_> = (0..15_000).collect();
            let sa = Sketch::new(&a, size, 69);
            let sb = Sketch::new(&b, size, 69);
            assert_eq!(sa.union(&sb).minima, Sketch::new(&all, size, 69).minima);
            assert_eq!(sa.jaccard(&sa), 1.0);
            assert!((sa.jaccard(&sb) - 1.0 / 3.0).abs() < 0.15);
            let empty = Sketch::new(&[], size, 69);
            assert_eq!(empty.jaccard(&sa), 0.0);
            assert_eq!(empty.jaccard(&empty), 1.0);
            assert_eq!(
                counts(&sa.minima, &sb.minima),
                counts_scalar(&sa.minima, &sb.minima)
            );
        }
    }
}
