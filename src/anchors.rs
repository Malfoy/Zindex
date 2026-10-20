//! Optional sequence anchors layered on ordinary k-mer/minimizer extraction.
//! Minimizer selection always precedes masking or randstrobe linking.

use crate::{canonical_encoded, encode_base, FeatureConfig, IndexMode, MinimizerScheme};
use anyhow::{ensure, Result};
use serde::{Deserialize, Serialize};
use std::{collections::VecDeque, str::FromStr};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SpacedSeed {
    width: u8,
    mask: u64,
}

impl FromStr for SpacedSeed {
    type Err = anyhow::Error;
    fn from_str(pattern: &str) -> Result<Self> {
        ensure!(
            !pattern.is_empty() && pattern.len() <= 31,
            "seed pattern width must be in 1..=31"
        );
        ensure!(
            pattern.bytes().all(|b| b == b'#' || b == b'_'),
            "seed pattern accepts only '#' and '_'"
        );
        ensure!(
            pattern.contains('#'),
            "seed pattern must contain at least one '#'"
        );
        ensure!(
            pattern.bytes().eq(pattern.bytes().rev()),
            "seed pattern must be symmetric under reversal for reverse-complement indexing"
        );
        let mut mask = 0;
        for b in pattern.bytes() {
            mask = (mask << 2) | if b == b'#' { 3 } else { 0 };
        }
        Ok(Self {
            width: pattern.len() as u8,
            mask,
        })
    }
}

impl std::fmt::Display for SpacedSeed {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in (0..self.width()).rev() {
            f.write_str(if self.mask & (3 << (2 * i)) != 0 {
                "#"
            } else {
                "_"
            })?;
        }
        Ok(())
    }
}

impl SpacedSeed {
    pub fn width(self) -> usize {
        self.width as usize
    }
    pub fn weight(self) -> usize {
        self.mask.count_ones() as usize / 2
    }

    /// Extract retained bases in their original order. BMI2 is optional.
    pub fn project(self, encoded: u64) -> u64 {
        #[cfg(target_arch = "x86_64")]
        if std::is_x86_feature_detected!("bmi2") {
            // SAFETY: checked CPU support, no memory accesses.
            return unsafe { std::arch::x86_64::_pext_u64(encoded, self.mask) };
        }
        self.project_portable(encoded)
    }

    fn project_portable(self, encoded: u64) -> u64 {
        let mut result = 0;
        for i in (0..self.width()).rev() {
            if self.mask & (3 << (2 * i)) != 0 {
                result = (result << 2) | ((encoded >> (2 * i)) & 3);
            }
        }
        result
    }

    pub fn canonical(self, encoded: u64) -> u64 {
        canonical_encoded(self.project(encoded), self.weight())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct RandstrobeConfig {
    /// Inclusive distance between the two strobe starts, in sequence bases.
    pub min_distance: usize,
    pub max_distance: usize,
}

impl RandstrobeConfig {
    pub fn validate(self) -> Result<()> {
        ensure!(
            self.min_distance > 0,
            "randstrobemer minimum distance must be greater than zero"
        );
        ensure!(
            self.min_distance <= self.max_distance,
            "randstrobemer minimum distance must not exceed maximum distance"
        );
        Ok(())
    }
}

pub fn print_settings(config: FeatureConfig) {
    println!(
        "seed_pattern\t{}",
        config
            .seed_pattern
            .map_or_else(|| "none".into(), |p| p.to_string())
    );
    println!("randstrobemers\t{}", config.randstrobemers.is_some());
    if let Some(strobes) = config.randstrobemers {
        println!(
            "strobe_min\t{}\nstrobe_max\t{}",
            strobes.min_distance, strobes.max_distance
        );
    }
}

#[derive(Clone, Copy)]
struct Strobe {
    position: usize,
    hash: u64,
}

/// XXH64 of one little-endian u64, seed zero, as used for canonical
/// strobes by strobealign. Standard XXH64 short-input rounds and avalanche.
fn strobe_hash(value: u64) -> u64 {
    let mut h = 0x27D4_EB2F_1656_67C5u64.wrapping_add(8)
        ^ value
            .wrapping_mul(0xC2B2_AE3D_27D4_EB4F)
            .rotate_left(31)
            .wrapping_mul(0x9E37_79B1_85EB_CA87);
    h = h
        .rotate_left(27)
        .wrapping_mul(0x9E37_79B1_85EB_CA87)
        .wrapping_add(0x85EB_CA77_C2B2_AE63);
    h ^= h >> 33;
    h = h.wrapping_mul(0xC2B2_AE3D_27D4_EB4F);
    h ^= h >> 29;
    h = h.wrapping_mul(0x1656_67B1_9E37_79F9);
    h ^ (h >> 32)
}

/// Versioned extraction identity. Selection follows strobealign's canonical
/// strobe XXH64 / masked-XOR popcount rule (q=255, first minimum wins).
/// Full-seed representation is the symmetric h1/2+h2/2 from Sahlin (2022),
/// not strobealign's newer partial-match/MCS index-specific bit layout.
pub const RANDSTROBE_SCHEME: &str = "strobealign-popcount8-symmetric64-downstream-v1";

fn reverse_encoded(mut key: u64, length: usize) -> u64 {
    let mut reverse = 0;
    for _ in 0..length {
        reverse = (reverse << 2) | ((key & 3) ^ 3);
        key >>= 2;
    }
    reverse
}

fn strobe(position: usize, forward: u64, length: usize, pattern: Option<SpacedSeed>) -> Strobe {
    let (forward, length) = pattern.map_or((forward, length), |p| (p.project(forward), p.weight()));
    let reverse = reverse_encoded(forward, length);
    Strobe {
        position,
        hash: strobe_hash(forward.min(reverse)),
    }
}

fn pair_key(a: Strobe, b: Strobe, _seed: u64) -> u64 {
    (a.hash / 2) + (b.hash / 2)
}

/// Buffer only the local candidate context, not all k-mers in a DNA run.
struct Linker<'a, F> {
    config: RandstrobeConfig,
    seed: u64,
    candidates: VecDeque<Strobe>,
    pending: usize,
    emit: &'a mut F,
}

impl<'a, F: FnMut(u64)> Linker<'a, F> {
    fn new(config: RandstrobeConfig, seed: u64, emit: &'a mut F) -> Self {
        Self {
            config,
            seed,
            candidates: VecDeque::new(),
            pending: 0,
            emit,
        }
    }
    fn push(&mut self, candidate: Strobe) {
        self.candidates.push_back(candidate);
        while self.pending < self.candidates.len()
            && candidate.position - self.candidates[self.pending].position
                > self.config.max_distance
        {
            self.finish_one();
        }
    }
    fn finish_one(&mut self) {
        let anchor = self.candidates[self.pending];
        // Exactly one reference seed per anchor. As in strobealign's iterator,
        // use the anchor itself when no downstream partner is available.
        // Reverse-complement partner selection belongs to QUERY extraction.
        let partner = self
            .best_partner(anchor, self.pending + 1..self.candidates.len())
            .unwrap_or(anchor);
        (self.emit)(pair_key(anchor, partner, self.seed));
        self.pending += 1;
        if self.pending < self.candidates.len() {
            let next = self.candidates[self.pending].position;
            while self.pending > 0 && next - self.candidates[0].position > self.config.max_distance
            {
                self.candidates.pop_front();
                self.pending -= 1;
            }
        }
    }
    fn best_partner(&self, anchor: Strobe, indices: impl Iterator<Item = usize>) -> Option<Strobe> {
        let mut best = None;
        for i in indices {
            let partner = self.candidates[i];
            let distance = anchor.position.abs_diff(partner.position);
            if distance > self.config.max_distance {
                break;
            }
            if distance < self.config.min_distance {
                continue;
            }
            let rank = ((anchor.hash ^ partner.hash) & 255).count_ones();
            if best.is_none_or(|(r, _)| rank < r) {
                best = Some((rank, partner));
                // Zero is optimal; stopping here preserves first-minimum ties.
                if rank == 0 {
                    break;
                }
            }
        }
        best.map(|(_, partner)| partner)
    }
    fn finish(&mut self) {
        while self.pending < self.candidates.len() {
            self.finish_one();
        }
    }
}

/// Preserve ambiguous separators, including non-ACGT bytes, under reversal.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| encode_base(b).map_or(b, |c| b"TGCA"[c as usize]))
        .collect()
}

pub(crate) fn scan<F: FnMut(u64)>(
    seq: &[u8],
    k: usize,
    seed: u64,
    config: FeatureConfig,
    emit: &mut F,
) {
    let mut ordinary = config;
    ordinary.seed_pattern = None;
    ordinary.randstrobemers = None;
    // `scan_features_uncompressed` already performed preprocessing.
    ordinary.homopolymer_compression = false;
    if config.randstrobemers.is_none() {
        let pattern = config.seed_pattern.expect("anchor extraction requested");
        crate::scan_plain_features_uncompressed(seq, k, seed, ordinary, &mut |key| {
            emit(pattern.canonical(key))
        });
        return;
    }
    let length = if config.mode == IndexMode::Minimizers {
        config.minimizer_size
    } else {
        k
    };
    for run in seq.split(|&b| encode_base(b).is_none()) {
        let mut linker = Linker::new(config.randstrobemers.unwrap(), seed, emit);
        if config.mode == IndexMode::Minimizers {
            let mut positions = Vec::new();
            match config.minimizer_scheme {
                MinimizerScheme::Simd => {
                    crate::scan_simd_minimizers_positions(run, k, length, seed, &mut |p, key| {
                        positions.push((p, key))
                    })
                }
                MinimizerScheme::OpenClosedMod => crate::open_closed_minimizers::scan_positions(
                    run,
                    k,
                    length,
                    seed,
                    &mut |p, key| positions.push((p, key)),
                ),
            }
            // Canonical tie handling may revisit positions; link each selected
            // physical position once, with distances independent of emission order.
            positions.sort_unstable_by_key(|&(p, _)| p);
            positions.dedup_by_key(|entry| entry.0);
            for (p, key) in positions {
                linker.push(if config.seed_pattern.is_some() {
                    strobe(p, key, length, config.seed_pattern)
                } else {
                    Strobe {
                        position: p,
                        hash: strobe_hash(key),
                    }
                });
            }
        } else if config.seed_pattern.is_none() {
            // The ordinary scanner rolls BOTH encodings in constant time,
            // rather than reversing each overlapping strobe from scratch.
            let mut position = 0;
            crate::scan_canonical_kmers(run, length, &mut |key| {
                linker.push(Strobe {
                    position,
                    hash: strobe_hash(key),
                });
                position += 1;
            });
        } else {
            let mask = (1u64 << (2 * length)) - 1;
            let mut forward = 0;
            for (end, &base) in run.iter().enumerate() {
                forward = ((forward << 2) | encode_base(base).unwrap() as u64) & mask;
                if end + 1 >= length {
                    linker.push(strobe(
                        end + 1 - length,
                        forward,
                        length,
                        config.seed_pattern,
                    ));
                }
            }
        }
        linker.finish();
    }
}

// Rank depends on only eight hash bits. For a dense candidate window, maintain
// its earliest member of each signature and intersect four-word bitsets with
// Hamming-distance classes, instead of rescanning every candidate per anchor.
fn signature_rank_masks() -> &'static [[[u64; 4]; 9]; 256] {
    static MASKS: std::sync::OnceLock<Box<[[[u64; 4]; 9]; 256]>> = std::sync::OnceLock::new();
    MASKS.get_or_init(|| {
        let mut masks = Box::new([[[0u64; 4]; 9]; 256]);
        for (query, row) in masks.iter_mut().enumerate() {
            for signature in 0..256 {
                row[(query as u8 ^ signature as u8).count_ones() as usize][signature / 64] |=
                    1u64 << (signature % 64);
            }
        }
        masks
    })
}

fn append_dense_pairs(
    hashes: &[u64],
    reverse: bool,
    config: RandstrobeConfig,
    keys: &mut Vec<u64>,
) {
    let n = hashes.len();
    let hash = |i: usize| hashes[if reverse { n - 1 - i } else { i }];
    if n <= config.min_distance {
        keys.extend((0..n).map(|i| {
            let h = hash(i) / 2;
            h + h
        }));
        return;
    }
    if n < 64 || config.max_distance - config.min_distance < 8 {
        let mut emit = |key| keys.push(key);
        let mut linker = Linker::new(config, 0, &mut emit);
        for i in 0..n {
            linker.push(Strobe {
                position: i,
                hash: hash(i),
            });
        }
        linker.finish();
        return;
    }
    let masks = signature_rank_masks();
    let mut heads = [usize::MAX; 256];
    let mut next = vec![usize::MAX; n];
    for i in (0..n).rev() {
        let signature = (hash(i) & 255) as usize;
        next[i] = heads[signature];
        heads[signature] = i;
    }
    heads.fill(usize::MAX);
    let mut present = [0u64; 4];
    for j in config.min_distance..n.min(config.max_distance.saturating_add(1)) {
        let signature = (hash(j) & 255) as usize;
        if heads[signature] == usize::MAX {
            heads[signature] = j;
            present[signature / 64] |= 1u64 << (signature % 64);
        }
    }
    for i in 0..n {
        let anchor = hash(i);
        let signature = (anchor & 255) as usize;
        let mut partner = heads[signature];
        if partner == usize::MAX && present.iter().any(|&word| word != 0) {
            for rank in 1..=8 {
                let hits: [u64; 4] =
                    std::array::from_fn(|word| masks[signature][rank][word] & present[word]);
                if hits.iter().any(|&word| word != 0) {
                    for (word, mut hits) in hits.into_iter().enumerate() {
                        while hits != 0 {
                            let signature = word * 64 + hits.trailing_zeros() as usize;
                            partner = partner.min(heads[signature]);
                            hits &= hits - 1;
                        }
                    }
                    break;
                }
            }
        }
        let partner = if partner == usize::MAX {
            anchor
        } else {
            hash(partner)
        };
        keys.push(anchor / 2 + partner / 2);
        // Expire the left endpoint, then insert the next right endpoint.
        if let Some(expired) = i.checked_add(config.min_distance).filter(|&j| j < n) {
            let signature = (hash(expired) & 255) as usize;
            debug_assert_eq!(heads[signature], expired);
            let upper = i.saturating_add(config.max_distance);
            heads[signature] = if next[expired] <= upper {
                next[expired]
            } else {
                usize::MAX
            };
            if heads[signature] == usize::MAX {
                present[signature / 64] &= !(1u64 << (signature % 64));
            }
        }
        if let Some(added) = i
            .checked_add(config.max_distance)
            .and_then(|j| j.checked_add(1))
            .filter(|&j| j < n)
        {
            let signature = (hash(added) & 255) as usize;
            if heads[signature] == usize::MAX {
                heads[signature] = added;
                present[signature / 64] |= 1u64 << (signature % 64);
            }
        }
    }
}

/// Bounded-record query helper. Reuse canonical strobe hashes across the two
/// independently selected orientations; reference construction remains one-way.
/// The caller bounds the read length before appending the two feature lists.
#[doc(hidden)]
pub fn append_dense_randstrobe_orientations(
    seq: &[u8],
    k: usize,
    _seed: u64,
    config: FeatureConfig,
    keys: &mut Vec<u64>,
) -> usize {
    assert!(config.mode == IndexMode::Kmers && config.seed_pattern.is_none());
    let distances = config.randstrobemers.expect("randstrobe configuration");
    let compressed;
    let seq = if config.homopolymer_compression {
        compressed = crate::homopolymer_compress(seq);
        &compressed
    } else {
        seq
    };
    let start = keys.len();
    let mut anchors = Vec::new();
    let mut runs = Vec::new();
    for run in seq.split(|&b| encode_base(b).is_none()) {
        let first = anchors.len();
        crate::scan_canonical_kmers(run, k, &mut |key| {
            anchors.push(strobe_hash(key));
        });
        let end = anchors.len();
        runs.push(first..end);
        append_dense_pairs(&anchors[first..end], false, distances, keys);
    }
    let forward_count = keys.len() - start;
    for range in runs.into_iter().rev() {
        append_dense_pairs(&anchors[range], true, distances, keys);
    }
    assert_eq!(keys.len() - start, forward_count * 2);
    forward_count
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encode_full_kmer;

    fn collect(seq: &[u8], k: usize, config: FeatureConfig) -> Vec<u64> {
        let mut keys = Vec::new();
        crate::scan_index_features(seq, k, 69, config, &mut |key| keys.push(key));
        keys
    }
    fn rc(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|&b| encode_base(b).map_or(b, |c| b"TGCA"[c as usize]))
            .collect()
    }
    fn dna(n: usize) -> Vec<u8> {
        (0..n)
            .map(|i| b"ACGT"[(crate::splitmix64(i as u64) & 3) as usize])
            .collect()
    }

    // Previous extraction path, retained independently of the optimized scan.
    fn legacy_strobes(seq: &[u8], k: usize, config: FeatureConfig) -> Vec<u64> {
        let length = if config.mode == IndexMode::Minimizers {
            config.minimizer_size
        } else {
            k
        };
        let mut result = Vec::new();
        let compressed;
        let seq = if config.homopolymer_compression {
            compressed = crate::homopolymer_compress(seq);
            &compressed
        } else {
            seq
        };
        for run in seq.split(|&b| encode_base(b).is_none()) {
            let mut emit = |key| result.push(key);
            let mut linker = Linker::new(config.randstrobemers.unwrap(), 69, &mut emit);
            if config.mode == IndexMode::Minimizers {
                let mut positions = Vec::new();
                match config.minimizer_scheme {
                    MinimizerScheme::Simd => {
                        crate::scan_simd_minimizers_positions(run, k, length, 69, &mut |p, _| {
                            positions.push(p)
                        })
                    }
                    MinimizerScheme::OpenClosedMod => {
                        crate::open_closed_minimizers::scan_positions(
                            run,
                            k,
                            length,
                            69,
                            &mut |p, _| positions.push(p),
                        )
                    }
                }
                positions.sort_unstable();
                positions.dedup();
                for p in positions {
                    linker.push(strobe(
                        p,
                        encode_full_kmer(&run[p..p + length]).unwrap(),
                        length,
                        config.seed_pattern,
                    ));
                }
            } else {
                let mask = (1u64 << (2 * length)) - 1;
                let mut forward = 0;
                for (end, &base) in run.iter().enumerate() {
                    forward = ((forward << 2) | encode_base(base).unwrap() as u64) & mask;
                    if end + 1 >= length {
                        linker.push(strobe(
                            end + 1 - length,
                            forward,
                            length,
                            config.seed_pattern,
                        ));
                    }
                }
            }
            linker.finish();
        }
        result
    }

    #[test]
    fn optimized_strobe_extraction_matches_previous_keys_in_order() {
        let random = dna(1500);
        let mut broken = random.clone();
        for p in (13..broken.len()).step_by(89) {
            broken[p] = b'N';
        }
        for seq in [random, broken, b"ACGT".repeat(375), b"A".repeat(1500)] {
            for seq in [&seq, &rc(&seq)] {
                for length in [1, 3, 7, 11, 21, 31] {
                    let mask: String = (0..length)
                        .map(|i| {
                            if i.min(length - 1 - i) % 3 == 1 {
                                '_'
                            } else {
                                '#'
                            }
                        })
                        .collect();
                    for mode in [IndexMode::Kmers, IndexMode::Minimizers] {
                        for scheme in [MinimizerScheme::Simd, MinimizerScheme::OpenClosedMod] {
                            for pattern in [None, Some(mask.parse().unwrap())] {
                                for homopolymer in [false, true] {
                                    let config = FeatureConfig::new(mode, length, 16)
                                        .with_minimizer_scheme(scheme)
                                        .with_seed_pattern(pattern)
                                        .with_homopolymer_compression(homopolymer)
                                        .with_randstrobemers(Some(RandstrobeConfig {
                                            min_distance: 15,
                                            max_distance: 50,
                                        }));
                                    let k = if mode == IndexMode::Kmers {
                                        length
                                    } else {
                                        length + 18
                                    };
                                    assert_eq!(collect(seq,k,config),legacy_strobes(seq,k,config),"{mode:?} {scheme:?} length={length} pattern={pattern:?} hpc={homopolymer}");
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn fused_dense_orientations_equal_independent_extraction_in_order() {
        let mut broken = dna(2000);
        for i in (3..broken.len()).step_by(37) {
            broken[i] = b'N';
        }
        for seq in [
            Vec::new(),
            dna(1),
            dna(10000),
            broken,
            b"A".repeat(3000),
            b"ACGT".repeat(700),
        ] {
            for k in [1, 3, 11, 21, 31] {
                for homopolymer in [false, true] {
                    for (min_distance, max_distance) in [(1, 1), (15, 50), (10000, 10000)] {
                        let config = FeatureConfig::new(IndexMode::Kmers, 11, 16)
                            .with_homopolymer_compression(homopolymer)
                            .with_randstrobemers(Some(RandstrobeConfig {
                                min_distance,
                                max_distance,
                            }));
                        let forward = collect(&seq, k, config);
                        let reverse = collect(&rc(&seq), k, config);
                        let mut actual = vec![123];
                        let n =
                            append_dense_randstrobe_orientations(&seq, k, 69, config, &mut actual);
                        assert_eq!(n, forward.len());
                        assert_eq!(actual[0], 123);
                        assert_eq!(&actual[1..1 + n], &forward);
                        assert_eq!(&actual[1 + n..], &reverse);
                    }
                }
            }
        }
    }

    #[test]
    fn signature_rank_masks_partition_every_possible_pair_exactly() {
        let masks = signature_rank_masks();
        for query in 0..256 {
            for signature in 0..256 {
                let expected = (query as u8 ^ signature as u8).count_ones() as usize;
                for rank in 0..=8 {
                    let present =
                        masks[query][rank][signature / 64] & (1u64 << (signature % 64)) != 0;
                    assert_eq!(present, rank == expected);
                }
            }
        }
    }

    #[test]
    fn signature_windows_match_scalar_partner_search_and_first_ties() {
        for n in [0, 1, 10, 63, 64, 65, 100, 1024] {
            for (min_distance, max_distance) in [
                (1, 1),
                (1, 9),
                (1, 50),
                (15, 50),
                (50, 50),
                (100, 200),
                (usize::MAX, usize::MAX),
            ] {
                let config = RandstrobeConfig {
                    min_distance,
                    max_distance,
                };
                for hashes in [
                    (0..n)
                        .map(|i| crate::splitmix64(i as u64))
                        .collect::<Vec<_>>(),
                    (0..n).map(|i| (i as u64) << 8).collect(),
                    (0..n)
                        .map(|i| ((i % 3) as u64) | ((i as u64) << 8))
                        .collect(),
                ] {
                    for reverse in [false, true] {
                        let input: Vec<_> = if reverse {
                            hashes.iter().rev().copied().collect()
                        } else {
                            hashes.clone()
                        };
                        let mut expected = Vec::new();
                        for i in 0..n {
                            let mut best = None;
                            for j in i + 1..n {
                                if j - i < min_distance || j - i > max_distance {
                                    continue;
                                }
                                let rank = ((input[i] ^ input[j]) & 255).count_ones();
                                if best.is_none_or(|(r, _)| rank < r) {
                                    best = Some((rank, j));
                                }
                            }
                            let partner = best.map_or(input[i], |(_, j)| input[j]);
                            expected.push(input[i] / 2 + partner / 2);
                        }
                        let mut actual = Vec::new();
                        append_dense_pairs(&hashes, reverse, config, &mut actual);
                        assert_eq!(
                            actual, expected,
                            "n={n} reverse={reverse} config={config:?}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn pattern_validation_and_projection() {
        for invalid in ["", "___", "#x#", "##_", "###############################_"] {
            assert!(invalid.parse::<SpacedSeed>().is_err(), "{invalid}");
        }
        let p: SpacedSeed = "##_##_##".parse().unwrap();
        assert_eq!(p.width(), 8);
        assert_eq!(p.weight(), 6);
        assert_eq!(p.to_string(), "##_##_##");
        assert_eq!(
            p.project(encode_full_kmer(b"ACGTAACG").unwrap()),
            encode_full_kmer(b"ACTACG").unwrap()
        );
        let a = b"ACATACCG";
        let b = b"ACGTAGCG";
        assert_eq!(
            p.canonical(encode_full_kmer(a).unwrap()),
            p.canonical(encode_full_kmer(b).unwrap())
        );
        assert_eq!(
            p.canonical(encode_full_kmer(a).unwrap()),
            p.canonical(encode_full_kmer(&rc(a)).unwrap())
        );
    }

    #[test]
    fn bit_extraction_matches_portable_projection_for_all_widths() {
        for width in 1..=31 {
            let mut pattern = vec![b'_'; width];
            for i in 0..width {
                if i == 0 || crate::splitmix64(i.min(width - 1 - i) as u64) & 1 != 0 {
                    pattern[i] = b'#';
                    pattern[width - 1 - i] = b'#';
                }
            }
            let p: SpacedSeed = std::str::from_utf8(&pattern).unwrap().parse().unwrap();
            for i in 0..100 {
                let key = crate::splitmix64(i);
                assert_eq!(p.project(key), p.project_portable(key));
            }
        }
    }

    #[test]
    fn randstrobes_link_selected_minimizers_not_unselected_mmers() {
        let seq = dna(31);
        for scheme in [MinimizerScheme::Simd, MinimizerScheme::OpenClosedMod] {
            let c = FeatureConfig::new(IndexMode::Minimizers, 9, 16).with_minimizer_scheme(scheme);
            assert_eq!(collect(&seq, 31, c).len(), 1);
            assert_eq!(
                collect(
                    &seq,
                    31,
                    c.with_randstrobemers(Some(RandstrobeConfig {
                        min_distance: 1,
                        max_distance: 30,
                    }))
                )
                .len(),
                1
            );
        }
        let c = FeatureConfig::new(IndexMode::Kmers, 9, 16).with_randstrobemers(Some(
            RandstrobeConfig {
                min_distance: 1,
                max_distance: 30,
            },
        ));
        assert!(!collect(&seq, 9, c).is_empty());
    }

    #[test]
    fn pair_keys_use_canonical_strobes_and_symmetric_paper_hash() {
        assert_eq!(strobe_hash(0x3332316f6c6c6568), 0x2d119fced5ee39ce);
        let a = strobe(0, encode_full_kmer(b"AC").unwrap(), 2, None);
        let b = strobe(5, encode_full_kmer(b"AG").unwrap(), 2, None);
        let reversed_a = strobe(0, encode_full_kmer(b"CT").unwrap(), 2, None);
        let reversed_b = strobe(5, encode_full_kmer(b"GT").unwrap(), 2, None);
        assert_eq!(pair_key(a, b, 69), pair_key(reversed_a, reversed_b, 69));
        let opposite = strobe(5, encode_full_kmer(b"CT").unwrap(), 2, None);
        assert_eq!(pair_key(a, b, 69), pair_key(a, opposite, 69));
        assert_eq!(pair_key(a, b, 69), pair_key(b, a, 69));
        assert_eq!(
            pair_key(a, b, 69),
            pair_key(a, Strobe { position: 100, ..b }, 69)
        );
    }

    #[test]
    fn all_retained_patterns_leave_default_keys_unchanged() {
        let seq = dna(4000);
        for mode in [IndexMode::Kmers, IndexMode::Minimizers] {
            for scheme in [MinimizerScheme::Simd, MinimizerScheme::OpenClosedMod] {
                let c = FeatureConfig::new(mode, 9, 16).with_minimizer_scheme(scheme);
                let pattern: SpacedSeed = "#########".parse().unwrap();
                assert_eq!(
                    collect(&seq, 9, c),
                    collect(&seq, 9, c.with_seed_pattern(Some(pattern)))
                );
            }
        }
    }

    #[test]
    fn minimizers_are_selected_before_masking() {
        for scheme in [MinimizerScheme::Simd, MinimizerScheme::OpenClosedMod] {
            let seq = dna(5000);
            let c = FeatureConfig::new(IndexMode::Minimizers, 9, 16).with_minimizer_scheme(scheme);
            let p: SpacedSeed = "##__#__##".parse().unwrap();
            let expected: Vec<_> = collect(&seq, 31, c)
                .into_iter()
                .map(|key| p.canonical(key))
                .collect();
            assert_eq!(expected, collect(&seq, 31, c.with_seed_pattern(Some(p))));
        }
    }

    #[test]
    fn reverse_complements_preserve_anchor_multisets() {
        for seq in [
            dna(3000),
            b"ACGTACGTACGTACGTACGTACGT".repeat(140),
            b"A".repeat(3000),
        ] {
            for mode in [IndexMode::Kmers, IndexMode::Minimizers] {
                for scheme in [MinimizerScheme::Simd, MinimizerScheme::OpenClosedMod] {
                    for pattern in [None, Some("##__#__##".parse().unwrap())] {
                        for strobes in [
                            None,
                            Some(RandstrobeConfig {
                                min_distance: 9,
                                max_distance: 35,
                            }),
                        ] {
                            let c = FeatureConfig::new(mode, 9, 16)
                                .with_minimizer_scheme(scheme)
                                .with_seed_pattern(pattern)
                                .with_randstrobemers(strobes);
                            let k = if mode == IndexMode::Kmers { 9 } else { 31 };
                            let mut a = collect(&seq, k, c);
                            let mut b = collect(&rc(&seq), k, c);
                            assert_eq!(a.len(), b.len());
                            if strobes.is_some() {
                                // Reference extraction is directional. Only the
                                // TWO query orientations form an invariant set.
                                a.extend(collect(&rc(&seq), k, c));
                                b.extend(collect(&seq, k, c));
                            }
                            a.sort_unstable();
                            b.sort_unstable();
                            assert_eq!(a, b, "{mode:?} {scheme:?} {pattern:?} {strobes:?}");
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn streaming_linker_matches_full_naive_context() {
        let seq = dna(2000);
        for distance in [(1, 1), (4, 13), (9, 100), (10000, 10000)] {
            let config = RandstrobeConfig {
                min_distance: distance.0,
                max_distance: distance.1,
            };
            for pattern in [None, Some("#_#_#_#".parse().unwrap())] {
                // Irregular spacing also exercises empty gaps between minimizers.
                let candidates: Vec<_> = (0..seq.len() - 7)
                    .filter(|&p| p % 71 < 10)
                    .map(|p| strobe(p, encode_full_kmer(&seq[p..p + 7]).unwrap(), 7, pattern))
                    .collect();
                let mut expected = Vec::new();
                for &anchor in &candidates {
                    let partner = candidates
                        .iter()
                        .copied()
                        .filter(|p| {
                            let distance = p.position.abs_diff(anchor.position);
                            p.position > anchor.position
                                && distance >= config.min_distance
                                && distance <= config.max_distance
                        })
                        .min_by_key(|p| ((anchor.hash ^ p.hash) & 255).count_ones())
                        .unwrap_or(anchor);
                    expected.push(pair_key(anchor, partner, 69));
                }
                let mut actual = Vec::new();
                let mut emit = |key| actual.push(key);
                let mut linker = Linker::new(config, 69, &mut emit);
                for candidate in candidates {
                    linker.push(candidate);
                }
                linker.finish();
                assert_eq!(actual, expected);
            }
        }
    }

    #[test]
    fn selector_matches_executed_upstream_strobealign_iterator_witnesses() {
        // Generated by check_strobealign_reference.py using the UNMODIFIED
        // upstream RandstrobeIterator at commit
        // 5dad538299f1479b00c6544acc3265ec6bc07555 and system libxxhash.
        // The aggregate covers every emitted key in order, including fallback.
        let seq = b"ACGTTGCACTGATCGATGGCATACCTAGGCTAACGTCGATCTTACGGAATCGCTAGCTGACCATGGTTACG";
        for (k, low, high, step, count, expected) in [
            (11, 4, 12, 1, 61, 11859818861789259928u64),
            (11, 2, 5, 3, 21, 7111390358099897965u64),
            (21, 7, 15, 1, 51, 7250723734144174666u64),
        ] {
            let mut keys = Vec::new();
            let mut emit = |key| keys.push(key);
            let mut linker = Linker::new(
                RandstrobeConfig {
                    min_distance: low * step,
                    max_distance: high * step,
                },
                69,
                &mut emit,
            );
            for p in (0..=seq.len() - k).step_by(step) {
                linker.push(strobe(
                    p,
                    encode_full_kmer(&seq[p..p + k]).unwrap(),
                    k,
                    None,
                ));
            }
            linker.finish();
            assert_eq!(keys.len(), count);
            let checksum = keys.iter().fold(0u64, |sum, &key| {
                sum.wrapping_mul(0x9e3779b185ebca87).wrapping_add(key)
            });
            assert_eq!(checksum, expected);
        }
    }

    #[test]
    fn reference_extraction_is_directional_and_never_doubles_anchor_density() {
        let seq = dna(1000);
        for mode in [IndexMode::Kmers, IndexMode::Minimizers] {
            let base = FeatureConfig::new(mode, 11, 16)
                .with_minimizer_scheme(MinimizerScheme::OpenClosedMod);
            let k = if mode == IndexMode::Kmers { 11 } else { 29 };
            let config = base.with_randstrobemers(Some(RandstrobeConfig {
                min_distance: 15,
                max_distance: 50,
            }));
            let mut forward = collect(&seq, k, config);
            let mut reverse = collect(&rc(&seq), k, config);
            assert_eq!(forward.len(), collect(&seq, k, base).len());
            assert_eq!(forward.len(), reverse.len());
            forward.sort_unstable();
            reverse.sort_unstable();
            assert_ne!(
                forward, reverse,
                "partner selection is not canonicalization"
            );
        }
    }

    #[test]
    fn masked_and_randstrobe_features_respect_breaks_and_compression() {
        let seq = dna(300);
        let mut input = seq.clone();
        input.extend_from_slice(b"NNN");
        input.extend_from_slice(&seq);
        for mode in [IndexMode::Kmers, IndexMode::Minimizers] {
            let c = FeatureConfig::new(mode, 9, 16)
                .with_minimizer_scheme(MinimizerScheme::OpenClosedMod)
                .with_seed_pattern(Some("##__#__##".parse().unwrap()))
                .with_randstrobemers(Some(RandstrobeConfig {
                    min_distance: 9,
                    max_distance: 40,
                }));
            let k = if mode == IndexMode::Kmers { 9 } else { 31 };
            assert_eq!(collect(&input, k, c), collect(&seq, k, c).repeat(2));
            let expanded: Vec<_> = input.iter().flat_map(|&b| [b, b, b]).collect();
            assert_eq!(
                collect(&expanded, k, c.with_homopolymer_compression(true)),
                collect(&crate::homopolymer_compress(&expanded), k, c)
            );
        }
    }

    #[test]
    fn distances_are_inclusive_and_unpaired_anchors_fall_back_to_self() {
        let c = FeatureConfig::new(IndexMode::Kmers, 3, 16).with_randstrobemers(Some(
            RandstrobeConfig {
                min_distance: 5,
                max_distance: 5,
            },
        ));
        assert!(collect(b"AC", 3, c).is_empty());
        assert_eq!(collect(b"ACGTACG", 3, c).len(), 5);
        assert_eq!(collect(b"ACGTACGT", 3, c).len(), 6);
        assert!(c
            .with_randstrobemers(Some(RandstrobeConfig {
                min_distance: 0,
                max_distance: 5
            }))
            .validate_anchors(3)
            .is_err());
        assert!(c
            .with_randstrobemers(Some(RandstrobeConfig {
                min_distance: 6,
                max_distance: 5
            }))
            .validate_anchors(3)
            .is_err());
    }
}
