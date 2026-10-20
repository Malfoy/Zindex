//! Canonical open-closed mod minimizers (Groot Koerkamp, Liu and Pibiri, 2025,
//! doi:10.1186/s13015-025-00270-0). `k` here is the *full window length*;
//! `m` is the representative length, called k in the paper.
//!
//! We use r=min(4,m), t=r+(m-r)%w, w=k-m+1. A t-mer anchor is preferred
//! if its smallest s-mer is in the middle (open), then at an end (closed),
//! otherwise it is a fallback. The representative starts at anchor_offset%w.
//!
//! Canonical adaptation: hash canonical t-/s-mers, choose s<=4 with t-s even,
//! and treat all occurrences of a tied inner minimum equally. Outer ties use
//! the orientation of the odd-length window (GT majority versus AC majority).
//! Since t == m (mod w), reversing an anchor maps the selected representative
//! to its reverse complement. This is window-local, including on repetitive
//! inputs, rather than merely canonicalizing a strand-dependent selection.

#[cfg(test)]
use std::collections::VecDeque;

use crate::{encode_base, mixsplit};

#[cfg(target_arch = "x86_64")]
mod avx2;

#[cfg(test)]
#[derive(Clone, Copy)]
struct Entry {
    position: usize,
    rank: (u8, u64),
}

#[cfg(test)]
struct Minimum {
    queue: VecDeque<Entry>,
    rightmost: bool,
}

#[cfg(test)]
impl Minimum {
    fn new(rightmost: bool) -> Self {
        Self {
            queue: VecDeque::new(),
            rightmost,
        }
    }

    fn push(&mut self, position: usize, rank: (u8, u64), start: usize) {
        while self.queue.front().is_some_and(|e| e.position < start) {
            self.queue.pop_front();
        }
        while self
            .queue
            .back()
            .is_some_and(|e| e.rank > rank || (self.rightmost && e.rank == rank))
        {
            self.queue.pop_back();
        }
        self.queue.push_back(Entry { position, rank });
    }

    fn get(&self) -> Entry {
        *self.queue.front().expect("nonempty sliding minimum")
    }
}

struct Rolling {
    forward: u64,
    reverse: u64,
    mask: u64,
    shift: usize,
}

impl Rolling {
    fn new(length: usize) -> Self {
        Self {
            forward: 0,
            reverse: 0,
            mask: (1u64 << (2 * length)) - 1,
            shift: 2 * (length - 1),
        }
    }
    fn push(&mut self, code: u8) -> u64 {
        self.forward = ((self.forward << 2) | code as u64) & self.mask;
        self.reverse = (self.reverse >> 2) | (((code ^ 3) as u64) << self.shift);
        self.forward.min(self.reverse)
    }
}

fn parameters(k: usize, m: usize) -> (usize, usize, usize) {
    let w = k - m + 1;
    let r = 4.min(m);
    let t = r + (m - r) % w;
    let mut s = 4.min(t);
    if (t - s) & 1 != 0 {
        s -= 1;
    }
    (w, t, s)
}

/// Extract in O(n) time with O(k) scratch plus fixed-size SIMD tile buffers,
/// splitting at ambiguous bases. AVX2 is selected automatically when supported.
/// Consecutive selections of the same position are emitted only once.
pub fn scan<F: FnMut(u64)>(seq: &[u8], k: usize, m: usize, seed: u64, emit: &mut F) {
    scan_positions(seq, k, m, seed, &mut |_, key| emit(key));
}

/// As `scan`, also reporting the selected representative's sequence offset.
pub fn scan_positions<F: FnMut(usize, u64)>(
    seq: &[u8],
    k: usize,
    m: usize,
    seed: u64,
    emit: &mut F,
) {
    if k == 0 || k & 1 == 0 || m == 0 || m > 31 || m > k {
        return;
    }
    let mut start = 0;
    while start < seq.len() {
        let end = start + first_invalid(&seq[start..]);
        if end - start >= k {
            scan_run(&seq[start..end], k, m, seed, &mut |position, key| {
                emit(start + position, key)
            });
        }
        start = end.saturating_add(1);
    }
}

fn first_invalid(seq: &[u8]) -> usize {
    #[cfg(target_arch = "x86_64")]
    if seq.len() >= 32 && std::is_x86_feature_detected!("avx2") {
        // SAFETY: runtime check. The vector scanner loads only complete blocks.
        return unsafe { avx2::first_invalid(seq) };
    }
    seq.iter()
        .position(|&base| encode_base(base).is_none())
        .unwrap_or(seq.len())
}

fn scan_run<F: FnMut(usize, u64)>(seq: &[u8], k: usize, m: usize, seed: u64, emit: &mut F) {
    if k == m {
        let mut rolling = Rolling::new(m);
        for (end, &base) in seq.iter().enumerate() {
            let key = rolling.push(encode_base(base).unwrap());
            if end + 1 >= m {
                emit(end + 1 - m, key);
            }
        }
        return;
    }
    #[cfg(target_arch = "x86_64")]
    if seq.len() - k + 1 >= k.max(48) && std::is_x86_feature_detected!("avx2") {
        // SAFETY: checked CPU support, and scan() supplies a valid DNA run.
        // Tests also call this with valid runs and supported parameters.
        return unsafe { avx2::scan_run(seq, k, m, seed, emit) };
    }
    scan_run_portable(seq, k, m, seed, emit);
}

/// Portable small-record path. The tiny s-mer domain lets us precompute its
/// exact 64-bit hash, without changing the ordering or introducing truncation.
fn scan_run_portable<F: FnMut(usize, u64)>(
    seq: &[u8],
    k: usize,
    m: usize,
    seed: u64,
    emit: &mut F,
) {
    let (w, t, s) = parameters(k, m);
    let inner_length = t - s + 1;
    let anchor_length = k - t + 1;
    let hashes = inner_hash_table(s, seed);
    let mut representative = Rolling::new(m);
    let t_mask = (1u64 << (2 * t)) - 1;
    let s_mask = (1u64 << (2 * s)) - 1;
    let ring_length = k.next_power_of_two();
    let mask = ring_length - 1;
    let mut keys = vec![0; ring_length];
    let mut inner_hashes = vec![0; ring_length];
    let mut minimum = ScalarWindow::new(inner_length);
    let mut anchors = ScalarAnchors::new(anchor_length);
    let mut gt_count = 0;
    let mut previous = None;
    for (end, &base) in seq.iter().enumerate() {
        let code = encode_base(base).unwrap();
        gt_count += usize::from(code >= 2);
        if end >= k {
            gt_count -= usize::from(encode_base(seq[end - k]).unwrap() >= 2);
        }
        let key = representative.push(code);
        if end + 1 >= m {
            keys[(end + 1 - m) & mask] = key;
        }
        if end + 1 < s {
            continue;
        }
        let h = hashes[(representative.forward & s_mask) as usize];
        let inner_pos = end + 1 - s;
        inner_hashes[inner_pos & mask] = h;
        let inner_min = minimum.push(h);
        if end + 1 < t {
            continue;
        }
        let pos = end + 1 - t;
        let class = if inner_hashes[(pos + (t - s) / 2) & mask] == inner_min {
            0
        } else if inner_hashes[pos & mask] == inner_min || h == inner_min {
            1
        } else {
            2
        };
        let anchor_key =
            (representative.forward & t_mask).min(representative.reverse >> (2 * (m - t)));
        let selected = anchors.push(ScalarRank {
            class,
            hash: mixsplit(anchor_key, seed),
            left: pos,
            right: pos,
        });
        if end + 1 < k {
            continue;
        }
        let start = end + 1 - k;
        let anchor = if gt_count > k / 2 {
            selected.right
        } else {
            selected.left
        };
        let pos = start + (anchor - start) % w;
        if previous != Some(pos) {
            emit(pos, keys[pos & mask]);
            previous = Some(pos);
        }
    }
}

fn inner_hash_table(s: usize, seed: u64) -> [u64; 256] {
    // Parser records usually share (s,seed). Copy the small cached table so
    // callbacks may safely invoke scan() again; no TLS borrow survives emit().
    thread_local! {
        static CACHE: std::cell::RefCell<(usize,u64,[u64;256])> = const {
            std::cell::RefCell::new((0,0,[0;256]))
        };
    }
    CACHE.with(|cache| {
        let mut cache = cache.borrow_mut();
        if cache.0 != s || cache.1 != seed {
            *cache = (s, seed, make_inner_hash_table(s, seed));
        }
        cache.2
    })
}

fn make_inner_hash_table(s: usize, seed: u64) -> [u64; 256] {
    let mut table = [0; 256];
    for (key, hash) in table[..1 << (2 * s)].iter_mut().enumerate() {
        *hash = mixsplit(
            crate::canonical_encoded(key as u64, s),
            seed ^ 0xC6BC_2796_92B5_CC83,
        );
    }
    table
}

struct ScalarWindow {
    ring: Vec<u64>,
    index: usize,
    prefix: u64,
}
impl ScalarWindow {
    fn new(length: usize) -> Self {
        Self {
            ring: vec![u64::MAX; length],
            index: 0,
            prefix: u64::MAX,
        }
    }
    fn push(&mut self, value: u64) -> u64 {
        self.ring[self.index] = value;
        self.prefix = self.prefix.min(value);
        self.index += 1;
        if self.index == self.ring.len() {
            self.index = 0;
            for i in (0..self.ring.len() - 1).rev() {
                self.ring[i] = self.ring[i].min(self.ring[i + 1]);
            }
            self.prefix = value;
        }
        self.prefix.min(self.ring[self.index])
    }
}

#[derive(Clone, Copy)]
struct ScalarRank {
    class: u8,
    hash: u64,
    left: usize,
    right: usize,
}
impl ScalarRank {
    const EMPTY: Self = Self {
        class: 3,
        hash: u64::MAX,
        left: 0,
        right: 0,
    };
    // a is temporally left of b: on equal ranks preserve both extremes.
    fn combine(a: Self, b: Self) -> Self {
        match (a.class, a.hash).cmp(&(b.class, b.hash)) {
            std::cmp::Ordering::Less => a,
            std::cmp::Ordering::Greater => b,
            std::cmp::Ordering::Equal => Self {
                right: b.right,
                ..a
            },
        }
    }
}
struct ScalarAnchors {
    ring: Vec<ScalarRank>,
    index: usize,
    prefix: ScalarRank,
}
impl ScalarAnchors {
    fn new(length: usize) -> Self {
        Self {
            ring: vec![ScalarRank::EMPTY; length],
            index: 0,
            prefix: ScalarRank::EMPTY,
        }
    }
    fn push(&mut self, value: ScalarRank) -> ScalarRank {
        self.ring[self.index] = value;
        self.prefix = ScalarRank::combine(self.prefix, value);
        self.index += 1;
        if self.index == self.ring.len() {
            self.index = 0;
            for i in (0..self.ring.len() - 1).rev() {
                self.ring[i] = ScalarRank::combine(self.ring[i], self.ring[i + 1]);
            }
            self.prefix = value;
        }
        ScalarRank::combine(self.ring[self.index], self.prefix)
    }
}

#[cfg(test)]
fn scan_run_scalar<F: FnMut(usize, u64)>(seq: &[u8], k: usize, m: usize, seed: u64, emit: &mut F) {
    let (w, t, s) = parameters(k, m);
    let inner_length = t - s + 1;
    let anchor_length = k - t + 1;
    let mut representative = Rolling::new(m);
    let mut anchor = Rolling::new(t);
    let mut inner = Rolling::new(s);
    let mut inner_min = Minimum::new(false);
    let mut left = Minimum::new(false);
    let mut right = Minimum::new(true);
    let mut inner_hashes = vec![0; inner_length];
    let mut representatives = vec![0; w];
    let mut gt_count = 0usize;
    let mut previous = None;
    for (end, &base) in seq.iter().enumerate() {
        let code = encode_base(base).unwrap();
        gt_count += usize::from(code >= 2);
        if end >= k {
            gt_count -= usize::from(encode_base(seq[end - k]).unwrap() >= 2);
        }
        let key = representative.push(code);
        let anchor_key = anchor.push(code);
        let inner_key = inner.push(code);
        if end + 1 >= m {
            representatives[(end + 1 - m) % w] = key;
        }
        if end + 1 < s {
            continue;
        }
        let inner_pos = end + 1 - s;
        let h = mixsplit(inner_key, seed ^ 0xC6BC_2796_92B5_CC83);
        inner_hashes[inner_pos % inner_length] = h;
        let inner_start = (end + 1).saturating_sub(t);
        inner_min.push(inner_pos, (0, h), inner_start);
        if end + 1 < t {
            continue;
        }
        let anchor_pos = end + 1 - t;
        let minimum = inner_min.get().rank.1;
        let middle = anchor_pos + (t - s) / 2;
        let class = if inner_hashes[middle % inner_length] == minimum {
            0
        } else if inner_hashes[anchor_pos % inner_length] == minimum || h == minimum {
            1
        } else {
            2
        };
        let start = (anchor_pos + 1).saturating_sub(anchor_length);
        let rank = (class, mixsplit(anchor_key, seed));
        left.push(anchor_pos, rank, start);
        right.push(anchor_pos, rank, start);
        if end + 1 < k {
            continue;
        }
        let chosen = if gt_count > k / 2 {
            right.get()
        } else {
            left.get()
        };
        let pos = start + (chosen.position - start) % w;
        if previous != Some(pos) {
            emit(pos, representatives[pos % w]);
            previous = Some(pos);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{canonical_encoded, encode_full_kmer, splitmix64};
    use std::collections::BTreeSet;

    fn reverse_complement(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|&b| b"TGCA"[encode_base(b).unwrap() as usize])
            .collect()
    }

    // Deliberately slow independent reference: evaluate every candidate in each
    // window, without rolling encodings, rings, or sliding queues.
    fn reference(seq: &[u8], k: usize, m: usize, seed: u64) -> Vec<(usize, u64)> {
        let (w, t, s) = parameters(k, m);
        let canonical =
            |slice: &[u8]| canonical_encoded(encode_full_kmer(slice).unwrap(), slice.len());
        let mut result = Vec::new();
        for (start, window) in seq.windows(k).enumerate() {
            let mut candidates = Vec::new();
            for (pos, anchor) in window.windows(t).enumerate() {
                let hashes: Vec<_> = anchor
                    .windows(s)
                    .map(|inner| mixsplit(canonical(inner), seed ^ 0xC6BC_2796_92B5_CC83))
                    .collect();
                let minimum = *hashes.iter().min().unwrap();
                let class = if hashes[(t - s) / 2] == minimum {
                    0
                } else if hashes[0] == minimum || *hashes.last().unwrap() == minimum {
                    1
                } else {
                    2
                };
                candidates.push(((class, mixsplit(canonical(anchor), seed)), pos));
            }
            let rank = candidates.iter().map(|&(rank, _)| rank).min().unwrap();
            let mut tied = candidates
                .iter()
                .filter(|&&(r, _)| r == rank)
                .map(|&(_, p)| p);
            let gt = window
                .iter()
                .filter(|&&b| encode_base(b).unwrap() >= 2)
                .count();
            let anchor = if gt > k / 2 {
                tied.next_back().unwrap()
            } else {
                tied.next().unwrap()
            };
            let pos = start + anchor % w;
            if result.last().map(|&(p, _)| p) != Some(pos) {
                result.push((pos, canonical(&seq[pos..pos + m])));
            }
        }
        result
    }

    #[test]
    fn rolling_matches_reference_and_reverse_complement_for_all_parameters() {
        for seq in [
            (0..190)
                .map(|i| b"ACGT"[(splitmix64(i) & 3) as usize])
                .collect::<Vec<_>>(),
            b"ACGT".repeat(48),
            b"A".repeat(190),
            b"AAAACCCCGGGGTTTT".repeat(12),
        ] {
            for k in (1..=65).step_by(2) {
                for m in 1..=31.min(k) {
                    for seed in [0, 69] {
                        let mut actual = Vec::new();
                        scan_run(&seq, k, m, seed, &mut |p, key| actual.push((p, key)));
                        assert_eq!(actual, reference(&seq, k, m, seed), "k={k},m={m}");
                        let mut portable = Vec::new();
                        scan_run_portable(&seq, k, m, seed, &mut |p, key| portable.push((p, key)));
                        assert_eq!(actual, portable, "portable k={k},m={m}");
                        let mut legacy = Vec::new();
                        scan_run_scalar(&seq, k, m, seed, &mut |p, key| legacy.push((p, key)));
                        assert_eq!(actual, legacy, "legacy k={k},m={m}");
                        let rc = reverse_complement(&seq);
                        let mut reversed = BTreeSet::new();
                        scan_run(&rc, k, m, seed, &mut |p, key| {
                            reversed.insert((seq.len() - m - p, key));
                        });
                        assert_eq!(
                            actual.iter().copied().collect::<BTreeSet<_>>(),
                            reversed,
                            "RC k={k},m={m}"
                        );
                        let indexed: BTreeSet<_> = actual.iter().map(|&(_, key)| key).collect();
                        for window in seq.windows(k) {
                            let mut query = Vec::new();
                            scan(window, k, m, seed, &mut |key| query.push(key));
                            assert_eq!(query.len(), 1);
                            assert!(indexed.contains(&query[0]));
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn ambiguities_split_runs_and_case_is_ignored() {
        let mut whole = Vec::new();
        let mut split = Vec::new();
        scan(b"acgtACGTAANNNacgtacgtaa", 9, 5, 69, &mut |key| {
            whole.push(key)
        });
        for _ in 0..2 {
            scan(b"ACGTACGTAA", 9, 5, 69, &mut |key| split.push(key));
        }
        assert_eq!(whole, split);
    }

    #[test]
    fn vector_tiles_and_tail_match_legacy_positions_and_keys_exactly() {
        for seq in [
            (0..9007)
                .map(|i| b"ACGT"[(splitmix64(i) & 3) as usize])
                .collect::<Vec<_>>(),
            b"aCgT".repeat(2253),
            b"A".repeat(9007),
            b"AAAACCCCGGGGTTTT".repeat(563),
        ] {
            for (k, m) in [
                (1, 1),
                (3, 2),
                (31, 31),
                (33, 31),
                (49, 31),
                (63, 21),
                (65, 30),
                (127, 31),
                (255, 31),
            ] {
                for seed in [0, 69, u64::MAX] {
                    let mut actual = Vec::new();
                    let mut legacy = Vec::new();
                    let mut portable = Vec::new();
                    scan_run(&seq, k, m, seed, &mut |p, key| actual.push((p, key)));
                    scan_run_scalar(&seq, k, m, seed, &mut |p, key| legacy.push((p, key)));
                    scan_run_portable(&seq, k, m, seed, &mut |p, key| portable.push((p, key)));
                    assert_eq!(actual, legacy, "k={k},m={m},seed={seed}");
                    assert_eq!(actual, portable, "portable k={k},m={m},seed={seed}");
                }
            }
        }
        let seq: Vec<_> = (0..4500)
            .map(|i| b"ACGT"[(splitmix64(i) & 3) as usize])
            .collect();
        // Around dispatch, four-lane, 1024-window tile and scalar-tail cuts.
        for windows in [
            1, 47, 48, 49, 50, 51, 64, 127, 128, 129, 131, 132, 4095, 4096, 4097, 4223, 4224, 4225,
        ] {
            let seq = &seq[..windows + 48];
            let mut actual = Vec::new();
            let mut legacy = Vec::new();
            scan_run(seq, 49, 31, 69, &mut |p, key| actual.push((p, key)));
            scan_run_scalar(seq, 49, 31, 69, &mut |p, key| legacy.push((p, key)));
            assert_eq!(actual, legacy, "windows={windows}");
        }
    }

    #[test]
    fn simd_validation_accepts_exactly_the_scalar_dna_alphabet() {
        let mut seq = b"aCgT".repeat(33);
        for offset in [0, 1, 30, 31, 32, 33, 63, 64, 127, 131] {
            for byte in 0..=255u8 {
                let old = seq[offset];
                seq[offset] = byte;
                let expected = seq
                    .iter()
                    .position(|&base| encode_base(base).is_none())
                    .unwrap_or(seq.len());
                assert_eq!(first_invalid(&seq), expected, "offset={offset},byte={byte}");
                seq[offset] = old;
            }
        }
    }

    #[test]
    fn cached_hash_tables_allow_reentrant_callbacks_and_multiple_seeds() {
        let seq = b"ACGTACGTAAAACCCCGGGGTTTT".repeat(50);
        let mut expected = Vec::new();
        scan_run_scalar(&seq, 49, 31, 69, &mut |_, key| expected.push(key));
        let mut actual = Vec::new();
        scan(&seq, 49, 31, 69, &mut |key| {
            actual.push(key);
            scan(
                b"ACGTACGTAAAACCCCGGGGTTTTACGTACGT",
                31,
                21,
                u64::MAX,
                &mut |_| {},
            );
        });
        assert_eq!(actual, expected);
    }

    #[test]
    #[ignore = "release-mode short-record measurement"]
    fn short_record_throughput() {
        let seq: Vec<_> = (0..5_000_000)
            .map(|i| b"ACGT"[(splitmix64(i) & 3) as usize])
            .collect();
        for length in [49, 75, 100, 150, 250, 1000] {
            for name in ["legacy", "portable", "dispatch", "avx2-small"] {
                let start = std::time::Instant::now();
                let mut count = 0;
                let mut checksum = 0;
                for record in seq.chunks_exact(length) {
                    let mut emit = |_: usize, key| {
                        count += 1;
                        checksum ^= key;
                    };
                    match name {
                        "legacy" => scan_run_scalar(record, 49, 31, 69, &mut emit),
                        "portable" => scan_run_portable(record, 49, 31, 69, &mut emit),
                        #[cfg(target_arch = "x86_64")]
                        "avx2-small" if std::is_x86_feature_detected!("avx2") => unsafe {
                            avx2::scan_run_at_threshold(record, 49, 31, 69, 4, &mut emit);
                        },
                        _ => scan_run(record, 49, 31, 69, &mut emit),
                    }
                }
                eprintln!("length={length} path={name} features={count} seconds={:.4} checksum={checksum}",start.elapsed().as_secs_f64());
            }
        }
    }

    #[test]
    #[ignore = "release-mode repeated extraction benchmark, optionally with real genomes"]
    fn repeated_extraction_benchmark() {
        use helicase::input::FromSlice;
        use helicase::{FastxParser, HelicaseParser};
        use std::io::Read;
        let mut records = Vec::new();
        let label = if let Ok(fof) = std::env::var("ZOR_OCM_GENOMES_FOF") {
            for input in crate::read_fof(std::path::Path::new(&fof))
                .unwrap()
                .into_iter()
                .take(32)
            {
                // Cached collections may use .fasta symlinks to .fasta.gz.
                // Resolve their actual encoding for the benchmark input reader.
                let path = std::fs::canonicalize(&input.path).unwrap();
                let mut data = Vec::new();
                if crate::is_compressed_path(&path) {
                    deko::read::AnyDecoder::new(std::fs::File::open(path).unwrap())
                        .read_to_end(&mut data)
                        .unwrap();
                } else {
                    data = std::fs::read(path).unwrap();
                }
                let mut parser =
                    FastxParser::<{ crate::PARSER_CONFIG }>::from_slice(&data).unwrap();
                while parser.next().is_some() {
                    records.push(parser.get_dna_string().to_vec());
                }
            }
            "real-32"
        } else {
            records.push(
                (0..5_000_000)
                    .map(|i| b"ACGT"[(splitmix64(i) & 3) as usize])
                    .collect(),
            );
            "synthetic"
        };
        let bases: usize = records.iter().map(Vec::len).sum();
        let paths = ["legacy", "portable", "optimized"];
        for (k, m) in [(49, 31), (63, 21)] {
            let mut times: [Vec<f64>; 3] = std::array::from_fn(|_| Vec::new());
            let mut expected = None;
            for repetition in 0..5 {
                for step in 0..3 {
                    let index = (step + repetition) % 3;
                    let start = std::time::Instant::now();
                    let mut count = 0u64;
                    let mut checksum = 0u64;
                    for record in records.iter().filter(|r| r.len() >= k) {
                        let mut emit = |pos: usize, key| {
                            count += 1;
                            checksum = checksum
                                .wrapping_mul(0x9E37_79B9)
                                .wrapping_add(key ^ pos as u64);
                        };
                        match paths[index] {
                            "legacy" => scan_run_scalar(record, k, m, 69, &mut emit),
                            "portable" => scan_run_portable(record, k, m, 69, &mut emit),
                            _ => scan_run(record, k, m, 69, &mut emit),
                        }
                    }
                    times[index].push(start.elapsed().as_secs_f64());
                    assert_eq!(
                        *expected.get_or_insert((count, checksum)),
                        (count, checksum)
                    );
                }
            }
            for (name, mut values) in paths.into_iter().zip(times) {
                values.sort_by(f64::total_cmp);
                eprintln!("input={label} bases={bases} k={k} m={m} path={name} median_s={:.6} min_s={:.6} features={} ordered_checksum={}",values[2],values[0],expected.unwrap().0,expected.unwrap().1);
            }
        }
    }

    #[test]
    #[ignore = "release-mode throughput measurement"]
    fn throughput() {
        let seq: Vec<_> = (0..5_000_000)
            .map(|i| b"ACGT"[(splitmix64(i) & 3) as usize])
            .collect();
        for (k, m) in [(49, 31), (63, 21)] {
            for name in ["legacy", "portable", "avx2"] {
                let start = std::time::Instant::now();
                let mut count = 0;
                let mut checksum = 0;
                let mut emit = |_: usize, key| {
                    count += 1;
                    checksum ^= key;
                };
                match name {
                    "legacy" => scan_run_scalar(&seq, k, m, 69, &mut emit),
                    "portable" => scan_run_portable(&seq, k, m, 69, &mut emit),
                    _ => scan_run(&seq, k, m, 69, &mut emit),
                }
                eprintln!(
                    "path={name} k={k} m={m} features={count} seconds={:.4} checksum={checksum}",
                    start.elapsed().as_secs_f64()
                );
            }
            for scheme in [
                crate::MinimizerScheme::Simd,
                crate::MinimizerScheme::OpenClosedMod,
            ] {
                let config = crate::FeatureConfig::new(crate::IndexMode::Minimizers, m, 16)
                    .with_minimizer_scheme(scheme);
                let start = std::time::Instant::now();
                let mut count = 0;
                let mut checksum = 0;
                crate::scan_index_features(&seq, k, 69, config, &mut |key| {
                    count += 1;
                    checksum ^= key;
                });
                eprintln!("scheme={} k={k} m={m} features={count} density={:.5} seconds={:.4} checksum={checksum}",scheme.label(),count as f64/seq.len() as f64,start.elapsed().as_secs_f64());
            }
        }
    }
}
