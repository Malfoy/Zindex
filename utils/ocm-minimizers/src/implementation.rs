// Canonical open-closed mod minimizers (Groot Koerkamp, Liu and Pibiri, 2025,
// doi:10.1186/s13015-025-00270-0). `k` here is the *full window length*;
// `m` is the representative length, called k in the paper.
//
// We use r=min(4,m), t=r+(m-r)%w, w=k-m+1. A t-mer anchor is preferred
// if its smallest s-mer is in the middle (open), then at an end (closed),
// otherwise it is a fallback. The representative starts at anchor_offset%w.
//
// Canonical adaptation: hash canonical t-/s-mers, choose s<=4 with t-s even,
// and treat all occurrences of a tied inner minimum equally. Outer ties use
// the orientation of the odd-length window (GT majority versus AC majority).
// Since t == m (mod w), reversing an anchor maps the selected representative
// to its reverse complement. This is window-local, including on repetitive
// inputs, rather than merely canonicalizing a strand-dependent selection.

#[cfg(test)]
use std::collections::VecDeque;

use crate::{encode_base, mixsplit};

#[cfg(target_arch = "x86_64")]
#[path = "avx2.rs"]
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
mod shared_reference_tests {
    use super::*;

    #[test]
    fn standalone_shared_core_matches_portable_and_legacy_kernels() {
        for seq in [
            (0..4097)
                .map(|i| b"ACGT"[(crate::splitmix64(i) & 3) as usize])
                .collect::<Vec<_>>(),
            b"ACGT".repeat(1025),
            vec![b'A'; 4097],
        ] {
            for (span, seed_length) in [(1, 1), (7, 7), (19, 12), (49, 31), (63, 21)] {
                for seed in [0, 69, u64::MAX] {
                    let mut actual = Vec::new();
                    scan_positions(&seq, span, seed_length, seed, &mut |p, k| {
                        actual.push((p, k))
                    });
                    let mut portable = Vec::new();
                    scan_run_portable(&seq, span, seed_length, seed, &mut |p, k| {
                        portable.push((p, k))
                    });
                    let mut legacy = Vec::new();
                    scan_run_scalar(&seq, span, seed_length, seed, &mut |p, k| {
                        legacy.push((p, k))
                    });
                    assert_eq!(actual, portable);
                    assert_eq!(actual, legacy);
                }
            }
        }
    }
}
