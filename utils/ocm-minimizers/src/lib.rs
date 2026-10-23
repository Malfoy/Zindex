//! Shared, dependency-free canonical open-closed mod minimizer iterator.
//! `k` is the full window span; `m` is the selected representative length.
mod implementation;
pub use implementation::{scan, scan_positions};

#[inline(always)]
pub fn encode_base(base: u8) -> Option<u8> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

pub fn encode_full_kmer(seq: &[u8]) -> Option<u64> {
    let mut value = 0u64;
    for &base in seq {
        value = (value << 2) | encode_base(base)? as u64;
    }
    Some(value)
}

pub fn canonical_encoded(encoded: u64, k: usize) -> u64 {
    let mut value = encoded;
    let mut rc = 0u64;
    for _ in 0..k {
        let code = value & 0b11;
        rc = (rc << 2) | (code ^ 0b11);
        value >>= 2;
    }
    encoded.min(rc)
}

#[inline(always)]
pub fn mixsplit(key: u64, seed: u64) -> u64 {
    splitmix64(key.wrapping_add(seed))
}

#[inline(always)]
pub fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}
