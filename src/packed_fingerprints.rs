//! Fused XOR/equality screens for tightly packed two- and four-bit fingerprints.
//! No unpacked cells or intermediate XOR row are allocated. Only matching field
//! indices are emitted; padding (including an adjacent origin plane) is ignored.

/// Screen corresponding packed fields across selected slices.
#[inline]
pub fn matches_rows(
    rows: &[&[u8]],
    bits: usize,
    fingerprint: u64,
    count: usize,
    mut hit: impl FnMut(usize),
) {
    debug_assert!(!rows.is_empty());
    debug_assert!(fingerprint < (1 << bits));
    debug_assert!(rows.iter().all(|r| r.len() == (count * bits).div_ceil(8)));
    match bits {
        2 => screen::<2>(rows, fingerprint, count, &mut hit),
        4 => screen::<4>(rows, fingerprint, count, &mut hit),
        _ => unreachable!("packed screen requires two- or four-bit fingerprints"),
    }
}

#[inline]
fn screen<const B: usize>(rows: &[&[u8]], fp: u64, count: usize, hit: &mut impl FnMut(usize)) {
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx2") {
        // Full blocks only: the last byte may also contain origin bits.
        unsafe { screen_avx2::<B>(rows, fp, count, hit) };
        return;
    }
    scalar::<B>(rows, fp, count, 0, hit);
}

#[inline(always)]
fn ones<const B: usize>() -> u64 {
    u64::MAX / ((1u64 << B) - 1)
}

#[inline(always)]
fn equal_fields<const B: usize>(diff: u64) -> u64 {
    let folded = diff | (diff >> 1);
    let folded = if B == 4 {
        folded | (folded >> 2)
    } else {
        folded
    };
    !folded & ones::<B>()
}

#[inline(always)]
fn emit<const B: usize>(mut mask: u64, base: usize, hit: &mut impl FnMut(usize)) {
    while mask != 0 {
        hit(base + mask.trailing_zeros() as usize / B);
        mask &= mask - 1;
    }
}

fn scalar<const B: usize>(
    rows: &[&[u8]],
    fp: u64,
    count: usize,
    start: usize,
    hit: &mut impl FnMut(usize),
) {
    let target = ones::<B>() * fp;
    let mut base = start;
    while count - base >= 64 / B {
        let byte = base * B / 8;
        let mut value = 0;
        for row in rows {
            value ^= u64::from_le_bytes(row[byte..byte + 8].try_into().unwrap());
        }
        emit::<B>(equal_fields::<B>(value ^ target), base, hit);
        base += 64 / B;
    }
    if base != count {
        let byte = base * B / 8;
        let bytes = ((count - base) * B).div_ceil(8);
        let mut value = 0;
        for row in rows {
            let mut word = [0; 8];
            word[..bytes].copy_from_slice(&row[byte..byte + bytes]);
            value ^= u64::from_le_bytes(word);
        }
        let valid = u64::MAX >> (64 - (count - base) * B);
        emit::<B>(equal_fields::<B>(value ^ target) & valid, base, hit);
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn screen_avx2<const B: usize>(
    rows: &[&[u8]],
    fp: u64,
    count: usize,
    hit: &mut impl FnMut(usize),
) {
    use std::arch::x86_64::*;
    let target = _mm256_set1_epi64x((ones::<B>() * fp) as i64);
    let low = _mm256_set1_epi64x(ones::<B>() as i64);
    let mut base = 0;
    while count - base >= 256 / B {
        let byte = base * B / 8;
        let mut value = _mm256_setzero_si256();
        for row in rows {
            value = _mm256_xor_si256(value, _mm256_loadu_si256(row.as_ptr().add(byte).cast()));
        }
        let diff = _mm256_xor_si256(value, target);
        let folded = _mm256_or_si256(diff, _mm256_srli_epi64::<1>(diff));
        let folded = if B == 4 {
            _mm256_or_si256(folded, _mm256_srli_epi64::<2>(folded))
        } else {
            folded
        };
        let mask = _mm256_andnot_si256(folded, low);
        if _mm256_testz_si256(mask, mask) == 0 {
            emit::<B>(_mm256_extract_epi64::<0>(mask) as u64, base, hit);
            emit::<B>(_mm256_extract_epi64::<1>(mask) as u64, base + 64 / B, hit);
            emit::<B>(_mm256_extract_epi64::<2>(mask) as u64, base + 128 / B, hit);
            emit::<B>(_mm256_extract_epi64::<3>(mask) as u64, base + 192 / B, hit);
        }
        base += 256 / B;
    }
    scalar::<B>(rows, fp, count, base, hit);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_byte_value_and_unaligned_simd_load_matches_independent_fields() {
        for bits in [2usize, 4] {
            for byte in 0..=255u8 {
                for offset in [1, 3, 7, 15, 31] {
                    let data = vec![byte; offset + 65];
                    for count in [
                        1,
                        3,
                        256 / bits - 1,
                        256 / bits,
                        256 / bits + 1,
                        512 / bits + 1,
                    ] {
                        let row = &data[offset..offset + (count * bits).div_ceil(8)];
                        for fp in 0..1 << bits {
                            let expected: Vec<_> = (0..count)
                                .filter(|&i| {
                                    ((byte >> (i * bits % 8)) as u64 & ((1 << bits) - 1)) == fp
                                })
                                .collect();
                            let mut actual = Vec::new();
                            matches_rows(&[row], bits, fp, count, |i| actual.push(i));
                            assert_eq!(
                                actual, expected,
                                "bits={bits} byte={byte} offset={offset} count={count} fp={fp}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn all_hash_counts_preserve_xor_cancellation_and_ignore_tail_padding() {
        for bits in [2usize, 4] {
            for count in [1, 3, 63, 64, 65, 127, 128, 129, 257] {
                for hashes in 1..=32 {
                    for fp in 0..1 << bits {
                        let repeated = ones_for_test(bits) * fp;
                        let mut data = vec![0; (count * bits).div_ceil(8)];
                        for (i, byte) in data.iter_mut().enumerate() {
                            *byte = (repeated >> (8 * (i % 8))) as u8;
                        }
                        // Corrupt every unused bit of the last byte. These may
                        // be the origin plane, not merely zero-filled padding.
                        if count * bits % 8 != 0 {
                            *data.last_mut().unwrap() |= 0xff << (count * bits % 8);
                        }
                        let rows = vec![data.as_slice(); hashes];
                        let target = if hashes % 2 == 0 { 0 } else { fp };
                        let mut actual = Vec::new();
                        matches_rows(&rows, bits, target, count, |i| actual.push(i));
                        assert_eq!(actual, (0..count).collect::<Vec<_>>());
                        actual.clear();
                        matches_rows(&rows, bits, target ^ 1, count, |i| actual.push(i));
                        assert!(actual.is_empty());
                    }
                }
            }
        }
    }

    fn ones_for_test(bits: usize) -> u64 {
        u64::MAX / ((1 << bits) - 1)
    }

    #[test]
    fn packed_simd_and_portable_screens_match_field_reference() {
        for bits in [2usize, 4] {
            for count in [
                0, 1, 2, 3, 7, 8, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129, 257,
            ] {
                for hashes in [1, 2, 3, 4, 8, 32] {
                    let data: Vec<Vec<u8>> = (0..hashes)
                        .map(|r| {
                            (0..(count * bits).div_ceil(8))
                                .map(|i| crate::splitmix64((r * 987 + i) as u64) as u8)
                                .collect()
                        })
                        .collect();
                    let rows: Vec<&[u8]> = data.iter().map(Vec::as_slice).collect();
                    for fp in 0..1u64 << bits {
                        let expected: Vec<usize> = (0..count)
                            .filter(|&i| {
                                rows.iter().fold(0u64, |v, r| {
                                    v ^ ((r[i * bits / 8] >> (i * bits % 8)) as u64
                                        & ((1 << bits) - 1))
                                }) == fp
                            })
                            .collect();
                        let mut actual = Vec::new();
                        matches_rows(&rows, bits, fp, count, |i| actual.push(i));
                        assert_eq!(
                            actual, expected,
                            "bits={bits}, count={count}, hashes={hashes}, fp={fp}"
                        );
                        actual.clear();
                        match bits {
                            2 => scalar::<2>(&rows, fp, count, 0, &mut |i| actual.push(i)),
                            4 => scalar::<4>(&rows, fp, count, 0, &mut |i| actual.push(i)),
                            _ => unreachable!(),
                        }
                        assert_eq!(actual, expected);
                    }
                }
            }
        }
    }
}
