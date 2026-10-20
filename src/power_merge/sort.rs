//! Allocation-free integer sorting for large extracted feature sets.
//! Small partitions retain std's comparison sort; full u64 order is unchanged.

pub(super) fn sort_unique(keys: &mut Vec<u64>) {
    radix(keys, 56);
    keys.dedup();
}

fn radix(keys: &mut [u64], shift: u32) {
    if keys.len() < 2048 {
        keys.sort_unstable();
        return;
    }
    let mut counts = [0usize; 256];
    for &key in keys.iter() {
        counts[((key >> shift) & 255) as usize] += 1;
    }
    let mut starts = [0usize; 257];
    for i in 0..256 {
        starts[i + 1] = starts[i] + counts[i];
    }
    let mut next: [usize; 256] = starts[..256].try_into().unwrap();
    for bucket in 0..256 {
        while next[bucket] < starts[bucket + 1] {
            let dest = ((keys[next[bucket]] >> shift) & 255) as usize;
            if dest == bucket {
                next[bucket] += 1;
            } else {
                keys.swap(next[bucket], next[dest]);
                next[dest] += 1;
            }
        }
    }
    if shift != 0 {
        for bucket in 0..256 {
            radix(&mut keys[starts[bucket]..starts[bucket + 1]], shift - 8);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn integer_sort_matches_std_for_duplicates_prefixes_and_extremes() {
        for size in [0, 1, 2047, 2048, 10_000, 100_000] {
            for mode in 0..5 {
                let mut keys: Vec<_> = (0..size)
                    .map(|i| match mode {
                        0 => crate::splitmix64(i as u64),
                        1 => crate::splitmix64(i as u64) % 17,
                        2 => u64::MAX - i as u64,
                        3 => 0xabcd_0000_0000_0000 | i as u64,
                        _ => u64::MAX,
                    })
                    .collect();
                let mut expected = keys.clone();
                expected.sort_unstable();
                expected.dedup();
                sort_unique(&mut keys);
                assert_eq!(keys, expected);
            }
        }
    }
}
