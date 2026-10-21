//! One open-addressed directory of 32-bit hashes with packed 32-bit dataset postings.
//! Different abandoned keys with the same hash intentionally share their origins.
use anyhow::{ensure, Context, Result};

const HASH_SEED: u64 = 0xD6E8_FEB8_6659_FD93;

#[inline]
pub(super) fn hash(key: u64, seed: u64) -> u32 {
    (crate::splitmix64(key ^ seed ^ HASH_SEED) >> 32) as u32
}

const MULTIPLE: u32 = 1 << 31;
const EMPTY: u32 = u32::MAX;
const WIDE: u32 = EMPTY - 1;

#[derive(Clone, Copy)]
struct Entry {
    hash: u32,
    // A normal document ID is inline. The high bit selects a multiple-origin
    // descriptor; EMPTY is a slot marker, never a reserved hash or origin ID.
    value: u32,
}
const EMPTY_ENTRY: Entry = Entry {
    hash: 0,
    value: EMPTY,
};

#[derive(Default)]
pub(super) struct OriginMap {
    directory: Vec<Entry>,
    multiple: Vec<std::ops::Range<usize>>,
    // Sorted overflow descriptors preserve the existing u32 key-count limit.
    wide: Vec<(u32, std::ops::Range<usize>)>,
    origins: Vec<u32>,
    hashes: usize,
    associations: usize,
    last_hash: Option<u32>,
    finalized: bool,
}

impl OriginMap {
    #[inline]
    fn first_slot(key: u32, capacity: usize) -> usize {
        // Recovery records arrive in hash order. Taking high hash bits makes
        // insertion a sequential pass instead of revisiting the entire table
        // once per high-bit prefix. Hashes are already uniformly mixed.
        let bits = capacity.trailing_zeros();
        if bits <= 32 {
            (key >> (32 - bits)) as usize
        } else {
            (key as usize) << (bits - 32)
        }
    }

    #[inline]
    pub(super) fn prefetch(&self, key: u32) {
        #[cfg(target_arch = "x86_64")]
        if !self.directory.is_empty() {
            let slot = Self::first_slot(key, self.directory.len());
            // Only the first directory bucket is known without doing the probe.
            // Actual collision resolution/postings remain unchanged in get().
            unsafe {
                std::arch::x86_64::_mm_prefetch(
                    self.directory.as_ptr().add(slot).cast(),
                    std::arch::x86_64::_MM_HINT_T0,
                );
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        let _ = key;
    }

    pub(super) fn with_capacity(keys: usize, origins: usize) -> Result<Self> {
        ensure!(keys <= u32::MAX as usize, "too many abandoned hashes");
        let capacity = if keys == 0 { 0 } else { Self::capacity(keys)? };
        let extra = origins.saturating_sub(keys);
        Ok(Self {
            directory: vec![EMPTY_ENTRY; capacity],
            multiple: Vec::with_capacity(extra),
            // At most twice the excess associations belong to non-singletons.
            // Do not reserve a posting for every inline singleton.
            origins: Vec::with_capacity(extra.saturating_mul(2).min(origins)),
            ..Default::default()
        })
    }

    fn capacity(keys: usize) -> Result<usize> {
        keys.checked_mul(4)
            .and_then(|n| n.div_ceil(3).max(4).checked_next_power_of_two())
            .context("abandoned directory size overflow")
    }

    fn insert(directory: &mut [Entry], entry: Entry) {
        let mask = directory.len() - 1;
        let mut slot = Self::first_slot(entry.hash, directory.len());
        while directory[slot].value != EMPTY {
            slot = (slot + 1) & mask;
        }
        directory[slot] = entry;
    }

    fn reserve_one(&mut self) -> Result<()> {
        if (self.hashes + 1) * 4 <= self.directory.len() * 3 {
            return Ok(());
        }
        let capacity = Self::capacity(self.hashes + 1)?;
        let mut directory = vec![EMPTY_ENTRY; capacity];
        for &entry in &self.directory {
            if entry.value != EMPTY {
                Self::insert(&mut directory, entry);
            }
        }
        self.directory = directory;
        Ok(())
    }
    fn add_multiple(&mut self, key: u32, origins: &[u32], limit: usize) -> u32 {
        let start = self.origins.len();
        self.origins.extend_from_slice(origins);
        let range = start..self.origins.len();
        if self.multiple.len() < limit {
            let id = self.multiple.len() as u32;
            self.multiple.push(range);
            MULTIPLE | id
        } else {
            self.wide.push((key, range));
            WIDE
        }
    }

    pub(super) fn push(&mut self, key: u32, origins: &[u32]) -> Result<()> {
        ensure!(!self.finalized, "abandoned directory already finalized");
        ensure!(
            self.last_hash.is_none_or(|previous| previous < key),
            "unsorted/duplicate abandoned hashes"
        );
        ensure!(
            !origins.is_empty() && origins.windows(2).all(|pair| pair[0] < pair[1]),
            "empty/unsorted/duplicate abandoned origins"
        );
        ensure!(self.hashes < u32::MAX as usize, "too many abandoned hashes");
        self.reserve_one()?;
        let value = if origins.len() == 1 && origins[0] < MULTIPLE {
            origins[0]
        } else {
            self.add_multiple(key, origins, (MULTIPLE - 2) as usize)
        };
        Self::insert(&mut self.directory, Entry { hash: key, value });
        self.hashes += 1;
        self.associations = self
            .associations
            .checked_add(origins.len())
            .context("abandoned origins overflow")?;
        self.last_hash = Some(key);
        Ok(())
    }

    pub(super) fn finish(&mut self) -> Result<()> {
        ensure!(!self.finalized, "abandoned directory already finalized");
        if self.hashes == 0 {
            return Ok(());
        }
        // Discard posting/descriptor growth slack, including legacy conversion
        // whose globally distinct hash count was not known in advance.
        self.multiple.shrink_to_fit();
        self.wide.shrink_to_fit();
        self.origins.shrink_to_fit();
        self.finalized = true;
        Ok(())
    }

    #[inline]
    pub(super) fn get(&self, key: u32) -> &[u32] {
        if self.directory.is_empty() {
            return &[];
        }
        let mask = self.directory.len() - 1;
        let mut slot = Self::first_slot(key, self.directory.len());
        loop {
            let entry = &self.directory[slot];
            if entry.value == EMPTY {
                return &[];
            }
            if entry.hash == key {
                return if entry.value & MULTIPLE == 0 {
                    std::slice::from_ref(&entry.value)
                } else if entry.value == WIDE {
                    let id = self
                        .wide
                        .binary_search_by_key(&key, |entry| entry.0)
                        .unwrap();
                    &self.origins[self.wide[id].1.clone()]
                } else {
                    &self.origins[self.multiple[(entry.value & !MULTIPLE) as usize].clone()]
                };
            }
            slot = (slot + 1) & mask;
        }
    }

    pub(super) fn len(&self) -> usize {
        self.hashes
    }
    pub(super) fn origins_len(&self) -> usize {
        self.associations
    }
    pub(super) fn memory_bytes(&self) -> usize {
        self.directory.capacity() * std::mem::size_of::<Entry>()
            + self.multiple.capacity() * std::mem::size_of::<std::ops::Range<usize>>()
            + self.wide.capacity() * std::mem::size_of::<(u32, std::ops::Range<usize>)>()
            + self.origins.capacity() * 4
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_map_and_zero_and_max_hashes() {
        let mut map = OriginMap::with_capacity(2, 3).unwrap();
        map.finish().unwrap();
        assert!(map.get(0).is_empty());
        map.push(0, &[0, u32::MAX]).unwrap();
        map.push(u32::MAX, &[19]).unwrap();
        map.finish().unwrap();
        assert_eq!(map.len(), 2);
        assert_eq!(map.origins_len(), 3);
        assert!(map.memory_bytes() >= 2 * 4 + 3 * 8 + 3 * 4);
        assert_eq!(map.get(0), [0, u32::MAX]);
        assert_eq!(map.get(u32::MAX), [19]);
        assert!(map.get(1).is_empty());
    }

    #[test]
    fn probing_collision_is_not_a_key_match() {
        let mut map = OriginMap::default();
        for key in [0, 8, 16, 24] {
            map.push(key, &[key + 1]).unwrap();
        }
        map.finish().unwrap();
        for key in [0, 8, 16, 24] {
            assert_eq!(map.get(key), [key + 1]);
        }
        assert!(map.get(32).is_empty());
    }

    #[test]
    fn validates_postings_and_record_order() {
        let mut map = OriginMap::default();
        assert!(map.push(1, &[]).is_err());
        assert!(map.push(1, &[2, 2]).is_err());
        assert!(map.push(1, &[2, 1]).is_err());
        map.push(1, &[1, 2]).unwrap();
        assert!(map.push(1, &[3]).is_err());
        map.finish().unwrap();
        assert!(map.push(2, &[3]).is_err());
    }

    #[test]
    fn wide_descriptor_escape_does_not_limit_origin_ids_or_hashes() {
        let mut map = OriginMap::with_capacity(2, 4).unwrap();
        for (hash, origins) in [(0, &[0, u32::MAX][..]), (u32::MAX, &[1, 2][..])] {
            let value = map.add_multiple(hash, origins, 0);
            OriginMap::insert(&mut map.directory, Entry { hash, value });
            assert_eq!(map.get(hash), origins);
        }
        assert!(map.get(1).is_empty());
    }

    #[test]
    fn inline_and_multiple_origins_match_a_reference_through_growth() {
        for reserved in [false, true] {
            let mut hashes: Vec<_> = (0..100_000u64)
                .map(|i| crate::splitmix64(i) as u32)
                .collect();
            hashes.extend([0, u32::MAX]);
            hashes.sort_unstable();
            hashes.dedup();
            let mut map = if reserved {
                OriginMap::with_capacity(hashes.len(), hashes.len() * 2).unwrap()
            } else {
                OriginMap::default()
            };
            let mut expected = std::collections::BTreeMap::new();
            for (i, hash) in hashes.into_iter().enumerate() {
                let values = if i % 127 == 0 {
                    vec![u32::MAX]
                } else if i % 25 == 0 {
                    vec![0, (i as u32 % 16_000) + 1, u32::MAX]
                } else {
                    vec![i as u32 % 16_384]
                };
                map.push(hash, &values).unwrap();
                expected.insert(hash, values);
            }
            map.finish().unwrap();
            assert_eq!(map.len(), expected.len());
            assert_eq!(
                map.origins_len(),
                expected.values().map(Vec::len).sum::<usize>()
            );
            for (&hash, values) in &expected {
                assert_eq!(map.get(hash), values);
            }
            for i in 0..1000 {
                let hash = crate::splitmix64(i + 100_000) as u32;
                assert_eq!(
                    map.get(hash),
                    expected.get(&hash).map_or(&[][..], Vec::as_slice)
                );
            }
        }
    }

    #[test]
    fn hash_uses_the_entire_input_key_and_seed() {
        assert_ne!(hash(1, 69), hash(1 + (1 << 32), 69));
        assert_ne!(hash(1, 69), hash(1, 70));
    }
}
