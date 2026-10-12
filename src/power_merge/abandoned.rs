//! One open-addressed directory of 32-bit hashes with packed 32-bit dataset postings.
//! Different abandoned keys with the same hash intentionally share their origins.
use anyhow::{ensure, Context, Result};

const HASH_SEED: u64 = 0xD6E8_FEB8_6659_FD93;

#[inline]
pub(super) fn hash(key: u64, seed: u64) -> u32 {
    (crate::splitmix64(key ^ seed ^ HASH_SEED) >> 32) as u32
}

#[derive(Default)]
pub(super) struct OriginMap {
    keys: Vec<u32>,
    offsets: Vec<u64>,
    origins: Vec<u32>,
    // Zero means empty; otherwise this is a dense record index plus one.
    // No per-key allocation, and no 64-bit abandoned key is retained.
    directory: Vec<u32>,
}

impl OriginMap {
    pub(super) fn with_capacity(keys: usize, origins: usize) -> Result<Self> {
        ensure!(keys <= u32::MAX as usize, "too many abandoned hashes");
        Ok(Self {
            keys: Vec::with_capacity(keys),
            offsets: Vec::with_capacity(keys.checked_add(1).context("abandoned offsets overflow")?),
            origins: Vec::with_capacity(origins),
            directory: Vec::new(),
        })
    }
    pub(super) fn push(&mut self, key: u32, origins: &[u32]) -> Result<()> {
        ensure!(
            self.directory.is_empty(),
            "abandoned directory already finalized"
        );
        ensure!(
            self.keys.last().is_none_or(|&previous| previous < key),
            "unsorted/duplicate abandoned hashes"
        );
        ensure!(
            !origins.is_empty() && origins.windows(2).all(|pair| pair[0] < pair[1]),
            "empty/unsorted/duplicate abandoned origins"
        );
        ensure!(
            self.keys.len() < u32::MAX as usize,
            "too many abandoned hashes"
        );
        if self.offsets.is_empty() {
            self.offsets.push(0);
        }
        self.keys.push(key);
        self.origins.extend_from_slice(origins);
        self.offsets.push(self.origins.len() as u64);
        Ok(())
    }

    pub(super) fn finish(&mut self) -> Result<()> {
        ensure!(
            self.directory.is_empty(),
            "abandoned directory already finalized"
        );
        if self.keys.is_empty() {
            return Ok(());
        }
        // Legacy conversion does not know the number of globally distinct hashes
        // in advance. Discard vector growth slack before creating the directory.
        self.keys.shrink_to_fit();
        self.offsets.shrink_to_fit();
        self.origins.shrink_to_fit();
        let capacity = self
            .keys
            .len()
            .checked_mul(4)
            .and_then(|n| n.div_ceil(3).checked_next_power_of_two())
            .context("abandoned directory size overflow")?;
        self.directory = vec![0; capacity];
        let mask = capacity - 1;
        for (index, &key) in self.keys.iter().enumerate() {
            let mut slot = key as usize & mask;
            while self.directory[slot] != 0 {
                slot = (slot + 1) & mask;
            }
            self.directory[slot] = index as u32 + 1;
        }
        Ok(())
    }

    #[inline]
    pub(super) fn get(&self, key: u32) -> &[u32] {
        if self.directory.is_empty() {
            return &[];
        }
        let mask = self.directory.len() - 1;
        let mut slot = key as usize & mask;
        loop {
            let entry = self.directory[slot];
            if entry == 0 {
                return &[];
            }
            let index = entry as usize - 1;
            if self.keys[index] == key {
                return &self.origins
                    [self.offsets[index] as usize..self.offsets[index + 1] as usize];
            }
            slot = (slot + 1) & mask;
        }
    }

    pub(super) fn len(&self) -> usize {
        self.keys.len()
    }
    pub(super) fn origins_len(&self) -> usize {
        self.origins.len()
    }
    pub(super) fn memory_bytes(&self) -> usize {
        self.keys.capacity() * 4
            + self.offsets.capacity() * 8
            + self.origins.capacity() * 4
            + self.directory.capacity() * 4
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
    fn hash_uses_the_entire_input_key_and_seed() {
        assert_ne!(hash(1, 69), hash(1 + (1 << 32), 69));
        assert_ne!(hash(1, 69), hash(1, 70));
    }
}
