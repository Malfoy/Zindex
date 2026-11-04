use super::*;

/// Sorted recovery records, decoded one hash and its dataset postings at a time.
pub(super) struct AbandonedReader<'a> {
    reader: BufReader<Box<dyn Read>>,
    manifest: &'a Manifest,
    delta_encoded: bool,
    record: usize,
    previous_hash: u32,
    associations: usize,
    pub(super) origins: Vec<u32>,
    bytes: Vec<u8>,
}

impl<'a> AbandonedReader<'a> {
    pub(super) fn new(directory: &Path, manifest: &'a Manifest) -> Result<Self> {
        ensure!(
            manifest.abandoned_file.as_deref() == Some(ABANDONED_FILE),
            "invalid abandoned map filename"
        );
        let mut reader = BufReader::with_capacity(
            256 * 1024,
            open_index_reader(&directory.join(ABANDONED_FILE))?,
        );
        let mut magic = [0; 8];
        reader.read_exact(&mut magic)?;
        ensure!(
            &magic == ABANDONED_MAGIC || &magic == ABANDONED_DELTA_MAGIC,
            "invalid abandoned map magic"
        );
        Ok(Self {
            reader,
            manifest,
            delta_encoded: &magic == ABANDONED_DELTA_MAGIC,
            record: 0,
            previous_hash: 0,
            associations: 0,
            origins: Vec::new(),
            bytes: vec![0; 64 * 1024],
        })
    }

    pub(super) fn advance(&mut self) -> Result<Option<u32>> {
        if self.record == self.manifest.abandoned_hashes {
            ensure!(
                self.associations == self.manifest.abandoned_origins,
                "abandoned origin count differs from manifest"
            );
            let mut extra = [0];
            ensure!(
                self.reader.read(&mut extra)? == 0,
                "unexpected abandoned map data"
            );
            return Ok(None);
        }
        let hash = if self.delta_encoded {
            let delta = read_var_u32(&mut self.reader)?;
            ensure!(
                self.record == 0 || delta != 0,
                "unsorted/duplicate abandoned hashes"
            );
            self.previous_hash
                .checked_add(delta)
                .context("abandoned hash delta overflow")?
        } else {
            read_u32(&mut self.reader)?
        };
        ensure!(
            self.record == 0 || self.previous_hash < hash,
            "unsorted/duplicate abandoned hashes"
        );
        self.previous_hash = hash;
        let count = if self.delta_encoded {
            read_var_u32(&mut self.reader)?
        } else {
            read_u32(&mut self.reader)?
        } as usize;
        ensure!(
            count > 0 && count <= self.manifest.colors.len(),
            "invalid abandoned origin count"
        );
        ensure!(
            count
                <= self
                    .manifest
                    .abandoned_origins
                    .saturating_sub(self.associations),
            "abandoned origins exceed manifest"
        );
        self.associations += count;
        self.origins.clear();
        if self.delta_encoded {
            let mut previous_origin = 0u32;
            for n in 0..count {
                let delta = read_var_u32(&mut self.reader)?;
                ensure!(n == 0 || delta != 0, "duplicate abandoned origins");
                let origin = previous_origin
                    .checked_add(delta)
                    .context("abandoned origin delta overflow")?;
                ensure!(
                    (origin as usize) < self.manifest.colors.len(),
                    "invalid abandoned dataset ID"
                );
                self.origins.push(origin);
                previous_origin = origin;
            }
        } else {
            let mut remaining = count;
            while remaining != 0 {
                let n = remaining.min(self.bytes.len() / 4);
                self.reader.read_exact(&mut self.bytes[..n * 4])?;
                for word in self.bytes[..n * 4].chunks_exact(4) {
                    let origin = u32::from_le_bytes(word.try_into().unwrap());
                    ensure!(
                        (origin as usize) < self.manifest.colors.len(),
                        "invalid abandoned dataset ID"
                    );
                    self.origins.push(origin);
                }
                remaining -= n;
            }
        }
        ensure!(
            self.origins.windows(2).all(|p| p[0] < p[1]),
            "unsorted/duplicate abandoned origins"
        );
        self.record += 1;
        Ok(Some(hash))
    }
}
