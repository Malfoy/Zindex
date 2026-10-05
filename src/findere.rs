//! Positional findere reconstruction shared by both index backends.
//! The backend supplies per-dataset s-mer hits; consecutive hits must belong
//! to the same dataset. Never combine scores across records or ambiguous bases.
use crate::{encode_base, is_compressed_path, scan_canonical_kmers, FastxStats, PARSER_CONFIG};
use anyhow::{ensure, Context, Result};
use deko::read::AnyDecoder;
use helicase::input::{FromMmap, FromSlice};
use helicase::{FastxParser, HelicaseParser};
use std::{fs::File, io::Read, path::Path};

pub const DEFAULT_Z: usize = 10;

pub fn indexed_length(k: usize, z: usize) -> Result<usize> {
    ensure!((1..=255).contains(&k), "findere query k must be in 1..=255");
    ensure!(z < k, "--findere-z must be smaller than query k");
    let s = k - z;
    ensure!(
        (1..=31).contains(&s),
        "findere indexed length k-z must be in 1..=31"
    );
    Ok(s)
}

#[derive(Debug)]
pub struct QueryResult {
    pub stats: FastxStats,
    /// Number of valid query k-mer positions, not the number of s-mer probes.
    pub total_kmers: u64,
    pub smer_queries: u64,
    pub scores: Vec<u64>,
}

/// Streaming reconstruction with O(number of datasets) scratch memory. Each
/// s-mer position is queried once, even when shared by z+1 adjacent k-mers.
pub struct Scorer {
    k: usize,
    s: usize,
    window: usize,
    hits: Vec<u64>,
    runs: Vec<usize>,
    pub result: QueryResult,
}

impl Scorer {
    pub fn new(k: usize, z: usize, colors: usize) -> Result<Self> {
        let s = indexed_length(k, z)?;
        Ok(Self {
            k,
            s,
            window: z + 1,
            hits: vec![0; colors],
            runs: vec![0; colors],
            result: QueryResult {
                stats: FastxStats::default(),
                total_kmers: 0,
                smer_queries: 0,
                scores: vec![0; colors],
            },
        })
    }

    pub fn sequence(&mut self, seq: &[u8], lookup: &mut impl FnMut(u64, &mut [u64])) {
        self.result.stats.bases += seq.len() as u64;
        self.result.stats.chunks += 1;
        for seq in seq.split(|&b| encode_base(b).is_none()) {
            if seq.len() < self.k {
                continue;
            }
            self.runs.fill(0);
            self.result.total_kmers += (seq.len() - self.k + 1) as u64;
            scan_canonical_kmers(seq, self.s, &mut |key| {
                self.hits.fill(0);
                lookup(key, &mut self.hits);
                self.result.smer_queries += 1;
                for ((run, &hit), score) in self
                    .runs
                    .iter_mut()
                    .zip(&self.hits)
                    .zip(&mut self.result.scores)
                {
                    *run = if hit == 0 {
                        0
                    } else {
                        (*run + 1).min(self.window)
                    };
                    if *run == self.window {
                        *score += 1;
                    }
                }
            });
        }
    }
}

pub fn query_file(
    path: &Path,
    k: usize,
    z: usize,
    colors: usize,
    mut lookup: impl FnMut(u64, &mut [u64]),
) -> Result<QueryResult> {
    let mut scorer = Scorer::new(k, z, colors)?;
    if is_compressed_path(path) {
        let mut data = Vec::new();
        AnyDecoder::new(File::open(path)?).read_to_end(&mut data)?;
        ensure!(!data.is_empty(), "{} is empty", path.display());
        let mut parser = FastxParser::<PARSER_CONFIG>::from_slice(&data)
            .with_context(|| format!("parsing {}", path.display()))?;
        while parser.next().is_some() {
            scorer.sequence(parser.get_dna_string(), &mut lookup);
        }
    } else {
        let mut parser = FastxParser::<PARSER_CONFIG>::from_file_mmap(path)
            .with_context(|| format!("parsing {}", path.display()))?;
        while parser.next().is_some() {
            scorer.sequence(parser.get_dna_string(), &mut lookup);
        }
    }
    Ok(scorer.result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn validates_lengths() {
        assert_eq!(indexed_length(31, DEFAULT_Z).unwrap(), 21);
        assert_eq!(indexed_length(41, 10).unwrap(), 31);
        assert_eq!(indexed_length(1, 0).unwrap(), 1);
        for (k, z) in [
            (0, 0),
            (10, 10),
            (31, 32),
            (42, 10),
            (256, 225),
            (31, usize::MAX),
        ] {
            assert!(indexed_length(k, z).is_err());
        }
    }

    #[test]
    fn hits_must_be_consecutive_and_in_the_same_color() {
        let mut scorer = Scorer::new(5, 2, 2).unwrap();
        let mut pos = 0;
        scorer.sequence(b"ACGTACGT", &mut |_, hits| {
            // Six s-mers: color 0 has runs of 2 and 3, color 1 only a run of 1.
            hits[if pos == 2 { 1 } else { 0 }] = 1;
            pos += 1;
        });
        assert_eq!(scorer.result.total_kmers, 4);
        assert_eq!(scorer.result.smer_queries, 6);
        assert_eq!(scorer.result.scores, [1, 0]);
    }

    #[test]
    fn boundaries_short_records_and_ambiguous_bases() {
        let mut scorer = Scorer::new(5, 2, 1).unwrap();
        let mut positive = |_: u64, hits: &mut [u64]| hits[0] = 1;
        scorer.sequence(b"ACGT", &mut positive);
        scorer.sequence(b"ACGT", &mut positive);
        scorer.sequence(b"ACGTNACGT", &mut positive);
        assert_eq!(scorer.result.total_kmers, 0);
        assert_eq!(scorer.result.smer_queries, 0);
        scorer.sequence(b"ACGTANACGTA", &mut positive);
        assert_eq!(scorer.result.total_kmers, 2);
        assert_eq!(scorer.result.scores, [2]);
        assert_eq!(scorer.result.smer_queries, 6);
    }

    #[test]
    fn matches_explicit_window_intersection_including_reverse_complement() {
        let reference = b"ACGTTAGCGGATCACGTAC";
        let queries: [&[u8]; 3] = [reference, b"GTACGTGATCCGCTAACGT", b"AAACGTTAGCGGTTCACGTAC"];
        for z in [0, 1, 3, 4] {
            let k = 7;
            let s = k - z;
            let mut set = HashSet::new();
            scan_canonical_kmers(reference, s, &mut |key| {
                set.insert(key);
            });
            for query in queries {
                let mut scorer = Scorer::new(k, z, 1).unwrap();
                scorer.sequence(query, &mut |key, hits| {
                    hits[0] = u64::from(set.contains(&key))
                });
                let expected = query
                    .windows(k)
                    .filter(|window| {
                        let mut all = true;
                        scan_canonical_kmers(window, s, &mut |key| all &= set.contains(&key));
                        all
                    })
                    .count() as u64;
                assert_eq!(scorer.result.scores, [expected]);
                assert_eq!(scorer.result.total_kmers, (query.len() - k + 1) as u64);
            }
        }
    }
}
