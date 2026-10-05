//! Positional findere reconstruction shared by both index backends.
//! The backend supplies per-dataset s-mer hits; consecutive hits must belong
//! to the same dataset. Never combine scores across records or ambiguous bases.
use crate::{encode_base, is_compressed_path, scan_canonical_kmers, FastxStats, PARSER_CONFIG};
use anyhow::{ensure, Context, Result};
use deko::read::AnyDecoder;
use helicase::input::{FromMmap, FromSlice};
use helicase::{FastxParser, HelicaseParser};
use rayon::prelude::*;
use std::{fs::File, io::Read, ops::Range, path::Path};

pub const DEFAULT_Z: usize = 10;
// Fixed partitioning makes matches AND actual probe counts independent of the
// worker count. Each task owns these k-mer starts plus k-1 trailing bases.
const TASK_KMERS: usize = 16_384;
const BATCH_BASES: usize = 1_048_576;

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
    /// Backend probes after optional union-filter rejection.
    pub backend_queries: u64,
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
                backend_queries: 0,
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

struct Batch {
    bases: Vec<u8>,
    ranges: Vec<Range<usize>>,
}

impl Batch {
    fn new() -> Self {
        Self {
            bases: Vec::with_capacity(BATCH_BASES + TASK_KMERS + 254),
            ranges: Vec::new(),
        }
    }
}

struct Worker<State> {
    scorer: Scorer,
    state: State,
}

fn flush_batch<State: Send>(
    batch: &mut Batch,
    workers: &mut Vec<Worker<State>>,
    k: usize,
    z: usize,
    colors: usize,
    init: &impl Fn() -> State,
    lookup: &(impl Fn(&mut State, u64, &mut [u64]) -> bool + Sync),
) {
    if batch.ranges.is_empty() {
        return;
    }
    let count = rayon::current_num_threads().min(batch.ranges.len());
    while workers.len() < count {
        workers.push(Worker {
            scorer: Scorer::new(k, z, colors).expect("validated findere lengths"),
            state: init(),
        });
    }
    let per_worker = batch.ranges.len().div_ceil(count);
    // Mutable state, decoding scratch, and counters belong to one worker. No
    // shared score-vector updates or atomics occur in the lookup hot path.
    workers[..count]
        .par_iter_mut()
        .zip(batch.ranges.par_chunks(per_worker))
        .for_each(|(worker, ranges)| {
            let mut backend_queries = 0;
            for range in ranges {
                worker
                    .scorer
                    .sequence(&batch.bases[range.clone()], &mut |key, hits| {
                        backend_queries += u64::from(lookup(&mut worker.state, key, hits));
                    });
            }
            worker.scorer.result.backend_queries += backend_queries;
        });
    batch.bases.clear();
    batch.ranges.clear();
}

fn query_parser<State: Send>(
    mut parser: FastxParser<'_, PARSER_CONFIG>,
    k: usize,
    z: usize,
    colors: usize,
    init: impl Fn() -> State,
    lookup: impl Fn(&mut State, u64, &mut [u64]) -> bool + Sync,
) -> Result<QueryResult> {
    let mut result = Scorer::new(k, z, colors)?.result;
    let mut workers = Vec::new();
    let mut batch = Batch::new();
    while parser.next().is_some() {
        let seq = parser.get_dna_string();
        result.stats.bases += seq.len() as u64;
        result.stats.chunks += 1;
        if seq.len() < k {
            continue;
        }
        let windows = seq.len() - k + 1;
        for start in (0..windows).step_by(TASK_KMERS) {
            let owned = TASK_KMERS.min(windows - start);
            let offset = batch.bases.len();
            // Overlap supplies context only: this sequence has exactly `owned`
            // k-mer windows, so no result is duplicated or omitted at a cut.
            batch
                .bases
                .extend_from_slice(&seq[start..start + owned + k - 1]);
            batch.ranges.push(offset..batch.bases.len());
            if batch.bases.len() >= BATCH_BASES {
                flush_batch(&mut batch, &mut workers, k, z, colors, &init, &lookup);
            }
        }
    }
    flush_batch(&mut batch, &mut workers, k, z, colors, &init, &lookup);
    for worker in workers {
        let local = worker.scorer.result;
        result.total_kmers += local.total_kmers;
        result.smer_queries += local.smer_queries;
        result.backend_queries += local.backend_queries;
        for (sum, value) in result.scores.iter_mut().zip(local.scores) {
            *sum += value;
        }
    }
    Ok(result)
}

/// Query on the caller's Rayon pool (the tools configure it with --threads).
/// Parsing/decompression feeds bounded batches; both many short records and a
/// single long record are distributed over workers. Scratch is reused across
/// batches. A lookup returns true if it actually probed the backend, false if
/// an optional prefilter rejected the key; hit entries start at zero either way.
pub fn query_file<State: Send>(
    path: &Path,
    k: usize,
    z: usize,
    colors: usize,
    init: impl Fn() -> State,
    lookup: impl Fn(&mut State, u64, &mut [u64]) -> bool + Sync,
) -> Result<QueryResult> {
    indexed_length(k, z)?;
    if is_compressed_path(path) {
        let mut data = Vec::new();
        AnyDecoder::new(File::open(path)?).read_to_end(&mut data)?;
        ensure!(!data.is_empty(), "{} is empty", path.display());
        let parser = FastxParser::<PARSER_CONFIG>::from_slice(&data)
            .with_context(|| format!("parsing {}", path.display()))?;
        query_parser(parser, k, z, colors, init, lookup)
    } else {
        let parser = FastxParser::<PARSER_CONFIG>::from_file_mmap(path)
            .with_context(|| format!("parsing {}", path.display()))?;
        query_parser(parser, k, z, colors, init, lookup)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;
    use std::sync::atomic::{AtomicUsize, Ordering};

    fn synthetic_hits(key: u64, hits: &mut [u64]) -> bool {
        hits[0] = 1;
        hits[1] = u64::from(key % 127 != 0);
        hits[2] = u64::from(key % 7 != 0);
        true
    }

    #[test]
    fn parallel_long_sequence_matches_sequential_reference_across_batches() {
        let mut seq: Vec<_> = (0..BATCH_BASES + 2 * TASK_KMERS + 503)
            .map(|i| b"ACGT"[(crate::splitmix64(i as u64) & 3) as usize])
            .collect();
        // Breaks around task cuts must not carry consecutive-hit state across N.
        seq[TASK_KMERS - 1] = b'N';
        seq[TASK_KMERS + 1] = b'N';
        seq[BATCH_BASES - 1] = b'N';
        let mut fasta = b">long\n".to_vec();
        fasta.extend_from_slice(&seq);
        fasta.extend_from_slice(
            b"\n>short\nACGT\n>second\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
        );
        for (k, z) in [(31, 10), (41, 10), (31, 0), (255, 224)] {
            let mut reference = Scorer::new(k, z, 3).unwrap();
            let mut parser = FastxParser::<PARSER_CONFIG>::from_slice(&fasta).unwrap();
            while parser.next().is_some() {
                reference.sequence(parser.get_dna_string(), &mut |key, hits| {
                    synthetic_hits(key, hits);
                });
            }
            let mut expected_probes = None;
            for threads in [1, 2, 4] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                let result = pool.install(|| {
                    query_parser(
                        FastxParser::<PARSER_CONFIG>::from_slice(&fasta).unwrap(),
                        k,
                        z,
                        3,
                        || (),
                        |_, key, hits| synthetic_hits(key, hits),
                    )
                    .unwrap()
                });
                assert_eq!(
                    result.scores, reference.result.scores,
                    "k={k},z={z},threads={threads}"
                );
                assert_eq!(result.total_kmers, reference.result.total_kmers);
                assert_eq!(result.stats.bases, reference.result.stats.bases);
                assert_eq!(result.stats.chunks, reference.result.stats.chunks);
                assert_eq!(result.backend_queries, result.smer_queries);
                assert!(result.smer_queries >= reference.result.smer_queries);
                if z == 0 {
                    assert_eq!(result.smer_queries, reference.result.smer_queries);
                }
                assert_eq!(
                    *expected_probes.get_or_insert(result.smer_queries),
                    result.smer_queries
                );
            }
        }
    }

    #[test]
    fn multiple_workers_execute_both_long_record_and_short_record_queries() {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap();
        for fasta in [
            format!(">long\n{}\n", "ACGT".repeat(TASK_KMERS * 4)),
            format!(">short\n{}\n", "ACGT".repeat(32)).repeat(2000),
        ] {
            let seen = AtomicUsize::new(0);
            let calls = AtomicUsize::new(0);
            let result = pool.install(|| {
                query_parser(
                    FastxParser::<PARSER_CONFIG>::from_slice(fasta.as_bytes()).unwrap(),
                    31,
                    10,
                    1,
                    || 0usize,
                    |local, _, hits| {
                        *local += 1;
                        seen.fetch_or(
                            1 << rayon::current_thread_index().expect("Rayon worker"),
                            Ordering::Relaxed,
                        );
                        calls.fetch_add(1, Ordering::Relaxed);
                            let probed = *local % 2 == 0;
                            hits[0] = u64::from(probed);
                            probed
                    },
                )
                .unwrap()
            });
            assert!(
                seen.load(Ordering::Relaxed).count_ones() > 1,
                "lookup ran on only one worker"
            );
            assert_eq!(result.scores[0], 0);
            assert_eq!(result.smer_queries, calls.load(Ordering::Relaxed) as u64);
            assert!(result.backend_queries > 0 && result.backend_queries < result.smer_queries);
        }
    }

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
