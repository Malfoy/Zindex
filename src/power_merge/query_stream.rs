//! Extract seeds once, scan the auxiliary file once, and keep results in RAM
//! while loading one independently compressed subindex at a time.
use super::*;

// Ordinary dense/minimizer queries do not need a global key->ID table or a vector of
// hits for every key/subindex pair. Keep extracted reads and query each loaded
// subindex with the same exact threshold-pruning worker as full-RAM mode.
struct QueryRecovery {
    entries: Vec<(u32, Vec<u32>)>,
    scans: usize,
    hashes: usize,
}

impl QueryRecovery {
    fn read(directory: &Path, m: &Manifest, reads: &[PreparedRead], recover: bool) -> Result<Self> {
        let mut result = Self {
            entries: Vec::new(),
            scans: 0,
            hashes: 0,
        };
        if !recover || !m.abandoned_recovery {
            return Ok(result);
        }
        let mut hashes: Vec<u32> = reads
            .iter()
            .flat_map(|r| &r.forward)
            .map(|&key| abandoned::hash(key, m.seed))
            .collect();
        hashes.par_sort_unstable();
        hashes.dedup();
        result.hashes = hashes.len();
        let mut at = 0;
        // Match the complete 32-bit recovery equivalence class, even when two
        // different queried 64-bit keys collide or route to different parts.
        visit_abandoned_file(directory, m, |hash, origins| {
            while at < hashes.len() && hashes[at] < hash {
                at += 1;
            }
            if at < hashes.len() && hashes[at] == hash {
                result.entries.push((hash, origins.to_vec()));
            }
            Ok(())
        })?;
        result.scans = 1;
        Ok(result)
    }

    fn selected(&self, remap: &[u32]) -> Result<Option<OriginMap>> {
        // These hashes arrive in increasing order. Growing the high-bit
        // directory while only a hash prefix has arrived creates long probe
        // clusters. Count in RAM first, then allocate the final directory once;
        // this does not rescan either the query input or the recovery file.
        let mut hashes = 0usize;
        let mut associations = 0usize;
        for (_, origins) in &self.entries {
            let count = origins
                .iter()
                .filter(|&&id| remap[id as usize] != u32::MAX)
                .count();
            if count != 0 {
                hashes += 1;
                associations = associations
                    .checked_add(count)
                    .context("selected recovery origins overflow")?;
            }
        }
        if hashes == 0 {
            return Ok(None);
        }
        let mut map = OriginMap::with_capacity(hashes, associations)?;
        let mut selected = Vec::new();
        for (hash, origins) in &self.entries {
            selected.clear();
            selected.extend(
                origins.iter().filter_map(|&id| {
                    (remap[id as usize] != u32::MAX).then_some(remap[id as usize])
                }),
            );
            if !selected.is_empty() {
                map.push(*hash, &selected)?;
            }
        }
        map.finish()?;
        debug_assert_eq!(map.len(), hashes);
        debug_assert_eq!(map.origins_len(), associations);
        Ok(Some(map))
    }
}

fn ordinary_features(index: &Index, f: FeatureConfig) -> bool {
    index.manifest.version >= 2
        && matches!(f.mode, IndexMode::Kmers | IndexMode::Minimizers)
        && f.randstrobemers.is_none()
}

// Counting-sort occurrences by partition, preserving read order within every
// partition. Storage is linear in occurrences+partitions, not reads*partitions.
fn route_reads(
    reads: &[PreparedRead],
    seed: u64,
    partitions: usize,
) -> (Vec<(usize, u64)>, Vec<usize>) {
    let mut offsets = vec![0usize; partitions + 1];
    for read in reads {
        for &key in &read.forward {
            offsets[partition_id(key, seed, partitions) + 1] += 1;
        }
    }
    for i in 1..offsets.len() {
        offsets[i] += offsets[i - 1];
    }
    let mut cursor = offsets[..partitions].to_vec();
    let mut routed = vec![(0, 0); offsets[partitions]];
    for (id, read) in reads.iter().enumerate() {
        for &key in &read.forward {
            let p = partition_id(key, seed, partitions);
            routed[cursor[p]] = (id, key);
            cursor[p] += 1;
        }
    }
    (routed, offsets)
}

fn install_bank(index: &mut Index, mut bank: Bank, remap: &[u32]) {
    assert!(index.banks.is_empty(), "one resident subindex at a time");
    for group in &mut bank.meta.groups {
        for origin in &mut group.members {
            *origin = remap[*origin] as usize;
        }
    }
    index.row_scratch_bytes = if bank.meta.fingerprint_first {
        0
    } else {
        bank.meta.row_bits.div_ceil(8)
    };
    index.banks.push(bank);
}

fn append_scores(
    out: &mut Vec<(u32, u64)>,
    scores: &RecordScores,
    documents: &[usize],
    total: u64,
    ratio: f64,
) {
    if ratio == 0.0 {
        out.extend(
            documents
                .iter()
                .enumerate()
                .map(|(local, &global)| (global as u32, scores.values[local])),
        );
    } else if total != 0 {
        for &local in &scores.touched {
            let count = scores.values[local];
            if count as f64 / total as f64 >= ratio {
                out.push((documents[local] as u32, count));
            }
        }
    }
}

fn query_ordinary(
    directory: &Path,
    query: &Path,
    ratio: f64,
    recover: bool,
    batch_size: usize,
    records: bool,
    format: RecordFormat,
    mut index: Index,
    f: FeatureConfig,
    start: Instant,
    metadata_seconds: f64,
) -> Result<()> {
    let preparation = Instant::now();
    let mut reads = if records {
        super::super::records::prepare_file(query, batch_size, |(name, seq)| {
            PreparedRead::new(name, seq, &index.manifest, f)
        })?
    } else {
        let mut keys = Vec::new();
        parse_features_from_file(query, index.manifest.k, index.manifest.seed, f, |key| {
            keys.push(key)
        })?;
        // Extraction has already preserved sequence/window boundaries. Splitting
        // feature occurrences is exact for aggregate dense/minimizer counts.
        keys.chunks(32 * 1024)
            .map(|chunk| PreparedRead {
                name: Vec::new(),
                forward: chunk.to_vec(),
                reverse: Vec::new(),
                runs: Vec::new(),
                total: chunk.len() as u64,
            })
            .collect()
    };
    let preparation_seconds = preparation.elapsed().as_secs_f64();
    let total: u64 = reads.iter().map(|r| r.total).sum();
    let mut groups = Vec::new();
    let mut begin = 0;
    while begin < index.manifest.banks.len() {
        let source = index.manifest.banks[begin].source_group();
        let mut end = begin + 1;
        while end < index.manifest.banks.len() && index.manifest.banks[end].source_group() == source
        {
            end += 1;
        }
        groups.push((begin, end));
        begin = end;
    }
    // Hide the largest load behind recovery, without ever materializing a
    // second subindex. Choose a partition actually used by a query, preserving
    // the existing ability to omit unused physical partition payloads.
    groups.sort_unstable_by_key(|&(a, _)| {
        Reverse(
            index.manifest.banks[a]
                .slots
                .saturating_mul(index.manifest.banks[a].row_bits.div_ceil(8)),
        )
    });
    let first_key = reads.iter().find_map(|r| r.forward.first()).copied();
    let first_id = first_key.map(|key| {
        let (a, b) = groups[0];
        a + partition_id(key, index.manifest.seed, b - a)
    });
    let acquire = Instant::now();
    let (first, recovery) = rayon::join(
        || -> Result<Option<(usize, Bank)>> {
            first_id
                .map(|id| load_modern_bank(directory, &index.manifest.banks[id]).map(|b| (id, b)))
                .transpose()
        },
        || QueryRecovery::read(directory, &index.manifest, &reads, recover),
    );
    let mut first = first?;
    let recovery = recovery?;
    let mut load_seconds = metadata_seconds + acquire.elapsed().as_secs_f64();
    let mut lookup_seconds = 0.0;
    let mut recovery_map_seconds = 0.0;
    let mut result = if records {
        vec![Vec::new(); reads.len()]
    } else {
        Vec::new()
    };
    let mut aggregate = vec![0u64; index.manifest.colors.len()];
    for (a, b) in groups {
        let documents = document_ids(&index.manifest.banks[a]);
        let mut remap = vec![u32::MAX; index.manifest.colors.len()];
        for (local, &global) in documents.iter().enumerate() {
            remap[global] = local as u32;
        }
        let recovery_start = Instant::now();
        index.abandoned = recovery.selected(&remap)?;
        recovery_map_seconds += recovery_start.elapsed().as_secs_f64();
        if b - a == 1 {
            if first_key.is_some() {
                let load = Instant::now();
                let bank = if first.as_ref().is_some_and(|(id, _)| *id == a) {
                    first.take().unwrap().1
                } else {
                    let bank = load_modern_bank(directory, &index.manifest.banks[a])?;
                    load_seconds += load.elapsed().as_secs_f64();
                    bank
                };
                install_bank(&mut index, bank, &remap);
                let lookup = Instant::now();
                let chunk = reads.len().div_ceil(rayon::current_num_threads()).max(1);
                if records {
                    reads
                        .par_chunks_mut(chunk)
                        .zip(result.par_chunks_mut(chunk))
                        .for_each(|(reads, outputs)| {
                            let mut worker = record_query::Worker::with_documents(documents.len());
                            for (read, out) in reads.iter_mut().zip(outputs) {
                                worker.query_prepared(&index, &mut read.forward, ratio);
                                append_scores(out, &worker.scores, &documents, read.total, ratio);
                            }
                        });
                } else {
                    let partials: Vec<_> = reads
                        .par_chunks_mut(chunk)
                        .map(|reads| {
                            let mut worker = record_query::Worker::with_documents(documents.len());
                            let mut sums = vec![0u64; documents.len()];
                            for read in reads {
                                worker.query_prepared(&index, &mut read.forward, 0.0);
                                for &id in &worker.scores.touched {
                                    sums[id] += worker.scores.values[id];
                                }
                            }
                            sums
                        })
                        .collect();
                    for sums in partials {
                        for (&global, count) in documents.iter().zip(sums) {
                            aggregate[global] += count;
                        }
                    }
                }
                lookup_seconds += lookup.elapsed().as_secs_f64();
                index.banks.clear();
            } else if records && ratio == 0.0 {
                for out in &mut result {
                    out.extend(documents.iter().map(|&global| (global as u32, 0)));
                }
            }
        } else {
            let (routed, offsets) = route_reads(&reads, index.manifest.seed, b - a);
            let mut partial = vec![Vec::<(usize, u64)>::new(); reads.len()];
            let mut remaining: Vec<_> = reads.iter().map(|r| r.total).collect();
            let mut parts: Vec<_> = (a..b).collect();
            if let Some((id, _)) = &first {
                if (a..b).contains(id) {
                    parts.swap(0, *id - a);
                }
            }
            let active_ratio = if records { ratio } else { 0.0 };
            for id in parts {
                let p = id - a;
                let keys = &routed[offsets[p]..offsets[p + 1]];
                if keys.is_empty() {
                    continue;
                }
                let load = Instant::now();
                let bank = if first.as_ref().is_some_and(|(loaded, _)| *loaded == id) {
                    first.take().unwrap().1
                } else {
                    let bank = load_modern_bank(directory, &index.manifest.banks[id])?;
                    load_seconds += load.elapsed().as_secs_f64();
                    bank
                };
                install_bank(&mut index, bank, &remap);
                let mut ranges = vec![0..0; reads.len()];
                let mut at = 0;
                while at < keys.len() {
                    let start = at;
                    let read = keys[at].0;
                    while at < keys.len() && keys[at].0 == read {
                        at += 1;
                    }
                    ranges[read] = start..at;
                }
                let lookup = Instant::now();
                partial
                    .par_iter_mut()
                    .zip(remaining.par_iter_mut())
                    .enumerate()
                    .for_each_init(
                        || record_query::Worker::with_documents(documents.len()),
                        |worker, (read, (previous, left))| {
                            let range = ranges[read].clone();
                            if range.is_empty() {
                                return;
                            }
                            let after = *left - range.len() as u64;
                            worker.query_partition(
                                &index,
                                keys[range].iter().map(|&(_, key)| key),
                                reads[read].total,
                                *left,
                                previous,
                                active_ratio,
                            );
                            previous.clear();
                            for &local in &worker.scores.touched {
                                let count = worker.scores.values[local];
                                if active_ratio == 0.0
                                    || (count + after) as f64 / reads[read].total as f64
                                        >= active_ratio
                                {
                                    previous.push((local, count));
                                }
                            }
                            *left = after;
                        },
                    );
                lookup_seconds += lookup.elapsed().as_secs_f64();
                index.banks.clear();
            }
            for (read, counts) in partial.into_iter().enumerate() {
                if records {
                    if ratio == 0.0 {
                        let start = result[read].len();
                        result[read].extend(documents.iter().map(|&global| (global as u32, 0)));
                        for (local, count) in counts {
                            result[read][start + local].1 = count;
                        }
                    } else if reads[read].total != 0 {
                        for (local, count) in counts {
                            if count as f64 / reads[read].total as f64 >= ratio {
                                result[read].push((documents[local] as u32, count));
                            }
                        }
                    }
                } else {
                    for (local, count) in counts {
                        aggregate[documents[local]] += count;
                    }
                }
            }
        }
        index.abandoned = None;
    }
    debug_assert!(first.is_none());
    let output = Instant::now();
    if records {
        let mut out = BufWriter::with_capacity(1024 * 1024, std::io::stdout());
        writeln!(out, "{}", format.header())?;
        for (read, hits) in reads.iter().zip(&mut result) {
            hits.sort_unstable_by_key(|hit| hit.0);
            format.write_record(
                &mut out,
                &read.name,
                read.total,
                hits.iter().map(|&(id, count)| (id as usize, count)),
                &index.manifest.colors,
            )?;
        }
        out.flush()?;
    } else {
        let classes = index
            .manifest
            .banks
            .iter()
            .map(BankMeta::source_group)
            .collect::<std::collections::BTreeSet<_>>()
            .len();
        print_aggregate(
            &index.manifest,
            query,
            &aggregate,
            total,
            total,
            total * classes as u64,
            ratio,
            start.elapsed().as_secs_f64() - load_seconds,
        );
    }
    let output_seconds = output.elapsed().as_secs_f64();
    eprintln!(
        "query_memory\tmode=subindex\tmax_loaded_subindexes=1\tsubindexes={}\tresult_storage=RAM",
        index.manifest.banks.len()
    );
    eprintln!("query_preparation\tinput_passes=1\tseed_extraction_passes=1\tabandoned_file_scans={}\tquery_features={total}\tqueried_recovery_hashes={}\tlookup=record-pruned",recovery.scans,recovery.hashes);
    eprintln!("query_phases\tprepare_s={preparation_seconds:.6}\trecovery_map_s={recovery_map_seconds:.6}\tlookup_s={lookup_seconds:.6}\toutput_s={output_seconds:.6}");
    eprintln!(
        "query_timings\tload_s={load_seconds:.6}\tquery_s={:.6}",
        start.elapsed().as_secs_f64() - load_seconds
    );
    Ok(())
}

struct PreparedRead {
    name: Vec<u8>,
    // Keys initially, then indices into the shared unique-key array.
    forward: Vec<u64>,
    reverse: Vec<u64>,
    runs: Vec<Range<usize>>,
    total: u64,
}

impl PreparedRead {
    fn new(name: &[u8], seq: &[u8], m: &Manifest, f: FeatureConfig) -> Result<Self> {
        // Reserve bounded short-read feature storage rather than repeatedly
        // growing it. Long or ambiguous records never reserve an unbounded
        // sequence-length estimate before their features are known.
        let capacity = if f.mode == IndexMode::Minimizers
            && f.randstrobemers.is_none()
            && !f.homopolymer_compression
        {
            seq.len()
                .checked_sub(m.k)
                .map_or(0, |n| (2 * (n + 1)).div_ceil(m.k - f.minimizer_size + 2))
                .min(256)
        } else if f.mode == IndexMode::Kmers && f.randstrobemers.is_none() {
            seq.len().checked_sub(m.k).map_or(0, |n| n + 1).min(1024)
        } else {
            0
        };
        let mut read = Self {
            name: name.to_vec(),
            forward: Vec::with_capacity(capacity),
            reverse: Vec::new(),
            runs: Vec::new(),
            total: 0,
        };
        if f.mode == IndexMode::Findere {
            let s = crate::findere::indexed_length(m.k, f.findere_z)?;
            for run in seq.split(|&b| crate::encode_base(b).is_none()) {
                let compressed;
                let run = if f.homopolymer_compression {
                    compressed = crate::homopolymer_compress(run);
                    compressed.as_slice()
                } else {
                    run
                };
                if run.len() < m.k {
                    continue;
                }
                let windows = run.len() - m.k + 1;
                read.total += windows as u64;
                // Preserve the existing parallel findere task boundaries and
                // logical probe counts, including their s-mer context overlap.
                for start in (0..windows).step_by(16_384) {
                    let owned = 16_384.min(windows - start);
                    let first = read.forward.len();
                    crate::scan_canonical_kmers(
                        &run[start..start + owned + m.k - 1],
                        s,
                        &mut |key| read.forward.push(key),
                    );
                    read.runs.push(first..read.forward.len());
                }
            }
        } else {
            crate::scan_index_features(seq, m.k, m.seed, f, &mut |key| read.forward.push(key));
            read.total = read.forward.len() as u64;
            if f.randstrobemers.is_some() {
                crate::scan_index_features(
                    &crate::reverse_complement(seq),
                    m.k,
                    m.seed,
                    f,
                    &mut |key| read.reverse.push(key),
                );
                ensure!(
                    read.forward.len() == read.reverse.len(),
                    "canonical anchor counts must agree across strands"
                );
            }
        }
        Ok(read)
    }
}

struct Prepared {
    reads: Vec<PreparedRead>,
    keys: Vec<u64>,
    abandoned: Vec<Vec<u32>>,
    total: u64,
    probes: u64,
}

impl Prepared {
    fn new(mut reads: Vec<PreparedRead>, f: FeatureConfig) -> Self {
        let mut keys: Vec<_> = reads
            .iter()
            .flat_map(|r| r.forward.iter().chain(&r.reverse).copied())
            .collect();
        super::super::sort::sort_unique(&mut keys);
        let ids: HashMap<_, _> = keys
            .iter()
            .enumerate()
            .map(|(id, &key)| (key, id as u64))
            .collect();
        let (mut total, mut probes) = (0, 0);
        for read in &mut reads {
            total += read.total;
            probes += (read.forward.len() + read.reverse.len()) as u64;
            for key in read.forward.iter_mut().chain(&mut read.reverse) {
                *key = ids[key];
            }
            if f.mode != IndexMode::Findere {
                read.forward.sort_unstable();
                read.reverse.sort_unstable();
            }
        }
        let abandoned = vec![Vec::new(); keys.len()];
        Self {
            reads,
            keys,
            abandoned,
            total,
            probes,
        }
    }

    fn recover(&mut self, directory: &Path, m: &Manifest, recover: bool) -> Result<usize> {
        if !recover || !m.abandoned_recovery || m.version == 1 {
            return Ok(0);
        }
        // One sequential merge join. Queried 64-bit keys sharing a 32-bit hash
        // intentionally receive the same union of abandoned-key origins.
        let mut hashes: Vec<_> = self
            .keys
            .iter()
            .enumerate()
            .map(|(id, &key)| (abandoned::hash(key, m.seed), id))
            .collect();
        hashes.sort_unstable_by_key(|entry| entry.0);
        let mut position = 0;
        visit_abandoned_file(directory, m, |hash, origins| {
            while position < hashes.len() && hashes[position].0 < hash {
                position += 1;
            }
            while position < hashes.len() && hashes[position].0 == hash {
                self.abandoned[hashes[position].1].extend_from_slice(origins);
                position += 1;
            }
            Ok(())
        })?;
        Ok(1)
    }
}

fn document_ids(meta: &BankMeta) -> Vec<usize> {
    let mut ids: Vec<_> = meta
        .groups
        .iter()
        .flat_map(|g| g.members.iter().copied())
        .collect();
    ids.sort_unstable();
    ids
}

fn route_queries(keys: &[u64], seed: u64, partitions: usize) -> (Vec<usize>, Vec<usize>) {
    if partitions == 1 {
        return ((0..keys.len()).collect(), vec![0, keys.len()]);
    }
    let mut offsets = vec![0usize; partitions + 1];
    for &key in keys {
        offsets[partition_id(key, seed, partitions) + 1] += 1;
    }
    for id in 1..offsets.len() {
        offsets[id] += offsets[id - 1];
    }
    let mut cursor = offsets[..partitions].to_vec();
    let mut order = vec![0; keys.len()];
    for (id, &key) in keys.iter().enumerate() {
        let part = partition_id(key, seed, partitions);
        order[cursor[part]] = id;
        cursor[part] += 1;
    }
    (order, offsets)
}

struct Scores {
    forward: RecordScores,
    reverse: RecordScores,
    last: Vec<usize>,
    runs: Vec<usize>,
}

impl Scores {
    fn new(n: usize, f: FeatureConfig) -> Self {
        Self {
            forward: RecordScores::new(n),
            reverse: RecordScores::new(if f.randstrobemers.is_some() { n } else { 0 }),
            last: if f.mode == IndexMode::Findere {
                vec![usize::MAX; n]
            } else {
                Vec::new()
            },
            runs: if f.mode == IndexMode::Findere {
                vec![0; n]
            } else {
                Vec::new()
            },
        }
    }
    fn strand(keys: &[u64], hits: &[Vec<u32>], scores: &mut RecordScores) {
        scores.clear();
        for run in keys.chunk_by(|a, b| a == b) {
            for &origin in &hits[run[0] as usize] {
                scores.add(origin as usize, run.len() as u64);
            }
        }
    }
    fn score(&mut self, read: &PreparedRead, hits: &[Vec<u32>], f: FeatureConfig) {
        if f.mode == IndexMode::Findere {
            self.forward.clear();
            for range in &read.runs {
                self.last.fill(usize::MAX);
                for (position, &key) in read.forward[range.clone()].iter().enumerate() {
                    for &origin in &hits[key as usize] {
                        let id = origin as usize;
                        self.runs[id] = if self.last[id].checked_add(1) == Some(position) {
                            self.runs[id] + 1
                        } else {
                            1
                        };
                        self.last[id] = position;
                        if self.runs[id] > f.findere_z {
                            self.forward.increment(id);
                        }
                    }
                }
            }
        } else {
            Self::strand(&read.forward, hits, &mut self.forward);
            if f.randstrobemers.is_some() {
                Self::strand(&read.reverse, hits, &mut self.reverse);
                for &id in &self.reverse.touched {
                    let count = self.reverse.values[id];
                    if count > self.forward.values[id] {
                        self.forward.add(id, count - self.forward.values[id]);
                    }
                }
            }
        }
    }
}

// Resolve each size class after all its partitions have contributed. This
// preserves findere windows crossing partitions and per-document strand maxima.
fn execute(
    directory: &Path,
    index: &mut Index,
    prepared: &Prepared,
    f: FeatureConfig,
    ratio: f64,
    records: bool,
    recover: bool,
) -> Result<(Vec<Vec<(u32, u64)>>, Vec<u64>, f64)> {
    let mut results = if records {
        vec![Vec::new(); prepared.reads.len()]
    } else {
        Vec::new()
    };
    let mut aggregate = vec![0u64; index.manifest.colors.len()];
    let mut load_seconds = 0.0;
    let mut first = 0;
    while first < index.manifest.banks.len() {
        let source = index.manifest.banks[first].source_group();
        let mut end = first + 1;
        while end < index.manifest.banks.len() && index.manifest.banks[end].source_group() == source
        {
            end += 1;
        }
        let documents = document_ids(&index.manifest.banks[first]);
        let mut remap = vec![u32::MAX; index.manifest.colors.len()];
        for (local, &global) in documents.iter().enumerate() {
            remap[global] = local as u32;
        }
        let mut hits: Vec<Vec<u32>> = prepared
            .abandoned
            .par_iter()
            .map(|origins| {
                origins
                    .iter()
                    .filter_map(|&origin| {
                        (remap[origin as usize] != u32::MAX).then_some(remap[origin as usize])
                    })
                    .collect()
            })
            .collect();
        let (order, offsets) = route_queries(&prepared.keys, index.manifest.seed, end - first);
        for id in first..end {
            let partition = id - first;
            let key_ids = &order[offsets[partition]..offsets[partition + 1]];
            if key_ids.is_empty() {
                continue;
            }
            let load = Instant::now();
            if index.manifest.version == 1 {
                let mut legacy = Index::load_selected(directory, recover, Some(id))?;
                index.banks = std::mem::take(&mut legacy.banks);
                index.abandoned = legacy.abandoned.take();
            } else {
                let mut bank = load_modern_bank(directory, &index.manifest.banks[id])?;
                for group in &mut bank.meta.groups {
                    for origin in &mut group.members {
                        *origin = remap[*origin] as usize;
                    }
                }
                index.banks.push(bank);
            }
            index.row_scratch_bytes = if index.banks[0].meta.fingerprint_first {
                0
            } else {
                index.banks[0].meta.row_bits.div_ceil(8)
            };
            load_seconds += load.elapsed().as_secs_f64();
            let current: &Index = index;
            let chunk = key_ids.len().div_ceil(rayon::current_num_threads()).max(1);
            let found: Vec<_> = key_ids
                .par_chunks(chunk)
                .map(|ids| {
                    let mut scores = RecordScores::new(documents.len());
                    let mut scratch = Vec::new();
                    ids.iter()
                        .map(|&key_id| {
                            let key = prepared.keys[key_id];
                            scores.clear();
                            current.lookup_into(key, &mut scores, &mut scratch);
                            (
                                key_id,
                                scores
                                    .touched
                                    .iter()
                                    .map(|&id| id as u32)
                                    .collect::<Vec<_>>(),
                            )
                        })
                        .collect::<Vec<_>>()
                })
                .collect();
            for chunk in found {
                for (key_id, origins) in chunk {
                    hits[key_id].extend(origins);
                }
            }
            index.banks.clear();
            index.abandoned = None;
        }
        hits.par_iter_mut().for_each(|origins| {
            origins.sort_unstable();
            origins.dedup();
        });
        if records {
            results.par_iter_mut().zip(&prepared.reads).for_each_init(
                || Scores::new(documents.len(), f),
                |scores, (out, read)| {
                    scores.score(read, &hits, f);
                    if ratio == 0.0 {
                        out.extend(
                            documents.iter().enumerate().map(|(local, &global)| {
                                (global as u32, scores.forward.values[local])
                            }),
                        );
                    } else if read.total != 0 {
                        for &local in &scores.forward.touched {
                            let count = scores.forward.values[local];
                            if count as f64 / read.total as f64 >= ratio {
                                out.push((documents[local] as u32, count));
                            }
                        }
                    }
                },
            );
        } else {
            let chunk = prepared
                .reads
                .len()
                .div_ceil(rayon::current_num_threads())
                .max(1);
            let partials: Vec<_> = prepared
                .reads
                .par_chunks(chunk)
                .map(|reads| {
                    let mut scores = Scores::new(documents.len(), f);
                    let mut sums = vec![0u64; documents.len()];
                    for read in reads {
                        scores.score(read, &hits, f);
                        for &id in &scores.forward.touched {
                            sums[id] += scores.forward.values[id];
                        }
                    }
                    sums
                })
                .collect();
            for partial in partials {
                for (&origin, count) in documents.iter().zip(partial) {
                    aggregate[origin] += count;
                }
            }
        }
        first = end;
    }
    Ok((results, aggregate, load_seconds))
}

pub(super) fn query_aggregate(
    directory: &Path,
    query: &Path,
    ratio: f64,
    recover: bool,
) -> Result<()> {
    let start = Instant::now();
    let mut index = Index::load_mode(directory, false, None, true)?;
    let mut load_seconds = start.elapsed().as_secs_f64();
    let f = index.manifest.features()?;
    if ordinary_features(&index, f) {
        return query_ordinary(
            directory,
            query,
            ratio,
            recover,
            4096,
            false,
            RecordFormat::Compact,
            index,
            f,
            start,
            load_seconds,
        );
    }
    let mut reads = Vec::new();
    if f.mode == IndexMode::Findere || f.randstrobemers.is_some() {
        let mut error = None;
        crate::parse_sequences_from_file(query, |seq| {
            match PreparedRead::new(&[], seq, &index.manifest, f) {
                Ok(read) => reads.push(read),
                Err(e) => error = Some(e),
            }
        })?;
        if let Some(error) = error {
            return Err(error);
        }
    } else {
        let mut forward = Vec::new();
        parse_features_from_file(query, index.manifest.k, index.manifest.seed, f, |key| {
            forward.push(key)
        })?;
        let total = forward.len() as u64;
        reads.push(PreparedRead {
            name: Vec::new(),
            forward,
            reverse: Vec::new(),
            runs: Vec::new(),
            total,
        });
    }
    let mut prepared = Prepared::new(reads, f);
    let load = Instant::now();
    let scans = prepared.recover(directory, &index.manifest, recover)?;
    load_seconds += load.elapsed().as_secs_f64();
    let (_, scores, loaded) = execute(directory, &mut index, &prepared, f, ratio, false, recover)?;
    load_seconds += loaded;
    let classes = index
        .manifest
        .banks
        .iter()
        .map(BankMeta::source_group)
        .collect::<std::collections::BTreeSet<_>>()
        .len();
    print_aggregate(
        &index.manifest,
        query,
        &scores,
        prepared.total,
        prepared.probes,
        prepared.probes * classes as u64,
        ratio,
        start.elapsed().as_secs_f64() - load_seconds,
    );
    report(
        &index,
        &prepared,
        scans,
        load_seconds,
        start.elapsed().as_secs_f64(),
    );
    Ok(())
}

pub(super) fn query_records(
    directory: &Path,
    query: &Path,
    ratio: f64,
    recover: bool,
    batch_size: usize,
    format: RecordFormat,
) -> Result<()> {
    let start = Instant::now();
    let mut index = Index::load_mode(directory, false, None, true)?;
    let mut load_seconds = start.elapsed().as_secs_f64();
    let f = index.manifest.features()?;
    if ordinary_features(&index, f) {
        return query_ordinary(
            directory,
            query,
            ratio,
            recover,
            batch_size,
            true,
            format,
            index,
            f,
            start,
            load_seconds,
        );
    }
    let reads = super::super::records::prepare_file(query, batch_size, |(name, seq)| {
        PreparedRead::new(name, seq, &index.manifest, f)
    })?;
    let mut prepared = Prepared::new(reads, f);
    let load = Instant::now();
    let scans = prepared.recover(directory, &index.manifest, recover)?;
    load_seconds += load.elapsed().as_secs_f64();
    let (mut results, _, loaded) =
        execute(directory, &mut index, &prepared, f, ratio, true, recover)?;
    load_seconds += loaded;
    let mut out = BufWriter::with_capacity(1024 * 1024, std::io::stdout());
    writeln!(out, "{}", format.header())?;
    for (read, hits) in prepared.reads.iter().zip(&mut results) {
        hits.sort_unstable_by_key(|hit| hit.0);
        format.write_record(
            &mut out,
            &read.name,
            read.total,
            hits.iter().map(|&(id, count)| (id as usize, count)),
            &index.manifest.colors,
        )?;
    }
    out.flush()?;
    report(
        &index,
        &prepared,
        scans,
        load_seconds,
        start.elapsed().as_secs_f64(),
    );
    Ok(())
}

fn report(index: &Index, prepared: &Prepared, scans: usize, load: f64, seconds: f64) {
    eprintln!(
        "query_memory\tmode=subindex\tmax_loaded_subindexes=1\tsubindexes={}\tresult_storage=RAM",
        index.manifest.banks.len()
    );
    eprintln!("query_preparation\tinput_passes=1\tseed_extraction_passes=1\tabandoned_file_scans={scans}\tunique_query_seeds={}", prepared.keys.len());
    eprintln!(
        "query_timings\tload_s={load:.6}\tquery_s={:.6}",
        seconds - load
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn occurrence_routing_preserves_read_ids_duplicates_and_empty_reads() {
        let reads: Vec<_> = [
            vec![83_843, 95_900, 123, 123],
            vec![],
            vec![0, u64::MAX, 95_900],
        ]
        .into_iter()
        .map(|forward| PreparedRead {
            total: forward.len() as u64,
            forward,
            name: Vec::new(),
            reverse: Vec::new(),
            runs: Vec::new(),
        })
        .collect();
        let mut expected: Vec<_> = reads
            .iter()
            .enumerate()
            .flat_map(|(id, read)| read.forward.iter().map(move |&key| (id, key)))
            .collect();
        expected.sort_unstable();
        for partitions in [1, 2, 8, 64] {
            let (routed, offsets) = route_reads(&reads, 69, partitions);
            let mut actual = routed.clone();
            actual.sort_unstable();
            assert_eq!(actual, expected);
            for p in 0..partitions {
                let keys = &routed[offsets[p]..offsets[p + 1]];
                assert!(keys.windows(2).all(|w| w[0].0 <= w[1].0));
                assert!(keys
                    .iter()
                    .all(|&(_, key)| partition_id(key, 69, partitions) == p));
            }
        }
    }

    #[test]
    fn selected_recovery_allocates_once_and_preserves_origin_unions() {
        let mut hashes: Vec<_> = (0..65_536u64)
            .map(|i| abandoned::hash(i, 69))
            .chain([0, u32::MAX])
            .collect();
        hashes.sort_unstable();
        hashes.dedup();
        let entries: Vec<_> = hashes
            .iter()
            .enumerate()
            .map(|(i, &hash)| {
                let origins = match i % 4 {
                    0 => vec![0, 2, 4],
                    1 => vec![1, 3, 5],
                    2 => vec![0, 1, 5],
                    _ => vec![5],
                };
                (hash, origins)
            })
            .collect();
        let recovery = QueryRecovery {
            entries,
            scans: 0,
            hashes: hashes.len(),
        };
        let map = recovery
            .selected(&[u32::MAX, 0, u32::MAX, 1, u32::MAX, 2])
            .unwrap()
            .unwrap();
        assert_eq!(map.directory_allocations(), 1);
        for (i, &hash) in hashes.iter().enumerate() {
            let expected: &[u32] = match i % 4 {
                0 => &[],
                1 => &[0, 1, 2],
                2 => &[0, 2],
                _ => &[2],
            };
            assert_eq!(map.get(hash), expected);
        }
        assert!(recovery.selected(&[u32::MAX; 6]).unwrap().is_none());
    }

    #[test]
    fn selected_recovery_empty_selection_has_no_directory() {
        let recovery = QueryRecovery {
            entries: Vec::new(),
            scans: 0,
            hashes: 0,
        };
        assert!(recovery.selected(&[]).unwrap().is_none());
    }

    #[test]
    fn direct_recovery_retains_colliding_query_keys_and_remaps_origins() {
        let reference = super::super::tests::union_index(123, &[0, 1]);
        let stamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let stage = Staging(std::env::temp_dir().join(format!(
            "zor-direct-recovery-{}-{stamp}",
            std::process::id()
        )));
        fs::create_dir(&stage.0).unwrap();
        let hash = abandoned::hash(83_843, 69);
        assert_eq!(hash, abandoned::hash(95_900, 69));
        write_abandoned_map(&stage.0, |emit| {
            emit(hash, &[0, 1])?;
            Ok((1, 2))
        })
        .unwrap();
        let reads = vec![PreparedRead {
            name: Vec::new(),
            forward: vec![83_843, 95_900, 123, 83_843],
            reverse: Vec::new(),
            runs: Vec::new(),
            total: 4,
        }];
        let recovered = QueryRecovery::read(&stage.0, &reference.manifest, &reads, true).unwrap();
        assert_eq!(recovered.scans, 1);
        assert_eq!(recovered.entries, vec![(hash, vec![0, 1])]);
        let map = recovered.selected(&[u32::MAX, 0]).unwrap().unwrap();
        assert_eq!(map.get(abandoned::hash(83_843, 69)), [0]);
        assert_eq!(map.get(abandoned::hash(95_900, 69)), [0]);
    }

    #[test]
    fn prepared_abandoned_matches_keep_both_distinct_keys_in_a_32_bit_collision() {
        let reference = super::super::tests::union_index(123, &[0, 1]);
        let stamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let stage = Staging(std::env::temp_dir().join(format!(
            "zor-query-collision-{}-{stamp}",
            std::process::id()
        )));
        fs::create_dir(&stage.0).unwrap();
        // Actual distinct 64-bit keys colliding under the stored 32-bit hash.
        let hash = abandoned::hash(83_843, 69);
        assert_eq!(hash, abandoned::hash(95_900, 69));
        assert_ne!(hash, abandoned::hash(123, 69));
        write_abandoned_map(&stage.0, |emit| {
            emit(hash, &[0, 1])?;
            Ok((1, 2))
        })
        .unwrap();
        let mut prepared = Prepared {
            reads: Vec::new(),
            keys: vec![83_843, 95_900, 123],
            abandoned: vec![Vec::new(); 3],
            total: 0,
            probes: 0,
        };
        assert_eq!(
            prepared
                .recover(&stage.0, &reference.manifest, true)
                .unwrap(),
            1
        );
        fs::remove_file(stage.0.join(ABANDONED_FILE)).unwrap();
        assert_eq!(prepared.abandoned, vec![vec![0, 1], vec![0, 1], vec![]]);
        assert_eq!(
            prepared
                .recover(&stage.0, &reference.manifest, false)
                .unwrap(),
            0
        );
    }
}
