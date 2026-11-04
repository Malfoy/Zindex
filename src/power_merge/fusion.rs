use super::*;

fn compatible(left: &Manifest, right: &Manifest) -> Result<()> {
    ensure!(
        left.version >= 2 && right.version >= 2,
        "fusion requires separate recovery files"
    );
    macro_rules! same {
        ($($field:ident),+ $(,)?) => {$(
            ensure!(left.$field == right.$field,
                "incompatible indexes: {} differs ({:?} versus {:?})",
                stringify!($field), left.$field, right.$field);
        )+};
    }
    same!(
        backend,
        k,
        index_mode,
        minimizer_size,
        minimizer_scheme,
        homopolymer_compression,
        seed_pattern,
        randstrobemers,
        randstrobe_scheme,
        modimizer_sampling,
        findere_z,
        hashes,
        seed,
        fingerprint_bits,
        segment_length,
        abandoned_recovery,
    );
    Ok(())
}

fn prepare_source_ids(manifest: &mut Manifest) -> Result<()> {
    if manifest.source_document_ids.is_empty() {
        manifest
            .source_document_ids
            .extend(0..manifest.colors.len());
    }
    let ids = &manifest.source_document_ids;
    let count = manifest.input_document_count.max(manifest.colors.len());
    ensure!(
        ids.len() == manifest.colors.len()
            && ids.iter().all(|&id| id < count)
            && ids.windows(2).all(|pair| pair[0] < pair[1]),
        "invalid source document IDs"
    );
    let mut seen = std::collections::HashSet::with_capacity(manifest.skipped_documents.len());
    for skipped in &manifest.skipped_documents {
        ensure!(
            skipped.source_document_id < count
                && ids.binary_search(&skipped.source_document_id).is_err()
                && seen.insert(skipped.source_document_id),
            "invalid skipped source document ID"
        );
    }
    Ok(())
}

fn validate_stream(reader: &mut impl Read, meta: &BankMeta, expected: usize) -> Result<()> {
    validate_bank_header(reader, meta, expected)?;
    let length = std::io::copy(
        &mut reader.by_ref().take(expected as u64),
        &mut std::io::sink(),
    )?;
    ensure!(length == expected as u64, "truncated bank payload");
    let mut extra = [0];
    ensure!(reader.read(&mut extra)? == 0, "unexpected bank data");
    Ok(())
}

struct CopyingReader<R, W> {
    input: R,
    output: W,
}

impl<R: Read, W: Write> Read for CopyingReader<R, W> {
    fn read(&mut self, bytes: &mut [u8]) -> std::io::Result<usize> {
        let n = self.input.read(bytes)?;
        self.output.write_all(&bytes[..n])?;
        Ok(n)
    }
}

fn copy_validated_stream(
    input: impl Read,
    output: &mut impl Write,
    meta: &BankMeta,
    expected: usize,
) -> Result<()> {
    let mut reader = BufReader::with_capacity(256 * 1024, CopyingReader { input, output });
    let mut magic = [0; 4];
    reader.read_exact(&mut magic)?;
    let mut reader = magic.as_slice().chain(reader);
    if magic == crate::ZSTD_FRAME_MAGIC {
        let mut decoder = zstd::stream::Decoder::with_buffer(&mut reader)?;
        validate_stream(&mut decoder, meta, expected)?;
    } else {
        validate_stream(&mut reader, meta, expected)?;
    }
    Ok(())
}

#[cfg(target_os = "linux")]
fn try_reflink(input: &File, output: &File) -> Result<bool> {
    use std::os::fd::AsRawFd;
    // FICLONE takes two live file descriptors and does not alter their offsets.
    let result = unsafe { libc::ioctl(output.as_raw_fd(), libc::FICLONE, input.as_raw_fd()) };
    if result == 0 {
        return Ok(true);
    }
    let error = std::io::Error::last_os_error();
    match error.raw_os_error() {
        Some(
            libc::EOPNOTSUPP
            | libc::ENOTTY
            | libc::EXDEV
            | libc::EINVAL
            | libc::ENOSYS
            | libc::EPERM
            | libc::EACCES,
        ) => Ok(false),
        _ => Err(error).context("cloning compressed subindex"),
    }
}

#[cfg(not(target_os = "linux"))]
fn try_reflink(_: &File, _: &File) -> Result<bool> {
    Ok(false)
}

fn copy_bank(path: &Path, meta: &BankMeta, destination: &Path) -> Result<()> {
    copy_bank_with(path, meta, destination, try_reflink)
}

fn copy_bank_with(
    path: &Path,
    meta: &BankMeta,
    destination: &Path,
    clone: impl FnOnce(&File, &File) -> Result<bool>,
) -> Result<()> {
    let expected = meta
        .slots
        .checked_mul(meta.row_bits.div_ceil(8))
        .context("bank size overflow")?;
    let input = File::open(path)?;
    let output = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(destination)?;
    let result = (|| {
        if clone(&input, &output)? {
            // Validate the independent filesystem snapshot, including its seek table.
            drop(output);
            bank_io::validate_table(destination, expected)?;
            validate_stream(&mut open_index_reader(destination)?, meta, expected)
        } else {
            bank_io::validate_table(path, expected)?;
            output.set_len(0)?;
            let mut writer = BufWriter::with_capacity(256 * 1024, output);
            copy_validated_stream(input, &mut writer, meta, expected)?;
            writer.flush()?;
            Ok(())
        }
    })();
    if result.is_err() {
        let _ = fs::remove_file(destination);
    }
    result
}

fn merge_recovery_many(
    inputs: &[PathBuf],
    manifests: &[Manifest],
    offsets: &[usize],
    stage: &Path,
) -> Result<(usize, usize)> {
    use std::{cmp::Reverse, collections::BinaryHeap};
    let mut readers = inputs
        .iter()
        .zip(manifests)
        .map(|(path, manifest)| AbandonedReader::new(path, manifest))
        .collect::<Result<Vec<_>>>()?;
    let mut heap = BinaryHeap::new();
    for (i, reader) in readers.iter_mut().enumerate() {
        if let Some(hash) = reader.advance()? {
            heap.push(Reverse((hash, i)));
        }
    }
    write_abandoned_map(stage, |emit| {
        let mut hashes = 0usize;
        let mut associations = 0usize;
        let mut origins = Vec::new();
        while let Some(Reverse((hash, _))) = heap.peek().copied() {
            origins.clear();
            // Equal hashes are visited in input order, preserving sorted dataset IDs.
            while let Some(Reverse((next, i))) = heap.peek().copied() {
                if next != hash {
                    break;
                }
                heap.pop();
                let offset =
                    u32::try_from(offsets[i]).context("dataset IDs must fit in 32 bits")?;
                for &id in &readers[i].origins {
                    origins.push(id.checked_add(offset).context("dataset ID overflow")?);
                }
                if let Some(next) = readers[i].advance()? {
                    heap.push(Reverse((next, i)));
                }
            }
            emit(hash, &origins)?;
            hashes = hashes
                .checked_add(1)
                .context("abandoned hash count overflow")?;
            associations = associations
                .checked_add(origins.len())
                .context("abandoned origin count overflow")?;
        }
        ensure!(hashes <= u32::MAX as usize, "too many abandoned hashes");
        Ok((hashes, associations))
    })
}

/// Pairwise compatibility wrapper for callers using the original API.
pub fn fuse_indexes(first: &Path, second: &Path, output: &Path) -> Result<()> {
    fuse_many_indexes(&[first.to_path_buf(), second.to_path_buf()], output)
}

/// Combine two or more collections without changing their compressed filters.
/// Inputs retain their order; recovery postings and document IDs are remapped once.
pub fn fuse_many_indexes(inputs: &[PathBuf], output: &Path) -> Result<()> {
    let start = Instant::now();
    ensure!(
        inputs.len() >= 2,
        "fusion requires at least two input indexes"
    );
    ensure!(
        !output.exists(),
        "output {} already exists, choose a new directory",
        output.display()
    );
    let parent = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let parent = fs::canonicalize(parent)?;
    let mut manifests = Vec::with_capacity(inputs.len());
    let mut offsets = Vec::with_capacity(inputs.len());
    let mut source_offsets = Vec::with_capacity(inputs.len());
    let mut dataset_count = 0usize;
    let mut input_document_count = 0usize;
    for input in inputs {
        ensure!(
            input.is_dir(),
            "input {} must be an index directory",
            input.display()
        );
        ensure!(
            !parent.starts_with(fs::canonicalize(input)?),
            "output must be outside the input index directories"
        );
        let mut manifest = Index::load_mode(input, false, None, true)
            .with_context(|| format!("reading input index {}", input.display()))?
            .manifest;
        if let Some(first) = manifests.first() {
            compatible(first, &manifest)?;
        }
        prepare_source_ids(&mut manifest)?;
        offsets.push(dataset_count);
        source_offsets.push(input_document_count);
        dataset_count = dataset_count
            .checked_add(manifest.colors.len())
            .context("dataset count overflow")?;
        input_document_count = input_document_count
            .checked_add(manifest.input_document_count.max(manifest.colors.len()))
            .context("source document count overflow")?;
        ensure!(
            dataset_count <= u32::MAX as usize,
            "dataset IDs must fit in 32 bits"
        );
        manifests.push(manifest);
    }
    let stamp = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
    let stage = Staging(parent.join(format!(".zor-fusing-{}-{stamp}", std::process::id())));
    fs::create_dir(&stage.0)?;
    let mut sources = Vec::new();
    let mut banks = Vec::new();
    let mut component = 0usize;
    for ((directory, manifest), &offset) in inputs.iter().zip(&manifests).zip(&offsets) {
        let mut previous = None;
        for original in &manifest.banks {
            let mut meta = original.clone();
            if previous != Some(meta.component) {
                component = component
                    .checked_add(1)
                    .context("component count overflow")?;
                previous = Some(meta.component);
            }
            sources.push(directory.join(&meta.file));
            meta.component = component;
            meta.file = meta.expected_filename();
            for group in &mut meta.groups {
                for id in &mut group.members {
                    *id = id.checked_add(offset).context("dataset ID overflow")?;
                }
            }
            banks.push(meta);
        }
    }
    // Every compressed subindex can be copied and validated independently.
    // Recovery streams are merged once, not through repeated pairwise passes.
    let (copied, recovery) = rayon::join(
        || {
            sources
                .par_iter()
                .zip(&banks)
                .try_for_each(|(source, meta)| {
                    copy_bank(source, meta, &stage.0.join(&meta.file))
                        .with_context(|| format!("copying subindex {}", source.display()))
                })
        },
        || {
            if manifests[0].abandoned_recovery {
                merge_recovery_many(inputs, &manifests, &offsets, &stage.0)
            } else {
                Ok((0, 0))
            }
        },
    );
    copied?;
    let recovery = recovery?;
    let mut iter = manifests.into_iter().zip(source_offsets);
    let (mut left, _) = iter.next().expect("at least two manifests");
    for (mut right, source_offset) in iter {
        for id in &mut right.source_document_ids {
            *id = id
                .checked_add(source_offset)
                .context("source document ID overflow")?;
        }
        for skipped in &mut right.skipped_documents {
            skipped.source_document_id = skipped
                .source_document_id
                .checked_add(source_offset)
                .context("source document ID overflow")?;
        }
        left.source_document_ids.extend(right.source_document_ids);
        left.colors.extend(right.colors);
        left.skipped_documents.extend(right.skipped_documents);
        macro_rules! sum {
            ($($field:ident),+ $(,)?) => {$(
                left.$field = left.$field.checked_add(right.$field).context(concat!(stringify!($field), " overflow"))?;
            )+};
        }
        sum!(
            planning_pair_evaluations,
            planning_possible_pair_evaluations,
            promoted_groups,
            singleton_logical_bits,
            predicted_logical_bits,
            planned_groups,
            rejected_groups,
            overflow_groups,
            exact_group_verifications,
            memory_guard_split_groups
        );
        left.sketch_seconds += right.sketch_seconds;
        left.plan_seconds += right.plan_seconds;
        left.max_subindex_size = left
            .max_subindex_size
            .zip(right.max_subindex_size)
            .map(|(a, b)| a.max(b));
    }
    left.banks = banks;
    left.version = 10;
    left.abandoned_hashes = recovery.0;
    left.abandoned_origins = recovery.1;
    left.input_document_count = input_document_count;
    left.max_subindexes = None;
    left.min_subindex_documents = 1;
    left.merge_memory_budget_bytes = None;
    left.construction_cell_storage = "copied-compressed-banks".into();
    left.column_order_method = None;
    left.construction_abandoned_storage = Some(
        if left.abandoned_recovery {
            "streaming-merge"
        } else {
            "omitted"
        }
        .into(),
    );
    left.build_seconds = start.elapsed().as_secs_f64();
    let mut writer = BufWriter::new(File::create(stage.0.join("manifest.json"))?);
    serde_json::to_writer_pretty(&mut writer, &left)?;
    writer.flush()?;
    drop(writer);
    if !left.skipped_documents.is_empty() {
        #[derive(Serialize)]
        struct Report<'a> {
            status: &'static str,
            input_documents: usize,
            indexed_documents: usize,
            skipped_document_count: usize,
            skipped_documents: &'a [SkippedDocument],
        }
        let report = Report {
            status: "fused",
            input_documents: left.input_document_count,
            indexed_documents: left.colors.len(),
            skipped_document_count: left.skipped_documents.len(),
            skipped_documents: &left.skipped_documents,
        };
        let mut writer = BufWriter::new(File::create(stage.0.join("skipped-documents.json"))?);
        serde_json::to_writer_pretty(&mut writer, &report)?;
        writer.flush()?;
    }
    let datasets = left.colors.len();
    let subindexes = left.banks.len();
    let status = left.collection_status();
    let skipped = left.skipped_documents.len();
    drop(left);
    // Reuse query validation before publishing the new manifest and bank mapping.
    Index::load_mode(&stage.0, false, None, true)?;
    ensure!(
        !output.exists(),
        "output appeared during fusion: {}",
        output.display()
    );
    fs::rename(&stage.0, output)?;
    eprintln!("collection\tstatus={status}\trequested={input_document_count}\tindexed={datasets}\tskipped={skipped}");
    eprintln!(
        "fused {} indexes into {} in {:.3}s: datasets={}, subindexes={}, worker_threads={}",
        inputs.len(),
        output.display(),
        start.elapsed().as_secs_f64(),
        datasets,
        subindexes,
        rayon::current_num_threads()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;

    fn stage() -> Staging {
        let stamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let stage = Staging(
            std::env::temp_dir().join(format!("zor-merge-bank-{}-{stamp}", std::process::id())),
        );
        fs::create_dir(&stage.0).unwrap();
        stage
    }

    fn meta(slots: usize) -> BankMeta {
        BankMeta {
            component: 0,
            file: format!("bank-{}.zorbank.zst", slots.ilog2()),
            slots,
            row_bits: 8,
            groups: vec![GroupMeta {
                members: vec![0],
                width: 8,
                offset: 0,
                union_cardinality: 0,
                fallback_keys: 0,
                abandoned_keys: Some(0),
            }],
            fingerprint_first: true,
            routing: None,
        }
    }

    struct ShortReads<'a> {
        bytes: &'a [u8],
        consumed: &'a Cell<usize>,
    }

    impl Read for ShortReads<'_> {
        fn read(&mut self, output: &mut [u8]) -> std::io::Result<usize> {
            let n = output.len().min(self.bytes.len()).min(3);
            output[..n].copy_from_slice(&self.bytes[..n]);
            self.bytes = &self.bytes[n..];
            self.consumed.set(self.consumed.get() + n);
            Ok(n)
        }
    }

    #[test]
    fn streaming_copy_consumes_each_byte_once_for_raw_and_concatenated_zstd_with_short_reads() {
        let meta = meta(256);
        let mut raw = Vec::new();
        raw.extend_from_slice(PLANE_MAGIC);
        for n in [meta.slots, meta.row_bits, meta.groups.len(), meta.slots] {
            raw.extend_from_slice(&(n as u64).to_le_bytes());
        }
        raw.extend((0..meta.slots).map(|n| n as u8));
        let mut encoded = zstd::stream::encode_all(&raw[..17], 0).unwrap();
        encoded.extend(zstd::stream::encode_all(&raw[17..], 0).unwrap());
        for bytes in [&raw, &encoded] {
            let consumed = Cell::new(0);
            let input = ShortReads {
                bytes,
                consumed: &consumed,
            };
            let mut copied = Vec::new();
            copy_validated_stream(input, &mut copied, &meta, meta.slots).unwrap();
            assert_eq!(consumed.get(), bytes.len());
            assert_eq!(&copied, bytes);
        }
    }

    #[test]
    fn streaming_copy_propagates_write_failures() {
        struct Broken;
        impl Write for Broken {
            fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
                Err(std::io::Error::other("injected write failure"))
            }
            fn flush(&mut self) -> std::io::Result<()> {
                Ok(())
            }
        }
        let meta = meta(256);
        let error =
            copy_validated_stream(&b"anything"[..], &mut Broken, &meta, meta.slots).unwrap_err();
        assert!(error.to_string().contains("injected write failure"));
    }

    #[test]
    fn copying_multiframe_banks_validates_payload_and_seek_table_without_reencoding() {
        let stage = stage();
        let meta = meta(1 << 23);
        let path = stage.0.join(&meta.file);
        let mut writer = bank_io::Writer::new(&path).unwrap();
        writer.write_all(PLANE_MAGIC).unwrap();
        for n in [meta.slots, meta.row_bits, meta.groups.len(), meta.slots] {
            writer.write_all(&(n as u64).to_le_bytes()).unwrap();
        }
        for n in 0..meta.slots / 65536 {
            // Force compressed input to span several copy buffers as well as Zstd frames.
            let block: Vec<_> = (0..65536)
                .map(|i| {
                    if n < 8 {
                        splitmix64((n * 65536 + i) as u64) as u8
                    } else {
                        n as u8
                    }
                })
                .collect();
            writer.write_all(&block).unwrap();
        }
        writer.finish().unwrap();
        let destination = stage.0.join("copy.zorbank.zst");
        copy_bank(&path, &meta, &destination).unwrap();
        let mut bytes = fs::read(&path).unwrap();
        assert_eq!(bytes, fs::read(destination).unwrap());
        let fallback = stage.0.join("fallback.zorbank.zst");
        copy_bank_with(&path, &meta, &fallback, |_, _| Ok(false)).unwrap();
        assert_eq!(bytes, fs::read(&fallback).unwrap());
        // Exercise snapshot validation even on filesystems without reflink support.
        let snapshot = stage.0.join("snapshot.zorbank.zst");
        copy_bank_with(&path, &meta, &snapshot, |mut input, mut output| {
            std::io::copy(&mut input, &mut output)?;
            Ok(true)
        })
        .unwrap();
        assert_eq!(bytes, fs::read(&snapshot).unwrap());
        let unchanged = fs::read(&fallback).unwrap();
        assert!(copy_bank(&path, &meta, &fallback).is_err());
        assert_eq!(unchanged, fs::read(&fallback).unwrap());
        let bad_snapshot = stage.0.join("bad-snapshot.zorbank.zst");
        assert!(
            copy_bank_with(&path, &meta, &bad_snapshot, |_, mut output| {
                output.write_all(b"invalid snapshot")?;
                Ok(true)
            })
            .is_err()
        );
        assert!(!bad_snapshot.exists());
        let source_bytes = bytes.clone();
        // Damage a frame checksum while keeping the independently checksummed seek table valid.
        let table_len = u64::from_le_bytes(bytes[bytes.len() - 8..].try_into().unwrap()) as usize;
        let table_start = bytes.len() - table_len;
        let frame_len = u64::from_le_bytes(
            bytes[table_start + 24..table_start + 32]
                .try_into()
                .unwrap(),
        ) as usize;
        bytes[frame_len - 1] ^= 1;
        fs::write(&path, &bytes).unwrap();
        let failed = stage.0.join("bad-frame.zorbank.zst");
        assert!(copy_bank_with(&path, &meta, &failed, |_, _| Ok(false)).is_err());
        assert!(!failed.exists());
        assert!(copy_bank(&path, &meta, &failed).is_err());
        assert!(!failed.exists());
        bytes = source_bytes;
        let checksum = bytes.len() - 24;
        bytes[checksum] ^= 1;
        fs::write(&path, bytes).unwrap();
        let missing = stage.0.join("failed.zorbank.zst");
        assert!(copy_bank(&path, &meta, &missing)
            .unwrap_err()
            .to_string()
            .contains("checksum mismatch"));
        assert!(!missing.exists());
    }
}
