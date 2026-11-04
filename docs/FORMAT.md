# Index format notes

## Zindex power-of-two merged directory

Zindex stores each index in a directory containing a manifest, compressed
subindexes and optional recovery data.
`manifest.json` contains version `8` (or `9` for hash-partitioned indexes and `10`
for combined indexes), feature/hash settings, dataset identities,
group members and widths, bank sizes and abandoned-key map metadata.
`input_document_count` records the original requested collection size,
`source_document_ids` maps indexed datasets to their original input positions,
and `skipped_documents` lists explicit memory exclusions. `info` derives
`collection_status` from these fields. An index with exclusions is `incomplete`.
Builds require all datasets unless skipping is explicitly enabled. Successful
skip reports use status `incomplete`, and failed all-skipped builds publish no index.
Version 4 adds `minimizer_scheme` (`simd` or `open-closed-mod`) and
`homopolymer_compression` (boolean). These govern both construction and queries.
Missing legacy fields default to SIMD selection and compression disabled.
Version 5 adds `seed_pattern` (a symmetric `#`/`_` string or null) and
`randstrobemers` (null or an object with `min_distance` and `max_distance`).
Missing legacy fields default to null. Distances are inclusive, in sequence
bases between strobe starts, after optional homopolymer compression. The bank
and abandoned-key layouts are unchanged.

Version 6 adds `randstrobe_scheme`, currently
`strobealign-popcount8-symmetric64-downstream-v1`. It uses one downstream
seed per anchor, canonical-strobe XXH64, masked-XOR popcount with q=255,
first-minimum tie breaking, self-pair fallback and the paper's symmetric
`h1/2+h2/2` full-seed key. Query orientations are scored independently and
combined by per-document maximum per record. Missing or unknown extraction
identities on randstrobemer indexes require rebuilding. Bank payloads and the
abandoned-key membership semantics are unchanged. Version 7 changes the lossless
recovery-file encoding, not bank payloads or feature extraction.
Version 8 adds `segment_length`: one fixed power-of-two segment size for every
filter/subindex, default 64 cells, independent of filter capacity. Array sizes
remain powers of two and must fit at least `hashes` complete segments. Queries
use this stored size. Versions 1 through 7 retain their historical size-dependent
position mapping. Changing a version number alone does not convert an index.

Version 9 adds hash-routed subindexes. `max_subindex_size` stores the build's
maximum uncompressed slice-payload bytes (default 8 GiB. A power of two).
Oversized finalized size classes are reconstructed as a power-of-two number
of partitions. Their bank metadata contains a `routing` object with
`source_slots`, `partition` (zero-based), and `partitions`. Each has
`slots * partitions == source_slots`. Both slot counts and the partition count
are powers of two. Unpartitioned banks omit `routing`. Partitions are ordered
by `(source_slots, partition)` and stored in
`bank-<source_exponent>-part-<partition>-of-<partitions>.zorbank.zst`.

The route is `splitmix64(key ^ seed ^ 0x34B61D8F29A7C053) & (partitions - 1)`.
Every key's H positions are computed using the selected bank's **local** slot
count. Keys are not replicated or queried in other partitions of the same size
class. Fingerprints and origin bits are unchanged. Each original document and
merged group appears once per partition of its size class. Cross-class duplicate
documents, missing partitions, inconsistent members and invalid geometries are
rejected. Group order may differ between partitions to improve compression.
Abandoned keys from every partition still enter the single global 32-bit-hash
multimap. Collisions retain the full union of original document IDs.

`--max-subindexes` and `min_subindex_documents` apply to original size classes
before hash partitioning. A size class can consequently contain multiple physical
bank files. Each bank's byte-padded slice payload fits `max_subindex_size`.
its 40-byte header, metadata and query-result memory are outside that cap.
The bank payload encoding itself is unchanged. Version-9 readers accept the
optional fixed segment geometry. Legacy library configurations without it keep
their historical position mapping. Readers before version 9 reject partitioned
directories, rather than accidentally querying every partition.

Version 10 combines existing index collections without rebuilding filters.
Each bank has a `component` ID identifying its original collection. Components
are flattened and assigned consecutive IDs starting at 1 during fusion, including
when an input is itself fused. Banks are ordered by
`(component, source_slots, partition)` and named
`component-<id>-bank-<exponent>[-part-<partition>-of-<partitions>].zorbank.zst`.
Within a component, the version-9 routing and dataset-mapping invariants apply.
Queries route independently within each `(component, source_slots)` group, so
collections with different partition geometries can coexist. Existing indexes
have implicit component 0 and retain their filenames and ordering.

Compressed bank bytes and their headers are copied unchanged. Linux fusion operations
attempt an independent copy-on-write reflink and validate the resulting
snapshot. If cloning is unavailable, one streaming pass validates decoded
payloads while copying the original compressed bytes. Both paths validate
seek-table geometry and checksums, bank headers, decoded lengths and payload
checksums when present. A failed copy removes its partial destination.
Independent banks are copied and validated in parallel using the requested worker
pool, alongside recovery merging. Each recovery input has a fixed-size buffer.
The combined manifest offsets each input's dataset IDs by the indexed dataset
count of all preceding inputs. A heap merges recovery files in hash order in one
pass, concatenating postings after the same ID remapping, including for colliding
hashes. Source document IDs and skipped-input reports use the original document
count of all preceding inputs as their offset. The combined index is self-contained. Readers before
version 10 reject it. Query and info validation reject invalid component order,
filenames and dataset assignments.

Each `bank-<exponent>.zorbank.zst` expands to magic `ZORPOW3\0`, followed by
four little-endian `u64` values: slots, bits per row, group count and data bytes.
Rows have `ceil(bits_per_row / 8)` bytes. Each row begins with all groups'
1-, 2-, 4-, 8- or 16-bit fingerprints, followed immediately by their packed origin bits. Singleton
groups have no origin bits. Group offsets in the manifest describe the original
complete-payload prefix. The origin-plane offset for group `i` is
`group_count * fingerprint_bits + group.offset - i * fingerprint_bits`.

`column_order_method` identifies construction ordering independently of the
`column_order` label. Current full/streaming builds record
`union-sketch-fingerprints-lexicographic-v1`: keep equal-width sections together,
sort by the existing final union-sketch fingerprint bytes, then break ties by
the first dataset ID. Only these 8-bit sketch bytes remain after planning.
The 32-bit union minima are freed. Rebuilt filters and memory/hash partitions
share their planned parent's key without new sketches or filter-cell reads.
The keys are construction-only and are not stored in the index. Naive builds
record `width-then-dataset-id-v1`. Fusion preserves input column layouts and
omits a single construction-order method for the combined collection.

Construction mode does not change the bank encoding or query layout. The default
two-pass builder records `construction_cell_storage` as
`two-pass-temporary-disk-packed-tiled-packing` (verified cells staged in files,
bit-packed at the group's payload width until bank compression). `--low-disk`
records `two-pass-memory-packed-no-scratch` (the same packed cells retained in
RAM). `--disk-backed` records `temporary-disk-u64-tiled-packing`. The two-pass
publication directory contains only final index files. First-pass projected
capacities may differ from exact final capacities. This does not change the
planned dataset membership partition.

New builds also record `construction_abandoned_storage`: `memory`,
`temporary-disk-sorted-runs`, or `omitted`. Missing legacy metadata means the
older in-RAM construction behavior. Both disk construction modes now stage
abandoned-key records as sorted temporary runs and use a bounded external merge.
the published `ZORABN2` file, 32-bit hashes, origin unions and query behavior are
unchanged. Temporary run records contain little-endian `(u32 hash, u32 group ID,
u64 origin mask)`, and are never part of the published index.

The optional manifest metadata `min_subindex_documents` records the build's
minimum original-document count per bank (new default 128, missing legacy
field defaults to 1). `promoted_groups` records filters rebuilt at a larger
power-of-two capacity to meet that constraint (missing legacy field defaults
to 0). Original documents are counted through `group.members`, not by the
number of merged groups. With fewer total documents than the requested minimum,
all groups occupy one bank. The planner selects a minimum-padded-bit-cost
contiguous partition of the exact slot-size buckets. Only promoted filters
require additional input reads and reconstruction. This changes neither the
bank encoding nor document partitions and creates no temporary key/cell files.
Fingerprint width and abandoned-key union semantics stay fixed, but changing
capacity may change which absent keys are false positives. Setting the minimum
to 1 preserves the old layout. Queries fuse row XOR and fingerprint comparison
in registers where supported. This is an execution change, not a format change.

`abandoned.zormap.zst` expands to magic `ZORABN2\0`, followed by records:
the delta from the previous 32-bit hash (initially zero), the origin count,
then the ascending origin IDs delta-coded from zero within each record.
Integers use unsigned base-128 variable-length encoding (at most five bytes
per u32). No hash or origin bits are discarded. The reader also accepts the
legacy `ZORABN1\0` stream: little-endian u32 hash, count and origin IDs.
This is one global multimap across all banks.
32-bit hash collisions union origin lists. It is omitted with `--drop-abandoned`.

The current reader also supports manifest versions 1 through 8 and their legacy
`ZORPOW1\0` interleaved banks. Version-1 per-group abandoned-key records are
converted to the global 32-bit-hash multimap on load. Older readers cannot read
new version-9 or version-10 directories. To add datasets, combine compatible
indexes with `zorindex fusion` or rebuild from the complete input list.

### Parallel loading of merged subindexes

New bank files concatenate independent Zstd frames of at most 4 MiB decoded
payload. Concatenating their decoded bytes gives the unchanged bank header and
matrix. The same compression level is used, and each frame has a Zstd content
checksum. Single-frame and raw legacy files remain supported.

An optional final Zstd skippable frame (magic `0x184D2A50`) indexes frame ranges.
Its payload is `ZORFRM1\0`, a little-endian u64 frame count, then four little-endian
u64 values per frame: compressed offset, compressed length, decoded offset,
decoded length. These are followed by the 64-bit FNV-1a checksum of the preceding
payload, `ZORSEEK1`, and the total payload length as a little-endian u64.
Offsets must cover the data contiguously, decoded coverage must match the bank
header, and the reader bounds frame and table allocations before decoding.

Streaming readers ignore the skippable table and concatenate the frames normally.
Parallel readers use disjoint output slices and bound compressed scratch to a
single frame per active worker. This compression framing does not change matrix
geometry, fingerprints, document IDs or membership semantics.
