# ZORINDEX

ZORINDEX contains two Rust command-line indexes for approximate colored
k-mer queries over many FASTA/FASTQ datasets:

- `zorindex`: Zindex, merged colored power-of-two ZOR banks by default;
  `--plain` selects the legacy interleaved per-color stack.
- `bloomindex`: Bindex, a matching interleaved stack of per-color Bloom filters.

Both binaries use the same input format, DNA feature extraction code, canonical
2-bit encoding, parallel parsing, and query reporting. This makes Bindex a
direct comparison target for Zindex when evaluating size, build time, and query
throughput.

The source tree is organized as:

```text
src/lib.rs              shared DNA parsing, feature extraction, hashing, FOF parsing
src/bin/zorindex.rs     Zindex CLI and pure ZOR stack implementation
src/bin/bloomindex.rs   Bindex CLI and interleaved Bloom stack implementation
src/power_merge/        Zindex merge planning, colored power-of-two banks and lookup
```

## Requirements

- Rust stable toolchain
- Enough RAM to hold the parsed feature sets during construction

## Build

```bash
cargo build --release
```

The release binaries are:

```text
target/release/zorindex
target/release/bloomindex
```

Both tools default to all available hardware threads. Override this with
`--threads`.

For copy-paste collaborator commands, see
[docs/COLLABORATOR.txt](docs/COLLABORATOR.txt). Build on the machine where you
will run the binaries: `.cargo/config.toml` enables `target-cpu=native`.

## Input Format

Both binaries take a file-of-files through `--fof`. Empty lines and lines
starting with `#` are ignored.

Single-column form:

```text
path/to/sample.fa.gz
path/to/other_sample.fna.zst
```

Named-color form:

```text
sample_a<TAB>path/to/sample_a.fa
sample_b<TAB>path/to/sample_b.fastq.gz
```

Whitespace-separated `name path` lines are also accepted. If no name is given,
the color name is derived from the file name.

Supported sequence input:

- FASTA
- FASTQ
- gzip-compressed FASTA/FASTQ
- zstd-compressed FASTA/FASTQ (`.zst`, `.zstd`)
- xz-compressed FASTA/FASTQ (`.xz`)
- bzip2-compressed FASTA/FASTQ (`.bz2`, `.bz`, `.bzip2`)

## Feature Modes

Both indexes support the same indexed feature modes:

```text
kmers       every canonical k-mer
minimizers  canonical minimizers from each k-window
modimizers  hash-sampled canonical k-mers
findere     all (k-z)-mers, with positional reconstruction of query k-mers
```

Full k-mer and modimizer modes currently require `k <= 31`, because canonical
DNA k-mers are encoded in a `u64`. Minimizer mode may use larger odd-length
windows because only the minimizer sequence is encoded; the default minimizer
length is 21. Both canonical minimizer schemes require odd `k`.

For Zindex, `--index-mode minimizers` defaults to **open-closed mod minimizers**
([Groot Koerkamp, Liu and Pibiri, 2025](https://doi.org/10.1186/s13015-025-00270-0)).
Use `--minimizer-scheme simd` to retain the standard canonical SIMD hash-minimizer
backend. Bindex still uses that SIMD backend. The overall feature-mode default
remains `kmers`, not `minimizers`.

The open-closed selector automatically uses AVX2 on supported x86-64 CPUs;
short records and other CPUs use its portable path. Four overlapping regions
are processed together on one thread with prefix/suffix sliding minima,
full-width 64-bit hashes, and bounded 1,024-window tile buffers. Exact selection,
strand-consistent ties and emitted order match the original scalar iterator;
the SIMD hash-minimizer scheme is still independently selectable.

The open-closed selector runs in linear time, with `w=k-m+1`,
`r=min(4,m)` and anchor length `t=r+(m-r)%w`, where `m` is `--minimizer-size`.
Open-syncmer anchors are preferred over closed-syncmer anchors, with ordinary
anchors as fallback; the selected representative offset is `anchor_offset%w`.
Our canonical adaptation hashes canonical anchors and inner s-mers, chooses the
largest `s<=min(4,t)` with `t-s` even, treats tied inner minima symmetrically, and
breaks outer ties using the odd window's GT/AC majority. This guarantees the
same canonical representatives on reverse-complemented sequences and matching
selection when an exact window is queried in isolation. It is not a byte-for-byte
implementation of the paper's directional tie rules.

Examples:

```bash
zorindex build --fof samples.fof --k 31 --output samples.zindex
zorindex build --fof samples.fof --k 63 --index-mode minimizers \
  --minimizer-size 21 --output samples.m21.zindex
zorindex build --fof samples.fof --k 63 --index-mode minimizers \
  --minimizer-size 21 --minimizer-scheme simd --output samples.simd.zindex
zorindex build --fof samples.fof --k 31 --index-mode modimizers \
  --modimizer-sampling 16 --output samples.mod16.zindex
```

The feature mode is stored in the index and reused by `query`. Legacy `--plain`
indexes also support `append`; merged directories currently require a rebuild.

### Optional homopolymer compression (Zindex)

Add `--homopolymer-compression` (alias `--homocompress`) to `build` to collapse
consecutive identical DNA bases before feature extraction, for any feature mode:

```bash
zorindex build --fof samples.fof --output samples.hpc.zindex \
  --k 49 --index-mode minimizers --minimizer-size 31 --homopolymer-compression
zorindex query --index samples.hpc.zindex --query reads.fa.gz
```

Compression is case-insensitive and never joins records or ambiguous-base runs.
The stored setting is reused automatically by queries and plain-index appends;
there is no separate query flag to accidentally mismatch it. `k` and minimizer
lengths apply to the compressed sequence, and query ratios count its features
(or its k-mer windows in findere mode), not original-base positions. Raw parsed
base counts are still reported. Compression is disabled by default.

New Zindex files/directories record both the minimizer scheme and compression
setting. Legacy indexes default to SIMD selection and no homopolymer compression.

### Optional spaced seeds and randstrobemers (Zindex)

For `kmers` or `minimizers`, `--seed-pattern` (alias `--spaced-seed`) retains
the bases marked `#` and discards those marked `_`. The pattern must contain
at least one `#`, be symmetric under reversal, and have width exactly `k` for
k-mers or `--minimizer-size` (`m`) for minimizers. **Ordinary minimizers are
selected first, then masked**; the pattern does not affect their selection.
The retained bases are packed and canonicalized, with optional BMI2 bit
extraction on supported x86-64 CPUs. An all-`#` pattern preserves ordinary keys.

```bash
zorindex build --fof samples.fof --output samples.spaced.zindex \
  --k 8 --seed-pattern '##_##_##'
zorindex build --fof samples.fof --output samples.spaced-min.zindex \
  --k 31 --index-mode minimizers --minimizer-size 9 --seed-pattern '##__#__##'
```

`--randstrobemers` (alias `--randstrobes`) links **two** strobes: all k-mer
positions in k-mer mode, or only already-selected minimizer positions in
minimizer mode. `k` or `m` specifies the length of **each strobe**, not the
combined key. `--strobe-min` and `--strobe-max` set inclusive distances between
their **start positions**, measured in bases, not in numbers of minimizers.
Their defaults are one and three strobe lengths, respectively. Minimum distance
must be positive and cannot exceed maximum distance; overlapping strobes are
allowed if explicitly requested by choosing a minimum below the strobe length.

```bash
zorindex build --fof samples.fof --output samples.rand.zindex \
  --k 21 --randstrobemers --strobe-min 25 --strobe-max 100
zorindex build --fof samples.fof --output samples.rand-min.zindex \
  --k 49 --index-mode minimizers --minimizer-size 31 \
  --randstrobemers --strobe-min 32 --strobe-max 100
zorindex query --index samples.rand-min.zindex --query reads.fa.gz
```

Construction emits **one downstream seed per selected anchor**, using
Strobealign's canonical-strobe XXH64 hashes and minimum masked-XOR popcount
selection (`q=255`, first/nearest candidate wins ties). An anchor without an
eligible downstream partner is paired with itself, as in Strobealign's
iterator. Full-seed keys use the symmetric `h1/2 + h2/2` representation from
[Sahlin (2022)](https://doi.org/10.1186/s13059-022-02831-7), not Strobealign's
newer partial-match/MCS-specific packed keys. Your all-k-mer or selected-
minimizer anchor set and inclusive base-distance settings remain unchanged;
this does not switch minimizer selection to syncmers.

Queries independently reselect partners in the forward and reverse-complement
directions and report the **larger count per document**, never the sum. Both
directions have one seed per selected anchor and the same denominator. For
multi-record queries the per-document maxima are summed **per record**, so
mixed-strand read files work correctly. Symmetric seed keys do not make partner
selection strand-invariant. Duplicate seed occurrences remain weighted. At
read boundaries, partner reselection or self-pair fallback may differ from a
long reference even without errors. The inter-strobe gap is not part of a key.

The earlier custom two-direction construction is incompatible: its indexes
must be rebuilt. New randstrobemer indexes carry a versioned extraction
identity, and old custom indexes are rejected rather than silently misqueried.

Both options may be combined: ordinary minimizer selection precedes masking,
and randstrobe linking then uses the masked strobes. Homopolymer compression,
if enabled, happens first, so all lengths and distances refer to compressed
bases. Neither option crosses record or ambiguous-base boundaries. The linker
buffers only the local distance window; minimizer mode additionally stores and
deduplicates selected positions. Settings are persisted and reused by queries
and plain-index appends, and `info` reports them. They are not supported with
findere or modimizers. Bindex does not expose these new options.

### Findere

Both tools accept `--index-mode findere`. Here `--k` is the **query** k-mer
length; the index stores all canonical `(k-z)`-mers from valid reference
stretches of at least `k` bases. `--findere-z` (alias `--z`) defaults to **10**.
For example, `--k 31` indexes 21-mers and requires 11 consecutive 21-mer hits
in the **same dataset** to report one 31-mer match.

```bash
zorindex build --fof samples.fof --k 31 --index-mode findere --output samples.findere.zindex
bloomindex build --fof samples.fof --k 31 --index-mode findere --output samples.findere.blmidx

# Override z; query automatically reuses the stored setting.
zorindex build --fof samples.fof --k 31 --index-mode findere --findere-z 3 --output samples.z3.zindex
zorindex query --index samples.findere.zindex --query reads.fa
bloomindex query --index samples.findere.blmidx --query reads.fa
```

Constraints are `0 <= z < k <= 255` and `1 <= k-z <= 31`. Setting `z=0`
recovers ordinary k-mer membership. Findere is a separate mode, not combined
with minimizer or modimizer sampling. `info` reports `findere_z` and `indexed_k`.

Query `matches`, `total_features`, and ratios count reconstructed **k-mer
positions**, not shorter-feature hits. Probe counters count actual shorter-feature
probes. Consecutive windows reuse those probes; reconstruction resets at each
FASTA/FASTQ record or ambiguous base and ignores stretches shorter than `k`.
Findere queries use the Rayon pool configured by `--threads` (all available
hardware threads by default). Bounded batches distribute both short records and
individual long sequences over workers, each with private decoding scratch and
score counters. Long-sequence tasks own disjoint k-mer windows and include `k-1`
trailing context bases: no matches are lost or duplicated at task boundaries.
Those boundaries repeat up to `z` shorter-feature probes; probe counters include
this work and are deterministic across thread counts. Parsing/decompression is
still serial. The additional findere negative-run skipping optimization is not
implemented.

Findere reduces some approximate-membership false positives but is not exact:
all constituent shorter features can exist without their full k-mer occurring.
Bloom's reported `estimated_fp_rate` remains the underlying shorter-feature
estimate, not the reconstructed k-mer false-positive rate. Pure ZOR can abandon
constraints; findere does not repair those false negatives and can propagate
one missing shorter feature to several query windows. Check `abandoned` when
choosing ZOR construction settings.

Findere uses feature-mode
tag 3 followed by a mode-specific `z` byte after the common feature settings;
older binaries reject this unknown mode. The stored mode and `z` survive append
and ZOR repack operations.

## Zindex

Zindex builds merged colored power-of-two banks by default, using streaming
planning, compression column ordering, sizing factor **1.0**, `most-deg2` cycle
breaking and `tie-scan=8`. Four hash locations are the default. These settings
do not need to be repeated on the command line:

```bash
zorindex build --fof samples.fof --output samples.zindex \
  --k 49 --index-mode minimizers --minimizer-size 31 --fingerprint-bits 16
zorindex query --index samples.zindex --query reads.fa.gz
zorindex info --index samples.zindex
```

`--output` must name a new directory, regardless of its extension. To use the
legacy single-stack layout instead, pass **`--plain`**. This builds one pure ZOR
filter per color and stores them in a slot-major interleaved layout:

```text
slot 0: color 0 fingerprint, color 1 fingerprint, ...
slot 1: color 0 fingerprint, color 1 fingerprint, ...
```

A query computes the ZOR positions once, XORs the selected rows, and compares
each color lane against the query fingerprint. The layout is intended to keep
the later SIMD/chunked XOR path straightforward.

### Build

```bash
zorindex build \
  --plain \
  --fof samples.fof \
  --k 31 \
  --fingerprint-bits 16 \
  --output samples.zoridx
```

Important options:

- `--hashes`: number of ZOR locations per key.
- `--fingerprint-bits`: `1`, `8` (default), or `16`. The 1-bit option is bit-packed and has approximately 50% false positives per absent indexed feature; sequence-level thresholds still apply.
- `--slot-scale`: multiplier applied to the largest color cardinality.
- `--slot-count`: explicit slot budget before segmented layout rounding.
- `--segment-length`: explicit power-of-two segment length.
- `--cycle-break`: heuristic used when the peel reaches a core.
- `--union-graph`: optional union ZOR filter used as a query prefilter.

### Append

```bash
zorindex append --index samples.zoridx --fof more_samples.fof \
  --output samples.extended.zoridx
```

Append reuses the original `k`, hash count, seed, layout, fingerprint width,
and feature mode. If the input index contains a union graph, append removes it
because the stored index does not contain all previous keys needed to update
the union graph exactly.

### Query

```bash
zorindex query --index samples.zoridx --query query.fa.gz --min-ratio 0.1
```

Output starts with a summary line, followed by one line per color passing the
ratio threshold:

```text
#query  query.fa.gz  k=31  index_mode=kmers  total_features=...  stack_probes=...  elapsed_s=...
color_name  hits  total_features  hit_ratio
```

Merged indexes also support one result per FASTA record:

```bash
zorindex query --index samples.pow2 --query reads.fa.gz --records --min-ratio 0.8
```

This emits `query_id`, `valid_kmers` (the indexed feature count), and ordered
numeric `document_id:count` pairs, matching the record-level benchmark TSV.
Lookup, thresholding and formatting run in parallel with persistent worker
counters; only matching documents are scanned when `--min-ratio` is positive.
At threshold zero, all documents, including zero counts, are still reported.
`--record-batch-size` defaults to 4096. The compatibility aliases
`--bench-records` and `--bench-batch` accept the original benchmark command line.
This mode accepts compressed FASTA and ordinary feature modes, not findere;
without `--records`, whole-file aggregate reporting is unchanged.

Record queries use a bounded parser/lookup/writer pipeline, recycling input and
output buffers. Ready input batches may be paired to amortize worker dispatch;
output remains in input-record and document order. Worker-local feature batches
hold at most 32,768 features and reuse complete lookup results for repeated keys,
including keys shared by different reads in that worker batch. Duplicate feature
occurrences retain their full multiplicity. Cached origin lists are capped at
262,144 32-bit IDs per worker and are discarded between feature batches; oversized
results are not retained. This is not a persistent index cache.

For positive `--min-ratio`, the total feature count is determined first. A filter
group is skipped only when **every** document in it has too few possible remaining
matches to pass. Contiguous retired groups are omitted from fingerprint reads;
1-bit fingerprints retire aligned 64-group blocks. Counts for every reported
document remain exact, including abandoned-key matches and 32-bit hash collisions.
Threshold-zero queries still report all counters, including zeros. Records longer
than the feature budget use a bounded-memory two-scan path rather than retaining
all their features. Whole-file aggregate queries deduplicate each feature batch
before partitioning it across workers, applying weighted counts.

Modern index loading decompresses independent banks in parallel, concurrently
with the global abandoned-key map; its postings use buffered bulk decoding.
Version-1 indexes retain legacy fallback conversion. These are execution-only
changes: no index-size increase, no new files or format changes, and no changes to
fingerprints, hash positions or abandoned-key collision handling. The positional
findere reconstruction path retains its existing per-feature lookup semantics.

### Metadata

```bash
zorindex info --index samples.zoridx
```

### Zindex: default merged power-of-two banks

```bash
zorindex build --fof samples.fof --k 31 --fingerprint-bits 16 --output samples.pow2
zorindex info --index samples.pow2
zorindex query --threads 4 --index samples.pow2 --query query.fa.gz
```

`--output` is a new directory. It contains `manifest.json` and one compressed
bank file per occupied power-of-two cell count. Each bank stores all groups
of that size in bit-packed slices, with a contiguous fingerprint plane followed
by origin bits. Queries XOR and screen fingerprints using AVX2 where available,
with a scalar fallback; origin bits are decoded only for matching fingerprints.

`--min-subindex-documents` defaults to **128 original documents**, counting all
members of each merged filter. Undersized slot-size buckets are combined by
promoting filters to a larger power of two, never by folding, dropping documents,
or changing the sketch merge plan. A dynamic program minimizes padded storage
among contiguous size-bucket partitions satisfying the minimum. Collections
with fewer than 128 documents use one subindex. Set the parameter to `1` to
retain the previous layout.

The constraint uses exact constructed sizes, not sketch estimates: only filters
needing promotion are rebuilt from their input documents. This may require an
additional input read and extra construction work, but introduces no temporary
key/cell disk files in the default mode. It trades larger filters for fewer
subindex probes. Fingerprint width, slot-scale (1.0), coloring and abandoned-key
semantics are unchanged; promoted filters can have different false-positive
identities because their hash positions change. `info` reports the original
document count of every final subindex and the number of promoted groups.
Groups are sorted by payload width, then by
samples of their assigned filter cells to cluster similar physical columns. Dataset
names and query output order stay the same as in the FOF; datasets are never split.

The planner sketches the actual indexed feature sets, so k-mers, minimizers,
modimizers, and findere use the same features for planning and construction.
Sketching uses vendored [Ragnar Groot Koerkamp's SimdSketch 0.3.0](https://github.com/RagnarGrootKoerkamp/simd-sketch)
with explicit feature-ID and union extensions (see `vendor/simd-sketch/NOTICE.md`).
Comparisons use **8-bit fingerprints**, AVX2 byte equality and correction for
accidental matches. The planner uses the dense-genome path: no empty-cell masks,
no occupancy checks during comparison, and all 256 fingerprint values available. AVX2
feature hashing/threshold collection reuses Ragnar's SIMD lane compaction.
Full **32-bit minima** are retained solely to support exact sketch unions:
truncated fingerprints cannot themselves be merged. At the default 4,096 buckets,
each sketch retains 4 KiB of fingerprints plus 16 KiB of merge state (20 KiB
total, versus the previous 32 KiB), with no mask allocation. Streaming additionally
packs the comparison bytes in one reusable, transposed arena (4 KiB per active
group at this setting, plus block padding and small position metadata).
This fast path assumes populated sketch buckets; similarity estimates for
tiny/sparse feature sets are not adjusted for unpopulated buckets.
The 8-bit sketch fingerprints are independent of the index fingerprint setting.
Unions update minima and fingerprints in place. Groups reaching the color cap
release their sketches immediately; all other sketches are freed when planning
finishes, before filter construction. Singleton cardinalities are counted exactly
during feature extraction; group cardinalities use collision-corrected Jaccard.
The new estimator can change merge decisions relative to the former 64-bit
sketches; planner pruning preserves exhaustive streaming decisions **with the
same new sketches**.

Planning modes are selectable with `--merge-planner`:

- `all-pairs`: repeatedly merge the pair with the largest predicted
  bit saving until no profitable pair remains. Candidate comparisons and heap
  storage can be quadratic in the number of datasets. Equal-saving ties use
  group IDs derived from input order, so its partition can also depend on order.
- `streaming` (default): process datasets in FOF order. Bucket current groups by
  slot exponent and origin count. Compare an incoming sketch only with buckets
  whose optimistic cost could save bits, excluding origin-limited buckets.
  Visit buckets in descending maximum possible saving, skipping those unable
  to beat the incumbent. Probe the earliest group first, stopping a bucket when
  that group achieves its upper bound. All comparison fingerprints are packed
  in one transposed arena: AVX2 compares the same fingerprint position across
  32 groups per vector, with independent vectors sharing query broadcasts. An exact,
  conservative cost bound determines the required number of matching bytes;
  blocks stop early when even all remaining bytes matching cannot suffice.
  Large scans use one coarse parallel dispatch, not repeated candidate batches
  or a freshly flattened candidate-ID vector. Intermediate scans use lazily
  created 8- or 16-worker pools (never exceeding the configured parallelism)
  below 16 or 32 MiB of comparison fingerprints; larger scans use the full
  configured pool. This avoids waking all workers for short scans. Auxiliary
  pools exist only during planning; extraction retains its configured parallelism.
  Removed columns are swap-compacted
  and freed blocks reused. Physical order may change, but logical group IDs
  retain exactly the original earliest-ID tie-breaking.
  Merge it into the group with the largest positive saving,
  or keep it separate. Update that group's sketch without revisiting old pairs.
  Pruning preserves the original streaming decisions and tie-breaking exactly;
  groups skipped for a small arrival remain eligible for larger future arrivals.
  There is no fixed size cap or merge reversion. Only current sketches are retained.
  Comparisons can still be quadratic in the worst case. Results depend on input order.
- `streaming-unpruned`: retain the original streaming scan over every current
  group, as a comparison baseline. It produces the same partition as `streaming`.

All modes use the same cost, including power-of-two rounding and payload width:
the fingerprint alone for a singleton, or fingerprint plus one membership bit
per dataset for a merged group. The fingerprint width is never reduced. These are
greedy optimizers, not globally optimal partition searches. Streaming is planning
only: all filters are constructed once after the batch's final plan is known.

Use `--merge-planner all-pairs` to opt into the all-pairs planner.

The streaming planner merges each dataset into the best existing group as it
arrives, so the input order matters: an order that keeps related genomes
adjacent gives smaller groups, fewer candidate evaluations and a smaller index.
`scripts/sort_fof_by_lineage.py LINEAGE_TSV IN.fof OUT.fof` orders a FOF by
phylum, class, order, family, genus and species from a lineage table (the
661k-bacteria `File2_taxid_lineage` layout); datasets without a lineage go last.
All modes default to `--merge-column-order compression`: within each bank's
equal-width sections, lexicographically sort samples from eight short windows
of the actual completed filter cells. This bounded-memory heuristic puts identical
sampled columns together without an all-pairs ordering search. It does not guarantee
smaller files on every collection. `--merge-column-order width` retains the previous
width/input-ID order for comparisons. Physical ordering never changes merges,
origin bit assignments, dataset identities, or query output order; offsets and member
IDs are stored in the manifest. Old manifests retain their legacy defaults.

Construction reuses the existing ZOR peel and cycle-breaking implementation.
Every retained equation's complete fingerprint/membership payload is verified.
By default, abandoned features are recovered through **one global multimap** shared
by all banks: a 32-bit hash of each abandoned feature maps to its originating
datasets (packed 32-bit IDs). Its open-addressed directory is queried once per
feature, not once per group. Hash collisions union the originating dataset lists;
there is no full-key collision check. Recovery adds to the filter's hits and never
overwrites them; each dataset is counted once per queried feature.
`build --drop-abandoned` omits both recovery storage and lookup,
accepting missing abandoned features. `query --ignore-abandoned` skips loading and
querying the map of an existing merged index. Recovery is enabled by default.
Version-1 and version-2 merged indexes remain readable. Version-1 fallback records are
converted into the global 32-bit-hash map on load, without rewriting the files.
New indexes use manifest version 5 (bank layout version 3) and require an updated reader.
**Merges are never
reverted**, even when exact union cardinalities, byte padding or fallback entries
erase the estimated savings. Actual capacity is computed from exact union keys
and remains a power of two; the planned membership partition is unchanged.
Construction defaults to **two input passes without temporary key files**:

1. Extract each dataset's exact cardinality and mergeable sketch, decide merges,
   and project the power-of-two bank layout. Only sketches, counts, feature
   signatures and the plan are retained; extracted keys are discarded.
2. Reread the inputs for each planned group, form its exact union and origin
   masks, and build its filter. Actual union cardinalities determine the final
   power-of-two sizes; estimates are not hard capacities and merges are never
   reverted. Changed indexed features between passes cause a clean error.

Completed filters are **bit-packed at their actual payload width in RAM**, not
stored as 64-bit cells. The default creates no temporary key/cell files; RAM must
fit the completed filters, planning metadata and active groups' keys/payloads/graphs.
The existing compression-aware column ordering is preserved using exactly the
same cell samples. Bank packing uses bounded row tiles and writes directly to
zstd; completed banks are released from RAM. The atomic publication directory
contains only the final index files. Both passes run as tasks on the
`--threads` rayon pool, which is the only concurrency setting. Parsing and
planning overlap: one pool thread inserts each parsed dataset into the
streaming planner, in input order, while the other threads parse the datasets
queued ahead of it, and the planner's own candidate scans use the same pool.
At most 1024 datasets are parsed ahead of the planner, so the slower side sets
the pace and memory stays bounded. Nothing is split by hand.

`--low-disk` explicitly selects the default scratch-free mode. To trade RAM for
scratch disk, use `--low-disk=false`: completed filters are staged as packed cell
files next to the output and removed after bank compression. `--disk-backed`
selects the previous single-extraction builder, using temporary `u64` feature
and cell caches followed by tiled bank packing. All modes produce identical bank files and the same global
abandoned-key map for unchanged inputs, settings and column order; only
construction metadata and timings differ.

Graph positions, adjacency, offsets and peel records use checked 32-bit indices
where possible, with a wide fallback. A mutable indexed degree heap avoids stale
priority entries while retaining the original tie-breaking and peel decisions.
Fingerprints, origin assignments, hash positions, false positives and query
counts are unchanged by construction mode.

Options: `--merge-sketch-size` (128 buckets), `--merge-max-colors` (48 datasets),
`--merge-min-saving` (0; a fractional safety margin), `--merge-plan-margin`
(0 by default: no inflation; an explicit positive value inflates the estimated
union only for merge planning). The manifest's `overflow_groups` counts groups
whose exact union required more slots than planned, before minimum-document
promotion. `--merge-disabled` selects
(power-of-two singleton banks for comparisons). Payloads currently fit at most
64 bits, including the 1-, 8- or 16-bit fingerprint. `--slot-count` acts as a minimum
per-group cell budget before rounding.

Query supports all feature modes and uses Rayon for ordinary and findere queries.
This format currently supports plain in-memory slices with zstd-compressed files;
append, repack, union prefilters, manual segment lengths, and split-large layouts
are rejected. A future append implementation must plan merges only within the
new update batch, never against old filters whose source documents may be gone.
`info` reports the planner, candidate evaluations, payload, padding, recovery-map
hashes/associations/RAM, fallback storage, metadata, and disk bytes
separately. The 1,024-genome before/after benchmark is summarized in
[docs/OPTIMIZATION_BENCHMARK.txt](docs/OPTIMIZATION_BENCHMARK.txt).

## Bindex

Bindex is the Bloom-filter comparison backend. It stores one Bloom filter per
color with bit-position-major interleaving:

```text
bit position 0: packed colors
bit position 1: packed colors
...
```

For a query, Bloom positions are computed once and the packed color words at
those positions are ANDed.

### Build

```bash
bloomindex build \
  --fof samples.fof \
  --k 31 \
  --false-positive-rate 0.00390625 \
  --hashes 6 \
  --output samples.blmidx
```

Sizing options:

- `--false-positive-rate`: target per-color Bloom false-positive rate.
- `--hashes`: Bloom hash count.
- `--bits-per-item`: override automatic bits-per-item sizing.
- `--bit-count`: explicit common bit count for every color.

Automatic sizing uses the largest color cardinality and applies it to all
colors. Appended colors larger than the original sizing target increase their
Bloom false-positive rate; append reports a warning in that case.

### Append, Query, Info

```bash
bloomindex append --index samples.blmidx --fof more_samples.fof \
  --output samples.extended.blmidx
bloomindex query --index samples.blmidx --query query.fa.gz --min-ratio 0.1
bloomindex info --index samples.blmidx
```

## Choosing Bloom Hash Counts

For a target false-positive rate `p` and `h` Bloom hashes, Bindex uses:

```text
bits_per_item = -h / ln(1 - p^(1/h))
```

The optimal Bloom hash count depends on the target `p`, memory budget, and
throughput target. For comparable experiments, fix the target false-positive
rate and sweep `--hashes`; the automatic sizing will adjust the bit count.

## Serialization

`build`, `append`, and `repack` always serialize an index to disk. Both Zindex
and Bindex write the serialized index stream through zstd by default. The
loaders detect zstd streams automatically and can still read raw legacy streams.

## Experimental Stack Compression

Zindex also has experimental per-row stack codecs behind `--stack-compression`.
This is separate from whole-index zstd serialization. The default stack mode is
`none`.

Available experimental modes:

```text
svb32-0124
lz4
lz4-lib
snappy
fastpfor256
lz4-hc
fastpfor-pack
```

Repack an existing index with an experimental stack codec:

```bash
zorindex repack --index samples.zoridx --output samples.lz4hc.zoridx \
  --stack-compression lz4-hc
```

## Testing

Run all unit tests:

```bash
cargo test --release
```

List test cases:

```bash
cargo test --release -- --list
```

The suite covers canonical encoding, reverse-complement equivalence, feature
modes, FOF parsing helpers, ZOR layouts, stack storage round trips, 8-bit and
16-bit query paths, Bloom sizing, and Bloom query behavior.

## Citation

Zindex uses ZOR filters:

Antoine Limasset. ZOR Filters: Fast and Smaller Than Fuse Filters. In 24th
International Symposium on Experimental Algorithms (SEA 2026), LIPIcs 371,
24:1-24:17, 2026. DOI:
[10.4230/LIPIcs.SEA.2026.24](https://doi.org/10.4230/LIPIcs.SEA.2026.24).

Compression applies to construction inputs and queries, including findere. Use the appropriate compression suffix; decoding errors are returned. Compressed construction inputs are decoded into memory before feature parsing.
