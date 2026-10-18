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
length is 21. The canonical SIMD minimizer backend requires odd `k`.

Examples:

```bash
zorindex build --fof samples.fof --k 31 --output samples.zindex
zorindex build --fof samples.fof --k 63 --index-mode minimizers \
  --minimizer-size 21 --output samples.m21.zindex
zorindex build --fof samples.fof --k 31 --index-mode modimizers \
  --modimizer-sampling 16 --output samples.mod16.zindex
```

The feature mode is stored in the index and reused by `query`. Legacy `--plain`
indexes also support `append`; merged directories currently require a rebuild.

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

Existing feature modes keep their serialized layout. Findere uses feature-mode
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
Groups are sorted by payload width, then by
samples of their assigned filter cells to cluster similar physical columns. Dataset
names and query output order stay the same as in the FOF; datasets are never split.

The planner sketches the actual indexed feature sets, so k-mers, minimizers,
modimizers, and findere use the same features for planning and construction.
Sketching follows the threshold-filtered mod-bucket method described by
[Ragnar Groot Koerkamp's SimdSketch](https://curiouscoding.nl/posts/simd-sketch/),
using full 64-bit minima to retain mergeability and AVX2 comparisons where
available. Singleton cardinalities are counted exactly during feature extraction;
group union cardinalities are estimated from sketch Jaccard values. This is an
implementation of the method for feature IDs, not a dependency on the upstream
sequence-sketching crate.

Planning modes are selectable with `--merge-planner`:

- `all-pairs`: repeatedly merge the pair with the largest predicted
  bit saving until no profitable pair remains. Candidate comparisons and heap
  storage can be quadratic in the number of datasets. Equal-saving ties use
  group IDs derived from input order, so its partition can also depend on order.
- `streaming` (default): process datasets in FOF order. Bucket current groups by
  slot exponent and origin count. Compare an incoming sketch only with buckets
  whose optimistic cost could save bits, excluding origin-limited buckets.
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
New indexes use manifest/bank layout version 3 and require an updated reader.
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

Completed filters are **bit-packed at their actual payload width**, not stored
as 64-bit cells, and staged by default as temporary cell files inside the
construction directory next to the output, so the scratch footprint is about the
uncompressed size of the final banks and RAM holds only planning metadata and
the active groups' keys/payloads/graphs. The existing compression-aware column
ordering is preserved using exactly the same cell samples. Bank packing reads
bounded row tiles from the staged cells and writes directly to zstd; each bank's
cell files are removed once it is written. The atomic publication directory
contains only the final index files. Both passes run as tasks on the
`--threads` rayon pool, which is the only concurrency setting. Parsing and
planning overlap: one pool thread inserts each parsed dataset into the
streaming planner, in input order, while the other threads parse the datasets
queued ahead of it, and the planner's own candidate scans use the same pool.
At most 1024 datasets are parsed ahead of the planner, so the slower side sets
the pace and memory stays bounded. Nothing is split by hand.

To build without any scratch disk, add `--low-disk`: the same packed filters are
then retained in RAM until bank compression, so RAM must also fit every
completed packed filter. `--disk-backed` retains the previous single-extraction
builder, using temporary `u64` feature and cell caches followed by tiled bank
packing. All modes produce identical bank files and the same global
abandoned-key map for unchanged inputs, settings and column order; only
construction metadata and timings differ.

Graph positions, adjacency, offsets and peel records use checked 32-bit indices
where possible, with a wide fallback. A mutable indexed degree heap avoids stale
priority entries while retaining the original tie-breaking and peel decisions.
Fingerprints, origin assignments, hash positions, false positives and query
counts are unchanged by construction mode.

Options: `--merge-sketch-size` (4096 buckets), `--merge-max-colors` (48 datasets),
`--merge-min-saving` (0; a fractional safety margin), `--merge-plan-margin`
(0.03: each merged group is planned as if its estimated union were 3% larger,
because the sketch estimate is noisy and an exact union that exceeds the planned
power of two doubles the group's slots at build time; the manifest's
`overflow_groups` counts the groups where that still happened), and `--merge-disabled`
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
