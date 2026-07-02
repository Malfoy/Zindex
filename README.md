# ZORINDEX

ZORINDEX contains two Rust command-line indexes for approximate colored
k-mer queries over many FASTA/FASTQ datasets:

- `zorindex`: Zindex, an interleaved stack of per-color pure ZOR filters.
- `bloomindex`: Bindex, a matching interleaved stack of per-color Bloom filters.

Both binaries use the same input format, DNA feature extraction code, canonical
2-bit encoding, parallel parsing, and query reporting. This makes Bindex a
direct comparison target for Zindex when evaluating size, build time, and query
throughput.

## Status

The implementation targets high-throughput experiments on large bacterial
collections. It is intentionally approximate:

- Zindex uses 8-bit or 16-bit fingerprints. A query hit is reported when the
  XOR of the selected cells equals the key fingerprint.
- Bindex uses a configurable Bloom false-positive rate and hash count.
- Zindex may abandon keys when a pure ZOR layer reaches an unpeelable core.
  Build output reports abandoned keys per color.

The source tree is organized as:

```text
src/lib.rs              shared DNA parsing, feature extraction, hashing, FOF parsing
src/bin/zorindex.rs     Zindex CLI and pure ZOR stack implementation
src/bin/bloomindex.rs   Bindex CLI and interleaved Bloom stack implementation
```

## Requirements

- Rust stable toolchain
- x86_64 or a platform supported by the selected compression crates
- Enough RAM to hold the parsed feature sets during construction

The release profile is tuned for local benchmark binaries:

- `opt-level = 3`
- `lto = "fat"`
- `codegen-units = 1`
- stripped symbols
- `panic = "abort"`
- `.cargo/config.toml` sets `-C target-cpu=native`

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
- zstd-compressed FASTA/FASTQ

Parsing uses `helicase`. Compressed files are decompressed through `deko`.

## Feature Modes

Both indexes support the same indexed feature modes:

```text
kmers       every canonical k-mer
minimizers  canonical minimizers from each k-window
modimizers  hash-sampled canonical k-mers
```

Full k-mer and modimizer modes currently require `k <= 31`, because canonical
DNA k-mers are encoded in a `u64`. Minimizer mode may use larger odd-length
windows because only the minimizer sequence is encoded; the default minimizer
length is 21. The canonical SIMD minimizer backend requires odd `k`.

Examples:

```bash
zorindex build --fof samples.fof --k 31 --output samples.zoridx
zorindex build --fof samples.fof --k 63 --index-mode minimizers \
  --minimizer-size 21 --output samples.m21.zoridx
zorindex build --fof samples.fof --k 31 --index-mode modimizers \
  --modimizer-sampling 16 --output samples.mod16.zoridx
```

The feature mode is stored in the index and reused by `query` and `append`.

## Zindex

Zindex builds one pure ZOR filter per color and stores the color filters in a
slot-major interleaved layout:

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
  --fof samples.fof \
  --k 31 \
  --hashes 4 \
  --fingerprint-bits 16 \
  --slot-scale 1.0 \
  --cycle-break most-deg2 \
  --tie-scan 8 \
  --output samples.zoridx
```

Important options:

- `--hashes`: number of ZOR locations per key.
- `--fingerprint-bits`: `8` or `16`.
- `--slot-scale`: multiplier applied to the largest color cardinality.
- `--slot-count`: explicit slot budget before segmented layout rounding.
- `--segment-length`: explicit power-of-two segment length.
- `--cycle-break`: heuristic used when the peel reaches a core.
- `--union-graph`: optional union ZOR filter used as a query prefilter.
- `--stack-compression`: row codec for the slot-major stack.

### Stack Compression

Available modes:

```text
none
svb32-0124
lz4
lz4-lib
snappy
fastpfor256
lz4-hc
fastpfor-pack
```

Compression is applied to complete slot rows. Querying compressed indexes keeps
the stack compressed and decodes only the rows touched by each feature. Use
`query --decode-stack` to decode to the plain layout before querying when that
is useful for a specific benchmark.

Repack an existing index:

```bash
zorindex repack --index samples.zoridx --output samples.lz4hc.zoridx \
  --stack-compression lz4-hc
```

### Append

```bash
zorindex append --index samples.zoridx --fof more_samples.fof \
  --output samples.extended.zoridx \
  --cycle-break most-deg2
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
modes, FOF parsing helpers, ZOR layouts, compressed-stack round trips, 8-bit
and 16-bit query paths, Bloom sizing, and Bloom query behavior.

## Citation

Zindex uses ZOR filters:

Antoine Limasset. ZOR Filters: Fast and Smaller Than Fuse Filters. In 24th
International Symposium on Experimental Algorithms (SEA 2026), LIPIcs 371,
24:1-24:17, 2026. DOI:
[10.4230/LIPIcs.SEA.2026.24](https://doi.org/10.4230/LIPIcs.SEA.2026.24).
