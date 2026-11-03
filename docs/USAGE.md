# Usage

Run commands from the repository root after `cargo build --release --locked`.
The bundled [bacterial example](../tests/data/bacteria/README.md) supplies the
input files used below.

## Build an index

### 1. Index all k-mers

```bash
./target/release/zorindex build --fof tests/data/bacteria/samples.fof \
  --output bacteria.zindex --k 31
```

Each FOF line is an input FASTA/FASTQ path, optionally preceded by a dataset
name and a tab. Relative paths use the working directory. The output must be
a new directory.

### 2. Index minimizers

```bash
./target/release/zorindex build --fof tests/data/bacteria/samples.fof \
  --output bacteria.min.zindex --k 31 --index-mode minimizers --window 19
```

Here `k=31` is the seed length and `w=19` counts consecutive candidate k-mers.
The selection window spans `k + w - 1 = 49` bases. Use `-k` and `-w` as short
forms. Minimizer mode requires `--window`.
Open-closed mod minimizers are the default selector. To choose the alternative
SIMD hash-minimizer selector, add `--minimizer-scheme simd`.

### 3. Use spaced seeds

```bash
./target/release/zorindex build --fof tests/data/bacteria/samples.fof \
  --output bacteria.spaced.zindex --k 21 --seed-pattern '#_#_#_#_#_#_#_#_#_#_#'
```

The pattern must be symmetric and its width must match `--k` in every mode.
A `#` retains a base and `_` skips one. In minimizer mode, select ordinary
minimizers first, then apply the pattern.

### 4. Use randstrobemers

```bash
./target/release/zorindex build --fof tests/data/bacteria/samples.fof \
  --output bacteria.strobe.zindex --k 21 --randstrobemers \
  --strobe-min 25 --strobe-max 100
```

Two strobes are linked. `k` is the length of each strobe. The distances are
between their starts. In minimizer mode, already-selected minimizers are linked.
Queries score both orientations and take the larger count per dataset and record.

### 5. Sample seeds by hash

```bash
./target/release/zorindex build --fof tests/data/bacteria/samples.fof \
  --output bacteria.sampled.zindex --k 31 --index-mode modimizers \
  --modimizer-sampling 16
```

This retains roughly one in 16 canonical k-mers, without a minimizer window.

### 6. Compress homopolymers

Add `--homopolymer-compression` to collapse repeated bases before extracting
seeds. The same setting is applied automatically to queries. Lengths and strobe
distances then refer to the compressed sequence.

### 7. Use findere

```bash
./target/release/zorindex build --fof tests/data/bacteria/samples.fof \
  --output bacteria.findere.zindex --k 31 --index-mode findere
```

Findere indexes `(k-z)`-mers and requires `z+1` consecutive matches in the
same dataset to report a k-mer match. `z` defaults to 10. Change it with
`--findere-z`. It is separate from sampling, spaced seeds and randstrobemers.

Dense k-mers, modimizers and minimizer seeds require lengths in 1..=31.
Minimizer windows count k-mers and must have an odd `k + w - 1` span of at most
255 bases. Findere requires
`0 <= z < k <= 255` and `1 <= k-z <= 31`.

## Query an index

### 1. Query one sequence

```bash
./target/release/zorindex query --index bacteria.zindex \
  --query tests/data/bacteria/query.fa
```

This reports aggregate counts over the query file:

```text
color    matches    ratio    unique_kmers    abandoned
```

`color` is the dataset name. `matches` counts query seed occurrences and
`ratio` is the matched fraction. `unique_kmers` is the historical column name
for the number of distinct indexed seeds in that dataset. `abandoned` reports
constraints abandoned during construction. Their keys are recovered by default.

### 2. Filter the results

```bash
./target/release/zorindex query --index bacteria.zindex \
  --query tests/data/bacteria/query.fa --min-ratio 0.8
```

Only datasets reaching the specified matched fraction are printed.

### 3. Query several records separately

```bash
./target/release/zorindex query --index bacteria.zindex \
  --query tests/data/bacteria/queries.fa --records --min-ratio 0.8
```

```text
query_id    valid_kmers    hits
read_1      100            0:100,2:93
```

Here 100 query seeds were extracted. Datasets 0 and 2 matched 100 and 93.
IDs index the manifest's `colors` array. `valid_kmers` counts selected seeds,
including in minimizer mode. It does not count all k-mer windows or covered
bases. Findere aggregate counts instead refer to reconstructed k-mer windows.

Per-record reporting accepts FASTA and is unavailable for findere. Without
`--records`, queries aggregate over the complete FASTA/FASTQ file.
Extraction settings, seed patterns and homopolymer compression are read from
the index. Approximate membership can produce matches for absent seeds.

### 4. Choose memory use and parallelism

Queries load one subindex at a time. Query parsing and seed extraction happen
once, and the abandoned-key file is scanned once. Results remain in RAM.
Use `--query-memory all` to load all subindexes together.

Both construction and queries use the available CPU threads. Supply
`--threads N` to choose a limit.

## Inspect and tune the index

```bash
./target/release/zorindex info --index bacteria.zindex
```

The default layout uses streaming planning and compression-aware column order.
Only estimated space-saving merges are accepted, accounting for colored cells
and power-of-two slot rounding. Merges are decided before constructing filters
and are not reverted afterward.

| Option | Purpose |
|---|---|
| `--fingerprint-bits 1,2,4,8,16` | Select fingerprint width, default 8 |
| `--hashes N` | Hash locations per filter, default 4 |
| `--merge-planner full` | Full pairwise planner for small collections (`all-pairs` is an alias) |
| `--merge-planner naive` | Ordered placement without dataset unions |
| `--merge-sketch-size N` | Sketch buckets, default 128 |
| `--max-subindexes N` | Cap size classes before physical hash partitioning |
| `--min-subindex-documents N` | Promote undersized classes, default 128 |
| `--max-subindex-size 8GiB` | Maximum uncompressed filter payload per subindex |
| `--low-disk` | Keep completed construction data in RAM instead of staging it |

The slot scale defaults to 1.0. The separate planning margin defaults to 0.03.
Sketch comparisons use 8-bit fingerprints, retaining full minima for merging.
Fixed 64-cell segments, no-heuristic cycle breaking and tie-scan 1 are defaults.

The subindex cap bounds filter payload, not total process RAM. Query seeds,
results, metadata and workers also use memory. Oversized documents can be
skipped by construction's memory guard. Inspect the skipped-document report.
Use `--skip-oversized-documents=false` to require all supplied inputs, and
`--max-ram-gb N` for the sampled process-RSS guard. Admission estimates and
sampled RSS checks do not replace an operating-system memory limit.

## Index maintenance

Merged indexes are rebuilt to incorporate new datasets. Legacy `--plain`
indexes also support `append` and ordinary `repack`.
`repack-subindexes` rewrites compression framing without changing filter cells:

```bash
./target/release/zorindex repack-subindexes --index bacteria.zindex \
  --output bacteria.reframed.zindex
```

See [format notes](FORMAT.md) and [troubleshooting](TROUBLESHOOTING.md).

## Advanced option reference

The README introduces the seed schemes, fingerprints, full/streaming planners
and normal queries. The following options provide more control.

| Build option | Effect |
|---|---|
| `--hashes N` | Hash locations per filter equation, 2..32 with default 4 |
| `--seed N` | Shared hash seed, default 69, stored for queries |
| `--merge-column-order compression\|width` | Compression-aware column order (default), or order by width and input ID |
| `--merge-sketch-size N` | Power-of-two sketch bucket count, default 128 |
| `--merge-max-colors N` | Maximum datasets sharing a filter, default 48, subject to the 64-bit payload limit |
| `--merge-min-saving F` | Minimum predicted fractional space saving for a merge, default 0 |
| `--merge-plan-margin F` | Inflate estimated unions for capacity planning, default 0.03 |
| `--merge-disabled` | Build singleton filters to measure layout without merging |
| `--merge-planner naive` | Place datasets in input order without unions |
| `--merge-planner streaming-unpruned` | Unpruned streaming baseline for comparisons |
| `--min-subindex-documents N` | Promote undersized size classes, default 128 |
| `--max-subindexes N` | Cap size classes before hash partitioning |
| `--max-subindex-size SIZE` | Uncompressed filter payload cap per physical subindex, default 8GiB |
| `--slot-scale F` | Filter-capacity scale, default 1.0 |
| `--slot-count N` | Override filter slot count |
| `--segment-length N` | Power-of-two segment size, default 64 for merged builds |
| `--cycle-break MODE`, `--tie-scan N` | Cycle-breaking strategy and tie scan count, defaults `no-heuristic` and 1 |
| `--drop-abandoned` | Omit recovery storage, allowing present seeds to be missed |
| `--disk-backed` | Stage feature, filter-cell and recovery data in temporary files |
| `--low-disk[=true\|false]` | Retain completed filters and recovery data in RAM, default false |
| `--build-memory-gib N` | Concurrent workspace admission budget, default 80% of available memory |
| `--max-ram-gb N` | Sampled process-RSS guard in decimal GB, alternative to `--build-memory-gib` |
| `--max-document-memory-gib N` | Per-document workspace admission limit |
| `--skip-oversized-documents[=true\|false]` | Skip and report oversized inputs (default true), or require every input |

The admission budgets and sampled RSS guard can reduce concurrent construction
work. A subindex payload cap alone does not cap process memory.

| Query option | Effect |
|---|---|
| `--record-batch-size N` | Parallel reporting batch size, default 4096 records |
| `--ignore-abandoned` | Skip recovery lookup, allowing present seeds to be missed |
| `--ignore-union` | Disable the optional union-graph precheck in legacy indexes |
| `--decode-stack` | Decode a legacy compressed stack before querying it |

### Legacy layouts and compatibility

`--plain` builds the legacy single-file layout. Its additional build controls
are `--union-graph` for a separate union precheck, `--stack-compression` for
slice codecs, and `--technical-layout flat|split-large` for splitting large
datasets into technical bins. `--target-technical-features` sets a split target,
`--max-split-parts` limits parts per dataset (default 64), and `--build-jobs`
limits concurrent dataset construction in this layout.

`append --index INDEX --fof INPUTS [--output OUTPUT]` adds datasets to a plain
index. Omitting `--output` replaces the input index. It accepts `--cycle-break`,
`--tie-scan` and `--stack-compression` for the new construction.

`repack --index INDEX --output OUTPUT --stack-compression CODEC` changes a plain
index's compression. `repack-subindexes --index INDEX --output OUTPUT` changes
merged-index compression framing. `info --index INDEX` inspects either layout.
All commands accept `--threads` and `--help`. The root command accepts `--version`.

`--power-of-two-merge` is a redundant compatibility flag for the default layout.
`--stack-block-rows` is deprecated and does not change the one-row slice unit.
`all-pairs` remains an alias for the `full` planner. Other aliases are
`--homocompress`, `--spaced-seed`, `--randstrobes`, `--strobe-window-min`,
`--strobe-window-max`, `--z`, `--bench-records` and `--bench-batch`.

### Migrating minimizer commands

Specify `--k 31 --window 19 --index-mode minimizers` for a 31-base seed in a
49-base span. `--minimizer-size` has been removed. The on-disk `k` field retains
its historical span meaning, with the seed length in `minimizer_size`. Existing
indexes need no conversion. `info` reports CLI-facing `k`, `window` and
`window_span` values for both old and new minimizer indexes.
