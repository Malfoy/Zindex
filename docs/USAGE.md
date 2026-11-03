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
  --output bacteria.min.zindex --k 49 --index-mode minimizers --minimizer-size 31
```

Here `k=49` is the selection-window span and `m=31` is the seed length.
Open-closed mod minimizers are the default selector. To choose the alternative
SIMD hash-minimizer selector, add `--minimizer-scheme simd`.

### 3. Use spaced seeds

```bash
./target/release/zorindex build --fof tests/data/bacteria/samples.fof \
  --output bacteria.spaced.zindex --k 21 --seed-pattern '#_#_#_#_#_#_#_#_#_#_#'
```

The pattern must be symmetric and match `k`, or `m` in minimizer mode.
A `#` retains a base and `_` skips one. In minimizer mode, select ordinary
minimizers first, then apply the pattern.

### 4. Use randstrobemers

```bash
./target/release/zorindex build --fof tests/data/bacteria/samples.fof \
  --output bacteria.strobe.zindex --k 21 --randstrobemers \
  --strobe-min 25 --strobe-max 100
```

Two strobes are linked. `k` is the length of each strobe; the distances are
between their starts. In minimizer mode, `m` specifies each strobe's length and
already-selected minimizers are linked. Queries score both orientations and
take the larger count per dataset and record.

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
same dataset to report a k-mer match. `z` defaults to 10; change it with
`--findere-z`. It is separate from sampling, spaced seeds and randstrobemers.

Dense k-mers and modimizers require lengths up to 31. Minimizer windows can be
odd lengths up to 255, with `m <= 31`. Findere requires
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
constraints abandoned during construction; their keys are recovered by default.

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

Here 100 query seeds were extracted; datasets 0 and 2 matched 100 and 93.
IDs index the manifest's `colors` array. `valid_kmers` counts selected seeds,
including in minimizer mode. It does not count all k-mer windows or covered
bases. Findere aggregate counts instead refer to reconstructed k-mer windows.

Per-record reporting accepts FASTA and is unavailable for findere. Without
`--records`, queries aggregate over the complete FASTA/FASTQ file.
Extraction settings, seed patterns and homopolymer compression are read from
the index. Approximate membership can produce matches for absent seeds.

### 4. Choose memory use and parallelism

Queries load one subindex at a time. Query parsing and seed extraction happen
once, and the abandoned-key file is scanned once; results remain in RAM.
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
| `--fingerprint-bits 1,2,4,8,16` | Select fingerprint width; default 8 |
| `--hashes N` | Hash locations per filter; default 4 |
| `--merge-planner all-pairs` | Full pairwise planner instead of streaming |
| `--merge-planner naive` | Ordered placement without dataset unions |
| `--merge-sketch-size N` | Sketch buckets; default 128 |
| `--max-subindexes N` | Cap size classes before physical hash partitioning |
| `--min-subindex-documents N` | Promote undersized classes; default 128 |
| `--max-subindex-size 8GiB` | Maximum uncompressed filter payload per subindex |
| `--low-disk` | Keep completed construction data in RAM instead of staging it |

The slot scale defaults to 1.0; the separate planning margin defaults to 0.03.
Sketch comparisons use 8-bit fingerprints, retaining full minima for merging.
Fixed 64-cell segments, no-heuristic cycle breaking and tie-scan 1 are defaults.

The subindex cap bounds filter payload, not total process RAM. Query seeds,
results, metadata and workers also use memory. Oversized documents can be
skipped by construction's memory guard; inspect the skipped-document report.
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
