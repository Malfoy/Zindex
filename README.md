# Zindex

Zindex is a resource-efficient approximate index for searching sequence
collections, in the same family as [Raptor](https://github.com/seqan/raptor),
[COBS](https://github.com/bingmann/cobs) and
[kmindex](https://github.com/tlemane/kmindex). It reports the number of query
seed matches in each indexed dataset.

Zindex combines [ZOR filters](https://doi.org/10.4230/LIPIcs.SEA.2026.24) with
color-preserving merging, allowing similar datasets to share a filter while
retaining their individual membership. Merge planning uses Ragnar Groot
Koerkamp's [SimdSketch](https://github.com/RagnarGrootKoerkamp/simd-sketch),
adapted to sketch indexed seeds and merge their sketches. Filters are grouped
into power-of-two subindexes with interleaved cells for SIMD queries.
Seed extraction includes an optimized implementation of the
[open-closed mod-minimizer algorithm](https://doi.org/10.1186/s13015-025-00270-0)
by Ragnar Groot Koerkamp, Daniel Liu and Giulio Ermanno Pibiri. FASTA/FASTQ
parsing uses [Helicase](https://doi.org/10.64898/2026.03.19.712912) by Igor
Martayan, Loup Lobet, Camille Marchet and Charles Paperman.

## Install

```bash
git clone https://github.com/Malfoy/Zindex.git
cd Zindex
cargo build --release --locked
```

## Build an index

[Ten real bacterial genomes](tests/data/bacteria/README.md), compressed as
FASTA.gz, are bundled with queries. Run these examples from the repository root.

Start by indexing every 31-mer:

```bash
./target/release/zorindex build --fof tests/data/bacteria/samples.fof \
  --output bacteria.zindex --k 31
```

`--fof` lists one FASTA/FASTQ path per dataset, optionally preceded by a name and
a tab. Relative paths use the working directory. Gzip, zstd, xz and bzip2 inputs
are supported. The output must be a new directory.

### Fingerprint size

`--fingerprint-bits` controls the space/accuracy tradeoff. Supported widths are
1, 2, 4, 8 and 16 bits, with 8 as the default. Wider fingerprints use more space
and reduce accidental seed matches. The filter fingerprint collision probability
is approximately `2^-b` for `b` bits, about 0.39% at 8 bits and 0.0015% at 16 bits.
This is a per-seed filter probability, rather than a whole-query error rate.

To use 16-bit fingerprints:

```bash
./target/release/zorindex build --fof tests/data/bacteria/samples.fof \
  --output bacteria.fp16.zindex --k 31 --fingerprint-bits 16
```

### Full or streaming construction

For a small collection, use `--merge-planner full` to consider all candidate
pairs and repeatedly choose the merge with the largest estimated space saving:

```bash
./target/release/zorindex build --fof tests/data/bacteria/samples.fof \
  --output bacteria.full.zindex --k 31 --merge-planner full
```

For large collections, use `--merge-planner streaming`, the default. It makes
merge decisions as datasets arrive, reducing planning work and memory. The
choice depends on collection size and similarity, with no fixed dataset-count
cutoff. Both planners preserve individual dataset identities.

### Minimizer sampling

For a smaller index, select 31-base seeds from windows of 19 consecutive k-mers:

```bash
./target/release/zorindex build --fof tests/data/bacteria/samples.fof \
  --output bacteria.min.zindex --k 31 --index-mode minimizers --window 19
```

`-k` / `--k` is the seed length. `-w` / `--window` is the number of candidate
k-mers, so each window spans `k + w - 1` bases, 49 here. Minimizer mode requires
`--window`, with `1 <= k <= 31` and an odd total span of at most 255 bases.
The default selector is open-closed mod. Add `--minimizer-scheme simd` to select
standard SIMD hash minimizers instead.

## Query an index

Query the simple 31-mer index with a 1-kb fragment from the first genome:

```bash
./target/release/zorindex query --index bacteria.zindex \
  --query tests/data/bacteria/query.fa
```

The result gives each dataset's match count and matched fraction, aggregated
over the query file. Membership is approximate, so absent seeds can also match.

To report only datasets matching at least 80% of query seeds:

```bash
./target/release/zorindex query --index bacteria.zindex \
  --query tests/data/bacteria/query.fa --min-ratio 0.8
```

To report each FASTA query separately, add `--records`. The bundled queries
include a fragment from each genome, its reverse complement and a negative
control:

```bash
./target/release/zorindex query --index bacteria.zindex \
  --query tests/data/bacteria/queries.fa --records --min-ratio 0.8
```

The columns are `query_id`, `valid_kmers` (selected seed count) and `hits`.
A hit such as `0:970` means 970 matches in dataset 0. Dataset IDs refer to the
index manifest's `colors` array. In minimizer mode the counts refer to sampled
seeds, rather than all k-mer positions or covered bases.

To query an alternative index, change only `--index`. Seed settings, patterns
and homopolymer compression are read from the index:

```bash
./target/release/zorindex query --index bacteria.min.zindex \
  --query tests/data/bacteria/queries.fa --records --min-ratio 0.8
```

Construction and queries use all available CPU threads by default. Add
`--threads N` to limit them. Queries load one subindex at a time. Use
`--query-memory all` to load all subindexes together when memory permits.

Inspect an index with:

```bash
./target/release/zorindex info --index bacteria.zindex
```

## Spaced seeds

A spaced seed retains the bases marked `#` and skips those marked `_`.
The pattern must be symmetric under reversal and its width must equal `--k`,
including in minimizer mode. Patterns can span at most 31 bases. For example,
this 21-base pattern retains 11 positions:

```bash
./target/release/zorindex build --fof tests/data/bacteria/samples.fof \
  --output bacteria.spaced.zindex --k 21 --seed-pattern '#_#_#_#_#_#_#_#_#_#_#'
```

Choose a mask for your read length, retained-base count and expected number of
substitutions. Jens Zentgraf and Sven Rahmann's
[Design of Worst-Case-Optimal Spaced Seeds (2025)](https://doi.org/10.4230/LIPIcs.WABI.2025.22)
optimizes the minimum number of hits or covered positions over substitution
placements. Their [optimization software](https://gitlab.com/rahmannlab/seed-optimization)
can help choose symmetric masks within Zindex's 31-base span limit. A wider
span can improve substitution tolerance but is more vulnerable to indels.

Moein Karami, Jens Zentgraf and Sven Rahmann's
[2026 follow-up](https://doi.org/10.4230/LIPIcs.WABI.2026.21) optimizes filters
using both hit counts and coverage. Zindex reports hit counts and does not
implement that combined criterion. The papers' dense-seed guarantees do not
automatically carry over to minimizer sampling.

## Randstrobemers

[Randstrobemers](https://doi.org/10.1101/gr.275648.121), introduced by Kristoffer
Sahlin, link separated sequence pieces to improve tolerance to indels. Here
`--k` is the length of each strobe, and distances are measured between their
start positions:

```bash
./target/release/zorindex build --fof tests/data/bacteria/samples.fof \
  --output bacteria.strobe.zindex --k 21 --randstrobemers \
  --strobe-min 25 --strobe-max 100
```

Zindex links two strobes using the selection rule from
[strobealign](https://doi.org/10.1186/s13059-022-02831-7). The default distance
range is `k` to `3*k`. With minimizer sampling, it links already-selected seeds.
Queries score both orientations and take the larger count per dataset and
record. Query commands are the same as above.

## Additional modes and options

- `--homopolymer-compression` collapses repeated bases before seed extraction.
- `--index-mode modimizers --modimizer-sampling 16` keeps roughly one in 16
  canonical k-mers by hash.
- `--index-mode findere --k 31 --findere-z 10` indexes 21-mers and reconstructs
  31-mer matches from consecutive hits. Findere uses aggregate queries.

The [usage guide](docs/USAGE.md) covers construction memory limits, subindex
tuning, recovery options and legacy index maintenance. The
[standalone minimizer utility](utils/minimizer-fasta/README.md) exports selected
seeds. Run `zorindex build --help` or `zorindex query --help` for all flags.

The optional example checker builds k-mer and minimizer indexes and verifies
their source matches:

```bash
bash utils/scripts/demo.sh
```

## Citation

TBD.

[MIT license](LICENSE).
