# Zindex

Zindex efficiently reports the number of matches between a query sequence and
every dataset in an indexed collection.

Zindex combines compact ZOR filters with color-preserving merging: similar
datasets share a filter while retaining their individual membership. A
sketch-based planner chooses space-saving merges and groups filters into
power-of-two subindexes, whose interleaved cells support SIMD queries.
An optimized implementation of low-density open-closed mod minimizers reduces
the number of indexed seeds; ordinary k-mers, spaced seeds and randstrobemers
are also supported.

## Install

```bash
git clone https://github.com/Malfoy/Zindex.git
cd Zindex
cargo build --release --locked
```

## Usage

[Ten real bacterial genomes](tests/data/bacteria/README.md), compressed as
FASTA.gz, are bundled with query sequences. Run the examples from the repository
root; no downloads are needed.

### Build an index

Start by indexing every 31-mer:

```bash
./target/release/zorindex build --fof tests/data/bacteria/samples.fof \
  --output bacteria.zindex --k 31
```

The file passed to `--fof` lists one dataset per line. The output is a new index
directory.

For a smaller index, select 31-base minimizers from 49-base windows:

```bash
./target/release/zorindex build --fof tests/data/bacteria/samples.fof \
  --output bacteria.min.zindex --k 49 --index-mode minimizers --minimizer-size 31
```

This uses low-density open-closed mod minimizers by default.

To use spaced seeds, supply a symmetric pattern (`#` retains a base, `_` skips it):

```bash
./target/release/zorindex build --fof tests/data/bacteria/samples.fof \
  --output bacteria.spaced.zindex --k 21 --seed-pattern '#_#_#_#_#_#_#_#_#_#_#'
```

To use randstrobemers, link two 21-base strobes with start positions 25–100 bases
apart:

```bash
./target/release/zorindex build --fof tests/data/bacteria/samples.fof \
  --output bacteria.strobe.zindex --k 21 --randstrobemers \
  --strobe-min 25 --strobe-max 100
```

### Query an index

Query the simple 31-mer index with a 1-kb fragment from the first genome:

```bash
./target/release/zorindex query --index bacteria.zindex \
  --query tests/data/bacteria/query.fa
```

The result gives each dataset's match count and the fraction of query seeds
matched. Membership is approximate, so absent seeds can also match.

To report only datasets matching at least 80% of the query seeds:

```bash
./target/release/zorindex query --index bacteria.zindex \
  --query tests/data/bacteria/query.fa --min-ratio 0.8
```

To report each FASTA query separately, use `--records`. The bundled queries
include a fragment from each genome, its reverse complement and a negative
control:

```bash
./target/release/zorindex query --index bacteria.zindex \
  --query tests/data/bacteria/queries.fa --records --min-ratio 0.8
```

The output contains `query_id`, `valid_kmers` (selected seed count) and `hits`.
A hit such as `0:970` means 970 matches in dataset 0; dataset IDs follow the
input list's order.

To query an alternative index, change only `--index`; seed settings are read
automatically from the index:

```bash
./target/release/zorindex query --index bacteria.min.zindex \
  --query tests/data/bacteria/queries.fa --records --min-ratio 0.8
```

The optional checker builds the k-mer and minimizer indexes and verifies that
every positive query matches its source genome:

```bash
bash utils/scripts/demo.sh
```

### Use your own data

Create a file with one FASTA/FASTQ path per line and pass it through `--fof`.
Dataset names can precede paths, separated by a tab. Gzip, zstd, xz and bzip2
inputs are supported.

See the [usage guide](docs/USAGE.md) for additional seed schemes, query options
and memory settings, or the [standalone minimizer utility](utils/minimizer-fasta/README.md)
for exporting selected seeds.

## Citation

TBD.

[MIT license](LICENSE).
