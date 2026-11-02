# Zindex

[![CI](https://github.com/Malfoy/Zindex/actions/workflows/ci.yml/badge.svg)](https://github.com/Malfoy/Zindex/actions/workflows/ci.yml)

Zindex searches sequence collections using compact approximate membership
indexes. It supports minimizers, k-mers, spaced seeds and randstrobemers, and
reports the datasets matching each query. Similar datasets can share a filter
while retaining their individual membership.

## Install

Use Rust through rustup (the repository selects Rust 1.97.0) and a C/C++ toolchain.
On Debian/Ubuntu, native build dependencies are `build-essential pkg-config`.
Python 3 and Bash are needed for the examples. Linux x86-64 is tested in CI.

```bash
git clone https://github.com/Malfoy/Zindex.git
cd Zindex
cargo build --release --locked
```

Build on the machine where you will run the tools: native CPU instructions are
enabled by default. See [troubleshooting](docs/TROUBLESHOOTING.md) for portable
builds, compiler details and installation problems.

## Try the bundled example

```bash
bash utils/scripts/demo.sh
```

No downloads are required. The [tiny dataset](tests/data/tiny/README.md) contains
three 8-kb references and six 1-kb queries, with expected answers for shared
sequences, dataset-specific sequences, reverse complements and a negative read.
The script checks all three tools and ends with:

```text
PASS: example queries and all three tools
```

Logs and results are retained in a new temporary directory whose path is printed.
For an optional example using **three real bacterial genomes**, see the
[checksum-verified download and query instructions](tests/data/bacteria/README.md).

## Index your own collection

Create `samples.fof` with one FASTA/FASTQ path per line, or a dataset name and
path separated by a tab. Relative paths use the command's working directory.
Gzip, zstd, xz and bzip2 inputs are supported; use the corresponding suffix.

```bash
./target/release/zorindex build --fof samples.fof --output samples.zindex \
  --k 49 --index-mode minimizers --minimizer-size 31 --fingerprint-bits 16
./target/release/zorindex query --index samples.zindex --query reads.fa.gz \
  --records --min-ratio 0.8
./target/release/zorindex info --index samples.zindex
```

The index output must be a new directory. Query extraction settings are reused
automatically from the index. Inspect `info` and any skipped-document report to
confirm which inputs were included.

## Read the results

`--records` reports each query FASTA record separately:

```text
query_id    valid_kmers    hits
read_1      100            0:100,2:93
```

This read has 100 selected seed occurrences; datasets 0 and 2 match 100 and 93.
Dataset IDs refer to the manifest's `colors` array. In minimizer mode, counts
refer to selected seeds, not all k-mer windows or covered bases. Membership is
approximate, so absent seeds can match.

Without `--records`, results aggregate over the entire FASTA/FASTQ file.
Per-record reporting is available for merged indexes except findere.
See [output formats and modes](docs/USAGE.md#query) for details.

## Documentation

- [Full usage, defaults, memory controls and seed options](docs/USAGE.md)
- [Standalone minimizer FASTA generator](utils/minimizer-fasta/README.md)
- [Bindex comparison commands](docs/COLLABORATOR.txt)
- [Benchmarking](docs/BENCHMARKING.md) and [manuscript artifact status](docs/REPRODUCIBILITY.md)
- [Formats](docs/FORMAT.md), [changelog](docs/CHANGELOG.md), and [contributing](.github/CONTRIBUTING.md)

The build produces `zorindex`, `bloomindex` (the Bloom-filter comparison tool),
and `minimizer-fasta` under `target/release/`. Zindex defaults to streaming merge
planning; minimizer mode uses open-closed mod selection. Other seed modes and
legacy layouts are described in the usage guide.

## Tests

```bash
cargo test --workspace --release --locked
```

The Rust suite has 389 unit/integration tests plus a doctest, covering extraction,
filter layouts, compression, query modes, parallelism and memory controls.
CI runs native and portable x86-64 release tests, the bundled demonstration,
and a tiny-input check of the benchmark recorder.

## Citation and license

See [CITATION.cff](CITATION.cff) for software citation metadata and record the
release or Git revision used. Zindex builds on Antoine Limasset's
[ZOR Filters: Fast and Smaller Than Fuse Filters (SEA 2026)](https://doi.org/10.4230/LIPIcs.SEA.2026.24).

[MIT license](LICENSE); vendored dependencies retain their licenses and notices.
