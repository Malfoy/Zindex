# Zindex

[![CI](https://github.com/Malfoy/Zindex/actions/workflows/ci.yml/badge.svg)](https://github.com/Malfoy/Zindex/actions/workflows/ci.yml)

Zindex searches collections of sequence datasets using approximate membership
filters. It indexes minimizers, ordinary k-mers, spaced seeds and randstrobemers,
and reports which datasets contain the query's seeds. Similar datasets can share
a ZOR filter while preserving their individual membership through colored cells;
interleaved slices allow queries across many datasets with SIMD bitwise operations.

The repository provides three Rust command-line tools:

| Binary | Purpose |
|---|---|
| `zorindex` | Zindex: merged, interleaved ZOR filters; the main collection index |
| `bloomindex` | Bindex: interleaved Bloom filters for comparison |
| `minimizer-fasta` | Independent, parallel minimizer FASTA generator |

## Build

Linux x86-64 is tested in CI with both native CPU instructions and a portable
baseline. AVX2 accelerates supported kernels; portable fallbacks are available.
Other platforms are not currently covered by CI.

Requirements: Git, Rust 1.89 or newer with Cargo, and a C/C++ toolchain.
The repository pins Rust 1.97.0 in `rust-toolchain.toml` for reproducible builds;
rustup selects it automatically. CI uses the same version. This also avoids an
LLVM code-generation failure observed with Rust 1.99.1 when compiling the
Bzip2 dependency for an AVX-512-capable native target.
On Debian/Ubuntu, install native build tools with
`sudo apt install build-essential pkg-config`. Python 3 and Bash are additionally
needed for the demonstration and benchmark runner.

```bash
git clone https://github.com/Malfoy/Zindex.git
cd Zindex
cargo build --release --locked
./target/release/zorindex --version
```

All three binaries are under `target/release/`. Commands below run from the
repository root; no separate installation or PATH modification is required.
The default build uses `target-cpu=native`: rebuild on the machine that will run
it. For a portable x86-64 build, use
`RUSTFLAGS='-C target-cpu=x86-64' cargo build --release --locked`.
Tools use all CPUs available to the process unless `--threads N` is supplied.

`--version` includes the package version, Git revision and working-tree status.
Source archives without Git metadata report that the revision is unavailable.
See [CHANGELOG.md](CHANGELOG.md) for unreleased changes; the source in a tagged
release can differ from `main`.

## Run the small demonstration

No downloads are needed after building:

```bash
bash scripts/demo.sh
```

This builds and queries Zindex and Bindex, runs the standalone minimizer utility,
and checks the expected source datasets, reverse-complement matches and negative
read. It writes its inputs' results and logs to a new temporary directory and
prints that directory. To choose the destination, use
`bash scripts/demo.sh /path/to/new-demo-directory`; an existing path is rejected.

The versioned [example inputs](examples/tiny/README.md) include the FOF,
FASTA records and expected document hits. This is a functional check, not a
performance or false-positive-rate benchmark.

## Index your collection

Create a file-of-files (FOF), with one input path per line or a dataset name and
path separated by a real tab:

```text
sample_a<TAB>/path/to/sample_a.fna.gz
sample_b<TAB>/path/to/sample_b.fasta.xz
```

Relative input paths use the command's working directory, not the FOF's directory.
FASTA and FASTQ are supported, including gzip, zstd, xz and bzip2 compression.
For indexer inputs, use the corresponding compression suffix. Blank FOF lines
and lines beginning with `#` are ignored.

A typical minimizer index and per-read query are:

```bash
./target/release/zorindex build --fof samples.fof --output samples.zindex \
  --k 49 --index-mode minimizers --minimizer-size 31 --fingerprint-bits 16
./target/release/zorindex query --index samples.zindex --query reads.fa.gz \
  --records --min-ratio 0.8
./target/release/zorindex info --index samples.zindex
```

The output index must be a new directory. The extraction settings are stored in
it and automatically reused by queries. Per-record reporting accepts FASTA and
is available for merged indexes except findere. Without `--records`, results
aggregate over the whole FASTA/FASTQ query file.

## Understand the output

Aggregate output starts with a `#query` summary containing `total_features`,
then a tab-separated table:

```text
color   matches   ratio   unique_kmers   abandoned
```

Here `color` is the dataset name, `matches` counts matched query seed occurrences,
and `ratio = matches / total_features`. The historical name `unique_kmers`
means the number of distinct indexed features in that dataset, including in
minimizer and other seed modes. Merged-index `abandoned` is a construction
diagnostic; recovery is enabled by default.

With `--records`, output is:

```text
query_id    valid_kmers    hits
read_1      100            0:100,2:93
```

`valid_kmers` is the number of query seed occurrences despite its historical
name; `hits` contains zero-based `document_id:matched_count` pairs that pass
the requested ratio. Document IDs index the manifest's `colors` array.
Empty hits are an empty third field. `info` also reports any skipped build inputs.

Minimizer counts refer to selected seeds, not all k-mer windows or percentage
of covered bases. Findere aggregate counts refer to reconstructed k-mer windows.
Randstrobemer queries take the larger forward/reverse count per dataset and
record; partner selection at read boundaries can reduce matches even for an
exact reference substring. Counts retain repeated seed occurrences.

Membership is approximate: absent seeds can match. Fingerprint false-positive
probabilities are approximately `2^-f` for an `f`-bit ZOR fingerprint; the
shared 32-bit abandoned-key hash map can also add matches through collisions.
A sequence-level ratio threshold is not a calibrated read classification
probability. Dropping abandoned-key recovery can additionally lose present seeds.

## Modes and defaults

| Setting | Default or supported values |
|---|---|
| Indexed features | `kmers`; opt into `minimizers`, `modimizers` or `findere` |
| Zindex minimizer selector | `open-closed-mod`; `simd` is also available |
| Bindex minimizer selector | `simd` |
| Zindex fingerprint width | 8 bits; choices: 1, 2, 4, 8, 16 |
| ZOR hash locations | 4 |
| Layout planner | `streaming`; `all-pairs`/`full`, `streaming-unpruned`, `naive` are available |
| Sketch | 128 buckets; 8-bit comparison fingerprints and 32-bit union minima |
| Sizing / planning margin | `--slot-scale 1.0` / `--merge-plan-margin 0.03` |
| Cycle breaking | `no-heuristic`, `--tie-scan 1`, fixed 64-cell segments |
| Column ordering | `compression` |
| Minimum documents per size class | 128; smaller collections use one class |
| Maximum uncompressed subindex | 8 GiB; larger classes are hash-partitioned |
| Query loading | One subindex at a time; `--query-memory all` loads all |
| Completed construction data | Staged on disk; `--low-disk` keeps it in RAM |
| Abandoned-key recovery | Enabled; one shared map of 32-bit hashes |

The subindex byte cap limits the uncompressed filter payload, not total process
RAM: query seeds, results, recovery matches, metadata and worker buffers also
use memory. Oversized input documents can be skipped by the default build
guard; inspect `info` and the skipped-document report to confirm what was indexed.
The [usage guide](docs/USAGE.md) explains the memory controls and their limits.

Both indexers support odd minimizer windows up to 255 bases with selected seeds
up to 31 bases. Dense k-mers and modimizers use lengths up to 31; findere uses
`1 <= k-z <= 31`. Spaced seeds, randstrobemers and homopolymer compression are
Zindex options. The standalone generator defaults to SIMD minimizers; select
`--minimizer-scheme open-closed-mod` to match Zindex's default selector.

## Documentation and reproducibility

- [Detailed usage and implementation notes](docs/USAGE.md)
- [Collaborator commands and matched Bindex comparison](docs/COLLABORATOR.txt)
- [Benchmarking and reproducibility](docs/BENCHMARKING.md)
- [Manuscript experiment artifact status](docs/REPRODUCIBILITY.md)
- [Index formats and compatibility](docs/FORMAT.md)
- [Standalone minimizer utility](utils/minimizer-fasta/README.md)
- [Contributing and reporting problems](CONTRIBUTING.md)

The default merged directory supports `build`, `query`, `info`, and
`repack-subindexes` for lossless compression reframing. `append` and ordinary
`repack` apply to legacy `--plain` indexes. Existing indexes keep their stored
extraction settings; rebuilding applies new layout choices.

## Tests

```bash
cargo fmt --all -- --check
cargo test --workspace --release --locked
```

CI runs release tests for native and portable x86-64 builds, plus the small
end-to-end demonstration. Performance benchmarks are separate from correctness
tests and are not run by the demonstration.

## Citation and license

Use [CITATION.cff](CITATION.cff) to cite the software with the exact revision or
release used. The underlying filter method is:

Antoine Limasset. *ZOR Filters: Fast and Smaller Than Fuse Filters.*
SEA 2026, LIPIcs 371, 24:1–24:17.
[doi:10.4230/LIPIcs.SEA.2026.24](https://doi.org/10.4230/LIPIcs.SEA.2026.24).

The project is distributed under the [MIT license](LICENSE). Vendored Helicase
and SimdSketch include their upstream licenses and provenance notices.
