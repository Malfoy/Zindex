# Reproducing results

The small demonstration and the measurement runner are available in this
repository. A complete public artifact bundle for the manuscript's RefSeq and
AllTheBacteria experiments has **not yet been assembled**. The demo and current
defaults must not be described as reproducing those large-collection results.

## Available material

| Material | Location | What it establishes |
|---|---|---|
| Fixed tiny inputs and expected hits | `tests/data/tiny/`, `utils/scripts/demo.sh` | Installation and functional behavior |
| Correctness tests | `tests/`, workspace unit tests | Covered extraction, storage, query and CLI combinations |
| One-run measurement runner | `utils/scripts/benchmark_zindex.py` | Reproducible commands, timing boundaries and artifact capture on a supplied collection |
| Historical optimization summaries | `docs/OPTIMIZATION_BENCHMARK.txt`, `docs/GENERAL_CASE_OPTIMIZATIONS.md` | Recorded measurements of their explicitly identified historical versions |
| Small bacterial accession list | `docs/benchmarks/accessions_128_bacterial.txt` | A list of 128 accessions, not the manuscript collections |
| Optional three-genome example | `tests/data/bacteria/`, `utils/scripts/prepare_bacteria_example.py` | Verified downloads and source-mapped forward/reverse queries |

Historical summaries refer to local, git-ignored experiment directories. Those
paths identify provenance; they are not downloadable artifacts. Their frozen
settings can differ from current defaults, including sketch size, sizing margin,
hash count, segment geometry and whether query time includes loading.

## Material still needed for manuscript-scale replication

For each RefSeq/AllTheBacteria result and seed/error-rate experiment, the public
bundle needs:

1. A named experiment and the exact code commit or immutable source archive,
   Cargo.lock, compiled binary checksum, compiler and CPU configuration.
2. The collection snapshot date, accession versions, retrieval URLs or commands,
   final ordered FOF and file checksums. Repository names and genome counts alone
   do not specify a reproducible input collection. Removed/unavailable accessions
   must be recorded rather than silently replaced.
3. Exact build and query commands, all non-default settings, observed included
   document counts and any skipped-document reports.
4. Query sequences or generation scripts, source-document truth, read lengths,
   random seeds and complete substitution/insertion/deletion models. Platform
   labels alone do not specify a simulator distribution.
5. Comparator versions, commands and equivalent extraction/threshold settings;
   native Zindex OCM and Bindex SIMD defaults do not select identical seeds.
6. Raw timing/RSS reports, index manifests and sizes, error exits and interrupted
   points, plus machine-readable plot data and plotting scripts.

Large input files can be distributed through a stable external archive; the
repository should provide the archive identifiers, checksums and retrieval
instructions. No archive URL, DOI, paper citation or missing historical command
is invented here. For publication replication, identify an exact commit or an
immutable tagged release together with its experiment artifacts; main can change.

## Measurement conventions

The [benchmarking guide](BENCHMARKING.md) describes the supplied runner. Each
configuration runs once: one construction and one positive/negative query
process. Query wall time includes process startup, index loading/decompression,
seed extraction, lookup, scoring, output and shutdown. Internal query timers can
exclude loading and must be labeled separately. The runner makes no claim that
filesystem caches are cold and does not drop system caches.

Random negative reads measure throughput; they do not independently calibrate
per-feature false-positive rates. Positive source recall requires a source-ID
truth file and validation of every expected source hit. The generic runner
preserves outputs but cannot infer that truth from arbitrary read names.

Do not rerun an old table using current defaults and present it as the same
experiment. Preserve unsuccessful points and their failure reasons. Benchmark
outputs and genomes belong outside the versioned source tree.
