# Benchmarking

Start with `bash utils/scripts/demo.sh` to verify installation. Keep real datasets,
indexes and measurement directories outside the versioned source tree.

## One-run measurement runner

On Linux with Python 3 and GNU `/usr/bin/time`:

```bash
cargo build --release --locked
python3 utils/scripts/benchmark_zindex.py \
  --fof /data/samples.fof --positive /data/positive.fa.gz \
  --negative /data/negative.fa.gz --threads 32 \
  --output /results/zindex-one-run
```

The output path must not already exist. The runner freezes a copy of the binary
and the ordered FOF with absolute paths, then performs **one build and one query
for each workload**. It never retries or repeats a point. Failures retain their
commands, logs and status instead of being silently omitted. Inputs are supplied
by the caller. No genomes are downloaded.

Defaults are OCM minimizers k31/w19, h4/fp8, seed 69, 128 sketch buckets,
slot-scale 1.0, planning margin 0.03, streaming planning, compression column
order, no-heuristic/tie-scan 1 and fixed 64-cell segments. The full build command
and actual manifest are retained. Parameter options include hash/fingerprint
width, planner, maximum subindexes, minimum documents, maximum subindex bytes,
and query memory mode. See `python3 utils/scripts/benchmark_zindex.py --help`.

The runner disables skipping oversized documents and checks the final document
count. A failed full-collection build is a failed point, not a smaller benchmark.
This does not turn memory admission estimates into hard memory guarantees.

Query reporting is per-FASTA-record at ratio 0.8 by default. Findere uses aggregate
reporting because its per-record mode is unsupported. Do not compare these two
reporting modes as though the scoring workloads were identical.

## Captured evidence and timing boundaries

The new result directory contains:

- `results.json`: host, CPU affinity, binary identity/version, settings,
  process wall time, user/system CPU time, peak RSS and index bytes.
- `documents.jsonl` and `inputs.fof`: ordered document identities and paths.
- `manifest.json`: the actual final index layout and extraction settings.
- Stage command JSON, complete stdout/stderr and GNU time reports.
- The frozen executable and the final index, which are retained for validation.

All times include startup, index loading/decompression, seed extraction, lookup,
scoring, output and shutdown. Internal `elapsed_s` counters have different
boundaries. Never substitute them for end-to-end time without labeling the
difference. RSS is GNU time's peak KiB converted to bytes. GB and GiB are distinct.

FOF, query files, executable and manifest receive SHA-256 identities. Genomes
use path/size/mtime identities by default. `--hash-inputs` additionally reads and
hashes their complete compressed bytes **before** measurement. `--hash-index`
hashes all final index files **after** the queries. These can be substantial
extra reads at repository scale, and are outside the measured stages.
Path/size/mtime records alone are not content-addressed input snapshots.

This runner does not drop system caches, enforce isolation from unrelated jobs,
or claim cold-cache measurements. Construction and earlier reads affect cache
state. Record that policy alongside any reported results. Run configurations
sequentially on an otherwise idle machine.

## Accuracy and historical experiments

The runner preserves query output. It cannot infer source-document truth from
arbitrary read names. Validate every positive against an explicit source-ID map,
and report skipped inputs, changed inputs, failures and incomplete points.
Random negatives are useful throughput workloads, not a standalone calibration
of per-feature false-positive probability.

Historical measurements retain their original settings and timing definitions.
See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for the current artifact inventory
and missing material needed to reproduce the manuscript's large collections.
The bundled example and this runner do not reproduce unpublished experiments.

## Legacy slice-compression experiments

`--stack-compression`, ordinary `repack` and `--decode-stack` concern the
legacy `--plain` format. Merged subindexes instead use Zstd-compressed slices.
To study a legacy slice codec:

```bash
./target/release/zorindex build --plain --fof samples.fof --k 31 \
  --index-mode minimizers --window 19 --output /results/base.zoridx
./target/release/zorindex repack --index /results/base.zoridx \
  --stack-compression lz4-hc --output /results/base.lz4hc.zoridx
```

Compare stored-query and decoded-query time separately, include loading in
end-to-end measurements, and retain repack time and both on-disk sizes.
