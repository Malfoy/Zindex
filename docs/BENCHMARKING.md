# Benchmarking

Keep benchmark datasets, generated indexes, and logs out of the source tree.
The `.gitignore` is configured for common FASTA/FASTQ files, generated index
files, paper PDFs, and external tool outputs.

## Recommended Build

```bash
cargo build --release
```

The repository release profile and `.cargo/config.toml` already enable local
native-code optimization. For controlled cross-machine comparisons, remove or
override `-C target-cpu=native`.

## Dataset File

Create a file-of-files:

```text
sample_0001	/path/to/genome_0001.fa.zst
sample_0002	/path/to/genome_0002.fa.zst
```

## Zindex Build

```bash
/usr/bin/time -v target/release/zorindex --threads "$(nproc)" build \
  --fof samples.fof \
  --k 63 \
  --index-mode minimizers \
  --minimizer-size 21 \
  --hashes 4 \
  --fingerprint-bits 8 \
  --cycle-break most-deg2 \
  --tie-scan 8 \
  --stack-compression none \
  --output results/samples.zoridx
```

Compressed-stack example:

```bash
/usr/bin/time -v target/release/zorindex --threads "$(nproc)" repack \
  --index results/samples.zoridx \
  --stack-compression lz4-hc \
  --output results/samples.lz4hc.zoridx
```

## Bindex Build

For a target false-positive rate of `1/256`:

```bash
/usr/bin/time -v target/release/bloomindex --threads "$(nproc)" build \
  --fof samples.fof \
  --k 31 \
  --false-positive-rate 0.00390625 \
  --hashes 6 \
  --output results/samples.blmidx
```

For `1/65536`:

```bash
/usr/bin/time -v target/release/bloomindex --threads "$(nproc)" build \
  --fof samples.fof \
  --k 31 \
  --false-positive-rate 0.0000152587890625 \
  --hashes 11 \
  --output results/samples.fp16.blmidx
```

Sweep `--hashes` for Bloom comparisons. Bindex will resize the Bloom bit count
from the target false-positive rate and selected hash count unless `--bit-count`
or `--bits-per-item` is supplied.

## Query Throughput

Use the same query file for every index:

```bash
/usr/bin/time -v target/release/zorindex --threads "$(nproc)" query \
  --index results/samples.zoridx \
  --query queries.fa.zst \
  --min-ratio 0 > results/zor.query.tsv

/usr/bin/time -v target/release/bloomindex --threads "$(nproc)" query \
  --index results/samples.blmidx \
  --query queries.fa.zst \
  --min-ratio 0 > results/bloom.query.tsv
```

Record at least:

- index size on disk
- build wall time
- build peak RSS
- query wall time
- query peak RSS
- reported feature count
- reported stack or Bloom probes

For compressed Zindex tests, compare:

- `query` against compressed storage
- `query --decode-stack` against a decoded in-memory stack
- `repack` time and output size
