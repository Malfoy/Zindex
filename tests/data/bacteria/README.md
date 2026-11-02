# Optional real bacterial example

This example downloads three complete RefSeq assemblies: E. coli K-12 MG1655,
E. coli O157:H7 Sakai, and B. subtilis 168. Together they contain about 14.5 Mb
of sequence and require only a few megabytes of compressed downloads.
Python 3.9 or newer and internet access are required. The regular synthetic
demonstration and CI tests do not download these genomes.

From the repository root, after `cargo build --release --locked`:

```bash
python3 utils/scripts/prepare_bacteria_example.py --output /tmp/zindex-bacteria
./target/release/zorindex build --fof /tmp/zindex-bacteria/samples.fof \
  --output /tmp/zindex-bacteria/index --k 49 --index-mode minimizers \
  --minimizer-size 31 --fingerprint-bits 16 --skip-oversized-documents=false
./target/release/zorindex query --index /tmp/zindex-bacteria/index \
  --query /tmp/zindex-bacteria/queries.fa --records --min-ratio 0.8 \
  > /tmp/zindex-bacteria/hits.tsv
python3 utils/scripts/prepare_bacteria_example.py --check-results /tmp/zindex-bacteria/hits.tsv
```

The download directory must be new. Each compressed genome is checked against
the NCBI MD5 value pinned in [assemblies.json](assemblies.json). Full gzip
decoding verifies its trailer. `sources.json` additionally records SHA-256
identities, URLs, retrieval time, byte/base counts and original document IDs.
`samples.fof` is only published after all files and query records validate.
Failed downloads retain diagnostics but do not publish a usable FOF.

For each assembly, the first contiguous 1,000-base ACGT run supplies a forward
query and its reverse complement. `source-truth.tsv` records the assembly,
contig, zero-based start and strand. A seventh query is a deterministic uniform
random control. The checker requires every positive's selected minimizers to
match its source dataset; additional hits in related genomes are allowed.
The random control is reported, not used as a false-positive-rate calibration.

Pinned sources:

- [E. coli K-12 MG1655, GCF_000005845.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/)
- [E. coli Sakai, GCF_000008865.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000008865.2/)
- [B. subtilis 168, GCF_000009045.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000009045.1/)

Each assembly's NCBI directory exposes `md5checksums.txt`; the pinned checksums
were verified there on 2026-11-02. If an upstream file is unavailable or its
checksum changes, the script fails instead of silently substituting an assembly.
This is a small functional example, not the manuscript benchmark collections.
