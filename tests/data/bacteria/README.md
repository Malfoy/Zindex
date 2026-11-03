# Bacterial example

Ten real RefSeq bacterial assemblies are bundled under `genomes/`.
They were selected from [the accession list](../../../docs/benchmarks/accessions_128_bacterial.txt)
with random seed 20261102. Accession versions, original NCBI URLs and checksums
are recorded in [assemblies.json](assemblies.json) and [sources.json](sources.json).

## Build and query

From the repository root:

```bash
./target/release/zorindex build --fof tests/data/bacteria/samples.fof \
  --output bacteria.zindex --k 31
./target/release/zorindex query --index bacteria.zindex \
  --query tests/data/bacteria/query.fa
```

`samples.fof` contains one named assembly per line. `query.fa` contains a
1-kb fragment from the first assembly.

## Use minimizers and multiple queries

```bash
./target/release/zorindex build --fof tests/data/bacteria/samples.fof \
  --output bacteria.min.zindex --k 31 --index-mode minimizers --window 19
./target/release/zorindex query --index bacteria.min.zindex \
  --query tests/data/bacteria/queries.fa --records --min-ratio 0.8 > matches.tsv
python3 utils/scripts/prepare_bacteria_example.py --check-results matches.tsv \
  --truth tests/data/bacteria/source-truth.tsv
```

`queries.fa` contains a 1-kb fragment from each genome, its reverse complement,
and a random negative control. `source-truth.tsv` records each positive's
source dataset, contig, offset and strand. Related assemblies can also match.
the checker verifies that every positive matches its source.

The optional checker runs both index modes in a new temporary directory:

```bash
bash utils/scripts/demo.sh
```

## Retrieve the original files again

```bash
python3 utils/scripts/prepare_bacteria_example.py --output /tmp/bacteria-originals
```

Downloads are checked against pinned NCBI checksums. A changed file is rejected.
