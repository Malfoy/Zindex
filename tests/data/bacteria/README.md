# Bacterial example

Ten real RefSeq bacterial assemblies are bundled under `genomes/`.
They were selected from [the accession list](../../../docs/benchmarks/accessions_128_bacterial.txt)
with random seed 20261102. Accession versions, original NCBI URLs and checksums
are recorded in [assemblies.json](assemblies.json) and [sources.json](sources.json).

## Run the example

The [README](../../../README.md#build-and-query) contains the build and query
commands. `query.fa` is a 1-kb fragment from the first assembly. `queries.fa`
contains a fragment from each genome, its reverse complement and a random
negative control. `source-truth.tsv` records each positive's source dataset,
contig, offset and strand. Related assemblies can also match.

`bash utils/scripts/demo.sh` verifies every positive against its source in both
k-mer and minimizer modes.

## Retrieve the original files again

```bash
python3 utils/scripts/prepare_bacteria_example.py --output /tmp/bacteria-originals
```

Downloads are checked against pinned NCBI checksums. A changed file is rejected.
