# Small functional example

Run `bash utils/scripts/demo.sh` from the repository root after
`cargo build --release --locked`. No genome downloads are required.

The three references each contain 8,192 synthetic DNA bases. `a` and `b` share
their first 4,096 bases and have independent suffixes; `c` is independent.
The checked-in bases are deterministic fixtures, not biological genomes.
They were generated with a 32-bit LCG (initial state 4271, multiplier 1664525,
increment 1013904223), selecting A/C/G/T from the top two bits after each update.
Generation order is a, b's independent suffix, c, then the negative query.

Queries are 1,000 bases: one from the shared prefix (offset 1000), one from
each reference's suffix (offset 6000), the reverse complement of `a_only`,
and an independent negative sequence. `expected-hits.tsv` gives the expected
dataset names at ratio 0.8 (`-` means no datasets). The demonstration also checks that each reported
positive has all its selected minimizers matched. These expectations are for
these fixed fixtures, not a guarantee that approximate indexes never return
false positives on other data.

`samples.fof` uses paths relative to the repository root. The demonstration uses
native open-closed mod minimizers for Zindex and SIMD minimizers for Bindex;
these select different seed sets. See the [comparison commands](../../../docs/COLLABORATOR.txt) for a comparison
that deliberately uses the same SIMD selector in both tools.

Output files are generated outside this example directory. The script preserves
its logs and outputs and refuses to reuse an existing destination.
