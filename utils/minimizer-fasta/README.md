# minimizer-fasta

Independent FASTA preprocessing utility from the Zindex repository. It emits
one canonical DNA seed per selected **position**, not one per distinct seed
value. Repeated values at different positions remain present.

The default selection is **canonical SimdMinimizers 2.3.1**, matching the original
`benchmark-simd-minimizers` / gzip benchmark generators. Select
`--minimizer-scheme open-closed-mod` (alias `ocm`) to use the same optimized
canonical open-closed mod iterator as native Zindex. Both tools share the
dependency-free `ocm-minimizers` crate; there is no separate copied iterator to
drift out of sync. Existing benchmark inputs and ongoing experiments are not
changed by installing it.

## Build

From the repository root:

```bash
cargo build --release -p minimizer-fasta
```

Only this utility and its dependencies are built; the index implementation is
not a dependency. The workspace release profile uses optimization level3,
fat LTO, one codegen unit, stripped symbols and panic=abort. The repository's
`.cargo/config.toml` sets `target-cpu=native`, enabling AVX2 or NEON on supported
hosts. Rebuild on each target machine: a native binary need not run on older CPUs.
The JSON summary reports `AVX2`, `NEON` or `scalar`; builds without either SIMD
feature use the dependency's portable scalar fallback.

Alternatively install the utility, explicitly keeping native CPU optimization:

```bash
RUSTFLAGS='-C target-cpu=native' cargo install --path utils/minimizer-fasta --locked
```

## Single file

```bash
./target/release/minimizer-fasta input.fna.gz output.fa.gz
./target/release/minimizer-fasta --k 21 --window 19 --threads 32 input.fa output.fa
./target/release/minimizer-fasta --minimizer-scheme open-closed-mod input.fna.xz output.fa.gz
./target/release/minimizer-fasta --minimizer-scheme open-closed-mod --seed 69 input.fna.zst output.fa
gzip -dc input.fa.gz | ./target/release/minimizer-fasta - - > output.fa
```

`--k` is the **seed length** (default31); `--window` / `-w` is the number of
candidate k-mers per window (default19). The total span is `k + window - 1`
(49 bases by default). The span must be odd for canonical strand tie-breaking.
Supported seed lengths are1..64 for `simd`, 1..31 for `open-closed-mod`, and
window counts1..65535. `--seed` defaults to69 and controls OCM's hash ordering;
the legacy SIMD selector keeps its original hashing regardless of this setting.
For example, utility `--k 31 --window 19 --minimizer-scheme open-closed-mod`
matches Zindex `--k 49 --minimizer-size 31 --minimizer-scheme open-closed-mod`.
No seed is selected
from an ACGT run shorter than the full window span. Upper/lowercase and sequence
whitespace are normalized; non-ACGT bases split runs. Record and ambiguity
boundaries are never crossed. FASTQ is intentionally rejected: this is a FASTA
utility, matching the original generator.

Headers retain the original generator's `>r1_p123` convention: one-based record
number and zero-based base position in the whitespace-normalized record.
ASCII lexical minimum of the seed and its reverse complement is emitted.
The `eligible_kmer_positions` statistic counts every valid k-mer position,
including positions in runs too short for a full selection window.

Input compression is detected from magic bytes: plain, gzip (including multiple
members), zstd, xz and bzip2 are supported. `-` supports stdin/stdout. Named
output files are **never overwritten**. An I/O error returns a nonzero status;
any partial output remains available for diagnosis and is not silently reused.

## Collections

```bash
./target/release/minimizer-fasta --fof samples.fof --output-dir minimizers --gzip
./target/release/minimizer-fasta --minimizer-scheme open-closed-mod \
  --fof samples.fof --output-dir ocm-minimizers --gzip
```

The FOF accepts `path` or `name<TAB>path`; blank/comment lines are ignored.
Relative input paths are resolved against the current working directory.
Files are processed concurrently using the same affinity-sized Rayon pool as
within-file chunk processing, not a nested worker pool for every genome.

Outputs are `00000000.minimizers.fa.gz`, etc. (`.fa` without `--gzip`), avoiding
collisions between inputs with identical basenames. `minimizers/minimizers.fof`
retains input names/order and contains absolute paths to the completed files.
The manifest is only written after all conversions succeed. Existing output
files or a manifest cause an error instead of being overwritten or resumed.

## Performance and compression

- Defaults to all CPUs available to the process; `--threads` overrides this.
- SIMD selection uses Ragnar Groot Koerkamp and Igor Martayan's SimdMinimizers.
- OCM selection uses Zindex's four-lane AVX2 kernel, full-width 64-bit hashes,
  prefix/suffix sliding minima, and a portable fallback. OCM writes positions
  directly into the reusable deduplication bitset, avoiding a packed-DNA copy
  and intermediate position vector.
- FASTA parsing uses Helicase's SIMD block lexer, not its record-buffering API.
  DNA spans are normalized and counted in bulk, without retaining complete
  chromosomes. `vendor/helicase` contains upstream0.1.1 with only the lexer API
  made public; parsing algorithms are unchanged. Magic-based gzip/xz/zstd/bzip2
  decoding also applies to stdin and the library `generate` API. Decoder errors
  are returned rather than triggering a parser panic or silently accepting EOF.
- Long chromosomes are streamed in bounded chunks, with both left and right
  context so every selected position exactly matches full-record selection.
- A reusable bitset replaces comparison sorting for exact sorted position
  deduplication. Distinct positions are never collapsed by k-mer value.
- Worker-local packed DNA, position/bitset arrays and output buffers are reused.
  Canonicalization uses a stack buffer, not a heap allocation per seed.
- Dedicated parsing/writing threads and bounded channels overlap I/O and compute;
  full input records, full output files and collection-wide seed sets are not
  retained. `--chunk-bases` defaults to262144; context is added automatically.
- Output gzip is compressed in parallel with the accelerated zlib-rs backend,
  using standard **concatenated gzip members** in original chunk order.
  Consumers must read all members (e.g. `gzip -dc`, Python `gzip`, Rust
  `MultiGzDecoder`). A single-member-only `GzDecoder` will stop too early.
- `.gz` / `.gz.partial` enables gzip automatically; `--gzip` forces it.
  `--gzip-level` defaults to1 and accepts0..9. The original benchmark helper
  used level6; compressed bytes/sizes may differ, but decompressed FASTA content
  and selected positions are unchanged.
- A normal single-member gzip input still has one sequential decompressor per
  file; the utility does not claim parallel inflation of one DEFLATE stream.

JSON statistics go to stdout, or stderr when FASTA itself is written to stdout.
They include `minimizer_scheme` and `seed` alongside the previous fields.
Output is deterministic across thread counts for identical parameters. Changing
chunk size may change gzip member boundaries, but never decompressed FASTA.

Tests compare the chunked bitset implementation to the original full-run
SIMD+sort+dedup algorithm, including repeated ties, reverse complements,
ambiguous bases, short/empty records, tiny chunks and long unwrapped records.
They also compare OCM output against full-record iterator positions for several
window/seed lengths, hash seeds, strands, periodic inputs and chunk sizes, and
cover all supported compression formats, compressed stdin, collection ordering,
truncated streams, refusal to overwrite and pipeline cancellation on I/O errors.

```bash
cargo test --release -p minimizer-fasta
```
