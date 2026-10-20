# SimdSketch provenance and local changes

Vendored from Ragnar Groot Koerkamp's published simd-sketch 0.3.0 crate,
upstream git revision a01c46dfaa9b0dc83edf9c6c74f6b371bcb967db:
https://github.com/RagnarGrootKoerkamp/simd-sketch

The upstream crate declares the MIT license. Its library and SIMD compaction
implementation are preserved here, with these explicit changes:

- Feature-ID adapter in src/features.rs: sketch the actual indexed features,
  not a different raw sequence/k-mer universe. SplitMix64 hashes are folded to
  32-bit hashes for the upstream mod-bucket representation (upper 32 bits of
  the mixed feature ID, without reserving any hash or fingerprint value).
- Retain full 32-bit bucket minima alongside 8-bit fingerprints for exact
  unions; compare fingerprints only. The indexed-feature adapter assumes dense
  genomes: it has no empty-cell masks or sparse/empty similarity handling.
- Runtime-dispatched AVX2 collection, byte comparisons and in-place unions,
  plus scalar fallbacks. Dense byte comparisons have no occupancy checks and
  adjust Jaccard for accidental 8-bit matches by (256*p - 1)/255, clamped at zero.
- Clone/equality derives for planner regression tests; byte comparisons route
  through the new SIMD kernel. The feature adapter calls the dense kernel
  directly, without the general sequence API's sparse-sketch dispatch.
- Read-only fingerprint access and an integer-match-count similarity helper
  support ZORINDEX's packed transposed planner arena, preserving the original
  estimator's exact f32 rounding. The transposed SIMD scan lives in ZORINDEX,
  not in the upstream SimdSketch implementation.
- Correct the generic scalar compaction fallback to accept u32 output.
- Library-only dependency manifest (upstream CLI dependencies are omitted).

BinDash 2.0 was examined for b-bit comparisons and SIMD/popcount techniques:
https://doi.org/10.1101/2024.03.13.584875
https://github.com/zhaoxiaofei/bindash/blob/master/src/genome.hpp
No BinDash source is copied. For 8-bit values, AVX2 byte equality plus
movemask/popcount avoids transposing to eight separate bit planes. No BinDash
or Rust BinDash2 implementation is used. Full minima remain mergeable without
densification; similarity comparisons assume populated buckets. The upstream
general sequence-sketching API remains available separately from this adapter.
