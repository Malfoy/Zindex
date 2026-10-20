//! Ragnar Groot Koerkamp's SimdSketch with feature-ID and exact-union extensions.
//! Compare 8-bit fingerprints; keep 32-bit minima only while planning unions.
//! See vendor/simd-sketch/NOTICE.md for upstream provenance and local changes.
pub use simd_sketch::features::MergeableSketch as Sketch;
