# Changelog

## Unreleased

This section describes the local candidate for review. It is not a published
release, and no existing release tag has been moved.

- Pin the tested Rust 1.97.0 toolchain for local and CI builds, avoiding the
  Rust 1.99.1 native-target LLVM failure observed in the Bzip2 dependency.
  Make the findere worker-participation test independent of scheduler speed.

- Integrate Helicase 1.2.0, propagate parsing errors and retain seeds crossing
  CRLF-wrapped sequence lines. Preserve bounded streaming in minimizer-fasta.
- Default merged queries to one resident subindex, extracting queries and
  scanning abandoned keys once and retaining results in RAM.
- Add a default 8-GiB uncompressed subindex cap with hash-routed partitions;
  preserve fingerprint width and merged document membership.
- Add construction memory controls and explicit reports for skipped documents.
- Add independent Zstd frames and lossless `repack-subindexes` for parallel
  loading, plus additional planner and query execution improvements.
- Accept odd minimizer windows up to 255 in Bindex, including loading those
  indexes; dense k-mer/modimizer keys remain limited to 31 bases.
- Add a checked small demonstration, build-revision reporting, corrected
  documentation, and a one-run benchmark artifact recorder.

Subindex partitioning changes filter cells and approximate false-positive
identities. Memory guards can exclude oversized input documents; their reports
must be checked before treating an index as the complete requested collection.
Feature extraction identities and index format versions are described in
[FORMAT.md](docs/FORMAT.md). Package version numbers are separate from those
format versions. Build output reports the Git revision and dirty status so two
builds sharing package version `0.1.0` can be distinguished.

## v0.1.0 — 2026-10-23

Published tag at commit `9ca1789f2de7922893def6766066c4b4cef91727`.
This tag predates the Helicase 1.2.0 update and the local query-memory/subindex
changes described above. Its historical behavior is retained by the tag;
unreleased changes should not be attributed to that release.
