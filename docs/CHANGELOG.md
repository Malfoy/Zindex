# Changelog

## Unreleased

- Default per-record queries to named TSV rows with dataset IDs, names, counts
  and ratios. Keep compact output selectable with `--output-format compact`.
  Send aggregate query diagnostics to stderr.
- Require every input dataset by default. Keep memory guards active, make
  exclusions opt-in, and report collection completeness in build summaries
  and `info`.
- Consolidate usage in the README and remove fusion's `--other` option.

- Order equal-width columns by the retained union-sketch fingerprints, removing
  filter-cell sampling. Free planning minima, share ordering keys across rebuilt
  filters and partitions, and record the ordering method in the manifest.

- Remove the capacity-floor, column-order, dataset-cap and minimum-saving CLI
  controls. Select column order internally, derive the group limit from the
  64-bit payload width, and require strictly positive predicted merge savings.

- Remove the query-time recovery opt-out. Queries automatically use recovery
  whenever it is stored in the index.

- Add the paper workflow figure and a build-to-query introduction to the README.
  Split regular help (`-h`) from advanced help (`-ah`), keeping ZOR heuristics
  in advanced help. Remove the unused aggregate-output column.

- Add `fusion -i FIRST SECOND [MORE ...] -o OUTPUT` for compatible index collections.
  Preserve compressed filter bytes and source match counts, stream recovery
  postings with remapped dataset IDs, and support repeated fusions and different
  subindex partition geometries in manifest version 10.
  Use reflinks where available or combine validation and copying in one pass.
  Copy independent subindexes with the requested worker pool and merge all
  recovery streams once, avoiding intermediate pairwise indexes. Move bank
  metadata and reuse recovery buffers.

- Remove the single-file index backend, row-codec experiments, maintenance
  commands and the unpruned streaming CLI mode. Each remaining option has
  one long name and one short flag.

- Use `--k` for the window span in minimizer mode and `--minimizer-size`
  (`-m`) for the selected sequence length. Use `-i` for `--index-mode`.
  Group help into mandatory, important and advanced options. Simplify the README
  and show the interleaved Bloom versus ZOR figure. Existing indexes are readable.

This section describes development changes on main that are not part of a
tagged release. Existing published release tags have not been moved.

- Simplify the landing page and consolidate the repository into five main
  directories: src, tests, utils, docs and vendor. Move fixtures under tests/data,
  helper scripts under utils/scripts, and build support under .cargo.
- Bundle ten real bacterial assemblies and source-mapped example queries.

- Pin the tested Rust 1.97.0 toolchain for local and CI builds, avoiding the
  Rust 1.99.1 native-target LLVM failure observed in the Bzip2 dependency.
  Make the findere worker-participation test independent of scheduler speed.

- Integrate Helicase 1.2.0, propagate parsing errors and retain seeds crossing
  CRLF-wrapped sequence lines. Preserve bounded streaming in minimizer-fasta.
- Default merged queries to one resident subindex, extracting queries and
  scanning abandoned keys once and retaining results in RAM.
- Add a default 8-GiB uncompressed subindex cap with hash-routed partitions.
  preserve fingerprint width and merged document membership.
- Add construction memory controls and explicit reports for skipped documents.
- Add independent Zstd frames for parallel loading, plus additional planner
  and query execution improvements.
- Add a checked small demonstration, build-revision reporting, corrected
  documentation, and a one-run benchmark artifact recorder.

Subindex partitioning changes filter cells and approximate false-positive
identities. Memory guards can exclude oversized input documents. Their reports
must be checked before treating an index as the complete requested collection.
Feature extraction identities and index format versions are described in
[FORMAT.md](FORMAT.md). Package version numbers are separate from those
format versions. Build output reports the Git revision and dirty status so two
builds sharing package version `0.1.0` can be distinguished.

## v0.1.0 — 2026-10-23

Published tag at commit `9ca1789f2de7922893def6766066c4b4cef91727`.
This tag predates the Helicase 1.2.0 update and the local query-memory/subindex
changes described above. Its historical behavior is retained by the tag.
unreleased changes should not be attributed to that release.
