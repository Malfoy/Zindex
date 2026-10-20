# General-case execution optimizations

These changes preserve the streaming merge policy, deterministic tie order,
feature hashes, exact union keys and origin masks, filter equations, fingerprint
widths, slot capacities, bank encoding and abandoned-key map semantics.

## Streaming planning

The current incumbent is published when a scan finds a better candidate. Other
workers observe a revision counter and refresh their match-count rejection
threshold using the same rounded Jaccard estimator and power-of-two cost rule.
The mutex is used only for incumbent publication or an observed revision, not
for each fingerprint comparison. Each physical sketch block also retains a
conservative minimum logical group ID. Once the cost upper bound cannot win,
blocks can be skipped without losing an earlier-ID tie. Removed IDs may leave
stale-low hints; these reduce pruning efficiency but cannot change decisions.
`planning_pair_evaluations` counts visited candidate pairs, including partial
SIMD scans, rather than the candidates skipped by the stronger bounds.

## Extraction and exact unions

Large extracted u64 arrays use in-place byte-radix partitioning; partitions
smaller than 2,048 values use the standard comparison sort. Deduplication still
produces exactly the ascending feature set. No feature-sized radix scratch
array or temporary file is allocated.

The scratch-free builder accumulates sorted colored runs using balanced binary
merges rather than repeatedly merging into one growing left-hand union. Small
runs are combined first. Similar genomes would otherwise leave many nearly
identical runs, so retained entry counts are collapsed when they exceed twice
the largest run. This is a bound on retained entries, not a strict bound on
peak allocator bytes: vector capacity and active merge workspaces also exist.
Every merge ORs the original document bits and keeps ascending key order.

## Construction cells

Assignment and verification use u8/u16/u32/u64 cells according to the complete
fingerprint-plus-origin payload width. The existing u64 assignment remains as
a correctness reference and for legacy callers. Eight-bit singleton cells can
be transferred directly into the completed packed filter without repacking.
Other widths use the same final bit encoding. Only construction buffers change;
graph indices and the persisted bank/recovery formats are unchanged.

## Per-record query pipeline

Up to 16 distinct weighted keys have their bank positions prepared together.
Four keys are prefetched ahead, hinting only the start of independent row
streams rather than every cache line. Cache results, scores and threshold
pruning are committed strictly in the original key order. Prepared addresses
do not speculate about future candidate survival or mutate query state.
The bounded long-record fallback remains available.

## Validation and measurements

The dedicated experiment is
`experiments/general_case_optimizations_20261020-jfh4po/`.
It freezes both source trees and release binaries and uses the retained real
4,096-genome bacterial collection from `../benchmarkzindex`: native OCM k49/m31,
128 sketch cells, h8/fp8, slot scale 1.0, planning margin zero, recovery enabled,
and all 32 runtime-affinity threads. It uses one build per version and
one query process per version/class using the same physical index and
32,768 reads of 2,000 bases per class. Query wall time includes loading and I/O.
Historical measurements are not overwritten or relabeled. The initially chosen
repetition policy was cancelled by the user; the primary comparison reuses only
the first completed build of each version, without rebuilding them.

The single-run comparison measured construction wall time 161.06 -> 144.83 s
(-10.08%) and peak construction RAM 5.304 -> 5.035 GiB (-5.07%). Query acceleration
was not demonstrated: positive-query wall time was unchanged at 1.94 s and
negative-query wall time increased from 1.80 to 1.93 s. Planning was already
negligible and increased from 0.060 to 0.074 s despite fewer visited candidates.
All serialized filter/recovery payloads and query results remained identical.

Canonical manifests, bank bytes, abandoned-map bytes and query outputs must
agree. Full counters are checked at threshold zero, with recovery enabled and
disabled. Positive source-document coverage is checked separately; random
negative reads are a throughput workload, not an individual-feature FPR
calibration. The experiment retains commands, raw measurements, manifests,
payload identities and query archives, and removes validated test indexes.
