# Index Format Notes

This document summarizes the binary formats used by the two command-line
indexes. The exact reader and writer implementations are the source of truth.

Both formats are written as zstd streams by default. The magic bytes described
below are the first bytes of the decompressed payload. Readers detect zstd frame
magic automatically and also accept raw legacy payloads.

## Zindex (`.zoridx`)

Magic:

```text
ZORIDX1\0
```

Current format version: `10` (readers also accept versions 2, 3, 4, 5, 7, 8 and 9).
Version 10 changes randstrobemer extraction to downstream-only construction
and separate two-orientation querying. Version-9 custom-randstrobemer files
are rejected with a rebuild message; other version-9 indexes remain readable.

Header fields:

- version
- `k`
- ZOR hash count
- fingerprint width (`1`, `8` or `16`); single-bit fingerprints are packed eight per byte, with each stack row rounded up to a whole byte
- flags
- seed
- feature mode
- minimizer size
- modimizer sampling
- findere z byte, only for findere mode
- minimizer scheme byte: 0 = SIMD, 1 = open-closed mod (version 8)
- homopolymer-compression byte: 0 = disabled, 1 = enabled (version 8)
- seed pattern width byte, followed by that many ASCII `#`/`_` bytes;
  width 0 disables masking (version 9)
- randstrobemer enabled byte: 0 = disabled, 1 = enabled; when enabled,
  followed by inclusive minimum and maximum start distances as little-endian
  `u64` values (version 9)
- ZOR segmented layout
- color metadata
- technical layout and technical-bin metadata (version 7)
- slot-major stack
- optional union graph

For 1-bit fingerprints, a plain stack row occupies
`ceil(technical_bin_count / 8)` bytes and a standalone filter (including a
union graph) occupies `ceil(array_length / 8)` bytes. Bits are least-significant
first within each byte. Slice codecs operate on these packed bytes with a
one-byte codec element width; the header fingerprint width remains `1`.
Unused final-byte bits are ignored by queries. The 8- and 16-bit layouts are
unchanged.

Version-9 anchor settings apply to construction, queries and plain-index
appends. Pattern syntax, symmetry and width, and randstrobemer distances are
validated on load. Older indexes have neither masking nor randstrobemers;
older readers reject version 9 rather than silently extracting different keys.

The stack can be stored plain or compressed. Plain stacks store all rows as
contiguous little-endian fingerprint bytes (or packed bits at width 1). Compressed stacks store row offsets
and codec payload bytes. Each compressed row expands to:

```text
ceil(technical_bin_count * fingerprint_bits / 8)
```

The query path computes row positions once per feature. In plain mode it XORs
directly from the interleaved rows. In compressed mode it decodes the selected
rows into scratch buffers and XORs those rows.

## Experimental Zindex Stack Codecs

The experimental stored stack-codec byte maps to:

```text
0  none
1  svb32 legacy
2  svb32-0124
3  bitpack legacy
5  lz4
6  snappy
7  fastpfor256
8  lz4-lib
9  lz4-hc
10 fastpfor-pack
```

Legacy codecs remain readable for old local indexes but are hidden from the
CLI when they are not useful for current experiments.

## Bindex (`.blmidx`)

Magic:

```text
BLMIDX1\0
```

Current format version: `2` (version 1 remains readable).

Header fields:

- version
- `k`
- Bloom hash count
- seed
- feature mode
- minimizer size
- modimizer sampling
- target false-positive rate
- common bit count
- color metadata
- technical layout and technical-bin metadata (version 2)
- bitstack words

Bindex stores `bit_count * color_words` little-endian `u64` words, where:

```text
color_words = ceil(technical_bin_count / 64)
```

For each Bloom bit position, all colors are packed into consecutive `u64`
words. Querying ANDs the packed words from all Bloom positions and increments
scores for the remaining set bits.

The default flat layout has one technical bin per output color. The optional
split-large layout maps several technical bins to one output color; queries
count each color at most once per feature.

## Zindex power-of-two merged directory

Zindex builds a new merged directory by default, regardless of the output
extension. `--plain` selects the legacy single-file `.zoridx` layout;
`--power-of-two-merge` remains an optional compatibility flag.
`manifest.json` contains version `8`, feature/hash settings, dataset identities,
group members and widths, bank sizes and abandoned-key map metadata.
Version 4 adds `minimizer_scheme` (`simd` or `open-closed-mod`) and
`homopolymer_compression` (boolean). These govern both construction and queries.
Missing legacy fields default to SIMD selection and compression disabled.
Version 5 adds `seed_pattern` (a symmetric `#`/`_` string or null) and
`randstrobemers` (null or an object with `min_distance` and `max_distance`).
Missing legacy fields default to null. Distances are inclusive, in sequence
bases between strobe starts, after optional homopolymer compression. The bank
and abandoned-key layouts are unchanged.

Version 6 adds `randstrobe_scheme`, currently
`strobealign-popcount8-symmetric64-downstream-v1`. It uses one downstream
seed per anchor, canonical-strobe XXH64, masked-XOR popcount with q=255,
first-minimum tie breaking, self-pair fallback and the paper's symmetric
`h1/2+h2/2` full-seed key. Query orientations are scored independently and
combined by per-document maximum per record. Missing or unknown extraction
identities on randstrobemer indexes require rebuilding. Bank payloads and the
abandoned-key membership semantics are unchanged. Version 7 changes the lossless
recovery-file encoding, not bank payloads or feature extraction.
Version 8 adds `segment_length`: one fixed power-of-two segment size for every
filter/subindex, default 64 cells, independent of filter capacity. Array sizes
remain powers of two and must fit at least `hashes` complete segments. Queries
use this stored size. Versions 1 through 7 retain their historical size-dependent
position mapping; changing a version number alone does not convert an index.

Each `bank-<exponent>.zorbank.zst` expands to magic `ZORPOW3\0`, followed by
four little-endian `u64` values: slots, bits per row, group count and data bytes.
Rows have `ceil(bits_per_row / 8)` bytes. Each row begins with all groups'
1-, 8- or 16-bit fingerprints, followed immediately by their packed origin bits. Singleton
groups have no origin bits. Group offsets in the manifest describe the original
complete-payload prefix; the origin-plane offset for group `i` is
`group_count * fingerprint_bits + group.offset - i * fingerprint_bits`.

Construction mode does not change the bank encoding or query layout. The default
two-pass builder records `construction_cell_storage` as
`two-pass-memory-packed-no-scratch` (verified cells retained in RAM bit-packed
at the group's payload width until bank compression). `--low-disk=false` records
`two-pass-temporary-disk-packed-tiled-packing` (the same packed cells staged in
files); `--disk-backed` records `temporary-disk-u64-tiled-packing`. The two-pass
publication directory contains only final index files. First-pass projected
capacities may differ from exact final capacities; this does not change the
planned dataset membership partition.

The optional manifest metadata `min_subindex_documents` records the build's
minimum original-document count per bank (new default 128, missing legacy
field defaults to 1). `promoted_groups` records filters rebuilt at a larger
power-of-two capacity to meet that constraint (missing legacy field defaults
to 0). Original documents are counted through `group.members`, not by the
number of merged groups. With fewer total documents than the requested minimum,
all groups occupy one bank. The planner selects a minimum-padded-bit-cost
contiguous partition of the exact slot-size buckets; only promoted filters
require additional input reads and reconstruction. This changes neither the
bank encoding nor document partitions and creates no temporary key/cell files.
Fingerprint width and abandoned-key union semantics stay fixed, but changing
capacity may change which absent keys are false positives. Setting the minimum
to 1 preserves the old layout. Queries fuse row XOR and fingerprint comparison
in registers where supported; this is an execution change, not a format change.

`abandoned.zormap.zst` expands to magic `ZORABN2\0`, followed by records:
the delta from the previous 32-bit hash (initially zero), the origin count,
then the ascending origin IDs delta-coded from zero within each record.
Integers use unsigned base-128 variable-length encoding (at most five bytes
per u32). No hash or origin bits are discarded. The reader also accepts the
legacy `ZORABN1\0` stream: little-endian u32 hash, count and origin IDs.
This is one global multimap across all banks;
32-bit hash collisions union origin lists. It is omitted with `--drop-abandoned`.

The current reader also supports manifest versions 1 through 7 and their legacy
`ZORPOW1\0` interleaved banks. Version-1 per-group abandoned-key records are
converted to the global 32-bit-hash multimap on load. Older readers cannot read
new version-8 directories. Merged directories currently cannot be appended or
repacked; rebuilding requires the complete input list.
