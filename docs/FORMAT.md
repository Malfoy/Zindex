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

Current format version: `7` (readers also accept legacy versions).

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
`manifest.json` contains version `3`, feature/hash settings, dataset identities,
group members and widths, bank sizes and abandoned-key map metadata.

Each `bank-<exponent>.zorbank.zst` expands to magic `ZORPOW3\0`, followed by
four little-endian `u64` values: slots, bits per row, group count and data bytes.
Rows have `ceil(bits_per_row / 8)` bytes. Each row begins with all groups'
1-, 8- or 16-bit fingerprints, followed immediately by their packed origin bits. Singleton
groups have no origin bits. Group offsets in the manifest describe the original
complete-payload prefix; the origin-plane offset for group `i` is
`group_count * fingerprint_bits + group.offset - i * fingerprint_bits`.

`abandoned.zormap.zst` expands to magic `ZORABN1\0`, followed by records:
32-bit abandoned-key hash, 32-bit origin count, then that many 32-bit origin IDs.
All integers are little-endian. This is one global multimap across all banks;
32-bit hash collisions union origin lists. It is omitted with `--drop-abandoned`.

The current reader also supports manifest versions 1 and 2 and their legacy
`ZORPOW1\0` interleaved banks. Version-1 per-group abandoned-key records are
converted to the global 32-bit-hash multimap on load. Older readers cannot read
new version-3 directories. Merged directories currently cannot be appended or
repacked; rebuilding requires the complete input list.
