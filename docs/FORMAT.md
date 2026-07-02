# Index Format Notes

This document summarizes the binary formats used by the two command-line
indexes. The exact reader and writer implementations are the source of truth.

## Zindex (`.zoridx`)

Magic:

```text
ZORIDX1\0
```

Current format version: `6`.

Header fields:

- version
- `k`
- ZOR hash count
- fingerprint width (`8` or `16`)
- flags
- seed
- feature mode
- minimizer size
- modimizer sampling
- ZOR segmented layout
- color metadata
- slot-major stack
- optional union graph

The stack can be stored plain or compressed. Plain stacks store all rows as
contiguous little-endian fingerprint bytes. Compressed stacks store row offsets
and codec payload bytes. Each compressed row expands to:

```text
color_count * fingerprint_bytes
```

The query path computes row positions once per feature. In plain mode it XORs
directly from the interleaved rows. In compressed mode it decodes the selected
rows into scratch buffers and XORs those rows.

## Zindex Stack Codecs

The stored codec byte maps to:

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

Current format version: `1`.

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
- bitstack words

Bindex stores `bit_count * color_words` little-endian `u64` words, where:

```text
color_words = ceil(color_count / 64)
```

For each Bloom bit position, all colors are packed into consecutive `u64`
words. Querying ANDs the packed words from all Bloom positions and increments
scores for the remaining set bits.
