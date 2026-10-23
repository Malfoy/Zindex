Vendored from the crates.io `helicase` 0.1.1 source distribution.
Upstream: https://github.com/imartayan/helicase
Authors: Igor Martayan and Charles Paperman. License: MIT; upstream LICENSE
is included verbatim (copyright 2025 Igor Martayan).

Local change: `src/lib.rs` exposes the existing `lexer` module publicly so
minimizer-fasta can use SIMD FASTA tokenization with bounded streaming memory.
No lexer, parser, compression, or SIMD algorithm has been changed.
