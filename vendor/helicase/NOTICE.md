Vendored from the crates.io `helicase` 1.2.0 source distribution.
Upstream: https://github.com/imartayan/helicase
Authors: Igor Martayan and Charles Paperman. License: MIT; upstream LICENSE
is included verbatim. Upstream source revision:
6f06193dab6204a2951b2cfb61569828406fc2c2 (see .cargo_vcs_info.json).

Local API exposure only: `src/lib.rs` makes the existing `lexer` module public;
`Lexer::next_block` is public and a `current_block()` accessor exposes block bytes.
This lets minimizer-fasta consume the upstream SIMD classifier with bounded
streaming memory. No upstream parser, compression, or SIMD algorithm is changed.
