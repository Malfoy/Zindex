# Installation and troubleshooting

## Toolchain and native dependencies

Use Rust through rustup so `rust-toolchain.toml` selects the tested compiler,
currently Rust 1.97.0. The package's declared language minimum is 1.89; CI uses
the pinned version. A native-target LLVM code-generation failure was observed
with Rust 1.99.1 in the Bzip2 dependency on an AVX-512-capable runner. Reproduce
builds with the pinned compiler before diagnosing an application failure.

On Debian/Ubuntu, install `build-essential pkg-config`. Python 3 and Bash are
needed for the example scripts. Run all documented commands from the repository
root unless stated otherwise.

## Copying binaries between machines

`.cargo/config.toml` enables `target-cpu=native`. A binary compiled on a newer CPU
can fail with `Illegal instruction` on another CPU. Rebuild on the destination,
or build the portable x86-64 baseline:

```bash
RUSTFLAGS='-C target-cpu=x86-64' cargo build --release --locked
```

Portable fallback warnings are expected in that configuration. Native and
portable x86-64 are covered by CI; other platforms are not currently covered.

## Input paths, compression and outputs

FOF-relative paths use the command's working directory, not the FOF's directory.
Use `name<TAB>path` with an actual tab when names or paths contain spaces.
The indexers select compressed streaming input using its suffix; use `.gz`,
`.gzip`, `.zst`, `.zstd`, `.xz`, `.bz2`, `.bz` or `.bzip2` as appropriate.
The standalone minimizer utility detects compression by magic bytes.

Build outputs and demonstration destinations must be new paths. Reusing an
existing path is rejected rather than overwriting a previous result.

## Memory and large collections

Merged queries load one subindex at a time by default; `--query-memory all`
loads the entire index. The default 8-GiB construction cap limits uncompressed
filter payload per subindex, not total process RAM. Query seeds, results,
recovery matches, metadata and worker buffers use additional memory.

Construction stages completed filters and abandoned keys on disk by default.
`--low-disk` retains them in RAM. Inspect `info` and skipped-document reports
before treating an index as representing the full requested collection.
See [memory controls](../README.md#memory-and-collection-completeness) for
admission limits and the sampled RSS guard. Oversized datasets fail by default.
Use `--skip-oversized-documents` only when an incomplete collection is acceptable.

## Report a failure

Capture the exact command, `--version` output, complete stderr, OS/CPU and
available RAM/scratch disk. Version output includes the Git revision and dirty
status; source archives report when Git metadata is unavailable. Keep failed
demo logs at the directory printed by the script. See
[contributing](../.github/CONTRIBUTING.md) for reporting instructions.
