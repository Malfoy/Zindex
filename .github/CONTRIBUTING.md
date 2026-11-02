# Contributing and reporting problems

Before reporting a failure, include the tool's `--version` output, exact command,
complete error message, OS/CPU, Rust version and available RAM/scratch disk.
For an existing index, include `info` output and the manifest when it contains
no confidential dataset identities. Do not publish private sequences or paths.

Build and run the small example with the commands in [README.md](../README.md).
For a development check:

```bash
cargo fmt --all -- --check
cargo test --workspace --locked
```

CI runs the release suite with native and portable x86-64 code generation.
Use the pinned Rust toolchain from `rust-toolchain.toml` when reproducing CI;
the workflow installs that same compiler. Compiler upgrades should be validated
on both targets before updating the pin.
Regression tests should exercise observable behavior, including stored settings
and loading when they are affected. Performance experiments are separate:
use the [benchmarking guide](../docs/BENCHMARKING.md), record one run per
configuration/workload, and preserve raw output rather than overwriting it.

Keep generated indexes, genomes and private experiment directories out of Git.
The tiny example FASTAs are intentional exceptions to the data ignore rules.
Do not change historical benchmark numbers when updating an implementation.

Vendored code changes must retain upstream licenses and update the corresponding
NOTICE.md with the exact origin and local modifications.

## Preparing a release

Run the correctness suite and example, review the unreleased changelog and format
compatibility, and choose an immutable commit/tag for the submission. Record the
toolchain and binary hashes. Publish release notes stating the defaults and
compatibility changes; do not move an existing published tag. Before describing
manuscript experiments as reproducible, supply the artifacts listed in
[REPRODUCIBILITY.md](../docs/REPRODUCIBILITY.md).
