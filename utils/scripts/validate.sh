#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$ROOT"
CPU=${1:-all}
case "$CPU" in
  all) CPUS=(native x86-64) ;;
  native|x86-64) CPUS=("$CPU") ;;
  *) echo "Usage: bash utils/scripts/validate.sh [all|native|x86-64]" >&2; exit 2 ;;
esac
export CARGO_INCREMENTAL=0
export CARGO_BUILD_JOBS=${CARGO_BUILD_JOBS:-2}
export RUST_TEST_THREADS=${RUST_TEST_THREADS:-4}

cargo fmt --all --check
python3 -m unittest discover -s tests -p 'test_bacteria_example.py'
for CPU in "${CPUS[@]}"; do
  echo "Validating $CPU"
  export RUSTFLAGS="-C target-cpu=$CPU"
  export CARGO_TARGET_DIR="$ROOT/target/validation/$CPU"
  cargo test --workspace --release --locked
  cargo build --workspace --release --locked --bins
  RESULTS=$(mktemp -d "${RUNNER_TEMP:-${TMPDIR:-/tmp}}/zindex-validation-$CPU.XXXXXX")
  ZINDEX_BIN_DIR="$CARGO_TARGET_DIR/release" bash utils/scripts/demo.sh "$RESULTS/demo"
  python3 utils/scripts/check_readme.py --bin "$CARGO_TARGET_DIR/release/zorindex"
done
echo "PASS: validation completed"
