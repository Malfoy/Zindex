#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$ROOT"
BIN=${ZINDEX_BIN_DIR:-"$ROOT/target/release"}
DATA="$ROOT/tests/data/bacteria"
if [[ ! -x "$BIN/zorindex" ]]; then
  echo "FAIL: run cargo build --release --locked first." >&2
  exit 1
fi
if [[ $# -gt 1 ]]; then
  echo "FAIL: usage: bash utils/scripts/demo.sh [new-output-directory]" >&2
  exit 1
fi
if [[ $# -eq 1 ]]; then
  if [[ -e "$1" ]]; then echo "FAIL: output already exists: $1" >&2; exit 1; fi
  mkdir -- "$1"
  OUT=$(cd -- "$1" && pwd)
else
  OUT=$(mktemp -d "${TMPDIR:-/tmp}/zindex-example.XXXXXX")
fi
echo "Example results: $OUT"
trap 'code=$?; if [[ $code -ne 0 ]]; then
  echo "FAIL: example checks did not all pass. Logs: $OUT" >&2
  for log in "$OUT"/*.stderr; do
    if [[ -f "$log" ]]; then tail -20 "$log" >&2; fi
  done
fi' EXIT

for mode in kmers minimizers; do
  args=(--k 31)
  if [[ "$mode" == minimizers ]]; then
    args=(--k 49 --index-mode minimizers --minimizer-size 31)
  fi
  "$BIN/zorindex" build --fof "$DATA/samples.fof" --output "$OUT/$mode.zindex" \
    --threads 4 "${args[@]}" > "$OUT/$mode-build.stdout" 2> "$OUT/$mode-build.stderr"
  "$BIN/zorindex" query --index "$OUT/$mode.zindex" --query "$DATA/queries.fa" \
    --records --min-ratio 0.8 --threads 4 \
    > "$OUT/$mode-hits.tsv" 2> "$OUT/$mode-query.stderr"
  python3 utils/scripts/prepare_bacteria_example.py \
    --check-results "$OUT/$mode-hits.tsv" --truth "$DATA/source-truth.tsv" \
    > "$OUT/$mode-check.stdout" 2> "$OUT/$mode-check.stderr"
done
echo "PASS: all example checks passed (k-mers and minimizers)"
