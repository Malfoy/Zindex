#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$ROOT"
BIN=${ZINDEX_BIN_DIR:-"$ROOT/target/release"}
for tool in zorindex bloomindex minimizer-fasta; do
  if [[ ! -x "$BIN/$tool" ]]; then
    echo "Missing $BIN/$tool; run cargo build --release --locked first." >&2
    exit 1
  fi
done
command -v python3 >/dev/null || { echo "Python 3 is required for validation." >&2; exit 1; }
if [[ $# -gt 1 ]]; then
  echo "Usage: bash utils/scripts/demo.sh [new-output-directory]" >&2
  exit 1
fi
if [[ $# -eq 1 ]]; then
  if [[ -e "$1" ]]; then
    echo "Refusing to reuse existing path: $1" >&2
    exit 1
  fi
  mkdir -- "$1"
  OUT=$(cd -- "$1" && pwd)
else
  OUT=$(mktemp -d "${TMPDIR:-/tmp}/zindex-demo.XXXXXX")
fi
echo "Demonstration results: $OUT"
trap 'status=$?; if [[ $status -ne 0 ]]; then
  echo "Demonstration failed; logs retained in $OUT" >&2
  for log in "$OUT"/*.stderr; do
    if [[ -f "$log" ]]; then echo "--- $log" >&2; tail -30 "$log" >&2; fi
  done
fi' EXIT
"$BIN/zorindex" --version > "$OUT/version.txt"
"$BIN/zorindex" build --fof tests/data/tiny/samples.fof --output "$OUT/zindex" \
  --k 49 --index-mode minimizers --minimizer-size 31 --fingerprint-bits 16 \
  --threads 4 > "$OUT/zindex-build.stdout" 2> "$OUT/zindex-build.stderr"
"$BIN/zorindex" query --index "$OUT/zindex" --query tests/data/tiny/queries.fa \
  --records --min-ratio 0.8 --threads 4 > "$OUT/zindex-records.tsv" 2> "$OUT/zindex-query.stderr"
"$BIN/zorindex" info --index "$OUT/zindex" > "$OUT/zindex-info.txt"
"$BIN/bloomindex" build --fof tests/data/tiny/samples.fof --output "$OUT/bindex.blmidx" \
  --k 49 --index-mode minimizers --minimizer-size 31 --hashes 4 \
  --false-positive-rate 0.0000152587890625 --threads 4 \
  > "$OUT/bindex-build.stdout" 2> "$OUT/bindex-build.stderr"
"$BIN/bloomindex" query --index "$OUT/bindex.blmidx" --query tests/data/tiny/a-only.fa \
  --min-ratio 0.8 --threads 4 > "$OUT/bindex-query.tsv"
"$BIN/minimizer-fasta" --minimizer-scheme open-closed-mod --threads 4 \
  tests/data/tiny/a.fa "$OUT/minimizers.fa.gz" \
  > "$OUT/minimizer.stdout" 2> "$OUT/minimizer.stderr"
python3 - "$OUT" <<'PY'
import csv
import gzip
import json
import sys
from pathlib import Path

out = Path(sys.argv[1])
manifest = json.loads((out / "zindex/manifest.json").read_text())
names = [color["name"] for color in manifest["colors"]]
expected = {}
with open("tests/data/tiny/expected-hits.tsv", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        expected[row["query_id"]] = row["datasets"].split(",") if row["datasets"] != "-" else []
with (out / "zindex-records.tsv").open(newline="") as handle:
    records = list(csv.DictReader(handle, delimiter="\t"))
assert [row["query_id"] for row in records] == list(expected), "Missing/reordered queries"
for row in records:
    count = int(row["valid_kmers"])
    assert count > 0, "Fixture must contain selected seeds"
    hits = [pair.split(":") for pair in row["hits"].split(",") if pair]
    assert [names[int(i)] for i, _ in hits] == expected[row["query_id"]], row
    assert all(int(n) == count for _, n in hits), row
with (out / "bindex-query.tsv").open() as handle:
    rows = [line.rstrip("\n").split("\t") for line in handle if not line.startswith("#")]
assert len(rows) == 2 and rows[1][0] == "a", rows
assert float(rows[1][2]) == 1.0, rows
with gzip.open(out / "minimizers.fa.gz", "rt") as handle:
    assert handle.readline().startswith(">"), "Empty minimizer output"
print("PASS: example queries and all three tools")
PY
