#!/usr/bin/env python3
"""Download three pinned bacterial assemblies, or check their example query results."""
import argparse
import csv
import datetime
import gzip
import hashlib
import json
from pathlib import Path
import sys
from urllib import request

ROOT = Path(__file__).resolve().parents[2]
ASSEMBLIES = ROOT / "tests/data/bacteria/assemblies.json"
QUERY_LENGTH = 1000


def save(path, value):
    temporary = path.with_name(path.name + ".partial")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def download(entry, output):
    path = output / (entry["accession"] + ".fna.gz")
    partial = path.with_name(path.name + ".partial")
    md5 = hashlib.md5(usedforsecurity=False)
    sha256 = hashlib.sha256()
    try:
        req = request.Request(entry["url"], headers={"User-Agent": "Zindex-example/0.1"})
        with request.urlopen(req, timeout=60) as source, partial.open("xb") as dest:
            for block in iter(lambda: source.read(256 * 1024), b""):
                md5.update(block)
                sha256.update(block)
                dest.write(block)
        if md5.hexdigest() != entry["ncbi_md5"]:
            raise ValueError("NCBI checksum mismatch for " + entry["accession"])
        partial.replace(path)
    finally:
        if partial.exists():
            partial.unlink()
    return path, sha256.hexdigest()


def source_query(path):
    """Select the first 1000-base ACGT run; retain contig-relative coordinates."""
    selected = None
    run = []
    run_start = 0
    contig = None
    offset = 0
    bases = 0
    contigs = 0
    with gzip.open(path, "rt", encoding="ascii") as source:
        for line in source:
            if line.startswith(">"):
                contig = line[1:].split()[0]
                contigs += 1
                offset = 0
                run.clear()
                continue
            seq = "".join(line.split()).upper()
            if seq and contig is None:
                raise ValueError("Sequence before FASTA header in " + str(path))
            bases += len(seq)
            if selected is None:
                for i, base in enumerate(seq):
                    if base not in "ACGT":
                        run.clear()
                        continue
                    if not run:
                        run_start = offset + i
                    run.append(base)
                    if len(run) == QUERY_LENGTH:
                        selected = ("".join(run), contig, run_start)
                        break
            offset += len(seq)
    if selected is None:
        raise ValueError("No 1000-base ACGT run in " + str(path))
    return (*selected, bases, contigs)


def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def negative_query():
    # Explicit generator avoids dependence on Python's random selection methods.
    state = 20261102
    result = []
    for _ in range(QUERY_LENGTH):
        state = (1664525 * state + 1013904223) & 0xFFFFFFFF
        result.append("ACGT"[state >> 30])
    return "".join(result)


def prepare(output, manifest=ASSEMBLIES):
    entries = json.loads(manifest.read_text())
    if not entries or len({e["accession"] for e in entries}) != len(entries):
        raise ValueError("Assembly manifest must contain distinct accessions")
    output = output.resolve()
    output.mkdir(parents=True)  # Existing paths are never reused.
    metadata = {"status": "downloading", "retrieved_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                "assemblies": [], "query_length": QUERY_LENGTH, "negative_lcg_seed": 20261102}
    save(output / "sources.json", metadata)
    try:
        fof, queries, truth = [], [], []
        for document_id, entry in enumerate(entries):
            path, sha256 = download(entry, output)
            seq, contig, start, bases, contigs = source_query(path)
            metadata["assemblies"].append(dict(entry, document_id=document_id,
                path=str(path), sha256=sha256, compressed_bytes=path.stat().st_size,
                bases=bases, contigs=contigs))
            save(output / "sources.json", metadata)
            fof.append(entry["accession"] + "\t" + str(path) + "\n")
            for strand, dna in [("forward", seq), ("reverse", reverse_complement(seq))]:
                name = entry["accession"] + "_" + strand
                queries.append(">" + name + "\n" + dna + "\n")
                truth.append([name, str(document_id), entry["accession"], contig, str(start), strand])
        queries.append(">random_negative\n" + negative_query() + "\n")
        truth.append(["random_negative", "-", "-", "-", "-", "-"])
        (output / "queries.fa").write_text("".join(queries))
        with (output / "source-truth.tsv").open("w", newline="") as dest:
            writer = csv.writer(dest, delimiter="\t", lineterminator="\n")
            writer.writerow(["query_id", "document_id", "accession", "contig", "start_0based", "strand"])
            writer.writerows(truth)
        metadata["queries_sha256"] = hashlib.sha256((output / "queries.fa").read_bytes()).hexdigest()
        # Publish the FOF only after all downloads, gzip trailers and queries validate.
        (output / "samples.fof").write_text("".join(fof))
        metadata["status"] = "complete"
    except Exception as error:
        metadata["status"] = "failed"
        metadata["error"] = str(error)
        raise
    finally:
        save(output / "sources.json", metadata)
    print("Prepared " + str(output / "samples.fof"))


def check_results(path):
    with path.with_name("source-truth.tsv").open(newline="") as source:
        expected = list(csv.DictReader(source, delimiter="\t"))
    with path.open(newline="") as source:
        rows = list(csv.DictReader(source, delimiter="\t"))
    if [r["query_id"] for r in rows] != [r["query_id"] for r in expected]:
        raise ValueError("Missing, additional or reordered query records")
    positives = 0
    for row, truth in zip(rows, expected):
        if truth["document_id"] == "-":
            if row["hits"]:
                print("Random negative control has reported matches: " + row["hits"])
            continue
        hits = dict(pair.split(":") for pair in row["hits"].split(",") if pair)
        total = int(row["valid_kmers"])
        if total <= 0 or int(hits.get(truth["document_id"], "0")) != total:
            raise ValueError("Source minimizers missing from " + row["query_id"])
        positives += 1
    print("PASS: all " + str(positives) + " bacterial source queries match, including reverse complements")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--output", type=Path, help="Download into a new directory (never reused)")
    modes.add_argument("--check-results", type=Path, help="Check hits.tsv next to the generated source-truth.tsv")
    args = parser.parse_args()
    if args.output is not None:
        prepare(args.output)
    else:
        check_results(args.check_results)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as error:
        print("Error:", error, file=sys.stderr)
        sys.exit(1)
