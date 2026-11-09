#!/usr/bin/env python3
"""One build and one query per workload, with frozen settings and GNU time logs."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def identity(path, hash_bytes=False):
    stat = path.stat()
    result = {"path": str(path), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    if hash_bytes:
        result["sha256"] = digest(path)
    return result


def color_name(path):
    name = path.name
    for suffix in [".zstd", ".zst", ".gzip", ".gz", ".xz", ".bz2", ".bz", ".bzip2",
                   ".fasta", ".fastq", ".fna", ".fa", ".fq", ".fas"]:
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    return name or "color"


def documents(fof):
    result = []
    for line in fof.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t") if "\t" in line else line.split()
        if len(fields) >= 2:
            name, file = fields[0].strip(), fields[1].strip()
        else:
            file = line
            name = color_name(Path(file))
        path = Path(file).resolve(strict=True)  # FOF paths use the caller's cwd.
        if not name or not path.is_file():
            raise ValueError("Invalid FOF document: " + line)
        result.append((name, path))
    if not result:
        raise ValueError("FOF contains no documents")
    return result


def save(path, value):
    temporary = path.with_name(path.name + ".partial")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def optional_command(args):
    try:
        p = subprocess.run(args, capture_output=True, text=True, check=False)
        return p.stdout.strip() if p.returncode == 0 else None
    except OSError:
        return None


def measured(root, label, command):
    save(root / (label + ".command.json"), command)
    env = dict(os.environ, LC_ALL="C")
    start = time.monotonic()
    with (root / (label + ".stdout")).open("wb") as out, (root / (label + ".stderr")).open("wb") as err:
        p = subprocess.run(["/usr/bin/time", "-v", "-o", str(root / (label + ".time.txt")),
                            *command], stdout=out, stderr=err, env=env, check=False)
    metrics = {"exit_status": p.returncode, "wall_seconds": time.monotonic() - start}
    # GNU time reports KiB on Linux. Keep its complete raw report as evidence.
    for line in (root / (label + ".time.txt")).read_text().splitlines():
        if "Maximum resident set size (kbytes):" in line:
            metrics["peak_rss_bytes"] = int(line.rsplit(":", 1)[1]) * 1024
        for text, field in [("User time (seconds):", "user_seconds"),
                            ("System time (seconds):", "system_seconds")]:
            if text in line:
                metrics[field] = float(line.rsplit(":", 1)[1])
    save(root / (label + ".metrics.json"), metrics)
    if p.returncode:
        raise RuntimeError(label + " failed; see " + str(root / (label + ".stderr")))
    return metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fof", type=Path, required=True)
    parser.add_argument("--positive", type=Path, required=True, help="Positive query FASTA")
    parser.add_argument("--negative", type=Path, required=True, help="Negative query FASTA")
    parser.add_argument("--output", type=Path, required=True, help="New result directory; never reused")
    parser.add_argument("--binary", type=Path, default=Path("target/release/zorindex"))
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--k", type=int, default=49, help="K-mer length or minimizer window length in bases")
    parser.add_argument("--minimizer-size", type=int, default=31, help="Selected minimizer length in bases")
    parser.add_argument("--minimizer-scheme", choices=["open-closed-mod", "simd"], default="open-closed-mod")
    parser.add_argument("--index-mode", choices=["kmers", "minimizers", "modimizers", "findere"], default="minimizers")
    parser.add_argument("--hashes", type=int, default=4)
    parser.add_argument("--fingerprint-bits", type=int, choices=[1, 2, 4, 8, 16], default=8)
    parser.add_argument("--merge-planner", choices=["streaming", "full", "naive"], default="streaming")
    parser.add_argument("--query-memory", choices=["subindex", "all"], default="subindex")
    parser.add_argument("--max-subindex-size", default="8GiB")
    parser.add_argument("--max-subindexes", type=int)
    parser.add_argument("--min-subindex-documents", type=int, default=128)
    parser.add_argument("--merge-plan-margin", type=float, default=0.03)
    parser.add_argument("--min-ratio", type=float, default=0.8)
    parser.add_argument("--hash-inputs", action="store_true", help="Also SHA-256 every source sequence file (extra full input read)")
    parser.add_argument("--hash-index", action="store_true", help="SHA-256 index files after queries (extra full index read)")
    args = parser.parse_args()
    if args.threads < 1 or not 0 <= args.min_ratio <= 1:
        parser.error("threads must be positive; min-ratio must be in [0,1]")
    if platform.system() != "Linux" or not Path("/usr/bin/time").is_file():
        parser.error("This measurement runner requires Linux and GNU /usr/bin/time")
    if "gnu" not in (optional_command(["/usr/bin/time", "--version"]) or "").lower():
        parser.error("/usr/bin/time is not GNU time")
    binary = args.binary.resolve(strict=True)
    if not os.access(binary, os.X_OK):
        parser.error("Binary is not executable")
    args.fof = args.fof.resolve(strict=True)
    args.positive = args.positive.resolve(strict=True)
    args.negative = args.negative.resolve(strict=True)
    inputs = documents(args.fof)
    root = args.output.resolve()
    root.mkdir(parents=True)  # An existing output is an error; never retry or overwrite.
    result = {"status": "preparing", "runs_per_stage": 1,
              "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "cwd": str(Path.cwd()), "host": platform.uname()._asdict(),
              "cpu_affinity": sorted(os.sched_getaffinity(0)),
              "options": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
              "stages": {}}
    save(root / "results.json", result)
    try:
        frozen = root / "zorindex.binary"
        shutil.copy2(binary, frozen)
        result["binary"] = identity(frozen, True)
        result["version"] = subprocess.check_output([str(frozen), "--version"], text=True).strip()
        result["rustc"] = optional_command(["rustc", "--version"])
        result["fof"] = identity(args.fof, True)
        result["queries"] = {name: identity(path, True) for name, path in
                             [("positive", args.positive), ("negative", args.negative)]}
        # Preserve FOF order, names and resolved paths for rerunning elsewhere.
        with (root / "documents.jsonl").open("w") as out:
            for i, (name, path) in enumerate(inputs):
                out.write(json.dumps(dict(identity(path, args.hash_inputs), document_id=i, name=name)) + "\n")
        fof = root / "inputs.fof"
        fof.write_text("".join(name + "\t" + str(path) + "\n" for name, path in inputs))
        result["status"] = "running"
        save(root / "results.json", result)
        build = [str(frozen), "build", "--fof", str(fof), "--output", str(root / "index"),
                 "--threads", str(args.threads), "--k", str(args.k), "--index-mode", args.index_mode,
                 "--minimizer-scheme", args.minimizer_scheme,
                 "--hashes", str(args.hashes), "--fingerprint-bits", str(args.fingerprint_bits),
                 "--merge-planner", args.merge_planner,
                 "--merge-sketch-size", "128", "--seed", "69",
                 "--slot-scale", "1.0", "--merge-plan-margin", str(args.merge_plan_margin),
                 "--cycle-break", "no-heuristic", "--tie-scan", "1", "--segment-length", "64",
                 "--min-subindex-documents", str(args.min_subindex_documents),
                 "--max-subindex-size", args.max_subindex_size, "--skip-oversized-documents=false"]
        if args.index_mode == "minimizers":
            build += ["--minimizer-size", str(args.minimizer_size)]
        if args.max_subindexes is not None:
            build += ["--max-subindexes", str(args.max_subindexes)]
        result["stages"]["build"] = measured(root, "build", build)
        manifest = json.loads((root / "index/manifest.json").read_text())
        if len(manifest["colors"]) != len(inputs):
            raise RuntimeError("Constructed document count differs from the requested collection")
        shutil.copy2(root / "index/manifest.json", root / "manifest.json")
        result["index_bytes"] = sum(p.stat().st_size for p in (root / "index").rglob("*") if p.is_file())
        result["index_files"] = [identity(p) for p in sorted((root / "index").rglob("*")) if p.is_file()]
        save(root / "results.json", result)
        for label, query in [("positive", args.positive), ("negative", args.negative)]:
            command = [str(frozen), "query", "--index", str(root / "index"), "--query", str(query),
                       "--threads", str(args.threads), "--query-memory", args.query_memory,
                       "--min-ratio", str(args.min_ratio)]
            if args.index_mode != "findere":
                command += ["--records", "--output-format", "compact"]
            result["stages"][label] = measured(root, label, command)
            save(root / "results.json", result)
        if args.hash_index:
            result["index_files"] = [identity(p, True) for p in sorted((root / "index").rglob("*")) if p.is_file()]
        result["manifest_sha256"] = digest(root / "manifest.json")
        result["status"] = "complete"
    except KeyboardInterrupt:
        result["status"] = "interrupted"
        raise
    except (OSError, RuntimeError, ValueError, subprocess.SubprocessError) as error:
        result["status"] = "failed"
        result["error"] = str(error)
        raise
    finally:
        save(root / "results.json", result)
    print(root / "results.json")


if __name__ == "__main__":
    try:
        main()
    except (OSError, RuntimeError, ValueError, subprocess.SubprocessError) as error:
        print("Error:", error, file=sys.stderr)
        sys.exit(1)
