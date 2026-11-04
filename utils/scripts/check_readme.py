#!/usr/bin/env python3
"""Execute README usage examples against the bundled data in a temporary directory."""
import argparse
import os
from pathlib import Path
import re
import subprocess
import tempfile


def main():
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bin", type=Path, default=root / "target/release/zorindex")
    args = parser.parse_args()
    binary = args.bin.resolve(strict=True)
    readme = (root / "README.md").read_text()
    blocks = re.findall(r"^```bash\n(.*?)^```", readme, re.M | re.S)
    # Installation is covered by the CI build and the example checker has its
    # own CI step. Check their syntax without cloning or rebuilding here.
    for block in blocks:
        subprocess.run(["bash", "-n"], input=block, text=True, check=True)
    usage = readme.split("## Build and query\n", 1)[1]
    examples = re.findall(r"^```bash\n(.*?)^```", usage, re.M | re.S)
    with tempfile.TemporaryDirectory(prefix="zindex-readme-") as temporary:
        work = Path(temporary)
        (work / "tests").symlink_to(root / "tests", target_is_directory=True)
        (work / "target/release").mkdir(parents=True)
        (work / "target/release/zorindex").symlink_to(binary)
        # Supply compatible inputs for the general four-shard example. The
        # two-input example creates first.zindex before these links are used.
        for i in range(1, 5):
            (work / f"shard{i}.zindex").symlink_to("first.zindex", target_is_directory=True)
        env = dict(os.environ)
        env["TMPDIR"] = str(work)
        for number, block in enumerate(examples, 1):
            result = subprocess.run(
                ["bash", "-euo", "pipefail", "-c", block], cwd=work,
                env=env, text=True, capture_output=True,
            )
            if result.returncode:
                raise RuntimeError(
                    f"README example {number} failed:\n{block}\n"
                    + result.stdout + result.stderr
                )
    print(f"PASS: all {len(examples)} README usage examples")


if __name__ == "__main__":
    main()
