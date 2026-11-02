"""Offline checks for the optional genome downloader and query truth files."""
import csv
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("bacteria_example", ROOT / "utils/scripts/prepare_bacteria_example.py")
example = importlib.util.module_from_spec(spec)
spec.loader.exec_module(example)


class BacteriaExampleTests(unittest.TestCase):
    def fixture(self, root, sequence=None, bad_checksum=False):
        sequence = sequence or "NNN" + "ACGTTGCA" * 200
        payload = gzip.compress((">contig description\n" + sequence + "\n").encode(), mtime=0)
        digest = hashlib.md5(payload, usedforsecurity=False).hexdigest()
        entry = {"accession": "GCF_example.1", "organism": "test", "url": "https://example.invalid/genome.gz",
                 "ncbi_md5": "0" * 32 if bad_checksum else digest}
        manifest = root / "assemblies.json"
        manifest.write_text(json.dumps([entry]))
        return manifest, payload

    def test_prepares_verified_files_source_coordinates_and_checked_queries(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest, payload = self.fixture(root)
            output = root / "prepared"
            with mock.patch.object(example.request, "urlopen", return_value=io.BytesIO(payload)) as fetch:
                example.prepare(output, manifest)
            self.assertEqual(fetch.call_count, 1)
            metadata = json.loads((output / "sources.json").read_text())
            self.assertEqual(metadata["status"], "complete")
            self.assertEqual(metadata["assemblies"][0]["sha256"], hashlib.sha256(payload).hexdigest())
            self.assertEqual(metadata["assemblies"][0]["bases"], 1603)
            self.assertEqual(metadata["assemblies"][0]["contigs"], 1)
            self.assertIn(str(output / "GCF_example.1.fna.gz"), (output / "samples.fof").read_text())
            with (output / "source-truth.tsv").open() as source:
                truth = list(csv.DictReader(source, delimiter="\t"))
            self.assertEqual([r["start_0based"] for r in truth[:2]], ["3", "3"])
            sequences = (output / "queries.fa").read_text().splitlines()[1::2]
            self.assertEqual(len(sequences), 3)
            self.assertEqual(sequences[1], example.reverse_complement(sequences[0]))
            hits = output / "hits.tsv"
            hits.write_text("query_id\tvalid_kmers\thits\n" + "".join(
                r["query_id"] + "\t100\t" + ("0:100" if r["document_id"] != "-" else "") + "\n"
                for r in truth))
            example.check_results(hits)
            hits.write_text(hits.read_text().replace("0:100", "0:99", 1))
            with self.assertRaisesRegex(ValueError, "Source minimizers missing"):
                example.check_results(hits)

    def test_checksum_failure_removes_partial_and_does_not_publish_fof(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest, payload = self.fixture(root, bad_checksum=True)
            output = root / "prepared"
            with mock.patch.object(example.request, "urlopen", return_value=io.BytesIO(payload)):
                with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                    example.prepare(output, manifest)
            self.assertEqual(json.loads((output / "sources.json").read_text())["status"], "failed")
            self.assertFalse((output / "samples.fof").exists())
            self.assertFalse(list(output.glob("*.partial")))
            self.assertFalse(list(output.glob("*.fna.gz")))

    def test_short_or_ambiguous_input_does_not_publish_fof(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest, payload = self.fixture(root, sequence="AN" * 1200)
            output = root / "prepared"
            with mock.patch.object(example.request, "urlopen", return_value=io.BytesIO(payload)):
                with self.assertRaisesRegex(ValueError, "No 1000-base ACGT run"):
                    example.prepare(output, manifest)
            self.assertFalse((output / "samples.fof").exists())
            self.assertEqual(json.loads((output / "sources.json").read_text())["status"], "failed")

    def test_existing_destination_is_untouched_without_network_access(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest, _ = self.fixture(root)
            output = root / "prepared"
            output.mkdir()
            sentinel = output / "keep.txt"
            sentinel.write_text("keep")
            with mock.patch.object(example.request, "urlopen") as fetch:
                with self.assertRaises(FileExistsError):
                    example.prepare(output, manifest)
            fetch.assert_not_called()
            self.assertEqual(sentinel.read_text(), "keep")
            self.assertEqual(list(output.iterdir()), [sentinel])


if __name__ == "__main__":
    unittest.main()
