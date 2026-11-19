import ast
import csv
import gzip
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class RepositoryValidationTests(unittest.TestCase):
    def test_reference_data_are_readable(self):
        known = (ROOT / "reference_data/known_AML_ref_fusion_list.txt").read_text(
            encoding="utf-8"
        ).splitlines()
        self.assertIn("BCR--ABL1", known)

        pfam = ROOT / "reference_data/pfam_domains.bed.gz"
        with gzip.open(pfam, "rt", encoding="utf-8") as handle:
            columns = handle.readline().rstrip("\n").split("\t")
        self.assertEqual(
            columns,
            [
                "Transcript_id",
                "Pfam_id",
                "Domain_name_abbreviation",
                "Domain_name_full",
                "Start",
                "End",
            ],
        )

    def test_ingest_source_parses(self):
        source = ROOT / "llm_mining/ingest.py"
        ast.parse(source.read_text(encoding="utf-8"), filename=str(source))

    def test_ingest_runner_requires_api_key(self):
        result = subprocess.run(
            ["bash", str(ROOT / "llm_mining/run_ingest.sh")],
            env={"PATH": "/usr/bin:/bin"},
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("NVIDIA_API_KEY", result.stderr)

    def test_combine_fusions_writes_header_and_sample_names(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            study = root / "study"
            for sample, fusion in [("S1", "BCR--ABL1"), ("S2", "ETV6--RUNX1")]:
                sample_dir = study / sample
                sample_dir.mkdir(parents=True)
                input_file = sample_dir / f"{sample}.FusionInspector.fusions.abridged.tsv"
                input_file.write_text(
                    f"X.FusionName\tvalue\n{fusion}\t1\n", encoding="utf-8"
                )

            output = root / "combined.tsv"
            result = subprocess.run(
                [
                    "bash",
                    str(ROOT / "workflow/scripts/combine_fusions.sh"),
                    "-i",
                    str(study),
                    "-o",
                    str(output),
                    "-p",
                    "*.FusionInspector.fusions.abridged.tsv",
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)

            with output.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["sample_name"] for row in rows], ["S1", "S2"])
            self.assertEqual(
                [row["X.FusionName"] for row in rows],
                ["BCR--ABL1", "ETV6--RUNX1"],
            )

    def test_combine_does_not_replace_output_when_no_files_match(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            study = root / "study"
            study.mkdir()
            output = root / "combined.tsv"
            output.write_text("existing\n", encoding="utf-8")

            result = subprocess.run(
                [
                    "bash",
                    str(ROOT / "workflow/scripts/combine_fusions.sh"),
                    "-i",
                    str(study),
                    "-o",
                    str(output),
                    "-p",
                    "*.not-present",
                ],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(output.read_text(encoding="utf-8"), "existing\n")


if __name__ == "__main__":
    unittest.main()
