import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


class PredictorIntegrationTests(unittest.TestCase):
    def test_bundled_model_regression(self):
        root = Path(__file__).resolve().parents[1]
        source = root / "tests" / "data" / "smoke_profiles.csv"
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "predictions.tsv"
            metrics = Path(tmp) / "metrics"
            command = [
                sys.executable,
                "-m",
                "cdiff_clademl.cli.predict",
                str(source),
                "-o",
                str(output),
                "--metrics_prefix",
                str(metrics),
                # The audit environment has sklearn 1.9; official installs are
                # pinned to the model's exact sklearn 1.7.1.
                "--allow-version-mismatch",
            ]
            completed = subprocess.run(
                command,
                cwd=root,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)

            predictions = pd.read_csv(output, sep="\t")
            self.assertEqual(
                predictions["final_prediction"].tolist(),
                ["Clade1", "Clade2", "Clade3", "Clade4", "Clade5", "Cryptic"],
            )
            self.assertEqual(set(predictions["quality_flag"]), {"PASS"})
            self.assertTrue((predictions["n_unknown_alleles"] == 0).all())
            self.assertTrue(Path(f"{metrics}_by_st_summary.tsv").is_file())

            summary_path = output.with_name("predictions_summary.json")
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            self.assertEqual(summary["model_version"], "1.0")
            self.assertEqual(summary["n_samples"], 6)
            self.assertEqual(summary["output_paths"]["predictions_tsv"], str(output.resolve()))
            self.assertTrue(output.with_name("report.html").is_file())

    def test_bundled_model_automatic_outputs_include_html_report(self):
        root = Path(__file__).resolve().parents[1]
        source = root / "tests" / "data" / "smoke_profiles.csv"
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            command = [
                sys.executable,
                "-m",
                "cdiff_clademl.cli.predict",
                str(source),
                "--allow-version-mismatch",
            ]
            completed = subprocess.run(
                command,
                cwd=work,
                text=True,
                capture_output=True,
                check=False,
                env={**os.environ, "PYTHONPATH": str(root)},
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            directories = list(work.glob("smoke_profiles_cdiffclademl_*"))
            self.assertEqual(len(directories), 1)
            result = directories[0]
            self.assertTrue((result / "predictions.tsv").is_file())
            self.assertTrue((result / "predictions_summary.json").is_file())
            self.assertTrue((result / "report.html").is_file())

    def test_invalid_input_creates_no_automatic_result_directory(self):
        root = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            completed = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "cdiff_clademl.cli.predict",
                    "missing_profiles.tsv",
                    "--allow-version-mismatch",
                ],
                cwd=work,
                text=True,
                capture_output=True,
                check=False,
                env={**os.environ, "PYTHONPATH": str(root)},
            )

            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("Input table does not exist", completed.stderr)
            self.assertEqual(list(work.glob("*_cdiffclademl_*")), [])


if __name__ == "__main__":
    unittest.main()
