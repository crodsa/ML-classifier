import json
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
            self.assertEqual(summary["model_version"], "5.3.0")
            self.assertEqual(summary["n_samples"], 6)


if __name__ == "__main__":
    unittest.main()
