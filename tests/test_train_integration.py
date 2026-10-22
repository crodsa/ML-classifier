import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd
import sklearn

from cdiff_clademl.config import LOCI


class TrainerIntegrationTests(unittest.TestCase):
    def _write_training_table(self, path: Path) -> None:
        rows = []
        labels = ["Clade1", "Clade2", "Clade3", "Clade4", "Clade5", "Cryptic"]
        for allele, label in enumerate(labels, start=1):
            for replicate in range(1, 6):
                fields = [f"{label}_{replicate}", "cdifficile", f"{allele}{replicate}"]
                fields.extend(f"{locus}({allele})" for locus in LOCI)
                fields.append(f"mlst_clade({label})")
                rows.append(",".join(fields))
        path.write_text("\n".join(rows) + "\n", encoding="utf-8")

    def _write_distance_matrices(self, directory: Path) -> None:
        directory.mkdir()
        lines = ["allele_1\tallele_2\tdistance"]
        for first in range(1, 7):
            for second in range(1, 7):
                if first == second:
                    distance = 0.0
                elif 6 in {first, second}:
                    distance = 0.8
                else:
                    distance = 0.2
                lines.append(f"{first}\t{second}\t{distance}")
        text = "\n".join(lines) + "\n"
        for locus in LOCI:
            (directory / f"{locus}.tsv").write_text(text, encoding="utf-8")

    def test_train_and_reload_bundle(self):
        root = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            training = work / "train.csv"
            matrices = work / "distances"
            model = work / "new_model.joblib"
            metrics = work / "oof"
            self._write_training_table(training)
            self._write_distance_matrices(matrices)

            train_command = [
                sys.executable,
                "-m",
                "cdiff_clademl.cli.train",
                str(training),
                "--distance_dir",
                str(matrices),
                "--metrics_prefix",
                str(metrics),
                "--min_cryptic_recall",
                "1.0",
                "--threshold_grid_size",
                "5",
                "-o",
                str(model),
            ]
            trained = subprocess.run(
                train_command,
                cwd=root,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(trained.returncode, 0, trained.stderr)
            self.assertTrue(model.is_file())
            card = json.loads((work / "new_model_model_card.json").read_text())
            self.assertEqual(card["bundle_schema_version"], 1)
            self.assertEqual(card["software"]["scikit_learn"], sklearn.__version__)

            oof = pd.read_csv(f"{metrics}_oof_meta.tsv", sep="\t")
            self.assertEqual(len(oof), 30)
            self.assertEqual(int(oof["is_cryptic"].sum()), 5)
            self.assertEqual(oof["isolate"].nunique(), 30)
            self.assertTrue(Path(f"{metrics}_by_st_summary.tsv").is_file())

            predictions = work / "predictions.tsv"
            predicted = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "cdiff_clademl.cli.predict",
                    str(training),
                    "--model",
                    str(model),
                    "-o",
                    str(predictions),
                    "--report-html",
                    str(work / "report.html"),
                ],
                cwd=root,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(predicted.returncode, 0, predicted.stderr)
            output = pd.read_csv(predictions, sep="\t")
            self.assertEqual(output.final_prediction.tolist(), output.mlst_clade.tolist())
            report = (work / "report.html").read_text(encoding="utf-8")
            self.assertIn("Detailed calls — all input profiles", report)
            self.assertIn("review_recommended", report)
            self.assertNotIn("Priority review queue", report)


if __name__ == "__main__":
    unittest.main()
