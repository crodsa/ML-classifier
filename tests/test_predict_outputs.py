"""Output-path tests for the one-command prediction workflow."""

from __future__ import annotations

import tempfile
import unittest
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import patch

from cdiff_clademl.cli.predict import (
    default_output_directory,
    prepare_output_paths,
    resolve_output_paths,
)


class PredictionOutputTests(unittest.TestCase):
    def test_default_output_directory_is_derived_from_input_and_utc_time(self) -> None:
        directory = default_output_directory(
            "inputs/profiles.tsv",
            now=datetime(2026, 8, 9, 12, 30, 45, tzinfo=timezone.utc),
        )

        self.assertEqual(directory, Path("profiles_cdiffclademl_20260809T123045Z"))

    def test_automatic_paths_include_every_routine_artifact(self) -> None:
        predictions, summary, report, automatic = resolve_output_paths(
            "profiles.tsv",
            None,
            None,
            now=datetime(2026, 8, 9, 12, 30, 45, tzinfo=timezone.utc),
        )

        self.assertTrue(automatic)
        self.assertEqual(
            predictions, Path("profiles_cdiffclademl_20260809T123045Z/predictions.tsv")
        )
        self.assertEqual(
            summary, Path("profiles_cdiffclademl_20260809T123045Z/predictions_summary.json")
        )
        self.assertEqual(report, Path("profiles_cdiffclademl_20260809T123045Z/report.html"))

    def test_explicit_prediction_path_keeps_report_beside_it(self) -> None:
        predictions, summary, report, automatic = resolve_output_paths(
            "profiles.tsv", "results/run.tsv", None
        )

        self.assertFalse(automatic)
        self.assertEqual(predictions, Path("results/run.tsv"))
        self.assertEqual(summary, Path("results/run_summary.json"))
        self.assertEqual(report, Path("results/report.html"))

    def test_explicit_report_path_keeps_predictions_beside_it(self) -> None:
        predictions, summary, report, automatic = resolve_output_paths(
            "profiles.tsv", None, "results/review.html"
        )

        self.assertFalse(automatic)
        self.assertEqual(predictions, Path("results/predictions.tsv"))
        self.assertEqual(summary, Path("results/predictions_summary.json"))
        self.assertEqual(report, Path("results/review.html"))

    def test_automatic_directory_collision_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            run_directory = Path(temporary) / "existing_run"
            run_directory.mkdir()
            with patch(
                "cdiff_clademl.cli.predict.default_output_directory",
                return_value=run_directory,
            ):
                with self.assertRaisesRegex(ValueError, "already exists"):
                    prepare_output_paths("profiles.tsv", None, None)

    def test_explicit_existing_artifact_is_not_overwritten(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            prediction_path = Path(temporary) / "predictions.tsv"
            prediction_path.write_text("existing\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "Refusing to overwrite"):
                prepare_output_paths("profiles.tsv", str(prediction_path), None)


if __name__ == "__main__":
    unittest.main()
