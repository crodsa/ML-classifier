"""Regression checks for release version and model-rights metadata."""

from __future__ import annotations

import json
import re
import unittest
from pathlib import Path


class ReleaseMetadataTests(unittest.TestCase):
    def test_release_versions_agree_across_package_metadata(self) -> None:
        root = Path(__file__).resolve().parents[1]
        package_init = (root / "cdiff_clademl" / "__init__.py").read_text(encoding="utf-8")
        pyproject = (root / "pyproject.toml").read_text(encoding="utf-8")
        citation = (root / "CITATION.cff").read_text(encoding="utf-8")
        version_match = re.search(r'^VERSION = "([^"]+)"$', package_init, re.MULTILINE)
        self.assertIsNotNone(version_match)
        version = version_match.group(1)

        self.assertRegex(
            pyproject, re.compile(rf'^version = "{re.escape(version)}"$', re.MULTILINE)
        )
        self.assertRegex(citation, re.compile(rf"^version: {re.escape(version)}$", re.MULTILINE))

    def test_model_card_records_authorized_distribution_boundary(self) -> None:
        root = Path(__file__).resolve().parents[1]
        card = json.loads((root / "model_model_card.json").read_text(encoding="utf-8"))
        distribution = card["distribution"]

        self.assertEqual(distribution["model_and_distance_matrices_license"], "MIT")
        self.assertFalse(distribution["training_tables_distributed"])
        self.assertIn("repository owner", distribution["rights_statement"])


if __name__ == "__main__":
    unittest.main()
