import tempfile
import unittest
from pathlib import Path

import pandas as pd

from cdiff_clademl.report import write_html_report


class HtmlReportTests(unittest.TestCase):
    def test_report_renders_review_and_escapes_identifiers(self):
        out = pd.DataFrame(
            {
                "isolate": ["sample <1>"],
                "st": ["new_ST"],
                "final_prediction": ["Cryptic"],
                "final_reason": ["distance+low_classical_support"],
                "quality_flag": ["UNKNOWN_ALLELE"],
                "stage1_support_band": ["LOW"],
                "stage1_confidence": [0.23],
                "n_missing_loci": [0],
                "n_approx_loci": [0],
                "n_unknown_alleles": [1],
                "sum_min_dist_any_classical": [0.2],
                "nearest_profile_dist_any_classical": [0.3],
                "review_recommended": [True],
            }
        )
        summary = {
            "class_counts": {"Cryptic": 1},
            "quality_flags": {"UNKNOWN_ALLELE": 1},
            "stage1_support_bands": {"LOW": 1},
            "n_review_recommended": 1,
            "n_cryptic_gate": 1,
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "report.html"
            write_html_report(path, out, summary, "test-model", "abc")
            report = path.read_text(encoding="utf-8")
            self.assertIn("Priority review queue", report)
            self.assertIn("sample &lt;1&gt;", report)
            self.assertIn("not validated for clinical diagnosis", report)


if __name__ == "__main__":
    unittest.main()
