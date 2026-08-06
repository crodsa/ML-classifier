import tempfile
import unittest
from pathlib import Path

import pandas as pd

from cdiff_clademl.report import write_html_report


class HtmlReportTests(unittest.TestCase):
    def test_report_renders_review_and_escapes_identifiers(self):
        out = pd.DataFrame(
            {
                "isolate": ["sample <1>", "sample 2"],
                "st": ["new_ST", "4"],
                "final_prediction": ["Cryptic", "Clade4"],
                "final_reason": ["distance+low_classical_support", "classical_stage1"],
                "quality_flag": ["UNKNOWN_ALLELE", "PASS"],
                "stage1_support_band": ["LOW", "HIGH"],
                "confidence_band": ["LOW", "HIGH"],
                "stage1_confidence": [0.23, 0.94],
                "n_missing_loci": [0, 0],
                "n_approx_loci": [0, 0],
                "n_unknown_alleles": [1, 0],
                "sum_min_dist_any_classical": [0.2, 0.003],
                "nearest_profile_dist_any_classical": [0.3, 0.004],
                "review_recommended": [True, False],
                "allelic_profile_note": ["review <this>", "complete profile"],
            }
        )
        summary = {
            "class_counts": {"Cryptic": 1, "Clade4": 1},
            "quality_flags": {"UNKNOWN_ALLELE": 1, "PASS": 1},
            "stage1_support_bands": {"LOW": 1, "HIGH": 1},
            "n_review_recommended": 1,
            "n_cryptic_gate": 1,
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "report.html"
            write_html_report(path, out, summary, "test-model", "abc")
            report = path.read_text(encoding="utf-8")
            self.assertIn("Priority review queue", report)
            self.assertIn("Detailed calls — all input profiles", report)
            self.assertIn("sample &lt;1&gt;", report)
            self.assertIn("sample 2", report)
            self.assertIn("review &lt;this&gt;", report)
            self.assertIn("Clade4", report)
            self.assertIn("not validated for clinical diagnosis", report)
            isolate_header = report.find("<th>isolate</th>")
            final_header = report.find("<th>final_prediction</th>")
            support_header = report.find("<th>stage1_confidence</th>")
            self.assertLess(isolate_header, final_header)
            self.assertLess(final_header, support_header)


if __name__ == "__main__":
    unittest.main()
