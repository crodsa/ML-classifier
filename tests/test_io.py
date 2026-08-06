import tempfile
import unittest
from pathlib import Path

import numpy as np

from cdiff_clademl.io import (
    evaluation_group_ids,
    parse_allele_field,
    read_mlst_table,
    write_grouped_eval_files,
)


class MlstIoTests(unittest.TestCase):
    def _table(self, text: str, require_labels: bool = False):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "input.csv"
            path.write_text(text, encoding="utf-8")
            return read_mlst_table(str(path), require_labels=require_labels)

    def test_allele_parser(self):
        self.assertEqual(parse_allele_field("adk(12)"), ("12", "exact"))
        self.assertEqual(parse_allele_field("adk(~12)"), ("12", "approx"))
        self.assertEqual(parse_allele_field("N/A"), ("__MISSING__", "missing"))
        self.assertEqual(parse_allele_field("~"), ("__MISSING__", "missing"))

    def test_header_and_unlabelled_input(self):
        text = (
            "isolate,species,st,adk,atpA,dxr,glyA,recA,sodA,tpi\n"
            "x,cdifficile,1,adk(1),1,2,3,4,5,6\n"
        )
        df = self._table(text)
        self.assertEqual(len(df), 1)
        self.assertNotIn("mlst_clade", df)
        self.assertEqual(df.loc[0, "adk_token"], "E:1")

    def test_training_rejects_missing_or_unknown_label(self):
        row10 = "x,cdifficile,1,1,1,1,1,1,1,1\n"
        with self.assertRaisesRegex(ValueError, "11-column"):
            self._table(row10, require_labels=True)

        row11 = row10.rstrip() + ",mlst_clade(Clade9)\n"
        with self.assertRaisesRegex(ValueError, "Unknown mlst_clade"):
            self._table(row11, require_labels=True)

    def test_parenthesized_missing_label_is_treated_as_unlabelled(self):
        row11 = "x,cdifficile,1,1,1,1,1,1,1,1,mlst_clade(nan)\n"
        df = self._table(row11)
        self.assertEqual(df.loc[0, "mlst_clade"], "")

    def test_rejects_malformed_width(self):
        with self.assertRaisesRegex(ValueError, "10 or 11 columns"):
            self._table("x,cdifficile,1,1,1\n")

    def test_grouped_metrics_require_consistent_predictions(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = str(Path(tmp) / "metrics")
            with self.assertRaisesRegex(ValueError, "Inconsistent"):
                write_grouped_eval_files(
                    prefix,
                    np.array(["Clade1", "Clade1"]),
                    np.array(["Clade1", "Clade2"]),
                    np.array(["1", "1"]),
                    ["Clade1", "Clade2"],
                )

    def test_new_st_values_are_separate_evaluation_groups(self):
        groups, n_unassigned = evaluation_group_ids(
            np.array(["new_ST", "1000", "new_ST"]), np.array(["a", "b", "c"])
        )
        self.assertEqual(n_unassigned, 2)
        self.assertEqual(groups.nunique(), 3)


if __name__ == "__main__":
    unittest.main()
