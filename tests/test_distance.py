import tempfile
import unittest
from pathlib import Path

import pandas as pd

from cdiff_clademl.config import LOCI
from cdiff_clademl.distance import (
    distance_to_allele,
    load_distance_stores,
    stores_for_reference,
)


class DistanceTests(unittest.TestCase):
    def _classical_df(self, allele: str = "1") -> pd.DataFrame:
        return pd.DataFrame({f"{locus}_base": [allele] for locus in LOCI})

    def _write_matrices(self, directory: Path, distance: float = 0.1) -> None:
        text = (
            f"allele_1\tallele_2\tdistance\n1\t1\t0\n1\t2\t{distance}\n2\t1\t{distance}\n2\t2\t0\n"
        )
        for locus in LOCI:
            (directory / f"{locus}.tsv").write_text(text, encoding="utf-8")

    def test_load_and_fold_specific_references(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            self._write_matrices(directory)
            stores = load_distance_stores(str(directory), self._classical_df("1"))
            self.assertEqual(stores["adk"].matrix[("1", "1")], 0.0)
            self.assertAlmostEqual(
                distance_to_allele("adk", "1", "exact", "2", stores["adk"]),
                0.1,
            )
            fold = stores_for_reference(stores, self._classical_df("2"))
            self.assertEqual(fold["adk"].classical_alleles, {"2"})
            self.assertEqual(stores["adk"].classical_alleles, {"1"})

    def test_negative_distance_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            self._write_matrices(directory)
            (directory / "adk.tsv").write_text(
                "allele_1\tallele_2\tdistance\n1\t1\t-1\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "negative or non-finite"):
                load_distance_stores(str(directory), self._classical_df("1"))


if __name__ == "__main__":
    unittest.main()
