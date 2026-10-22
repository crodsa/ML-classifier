import unittest

import numpy as np
import pandas as pd

from cdiff_clademl.model import apply_cryptic_gate, threshold_scan_distance_stage1


class GateTests(unittest.TestCase):
    def test_threshold_search_separates_simple_cases(self):
        meta = pd.DataFrame(
            {
                "sum_min_dist_any_classical": [0.0, 0.01, 1.0, 1.2],
                "nearest_profile_dist_any_classical": [0.0, 0.01, 1.0, 1.2],
                "stage1_confidence": [0.95, 0.9, 0.3, 0.2],
                "stage1_margin": [0.9, 0.8, 0.1, 0.05],
                "n_missing_loci": [0, 0, 0, 0],
            }
        )
        truth = np.array([0, 0, 1, 1])
        rules = threshold_scan_distance_stage1(meta, truth, min_cryptic_recall=1.0, grid_size=5)
        np.testing.assert_array_equal(apply_cryptic_gate(meta, rules), truth)

    def test_threshold_search_validates_inputs(self):
        meta = pd.DataFrame(
            {
                "sum_min_dist_any_classical": [0.0],
                "nearest_profile_dist_any_classical": [0.0],
                "stage1_confidence": [1.0],
                "stage1_margin": [1.0],
                "n_missing_loci": [0],
            }
        )
        with self.assertRaisesRegex(ValueError, "both Cryptic and classical"):
            threshold_scan_distance_stage1(meta, np.array([0]), 0.9)


if __name__ == "__main__":
    unittest.main()
