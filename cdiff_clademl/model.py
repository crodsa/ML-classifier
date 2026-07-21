"""
Stage 1 pipeline construction, Cryptic gate threshold search, and gate application.
"""

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder

from .config import TOKEN_COLS
from .utils import vprint


def build_stage1_pipeline(n_jobs: int = -1) -> Pipeline:
    """Build the Stage 1 sklearn Pipeline (OHE tokenised alleles → RandomForest).

    The pipeline:
        pre: ColumnTransformer with OneHotEncoder on TOKEN_COLS (unknown categories
             handled silently at inference via handle_unknown='ignore').
        clf: RandomForestClassifier with 500 trees, balanced class weighting,
             and all CPU cores.

    Returns:
        Unfitted sklearn Pipeline ready for .fit() and .predict_proba().
    """
    pre = ColumnTransformer(
        transformers=[("cat", OneHotEncoder(handle_unknown="ignore"), TOKEN_COLS)],
        remainder="drop",
    )
    clf = RandomForestClassifier(
        n_estimators=500,
        random_state=42,
        n_jobs=n_jobs,
        class_weight="balanced_subsample",
    )
    return Pipeline([("pre", pre), ("clf", clf)])


def threshold_scan_distance_stage1(
    meta: pd.DataFrame,
    y_true_bin: np.ndarray,
    min_cryptic_recall: float,
    grid_size: int = 25,
    verbose: bool = False,
    tag: str = "",
) -> dict[str, float]:
    """Grid-search the best Cryptic gate thresholds on OOF meta-features.

    Searches over four independent threshold dimensions:
        sum_min_dist_any_classical   >= t1   (distance signal, OR)
        nearest_profile_dist_any_classical >= t2   (profile distance, OR)
        stage1_confidence            <= t3   (low Stage 1 confidence, OR)
        stage1_margin                <= t4   (low margin between top-2 classes)

    A sample is called Cryptic when (dist OR profile) AND (low_conf OR low_margin)
    AND n_missing_loci <= max_missing.

    Selection criterion: among configurations that meet min_cryptic_recall,
    minimize false positives first, then maximize balanced accuracy.

    Args:
        meta:               OOF meta-feature DataFrame.
        y_true_bin:         Binary array (1 = Cryptic, 0 = classical).
        min_cryptic_recall: Minimum acceptable Cryptic recall (e.g. 0.90).
        verbose:            Log search progress when True.
        tag:                Label for verbose messages.

    Returns:
        Dict with keys: sum_min_dist_any_threshold, nearest_profile_any_threshold,
        stage1_confidence_max, stage1_margin_max, max_missing.
    """
    if not 0 < min_cryptic_recall <= 1:
        raise ValueError("min_cryptic_recall must be in the interval (0, 1].")
    if not 5 <= grid_size <= 100:
        raise ValueError("grid_size must be between 5 and 100.")
    y_true_bin = np.asarray(y_true_bin)
    if len(y_true_bin) != len(meta) or not set(np.unique(y_true_bin)).issubset({0, 1}):
        raise ValueError("y_true_bin must be a binary vector aligned with meta.")
    if not (y_true_bin == 1).any() or not (y_true_bin == 0).any():
        raise ValueError("Threshold tuning requires both Cryptic and classical samples.")

    def quant_grid(vals: np.ndarray) -> np.ndarray:
        if len(vals) <= grid_size:
            return vals
        return np.unique(np.quantile(vals, np.linspace(0, 1, grid_size)))

    sum_grid = quant_grid(np.unique(meta["sum_min_dist_any_classical"].to_numpy()))
    near_grid = quant_grid(np.unique(meta["nearest_profile_dist_any_classical"].to_numpy()))
    conf_grid = quant_grid(np.unique(meta["stage1_confidence"].to_numpy()))
    margin_grid = quant_grid(np.unique(meta["stage1_margin"].to_numpy()))

    vprint(
        f"threshold_scan {tag}: grids "
        f"sum={len(sum_grid)} near={len(near_grid)} "
        f"conf={len(conf_grid)} margin={len(margin_grid)}",
        verbose,
    )

    n_cryptic = int((y_true_bin == 1).sum())
    n_classical = int((y_true_bin == 0).sum())
    best = None
    outer_total = max(1, len(sum_grid) * len(near_grid))
    outer_i = 0

    sum_arr = meta["sum_min_dist_any_classical"].to_numpy()
    near_arr = meta["nearest_profile_dist_any_classical"].to_numpy()
    conf_arr = meta["stage1_confidence"].to_numpy()
    margin_arr = meta["stage1_margin"].to_numpy()
    miss_arr = meta["n_missing_loci"].to_numpy()

    for t1 in sum_grid:
        dist1 = sum_arr >= t1
        for t2 in near_grid:
            outer_i += 1
            if verbose and (outer_i == 1 or outer_i % 50 == 0 or outer_i == outer_total):
                vprint(f"threshold_scan {tag}: outer {outer_i}/{outer_total}", verbose)
            dist = dist1 | (near_arr >= t2)
            if not dist.any():
                continue
            for t3 in conf_grid:
                low1 = conf_arr <= t3
                for t4 in margin_grid:
                    low = low1 | (margin_arr <= t4)
                    pred = dist & low & (miss_arr <= 1)
                    tp = int(((pred == 1) & (y_true_bin == 1)).sum())
                    fp = int(((pred == 1) & (y_true_bin == 0)).sum())
                    fn = int(((pred == 0) & (y_true_bin == 1)).sum())
                    recall = tp / max(1, n_cryptic)
                    if recall < min_cryptic_recall:
                        continue
                    precision = tp / max(1, tp + fp)
                    f1 = (
                        0.0
                        if precision + recall == 0
                        else 2 * precision * recall / (precision + recall)
                    )
                    tn = int(((pred == 0) & (y_true_bin == 0)).sum())
                    spec = tn / max(1, n_classical)
                    bal = 0.5 * (recall + spec)
                    key = (fp, fn, -bal, -f1, float(t1 + t2), float(t3 + t4))
                    if best is None or key < best[0]:
                        best = (
                            key,
                            {
                                "sum_min_dist_any_threshold": float(t1),
                                "nearest_profile_any_threshold": float(t2),
                                "stage1_confidence_max": float(t3),
                                "stage1_margin_max": float(t4),
                                "max_missing": 1,
                            },
                        )

    if best is None:
        raise ValueError(
            "No Cryptic-gate thresholds satisfy the requested minimum recall. "
            "Inspect missing loci and lower --min_cryptic_recall only with a "
            "documented scientific rationale."
        )

    vprint(f"threshold_scan {tag}: selected rules {best[1]}", verbose)
    return best[1]


def apply_v5_gate(meta: pd.DataFrame, rules: dict[str, float]) -> np.ndarray:
    """Apply the Cryptic gate rules to a meta-feature DataFrame.

    A sample is called Cryptic (1) when:
        (sum_min_dist >= t1 OR nearest_profile_dist >= t2)
        AND (stage1_confidence <= t3 OR stage1_margin <= t4)
        AND n_missing_loci <= max_missing

    Args:
        meta:  Meta-feature DataFrame (from build_meta_features).
        rules: Dict of thresholds from threshold_scan_distance_stage1.

    Returns:
        Integer array of shape (n_samples,) with 1 = Cryptic, 0 = classical.
    """
    dist = (
        meta["sum_min_dist_any_classical"].to_numpy() >= rules["sum_min_dist_any_threshold"]
    ) | (
        meta["nearest_profile_dist_any_classical"].to_numpy()
        >= rules["nearest_profile_any_threshold"]
    )
    low = (meta["stage1_confidence"].to_numpy() <= rules["stage1_confidence_max"]) | (
        meta["stage1_margin"].to_numpy() <= rules["stage1_margin_max"]
    )
    ok_missing = meta["n_missing_loci"].to_numpy() <= int(rules.get("max_missing", 1))
    return (dist & low & ok_missing).astype(int)
