"""
Stage 1 feature engineering: Shannon entropy, per-row meta-feature building,
and the build_meta_features() assembler used by both train and predict.
"""

import numpy as np
import pandas as pd

from .config import LOCI
from .distance import DistanceStore, best_profile_distance, min_distance_to_set
from .utils import vprint


def shannon_entropy(probs: np.ndarray) -> np.ndarray:
    """Compute Shannon entropy for each row of a probability matrix.

    Args:
        probs: (n_samples, n_classes) probability array.

    Returns:
        1-D array of length n_samples with per-row entropy values.
    """
    probs = np.clip(probs, 1e-12, 1.0)
    return -np.sum(probs * np.log(probs), axis=1)


def build_meta_features(
    query_df: pd.DataFrame,
    stage1_pred: np.ndarray,
    stage1_proba: np.ndarray,
    classical_ref_df: pd.DataFrame,
    stores: dict[str, DistanceStore],
    stage1_classes: list[str],
    verbose: bool = False,
    tag: str = "",
) -> pd.DataFrame:
    """Build Stage 2 meta-features from Stage 1 predictions and distance computations.

    Meta-features
    -------------
    stage1_predicted_clade:          Predicted class from Stage 1 RF.
    stage1_confidence:               Max class probability (confidence).
    stage1_margin:                   Gap between top-1 and top-2 probabilities.
    stage1_entropy:                  Shannon entropy of the probability vector.
    sum_min_dist_any_classical:      Sum across loci of min distance to any
                                     classical allele in the training set.
    nearest_profile_dist_any_classical: Min distance to any unique classical
                                     allele profile (vectorized; see distance.py).
    n_missing_loci:                  Count of loci with missing alleles.
    rare_locus_count_report:         Count of loci with non-zero distance to
                                     nearest classical allele.

    Args:
        query_df:         MLST DataFrame with {locus}_base/status columns.
        stage1_pred:      Stage 1 predicted class array (n_samples,).
        stage1_proba:     Stage 1 class probability matrix (n_samples, n_classes).
        classical_ref_df: Collapsed unique classical profiles from training.
        stores:           Dict of DistanceStore per locus.
        stage1_classes:   Ordered list of class names matching stage1_proba columns.
        verbose:          Log progress when True.
        tag:              Label for verbose log messages.

    Returns:
        DataFrame indexed like query_df with all meta-feature columns.
    """
    n = len(query_df)
    stage1_proba = np.asarray(stage1_proba, dtype=float)
    stage1_pred = np.asarray(stage1_pred)
    if len(stage1_pred) != n or stage1_proba.shape != (n, len(stage1_classes)):
        raise ValueError("Stage 1 prediction dimensions do not match the query table/classes.")
    if not np.isfinite(stage1_proba).all() or (stage1_proba < 0).any():
        raise ValueError("Stage 1 produced invalid probabilities.")
    if n and not np.allclose(stage1_proba.sum(axis=1), 1.0, atol=1e-6):
        raise ValueError("Stage 1 probabilities do not sum to one.")
    vprint(
        f"build_meta_features {tag}: n_query={n}, n_classical_ref_unique={len(classical_ref_df)}",
        verbose,
    )

    # ── Vectorized nearest-profile distance ──────────────────────────────────
    nearest_profile = best_profile_distance(
        query_df,
        classical_ref_df,
        stores,
        verbose=verbose,
        tag=tag,
    )

    # ── Per-row Stage 1 stats and per-locus distance features ────────────────
    sorted_proba = np.sort(stage1_proba, axis=1)[:, ::-1]  # descending
    confidences = sorted_proba[:, 0]
    margins = sorted_proba[:, 0] - (sorted_proba[:, 1] if stage1_proba.shape[1] > 1 else 0)
    entropies = shannon_entropy(stage1_proba)

    rows = []
    q_base_arr = {locus: query_df[f"{locus}_base"].to_numpy(dtype=str) for locus in LOCI}
    q_status_arr = {locus: query_df[f"{locus}_status"].to_numpy(dtype=str) for locus in LOCI}
    min_dist_cache: dict[str, dict[tuple[str, str], float]] = {locus: {} for locus in LOCI}

    for i in range(n):
        sum_min = 0.0
        n_missing = 0
        n_approx = 0
        n_unknown = 0
        rare_count = 0
        for locus in LOCI:
            base, status = q_base_arr[locus][i], q_status_arr[locus][i]
            if status == "missing" or base == "__MISSING__":
                n_missing += 1
            elif status == "approx":
                n_approx += 1
            if base != "__MISSING__" and base not in stores[locus].allele_set:
                n_unknown += 1
            locus_cache = min_dist_cache[locus]
            key = (base, status)
            d = locus_cache.get(key)
            if d is None:
                d = min_distance_to_set(
                    locus,
                    base,
                    status,
                    stores[locus].classical_alleles,
                    stores[locus],
                )
                locus_cache[key] = d
            sum_min += d
            if d > 0:
                rare_count += 1
        rows.append(
            {
                "stage1_predicted_clade": stage1_pred[i],
                "stage1_confidence": float(confidences[i]),
                "stage1_margin": float(margins[i]),
                "stage1_entropy": float(entropies[i]),
                "sum_min_dist_any_classical": float(sum_min),
                "n_missing_loci": int(n_missing),
                "n_approx_loci": int(n_approx),
                "n_unknown_alleles": int(n_unknown),
                "rare_locus_count_report": int(rare_count),
            }
        )
        if verbose and (i == 0 or (i + 1) % 250 == 0 or i == n - 1):
            vprint(f"build_meta_features {tag}: processed {i + 1}/{n}", verbose)

    out = pd.DataFrame(rows, index=query_df.index)
    out["nearest_profile_dist_any_classical"] = nearest_profile[
        "nearest_profile_dist_any_classical"
    ].values
    return out
