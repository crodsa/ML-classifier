"""
Training entry point for the two-stage C. difficile MLST clade classifier.
"""

import argparse
import json
import os
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import sklearn
from sklearn.metrics import average_precision_score
from sklearn.model_selection import GroupKFold, StratifiedGroupKFold

from cdiff_clademl import PROGRAM, VERSION
from cdiff_clademl.bundle import BUNDLE_SCHEMA_VERSION, sha256_file
from cdiff_clademl.config import ALL_CLADES, CLASSICAL_CLADES
from cdiff_clademl.distance import (
    collapse_unique_profiles,
    load_distance_stores,
    stores_for_reference,
)
from cdiff_clademl.features import build_meta_features
from cdiff_clademl.io import (
    read_mlst_table,
    summarize_multiclass,
    validate_training_table,
    write_eval_files,
    write_grouped_eval_files,
    write_model_card,
)
from cdiff_clademl.model import (
    apply_cryptic_gate,
    build_stage1_pipeline,
    threshold_scan_distance_stage1,
)
from cdiff_clademl.utils import setup_logging, vprint

# ==============================================================================
# Parallel fold worker (must be top-level for ProcessPoolExecutor pickling)
# ==============================================================================


def _fold_worker(args: tuple) -> pd.DataFrame:
    """Train Stage 1 on one CV fold and return OOF meta-features.

    This function is called in a worker process by ProcessPoolExecutor.
    All arguments must be picklable (DataFrames, dicts, dataclasses with
    primitive fields). Logging is suppressed inside workers to avoid
    interleaved stderr output.

    Args:
        args: Tuple of fold number, train/validation data, distance stores,
              RF worker count, and verbose flag.

    Returns:
        OOF meta-feature DataFrame for the validation split of this fold,
        with columns 'true_label' and 'is_cryptic' appended.
    """
    fold, tr, va, stores, rf_n_jobs, verbose = args
    tag = f"[fold {fold}/5]"
    stage1 = build_stage1_pipeline(n_jobs=rf_n_jobs)
    stage1.fit(tr, tr["mlst_clade"])
    probs = stage1.predict_proba(va)
    classes = list(stage1.named_steps["clf"].classes_)
    pred = stage1.predict(va)
    fold_stores = stores_for_reference(stores, tr)
    meta = build_meta_features(
        va,
        pred,
        probs,
        collapse_unique_profiles(tr),
        fold_stores,
        classes,
        verbose=verbose,
        tag=tag,
    )
    meta["true_label"] = va["mlst_clade"].values
    meta["is_cryptic"] = meta["true_label"].eq("Cryptic").astype(int)
    meta["st"] = va["st"].values
    meta["isolate"] = va["isolate"].values
    meta["cv_fold"] = fold
    return meta


# ==============================================================================
# Main
# ==============================================================================


def main() -> None:
    ap = argparse.ArgumentParser(
        prog="CdiffCladeML-train",
        description=f"{PROGRAM} v{VERSION} — Train the two-stage MLST clade classifier",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    ap.add_argument("training_table", help="11-column MLST table with mlst_clade labels.")
    ap.add_argument(
        "--distance_dir", required=True, help="Directory of per-locus distance TSVs ({locus}.tsv)."
    )
    ap.add_argument(
        "-o", "--output_model", required=True, help="Path for the output .joblib model bundle."
    )
    ap.add_argument(
        "--metrics_prefix", default=None, help="File-path prefix for OOF evaluation outputs."
    )
    ap.add_argument(
        "--min_cryptic_recall",
        type=float,
        default=0.90,
        help="Minimum Cryptic recall enforced during threshold search.",
    )
    ap.add_argument(
        "--threshold_grid_size",
        type=int,
        default=25,
        help="Quantile candidates per Cryptic-gate threshold dimension (5-100).",
    )
    ap.add_argument(
        "-j",
        "--n_jobs",
        type=int,
        default=1,
        help=(
            "Number of CV folds to run in parallel. "
            "Max useful value is 5 (one per fold). "
            "Note: each fold already uses n_jobs=-1 inside RF; "
            "set this to 1 on memory-constrained machines."
        ),
    )
    ap.add_argument("--verbose", action="store_true", help="Verbose logging (INFO level).")
    args = ap.parse_args()
    setup_logging(args.verbose)

    if args.n_jobs < 1:
        ap.error("--n_jobs must be at least 1.")
    if not 0 < args.min_cryptic_recall <= 1:
        ap.error("--min_cryptic_recall must be in (0, 1].")
    if not 5 <= args.threshold_grid_size <= 100:
        ap.error("--threshold_grid_size must be between 5 and 100.")

    t0 = time.time()
    vprint(f"{PROGRAM} v{VERSION} — training started", args.verbose)

    # ── Load data ─────────────────────────────────────────────────────────────
    vprint("Loading training data...", args.verbose)
    df = read_mlst_table(args.training_table, require_labels=True)
    validate_training_table(df, n_splits=5)
    classical_df = df[df["mlst_clade"].isin(CLASSICAL_CLADES)].copy()
    cryptic_df = df[df["mlst_clade"].eq("Cryptic")].copy()
    vprint(
        f"Loaded {len(df)} samples: {len(classical_df)} classical, {len(cryptic_df)} cryptic.",
        args.verbose,
    )

    stores = load_distance_stores(
        args.distance_dir,
        classical_df,
        verbose=args.verbose,
    )

    # ── Stage 1 OOF cross-validation ──────────────────────────────────────────
    sgkf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
    cryptic_gkf = GroupKFold(n_splits=5, shuffle=True, random_state=42)
    cryptic_splits = list(cryptic_gkf.split(cryptic_df, groups=cryptic_df["st"]))

    n_workers = min(args.n_jobs, 5)
    # Avoid multiplying all-core RF jobs by five process workers.
    fold_rf_n_jobs = -1 if n_workers == 1 else 1

    fold_args: list[tuple] = []
    classical_splits = sgkf.split(
        classical_df,
        classical_df["mlst_clade"],
        groups=classical_df["st"],
    )
    for fold, ((tr_idx, va_idx), (_, cryptic_va_idx)) in enumerate(
        zip(classical_splits, cryptic_splits, strict=True),
        start=1,
    ):
        tr = classical_df.iloc[tr_idx].copy()
        va_class = classical_df.iloc[va_idx].copy()
        va_cryptic = cryptic_df.iloc[cryptic_va_idx].copy()
        va = pd.concat([va_class, va_cryptic], axis=0).sort_index()
        fold_args.append((fold, tr, va, stores, fold_rf_n_jobs, args.verbose))

    vprint(
        f"Stage 1 OOF CV: 5 folds, {n_workers} parallel worker(s)...",
        args.verbose,
    )

    oof_meta_parts: list[pd.DataFrame] = [None] * 5  # preserve fold order

    if n_workers <= 1:
        for fa in fold_args:
            fold_num = fa[0]
            vprint(f"  Fold {fold_num}/5 starting...", args.verbose)
            meta = _fold_worker(fa)
            oof_meta_parts[fold_num - 1] = meta
            vprint(f"  Fold {fold_num}/5 complete (n_val={len(meta)}).", args.verbose)
    else:
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            future_to_fold = {pool.submit(_fold_worker, fa): fa[0] for fa in fold_args}
            done = 0
            for future in as_completed(future_to_fold):
                fold_num = future_to_fold[future]
                done += 1
                try:
                    meta = future.result()
                except Exception as exc:
                    print(
                        f"ERROR: fold {fold_num} worker raised: {exc}",
                        file=sys.stderr,
                    )
                    sys.exit(1)
                oof_meta_parts[fold_num - 1] = meta
                vprint(
                    f"  [{done}/5] Fold {fold_num} complete (n_val={len(meta)}).",
                    args.verbose,
                )

    oof_meta = pd.concat(oof_meta_parts, axis=0).sort_index()
    if len(oof_meta) != len(df) or not oof_meta.index.is_unique:
        raise RuntimeError(
            "Internal CV error: every training row must receive exactly one OOF prediction."
        )
    vprint(f"OOF meta-features assembled: {len(oof_meta)} rows.", args.verbose)

    # ── Stage 2: threshold search on OOF features ─────────────────────────────
    gate_y = oof_meta["is_cryptic"].to_numpy()
    ap_score = float(
        average_precision_score(
            gate_y,
            -oof_meta["stage1_confidence"].to_numpy()
            + oof_meta["sum_min_dist_any_classical"].to_numpy(),
        )
    )
    vprint(f"Stage 2 OOF proxy PR-AUC = {ap_score:.4f}", args.verbose)

    t_tune = time.time()
    rules = threshold_scan_distance_stage1(
        oof_meta,
        gate_y,
        args.min_cryptic_recall,
        grid_size=args.threshold_grid_size,
        verbose=args.verbose,
        tag="[OOF tuning]",
    )
    vprint(f"Threshold tuning: {time.time() - t_tune:.1f}s", args.verbose)

    # Reduced-bias gate diagnostic: each fold is predicted with thresholds
    # tuned on the other four OOF folds. This avoids direct resubstitution but
    # is not a fully nested CV of the complete two-stage system; external
    # validation remains the primary performance estimate. Final deployment
    # rules are tuned above on all OOF features and stored separately.
    crossfit_call = np.zeros(len(oof_meta), dtype=int)
    crossfit_rules = {}
    fold_values = oof_meta["cv_fold"].to_numpy()
    for held_out_fold in range(1, 6):
        tune_mask = fold_values != held_out_fold
        eval_mask = ~tune_mask
        fold_rules = threshold_scan_distance_stage1(
            oof_meta.loc[tune_mask],
            oof_meta.loc[tune_mask, "is_cryptic"].to_numpy(),
            args.min_cryptic_recall,
            grid_size=args.threshold_grid_size,
            verbose=args.verbose,
            tag=f"[gate cross-fit {held_out_fold}/5]",
        )
        crossfit_rules[str(held_out_fold)] = fold_rules
        crossfit_call[eval_mask] = apply_cryptic_gate(oof_meta.loc[eval_mask], fold_rules)

    final_pred_crossfit = np.where(
        crossfit_call == 1,
        "Cryptic",
        oof_meta["stage1_predicted_clade"].to_numpy(),
    )

    # In-sample application of final rules is retained only as a tuning
    # diagnostic and must not be reported as validation performance.
    final_rules_call = apply_cryptic_gate(oof_meta, rules)

    # ── OOF evaluation outputs ────────────────────────────────────────────────
    if args.metrics_prefix:
        import os as _os

        _prefix_dir = _os.path.dirname(args.metrics_prefix)
        if _prefix_dir:
            _os.makedirs(_prefix_dir, exist_ok=True)
        oof_out = oof_meta.copy()
        oof_out["cryptic_crossfit_call"] = crossfit_call
        oof_out["final_prediction_crossfit"] = final_pred_crossfit
        oof_out["cryptic_final_rules_call"] = final_rules_call
        oof_out.to_csv(f"{args.metrics_prefix}_oof_meta.tsv", sep="\t", index=True)
        write_eval_files(
            args.metrics_prefix,
            oof_meta["true_label"].to_numpy(),
            final_pred_crossfit,
            ALL_CLADES,
        )
        write_grouped_eval_files(
            f"{args.metrics_prefix}_by_st",
            oof_meta["true_label"].to_numpy(),
            final_pred_crossfit,
            oof_meta["st"].to_numpy(),
            ALL_CLADES,
        )
        vprint("OOF evaluation files written.", args.verbose)

    # ── Final Stage 1 fit on all classical samples ────────────────────────────
    vprint("Fitting final Stage 1 on full training set...", args.verbose)
    t_fit = time.time()
    stage1_final = build_stage1_pipeline()
    stage1_final.fit(classical_df, classical_df["mlst_clade"])
    vprint(f"Final Stage 1 fit: {time.time() - t_fit:.1f}s", args.verbose)

    # ── Assemble bundle and write outputs ──────────────────────────────────────
    crossfit_metrics = summarize_multiclass(
        oof_meta["true_label"].to_numpy(), final_pred_crossfit, labels=ALL_CLADES
    )
    one_per_st = pd.DataFrame(
        {
            "st": oof_meta["st"].to_numpy(),
            "true": oof_meta["true_label"].to_numpy(),
            "pred": final_pred_crossfit,
        }
    ).drop_duplicates("st")
    crossfit_metrics_by_st = summarize_multiclass(
        one_per_st["true"].to_numpy(),
        one_per_st["pred"].to_numpy(),
        labels=ALL_CLADES,
    )
    summary = {
        "version": VERSION,
        "bundle_schema_version": BUNDLE_SCHEMA_VERSION,
        "n_total": int(len(df)),
        "n_classical": int(len(classical_df)),
        "n_cryptic": int(len(cryptic_df)),
        "classical_distribution": {
            k: int(v) for k, v in classical_df["mlst_clade"].value_counts().to_dict().items()
        },
        "unique_st_distribution": {
            label: int(count) for label, count in df.groupby("mlst_clade")["st"].nunique().items()
        },
        "oof_pr_proxy": ap_score,
        "cross_fitted_oof_metrics_isolate_weighted": crossfit_metrics,
        "cross_fitted_oof_metrics_st_weighted": crossfit_metrics_by_st,
        "rules": rules,
        "crossfit_rules": crossfit_rules,
    }

    if args.metrics_prefix:
        with open(f"{args.metrics_prefix}_training_summary.json", "w") as fh:
            json.dump(summary, fh, indent=2)

    bundle = {
        "bundle_schema_version": BUNDLE_SCHEMA_VERSION,
        "stage1_model": stage1_final,
        "distance_rules": rules,
        "distance_stores": stores,
        "classical_profile_refs": collapse_unique_profiles(classical_df),
        "summary": summary,
    }
    vprint(f"Saving model bundle to {args.output_model}...", args.verbose)
    _model_dir = os.path.dirname(args.output_model)
    if _model_dir:
        os.makedirs(_model_dir, exist_ok=True)
    joblib.dump(bundle, args.output_model)

    # Model card — written alongside the .joblib file
    card_path = os.path.splitext(args.output_model)[0] + "_model_card.json"
    distance_hashes = {
        locus: sha256_file(Path(args.distance_dir) / f"{locus}.tsv") for locus in sorted(stores)
    }
    write_model_card(
        card_path,
        summary=summary,
        rules=rules,
        oof_metrics=crossfit_metrics,
        extra={
            "model_sha256": sha256_file(args.output_model),
            "bundle_schema_version": BUNDLE_SCHEMA_VERSION,
            "training_table": Path(args.training_table).name,
            "training_table_sha256": sha256_file(args.training_table),
            "distance_dir": Path(args.distance_dir).name,
            "distance_matrix_sha256": distance_hashes,
            "min_cryptic_recall": args.min_cryptic_recall,
            "threshold_grid_size": args.threshold_grid_size,
            "n_cv_folds": 5,
            "n_parallel_workers": n_workers,
            "cross_fitted_oof_metrics_st_weighted": crossfit_metrics_by_st,
            "software": {
                "python": platform.python_version(),
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "scikit_learn": sklearn.__version__,
                "joblib": joblib.__version__,
            },
            "total_train_sec": round(time.time() - t0, 1),
        },
    )

    elapsed = time.time() - t0
    vprint(f"\nDone. Total time: {elapsed:.1f}s", args.verbose)
    vprint(f"  Model   : {args.output_model}", args.verbose)
    vprint(f"  Card    : {card_path}", args.verbose)


if __name__ == "__main__":
    main()
