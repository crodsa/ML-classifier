"""
Prediction entry point for the v5 Hybrid C. difficile clade classifier.

Priority 3 — Structured output
--------------------------------
- Prediction output includes a 'confidence_band' column (HIGH / MEDIUM / LOW)
  derived from stage1_confidence, making it easier to triage samples for
  manual review without reading raw probability values.
- A per-run prediction summary JSON is written alongside the TSV output.
"""

import argparse
import json
import logging
import os
import sys
import types
import warnings
from datetime import datetime, timezone
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import sklearn
from sklearn.exceptions import InconsistentVersionWarning

from cdiff_clademl import PROGRAM, VERSION
from cdiff_clademl.bundle import (
    ModelBundleError,
    find_bundled_model,
    load_model_card,
    validate_model_bundle,
    verify_model_card_hash,
)
from cdiff_clademl.config import ALL_CLADES, CONFIDENCE_BAND_THRESHOLDS
from cdiff_clademl.features import build_meta_features
from cdiff_clademl.io import (
    read_mlst_table,
    write_eval_files,
    write_grouped_eval_files,
)
from cdiff_clademl.model import apply_v5_gate
from cdiff_clademl.utils import setup_logging, vprint


# ---------------------------------------------------------------------------
# Pickle compatibility shim: models serialised under the old package name
# 'cdiff_polifasico' can still be loaded after the rename to 'cdiff_clademl'.
# ---------------------------------------------------------------------------
def _register_compat_modules() -> None:
    import cdiff_clademl  # noqa: PLC0415

    _submodules = ["model", "features", "distance", "config", "utils", "io"]
    for _sub in _submodules:
        try:
            import importlib

            importlib.import_module(f"cdiff_clademl.{_sub}")
        except ImportError:
            pass
    _pkg = types.ModuleType("cdiff_polifasico")
    for _attr in dir(cdiff_clademl):
        setattr(_pkg, _attr, getattr(cdiff_clademl, _attr))
    sys.modules.setdefault("cdiff_polifasico", _pkg)
    for _sub in _submodules:
        _real = sys.modules.get(f"cdiff_clademl.{_sub}")
        if _real:
            sys.modules.setdefault(f"cdiff_polifasico.{_sub}", _real)


# ==============================================================================
# Confidence band helper (Priority 3)
# ==============================================================================


def assign_confidence_band(confidence_series: pd.Series) -> pd.Series:
    """Assign a human-readable confidence tier to each prediction.

    Tiers are derived from stage1_confidence (max Stage 1 class probability).
    They describe support for the classical Stage 1 call, not calibrated
    confidence in the final classical/Cryptic prediction.

    Thresholds are defined in config.CONFIDENCE_BAND_THRESHOLDS.

    Args:
        confidence_series: Series of float confidence values (0–1).

    Returns:
        Series of band strings ('HIGH', 'MEDIUM', 'LOW').
    """
    high_thresh = CONFIDENCE_BAND_THRESHOLDS["HIGH"]
    medium_thresh = CONFIDENCE_BAND_THRESHOLDS["MEDIUM"]

    def _band(c: float) -> str:
        if c >= high_thresh:
            return "HIGH"
        if c >= medium_thresh:
            return "MEDIUM"
        return "LOW"

    return confidence_series.map(_band)


# ==============================================================================
# Prediction summary (Priority 3)
# ==============================================================================


def write_prediction_summary(
    path: str,
    out_df: pd.DataFrame,
    model_version: str | None = None,
    model_sha256: str | None = None,
) -> None:
    """Write a JSON summary of prediction results for downstream automation.

    Includes class counts, confidence band distribution, and cryptic call rate.

    Args:
        path:          Output JSON file path.
        out_df:        Prediction output DataFrame (must contain final_prediction
                       and confidence_band columns).
        model_version: Optional version string from the loaded model bundle.
    """
    n = len(out_df)
    class_counts = out_df["final_prediction"].value_counts().to_dict()
    band_counts = out_df["stage1_support_band"].value_counts().to_dict()
    quality_counts = out_df["quality_flag"].value_counts().to_dict()
    n_cryptic = int(out_df["cryptic_distance_call"].sum())

    summary = {
        "tool": PROGRAM,
        "version": VERSION,
        "model_version": model_version,
        "model_sha256": model_sha256,
        "runtime": {
            "python": sys.version.split()[0],
            "scikit_learn": sklearn.__version__,
        },
        "predicted_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n_samples": int(n),
        "class_counts": class_counts,
        "stage1_support_bands": band_counts,
        "quality_flags": quality_counts,
        "n_review_recommended": int(out_df["review_recommended"].sum()),
        "n_cryptic_gate": n_cryptic,
        "cryptic_rate_pct": round(100.0 * n_cryptic / max(1, n), 2),
    }
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)


# ==============================================================================
# Main
# ==============================================================================


def main() -> None:
    ap = argparse.ArgumentParser(
        prog="CdiffCladeML-predict",
        description=f"{PROGRAM} v{VERSION} — Run v5 hybrid classifier on new samples",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    ap.add_argument("input_table", help="10- or 11-column MLST table (mlst_clade column optional).")
    ap.add_argument(
        "--model",
        default=None,
        help=(
            "Path to a trusted .joblib bundle produced by CdiffCladeML-train. "
            "Uses the model bundled with the official wheel when omitted."
        ),
    )
    ap.add_argument(
        "-o", "--output_predictions", required=True, help="Output TSV path for predictions."
    )
    ap.add_argument(
        "--metrics_prefix",
        default=None,
        help="If set and mlst_clade column is present, write eval files with this prefix.",
    )
    ap.add_argument("--verbose", action="store_true", help="Verbose logging (INFO level).")
    ap.add_argument(
        "--allow-version-mismatch",
        action="store_true",
        help="Allow loading a model serialized by another scikit-learn version (unsafe).",
    )
    args = ap.parse_args()
    setup_logging(args.verbose)

    vprint(f"{PROGRAM} v{VERSION} — prediction started", args.verbose)

    # ── Load model and data ───────────────────────────────────────────────────
    try:
        model_path = Path(args.model).expanduser() if args.model else find_bundled_model()
        if not model_path.is_file():
            raise ModelBundleError(f"Model file does not exist: {model_path}")
        card, card_path = load_model_card(model_path)
        model_sha256 = verify_model_card_hash(model_path, card)
        if card is None:
            logging.warning(
                "No model card found at %s; provenance and hash cannot be independently checked.",
                card_path,
            )
        recorded_sklearn = ((card or {}).get("serialization") or {}).get("scikit_learn") or (
            (card or {}).get("software") or {}
        ).get("scikit_learn")
        if (
            recorded_sklearn
            and recorded_sklearn != sklearn.__version__
            and not args.allow_version_mismatch
        ):
            raise ModelBundleError(
                f"Model requires scikit-learn {recorded_sklearn}, but runtime is "
                f"{sklearn.__version__}. Install the dependencies pinned by this release."
            )

        # joblib/pickle can execute code. Only models from trusted sources are safe.
        _register_compat_modules()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", InconsistentVersionWarning)
            bundle = joblib.load(model_path)
        mismatches = [
            item for item in caught if issubclass(item.category, InconsistentVersionWarning)
        ]
        if mismatches and not args.allow_version_mismatch:
            raise ModelBundleError(
                "The model was serialized with a different scikit-learn version. "
                "Install the dependencies pinned by this release; use "
                "--allow-version-mismatch only for expert diagnostics."
            )
        if mismatches:
            model_versions = sorted(
                {
                    str(getattr(item.message, "original_sklearn_version", "unknown"))
                    for item in mismatches
                }
            )
            logging.warning(
                "Expert override: loading sklearn model version(s) %s with runtime %s; "
                "%d compatibility warnings suppressed.",
                ", ".join(model_versions),
                sklearn.__version__,
                len(mismatches),
            )
        bundle = validate_model_bundle(bundle)
    except (EOFError, ImportError, OSError, TypeError, ModelBundleError, ValueError) as exc:
        ap.error(str(exc))

    summary = bundle.get("summary", {})
    model_version = (
        (card or {}).get("model_version")
        or (card or {}).get("version")
        or summary.get("model_version")
        or summary.get("version")
        or "legacy-unknown"
    )
    vprint(f"Using model {model_path} (SHA-256 {model_sha256[:12]}…)", args.verbose)

    try:
        df = read_mlst_table(args.input_table, require_labels=False)
    except (OSError, ValueError) as exc:
        ap.error(str(exc))
    vprint(f"Loaded {len(df)} rows for prediction.", args.verbose)

    stage1 = bundle["stage1_model"]
    probs = stage1.predict_proba(df)
    pred = stage1.predict(df)
    classes = list(stage1.named_steps["clf"].classes_)

    meta = build_meta_features(
        df,
        pred,
        probs,
        bundle["classical_profile_refs"],
        bundle["distance_stores"],
        classes,
        verbose=args.verbose,
        tag="[predict]",
    )

    cryptic_call = apply_v5_gate(meta, bundle["distance_rules"])
    final_pred = np.where(cryptic_call == 1, "Cryptic", meta["stage1_predicted_clade"].to_numpy())
    reason = np.where(
        cryptic_call == 1,
        "distance+low_classical_support",
        "classical_stage1",
    )

    # ── Assemble output DataFrame ─────────────────────────────────────────────
    out = df.copy()
    for col in meta.columns:
        out[col] = meta[col].values
    out["cryptic_distance_call"] = cryptic_call
    out["final_prediction"] = final_pred
    out["final_reason"] = reason

    # This is support for the classical Stage 1 model, not final-call confidence.
    out["stage1_support_band"] = assign_confidence_band(out["stage1_confidence"])
    out["confidence_band"] = out["stage1_support_band"]  # v5.3 compatibility alias
    out["quality_flag"] = np.select(
        [
            out["n_missing_loci"].gt(1),
            out["n_unknown_alleles"].gt(0),
            out["n_missing_loci"].gt(0) | out["n_approx_loci"].gt(0),
        ],
        ["INSUFFICIENT_LOCI", "UNKNOWN_ALLELE", "PARTIAL_PROFILE"],
        default="PASS",
    )
    out["review_recommended"] = out["stage1_support_band"].eq("LOW") | out["quality_flag"].ne(
        "PASS"
    )

    _out_dir = os.path.dirname(args.output_predictions)
    if _out_dir:
        os.makedirs(_out_dir, exist_ok=True)
    out.to_csv(args.output_predictions, sep="\t", index=False)
    vprint(f"Predictions written to {args.output_predictions}", args.verbose)

    # Prediction summary JSON (Priority 3)
    summary_path = os.path.splitext(args.output_predictions)[0] + "_summary.json"
    write_prediction_summary(
        summary_path,
        out,
        model_version=model_version,
        model_sha256=model_sha256,
    )
    vprint(f"Prediction summary written to {summary_path}", args.verbose)

    # ── Optional evaluation (requires ground truth) ────────────────────────────
    _NULL_LIKE = {"nan", "none", "na", "", "-"}
    _labelled = (
        df.get("mlst_clade", pd.Series(dtype=str)).apply(
            lambda v: str(v).strip().lower() not in _NULL_LIKE
        )
        if "mlst_clade" in df.columns
        else pd.Series(False, index=df.index)
    )
    if _labelled.any() and args.metrics_prefix:
        _mprefix_dir = os.path.dirname(args.metrics_prefix)
        if _mprefix_dir:
            os.makedirs(_mprefix_dir, exist_ok=True)
        write_eval_files(
            args.metrics_prefix,
            df.loc[_labelled, "mlst_clade"].to_numpy(),
            final_pred[_labelled],
            ALL_CLADES,
        )
        labelled_st = df.loc[_labelled, "st"].astype(str).str.strip()
        missing_st = labelled_st.str.casefold().isin({"", "-", "na", "n/a", "nan", "none"})
        if missing_st.any():
            logging.warning(
                "Skipping ST-weighted metrics because %d labelled row(s) lack ST.",
                int(missing_st.sum()),
            )
        else:
            write_grouped_eval_files(
                f"{args.metrics_prefix}_by_st",
                df.loc[_labelled, "mlst_clade"].to_numpy(),
                final_pred[_labelled],
                labelled_st.to_numpy(),
                ALL_CLADES,
            )
        vprint("Evaluation files written.", args.verbose)


if __name__ == "__main__":
    main()
