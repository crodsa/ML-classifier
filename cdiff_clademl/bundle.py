"""Model-bundle discovery, integrity checks, and structural validation."""

from __future__ import annotations

import hashlib
import json
import math
import sysconfig
from pathlib import Path
from typing import Any

from . import PROGRAM
from .config import BASE_COLS, CLASSICAL_CLADES, LOCI

BUNDLE_SCHEMA_VERSION = 1
MODEL_FILENAME = "model_v5.joblib"
MODEL_CARD_FILENAME = "model_v5_model_card.json"

_REQUIRED_BUNDLE_KEYS = {
    "stage1_model",
    "distance_rules",
    "distance_stores",
    "classical_profile_refs",
    "summary",
}
_REQUIRED_RULE_KEYS = {
    "sum_min_dist_any_threshold",
    "nearest_profile_any_threshold",
    "stage1_confidence_max",
    "stage1_margin_max",
}


class ModelBundleError(ValueError):
    """Raised when a serialized model bundle does not satisfy the contract."""


def sha256_file(path: str | Path) -> str:
    """Return a streaming SHA-256 digest for *path*."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sidecar_path(model_path: str | Path) -> Path:
    """Return the conventional model-card path next to a joblib file."""
    path = Path(model_path)
    return path.with_name(f"{path.stem}_model_card.json")


def load_model_card(model_path: str | Path) -> tuple[dict[str, Any] | None, Path]:
    """Load the adjacent model card when present."""
    path = sidecar_path(model_path)
    if not path.is_file():
        return None, path
    try:
        with path.open(encoding="utf-8") as handle:
            card = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise ModelBundleError(f"Cannot read model card {path}: {exc}") from exc
    if not isinstance(card, dict):
        raise ModelBundleError(f"Model card {path} must contain a JSON object.")
    return card, path


def verify_model_card_hash(
    model_path: str | Path,
    card: dict[str, Any] | None,
) -> str:
    """Verify the model digest recorded in a card and return the actual digest."""
    actual = sha256_file(model_path)
    expected = card.get("model_sha256") if card else None
    if expected and str(expected).lower() != actual:
        raise ModelBundleError(
            "Model integrity check failed: the joblib SHA-256 does not match its "
            "model card. Keep each model together with the card produced for it."
        )
    return actual


def find_bundled_model() -> Path:
    """Locate the model installed by the wheel, with a source-tree fallback."""
    installed = Path(sysconfig.get_path("data")) / "share" / PROGRAM / "models" / MODEL_FILENAME
    if installed.is_file():
        return installed

    source_tree = Path(__file__).resolve().parent.parent / MODEL_FILENAME
    if source_tree.is_file():
        return source_tree

    raise ModelBundleError(
        "No bundled model was found. Install the official wheel or pass "
        "--model /path/to/model.joblib."
    )


def _finite_rule(rules: dict[str, Any], key: str) -> None:
    try:
        value = float(rules[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise ModelBundleError(f"Invalid or missing model rule: {key}") from exc
    if not math.isfinite(value):
        raise ModelBundleError(f"Model rule {key} must be finite, got {value!r}.")


def validate_model_bundle(bundle: Any) -> dict[str, Any]:
    """Validate the runtime contract of a loaded joblib model bundle."""
    if not isinstance(bundle, dict):
        raise ModelBundleError("The model bundle must be a dictionary.")
    schema = bundle.get("bundle_schema_version", 0)
    if not isinstance(schema, int) or schema < 0 or schema > BUNDLE_SCHEMA_VERSION:
        raise ModelBundleError(
            f"Unsupported model bundle schema {schema!r}; runtime supports up to "
            f"{BUNDLE_SCHEMA_VERSION}."
        )

    missing = sorted(_REQUIRED_BUNDLE_KEYS - set(bundle))
    if missing:
        raise ModelBundleError(f"Model bundle is missing keys: {', '.join(missing)}")

    stage1 = bundle["stage1_model"]
    if not all(hasattr(stage1, method) for method in ("predict", "predict_proba")):
        raise ModelBundleError("stage1_model does not implement predict/predict_proba.")
    try:
        classes = [str(value) for value in stage1.named_steps["clf"].classes_]
    except (AttributeError, KeyError) as exc:
        raise ModelBundleError("stage1_model is not a fitted CdiffCladeML pipeline.") from exc
    if set(classes) != set(CLASSICAL_CLADES):
        raise ModelBundleError(
            f"stage1_model classes are incompatible: expected {CLASSICAL_CLADES}, got {classes}."
        )

    rules = bundle["distance_rules"]
    if not isinstance(rules, dict):
        raise ModelBundleError("distance_rules must be a dictionary.")
    missing_rules = sorted(_REQUIRED_RULE_KEYS - set(rules))
    if missing_rules:
        raise ModelBundleError(f"Model rules are missing: {', '.join(missing_rules)}")
    for key in _REQUIRED_RULE_KEYS:
        _finite_rule(rules, key)
    for key in ("sum_min_dist_any_threshold", "nearest_profile_any_threshold"):
        if float(rules[key]) < 0:
            raise ModelBundleError(f"Model distance rule {key} cannot be negative.")
    for key in ("stage1_confidence_max", "stage1_margin_max"):
        if not 0 <= float(rules[key]) <= 1:
            raise ModelBundleError(f"Model support rule {key} must be between 0 and 1.")
    max_missing = rules.get("max_missing", 1)
    if (
        isinstance(max_missing, bool)
        or not isinstance(max_missing, (int, float))
        or not float(max_missing).is_integer()
        or not 0 <= int(max_missing) <= len(LOCI)
    ):
        raise ModelBundleError(
            f"Model rule max_missing must be an integer between 0 and {len(LOCI)}."
        )

    stores = bundle["distance_stores"]
    if not isinstance(stores, dict) or set(stores) != set(LOCI):
        raise ModelBundleError(f"distance_stores must contain exactly these loci: {LOCI}.")
    for locus, store in stores.items():
        required = ("matrix", "allele_set", "classical_alleles", "stats")
        if not all(hasattr(store, attr) for attr in required):
            raise ModelBundleError(f"Distance store for {locus} is malformed.")
        for stat in ("mean_positive", "q95_positive"):
            try:
                value = float(store.stats[stat])
            except (KeyError, TypeError, ValueError) as exc:
                raise ModelBundleError(f"Distance store {locus} lacks valid {stat}.") from exc
            if not math.isfinite(value) or value < 0:
                raise ModelBundleError(f"Distance store {locus} has invalid {stat}={value!r}.")

    refs = bundle["classical_profile_refs"]
    if not hasattr(refs, "columns") or refs.empty:
        raise ModelBundleError("classical_profile_refs must be a non-empty table.")
    missing_columns = [column for column in BASE_COLS if column not in refs.columns]
    if missing_columns:
        raise ModelBundleError(
            f"classical_profile_refs lacks columns: {', '.join(missing_columns)}"
        )

    if not isinstance(bundle["summary"], dict):
        raise ModelBundleError("summary must be a dictionary.")
    return bundle
