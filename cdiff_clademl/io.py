"""
MLST table I/O, allele field parsing, evaluation file writing,
and model card output (Priority 3 structured output).
"""

import csv
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    classification_report,
    confusion_matrix,
    f1_score,
)

from . import PROGRAM, VERSION
from .config import ALL_CLADES, BASE_COLS, CLASSICAL_CLADES, LOCI

_MISSING_VALUES = {"", "-", "na", "n/a", "nan", "none", "null", "<na>", "?"}


# ==============================================================================
# Allele and clade field parsing
# ==============================================================================


def parse_allele_inner(s: str) -> tuple[str, str]:
    """Parse an inner allele string into (base, status).

    Args:
        s: Raw allele string after stripping parenthetical wrappers.

    Returns:
        (base, status) where status is 'missing', 'approx', or 'exact'.
    """
    s = str(s).strip()
    if s.casefold() in _MISSING_VALUES:
        return "__MISSING__", "missing"
    if s.startswith("~"):
        base = s[1:].strip()
        if not base:
            return "__MISSING__", "missing"
        return base, "approx"
    return s, "exact"


def extract_parenthetical(raw: str, prefix: str | None = None) -> str:
    """Extract the first parenthetical value, or strip a prefix.

    Args:
        raw:    Raw field value from the MLST table.
        prefix: Optional string prefix to strip if no parenthetical is found.

    Returns:
        Extracted or stripped string.
    """
    s = str(raw).strip()
    m = re.search(r"\((.*?)\)", s)
    if m:
        return m.group(1)
    if prefix and s.startswith(prefix):
        return s[len(prefix) :]
    return s


def parse_allele_field(raw: str) -> tuple[str, str]:
    """Parse a raw MLST allele field into (base, status).

    Args:
        raw: Raw cell value, possibly wrapped in parentheses or prefixed with '~'.

    Returns:
        (base, status) tuple.
    """
    return parse_allele_inner(extract_parenthetical(raw))


def parse_clade_field(raw: str) -> str:
    """Parse a raw clade field, stripping parenthetical wrappers.

    Args:
        raw: Raw cell value.

    Returns:
        Cleaned clade string.
    """
    value = str(raw).strip()
    if value.casefold() in _MISSING_VALUES:
        return ""
    inner = extract_parenthetical(value)
    return inner.strip() if inner else value


# ==============================================================================
# MLST table reading
# ==============================================================================


def _read_mlst_rows(path: str) -> tuple[list[list[str]], int]:
    """Read a comma- or tab-delimited MLST file without silently dropping rows."""
    source = Path(path)
    if not source.is_file():
        raise ValueError(f"Input table does not exist or is not a file: {path}")

    with source.open(encoding="utf-8-sig", newline="") as handle:
        text = handle.read()
    if not text.strip():
        raise ValueError(f"Input table is empty: {path}")

    sample = text[:32768]
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        delimiter = dialect.delimiter
    except csv.Error:
        first = next((line for line in text.splitlines() if line.strip()), "")
        delimiter = "\t" if first.count("\t") > first.count(",") else ","

    rows = [
        [cell.strip() for cell in row]
        for row in csv.reader(text.splitlines(), delimiter=delimiter)
        if row and any(cell.strip() for cell in row)
    ]
    if not rows:
        raise ValueError(f"Input table contains no data rows: {path}")

    expected_header = ["isolate", "species", "st", *(locus.casefold() for locus in LOCI)]
    first_lower = [cell.casefold() for cell in rows[0]]
    if first_lower in (expected_header, [*expected_header, "mlst_clade"]):
        rows = rows[1:]
        if not rows:
            raise ValueError(f"Input table contains a header but no data rows: {path}")

    invalid = [
        (line_no, len(row)) for line_no, row in enumerate(rows, start=1) if len(row) not in {10, 11}
    ]
    if invalid:
        details = ", ".join(f"row {line}: {width} columns" for line, width in invalid[:5])
        raise ValueError(f"Every row in {path} must have 10 or 11 columns; found {details}.")
    return rows, max(map(len, rows))


def read_mlst_table(path: str, require_labels: bool = False) -> pd.DataFrame:
    """Read a PubMLST-style MLST table and expand allele fields into feature columns.

    The table may have 10 columns (no clade label) or 11 (with clade label).
    Column order: isolate, species, st, adk, atpA, dxr, glyA, recA, sodA, tpi[, mlst_clade].

    For each locus, three feature columns are derived:
        {locus}_base:   allele number string (or '__MISSING__').
        {locus}_status: 'exact', 'approx', or 'missing'.
        {locus}_token:  'E:{base}', 'A:{base}', or '__MISSING__' for OHE.

    Args:
        path:           Path to the CSV/TSV table (auto-detected separator).
        require_labels: If True, raises ValueError when mlst_clade column is absent.

    Returns:
        Expanded DataFrame with all feature columns added.

    Raises:
        ValueError: On wrong column count or missing required label column.
    """
    cols10 = ["isolate", "species", "st", *LOCI]
    cols11 = cols10 + ["mlst_clade"]

    rows, n_cols = _read_mlst_rows(path)
    if n_cols == 11:
        rows = [row if len(row) == 11 else [*row, ""] for row in rows]
        raw = pd.DataFrame(rows, columns=cols11, dtype=str)
    else:
        raw = pd.DataFrame(rows, columns=cols10, dtype=str)

    if require_labels and n_cols != 11:
        raise ValueError("Training requires an 11-column table with mlst_clade labels.")

    df = raw.copy()
    for locus in LOCI:
        bases, statuses = zip(*df[locus].map(parse_allele_field), strict=True)
        df[f"{locus}_base"] = list(bases)
        df[f"{locus}_status"] = list(statuses)
        df[f"{locus}_token"] = [
            "__MISSING__" if st == "missing" else ("A:" + b if st == "approx" else "E:" + b)
            for b, st in zip(bases, statuses, strict=True)
        ]

    if "mlst_clade" in df.columns:
        df["mlst_clade"] = df["mlst_clade"].map(parse_clade_field)
        labelled = df["mlst_clade"].ne("")
        invalid_labels = sorted(set(df.loc[labelled, "mlst_clade"]) - set(ALL_CLADES))
        if invalid_labels:
            raise ValueError(
                "Unknown mlst_clade label(s): "
                f"{', '.join(invalid_labels)}. Expected one of {ALL_CLADES}."
            )
        if require_labels and not labelled.all():
            rows_missing = (df.index[~labelled] + 1).tolist()[:10]
            raise ValueError(f"Training rows without mlst_clade labels: {rows_missing}")

    return df


def validate_training_table(df: pd.DataFrame, n_splits: int = 5) -> None:
    """Validate biological and cross-validation invariants before training."""
    if "mlst_clade" not in df or df["mlst_clade"].eq("").any():
        raise ValueError("Every training row must have a valid mlst_clade label.")

    missing_classes = [label for label in ALL_CLADES if label not in set(df["mlst_clade"])]
    if missing_classes:
        raise ValueError(f"Training data is missing required classes: {missing_classes}")

    duplicated_isolates = df["isolate"].duplicated(keep=False)
    if duplicated_isolates.any():
        examples = sorted(df.loc[duplicated_isolates, "isolate"].unique())[:5]
        raise ValueError(f"Duplicate isolate identifiers in training data: {examples}")

    missing_st = df["st"].map(lambda value: str(value).strip().casefold() in _MISSING_VALUES)
    if missing_st.any():
        raise ValueError("Training requires a non-missing ST for grouped validation.")

    st_label_counts = df.groupby("st", dropna=False)["mlst_clade"].nunique()
    if (st_label_counts > 1).any():
        examples = st_label_counts[st_label_counts > 1].index.astype(str).tolist()[:5]
        raise ValueError(f"A single ST maps to multiple clades: {examples}")

    profile_label_counts = df.groupby(BASE_COLS, dropna=False)["mlst_clade"].nunique()
    if (profile_label_counts > 1).any():
        raise ValueError("One or more seven-locus profiles map to multiple clades.")

    required = [*CLASSICAL_CLADES, "Cryptic"]
    group_counts = df.groupby("mlst_clade")["st"].nunique()
    too_few = {
        label: int(group_counts.get(label, 0))
        for label in required
        if group_counts.get(label, 0) < n_splits
    }
    if too_few:
        raise ValueError(
            f"At least {n_splits} distinct ST groups per class are required; got {too_few}."
        )


# ==============================================================================
# Evaluation file writing
# ==============================================================================


def summarize_multiclass(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    labels: list[str] | None = None,
) -> dict[str, float]:
    """Compute aggregate multiclass classification metrics.

    Args:
        y_true: Ground-truth label array.
        y_pred: Predicted label array.

    Returns:
        Dict with accuracy, balanced_accuracy, macro_f1, weighted_f1.
    """
    return {
        "n_samples": int(len(y_true)),
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "macro_f1": float(
            f1_score(y_true, y_pred, labels=labels, average="macro", zero_division=0)
        ),
        "weighted_f1": float(
            f1_score(y_true, y_pred, labels=labels, average="weighted", zero_division=0)
        ),
    }


def write_eval_files(
    prefix: str,
    y_true: np.ndarray,
    y_pred: np.ndarray,
    labels: list[str],
) -> None:
    """Write the standard set of evaluation output files.

    Writes four files:
        {prefix}_summary.tsv:                 aggregate metrics row.
        {prefix}_per_class.tsv:               per-class precision/recall/F1.
        {prefix}_confusion.tsv:               confusion matrix.
        {prefix}_classification_report.json:  sklearn full report dict.

    Args:
        prefix: File path prefix (may include directory).
        y_true: Ground-truth labels.
        y_pred: Predicted labels.
        labels: Ordered label list for consistent matrix/report ordering.
    """
    summary = summarize_multiclass(y_true, y_pred, labels=labels)
    pd.DataFrame([summary]).to_csv(f"{prefix}_summary.tsv", sep="\t", index=False)

    per = classification_report(
        y_true,
        y_pred,
        labels=labels,
        output_dict=True,
        zero_division=0,
    )
    pd.DataFrame(per).T.to_csv(f"{prefix}_per_class.tsv", sep="\t")

    cm = confusion_matrix(y_true, y_pred, labels=labels)
    pd.DataFrame(cm, index=labels, columns=labels).to_csv(f"{prefix}_confusion.tsv", sep="\t")

    with open(f"{prefix}_classification_report.json", "w", encoding="utf-8") as fh:
        json.dump(per, fh, indent=2)


def write_grouped_eval_files(
    prefix: str,
    y_true: np.ndarray,
    y_pred: np.ndarray,
    groups: np.ndarray,
    labels: list[str],
) -> None:
    """Write one-vote-per-group metrics after checking group consistency."""
    grouped = pd.DataFrame({"group": groups, "true": y_true, "pred": y_pred})
    consistency = grouped.groupby("group", dropna=False).agg(
        n_true=("true", "nunique"), n_pred=("pred", "nunique")
    )
    if (consistency > 1).any(axis=None):
        bad = consistency[(consistency > 1).any(axis=1)].index.astype(str).tolist()[:5]
        raise ValueError(f"Inconsistent labels or predictions within ST group(s): {bad}")
    one_per_group = grouped.drop_duplicates("group")
    write_eval_files(
        prefix,
        one_per_group["true"].to_numpy(),
        one_per_group["pred"].to_numpy(),
        labels,
    )


# ==============================================================================
# Model card (Priority 3 — structured output)
# ==============================================================================


def write_model_card(
    path: str,
    summary: dict,
    rules: dict,
    oof_metrics: dict[str, float] | None = None,
    extra: dict | None = None,
) -> None:
    """Write a reproducibility-focused model card JSON.

    Contains version, training date, data counts, learned threshold rules,
    and OOF performance metrics. Designed to be committed alongside the
    .joblib model file so that any future user can understand what data and
    hyperparameters produced it.

    Args:
        path:        Output file path (e.g. 'model_card.json').
        summary:     Training summary dict (from the bundle).
        rules:       Learned Cryptic gate threshold rules.
        oof_metrics: Optional dict of OOF multiclass metrics.
        extra:       Optional additional key-value pairs to include.
    """
    card = {
        "tool": PROGRAM,
        "version": VERSION,
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "architecture": {
            "stage1": (
                "RandomForest (500 trees, balanced_subsample, OOF 5-fold "
                "StratifiedGroupKFold by ST)"
            ),
            "stage2": (
                "Rule-based Cryptic gate (distance + stage1 support thresholds; "
                "cross-fitted evaluation by ST)"
            ),
        },
        "training_data": {
            "n_total": summary.get("n_total"),
            "n_classical": summary.get("n_classical"),
            "n_cryptic": summary.get("n_cryptic"),
            "classical_distribution": summary.get("classical_distribution"),
            "unique_st_distribution": summary.get("unique_st_distribution"),
        },
        "learned_rules": rules,
        "oof_pr_proxy_tuning_diagnostic": summary.get("oof_pr_proxy"),
        "oof_validation_status": "cross_fitted_development_diagnostic_not_nested_cv",
        "cross_fitted_oof_metrics_isolate_weighted": oof_metrics or {},
    }
    if extra:
        card.update(extra)
    software = card.get("software", {})
    card.setdefault(
        "serialization",
        {
            "format": "joblib/pickle",
            "scikit_learn": software.get("scikit_learn"),
            "python": software.get("python"),
            "warning": "Load only model files from trusted sources.",
        },
    )

    with open(path, "w", encoding="utf-8") as fh:
        json.dump(card, fh, indent=2)
