# ruff: noqa: E501
"""Standalone, human-readable HTML reports for prediction runs.

The report is intentionally a static HTML file: it can be archived alongside a
TSV/JSON result, opened without a server, and reviewed without transmitting
sample identifiers to a third party.
"""

from __future__ import annotations

import html
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from . import PROGRAM, VERSION
from .config import ALL_CLADES
from .io import evaluation_group_ids, summarize_multiclass


def _esc(value: Any) -> str:
    """Return an HTML-escaped, display-safe value."""
    if pd.isna(value):
        return ""
    return html.escape(str(value))


def _table(headers: list[str], rows: list[list[Any]], css_class: str = "") -> str:
    """Render a compact HTML table without external JavaScript or assets."""
    head = "".join(f"<th>{_esc(header)}</th>" for header in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{_esc(cell)}</td>" for cell in row) + "</tr>" for row in rows
    )
    return f'<table class="{css_class}"><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>'


def _metric_rows(metrics: dict[str, float]) -> list[list[str]]:
    return [
        ["Profiles evaluated", str(int(metrics["n_samples"]))],
        ["Accuracy", f'{metrics["accuracy"]:.3f}'],
        ["Balanced accuracy", f'{metrics["balanced_accuracy"]:.3f}'],
        ["Macro F1", f'{metrics["macro_f1"]:.3f}'],
        ["Weighted F1", f'{metrics["weighted_f1"]:.3f}'],
    ]


def _evaluation_sections(out_df: pd.DataFrame) -> str:
    """Create performance sections only for rows with an actual reference label."""
    if "mlst_clade" not in out_df.columns:
        return ""
    labels = out_df["mlst_clade"].astype(str).str.strip()
    labelled = labels.ne("") & labels.str.casefold().ne("nan")
    if not labelled.any():
        return (
            "<section><h2>Evaluation</h2><p class=\"note\">No reference clade labels were supplied. "
            "This run provides classifications, not a new estimate of performance.</p></section>"
        )

    evaluation = out_df.loc[
        labelled, ["isolate", "st", "mlst_clade", "final_prediction"]
    ].copy()
    metrics = summarize_multiclass(
        evaluation["mlst_clade"].to_numpy(),
        evaluation["final_prediction"].to_numpy(),
        labels=ALL_CLADES,
    )
    observed_labels = sorted(evaluation["mlst_clade"].unique())
    observed_metrics = summarize_multiclass(
        evaluation["mlst_clade"].to_numpy(),
        evaluation["final_prediction"].to_numpy(),
        labels=observed_labels,
    )
    blocks = [
        "<section><h2>Evaluation against supplied labels</h2>",
        _table(["Metric", "Value"], _metric_rows(metrics), "metrics"),
        (
            '<p class="note">Reference coverage: '
            f"{len(observed_labels)}/{len(ALL_CLADES)} model classes. Macro F1 across observed "
            f"reference classes: {observed_metrics['macro_f1']:.3f}. The fixed six-class macro F1 "
            "above assigns zero to unobserved classes.</p>"
        ),
    ]
    groups, n_unassigned = evaluation_group_ids(evaluation["st"], evaluation["isolate"])
    grouped = evaluation.assign(_group=groups).drop_duplicates("_group")
    grouped_metrics = summarize_multiclass(
        grouped["mlst_clade"].to_numpy(),
        grouped["final_prediction"].to_numpy(),
        labels=ALL_CLADES,
    )
    note = (
        f" {n_unassigned} row(s) with unassigned ST were counted as separate groups."
        if n_unassigned
        else ""
    )
    blocks.extend(
        [
            "<h3>One vote per sequence type</h3>",
            _table(["Metric", "Value"], _metric_rows(grouped_metrics), "metrics"),
            f'<p class="note">{note.strip()}</p>' if note else "",
        ]
    )
    blocks.append("</section>")
    return "".join(blocks)


def write_html_report(
    path: str | Path,
    out_df: pd.DataFrame,
    summary: dict[str, Any],
    model_version: str,
    model_sha256: str | None,
) -> None:
    """Write a static report for review, interpretation and archiving.

    The report does not transform model outputs. It makes run-level quality,
    review priorities, and model provenance visible while retaining the TSV as
    the complete machine-readable result.
    """
    n = len(out_df)
    review = out_df.loc[out_df["review_recommended"]].copy()
    priority = {"INSUFFICIENT_LOCI": 0, "UNKNOWN_ALLELE": 1, "PARTIAL_PROFILE": 2, "PASS": 3}
    review["_priority"] = review["quality_flag"].map(priority).fillna(4)
    review["_cryptic"] = review["final_prediction"].eq("Cryptic").astype(int)
    review = review.sort_values(
        ["_priority", "_cryptic", "stage1_confidence"], ascending=[True, False, True]
    ).head(100)

    class_rows = [[label, int(summary.get("class_counts", {}).get(label, 0))] for label in ALL_CLADES]
    quality_rows = [[label, int(count)] for label, count in summary.get("quality_flags", {}).items()]
    support_rows = [[label, int(count)] for label, count in summary.get("stage1_support_bands", {}).items()]
    queue_columns = [
        "isolate",
        "st",
        "final_prediction",
        "final_reason",
        "quality_flag",
        "stage1_support_band",
        "n_missing_loci",
        "n_approx_loci",
        "n_unknown_alleles",
        "rare_locus_count_report",
        "sum_min_dist_any_classical",
        "nearest_profile_dist_any_classical",
    ]
    queue = review.reindex(columns=queue_columns)
    queue_rows = []
    for _, row in queue.iterrows():
        rendered = []
        for column in queue_columns:
            value = row[column]
            if isinstance(value, (float, np.floating)) and column.startswith(("sum_", "nearest_")):
                rendered.append(f"{value:.4f}")
            else:
                rendered.append(value)
        queue_rows.append(rendered)

    hash_display = (model_sha256 or "not recorded")
    provenance_rows = [
        ["Software", f"{PROGRAM} {VERSION}"],
        ["Model", model_version],
        ["Model SHA-256", hash_display],
        ["Profiles processed", n],
        ["Review recommended", int(summary.get("n_review_recommended", 0))],
        ["Cryptic gate calls", int(summary.get("n_cryptic_gate", 0))],
    ]

    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{PROGRAM} run report</title>
<style>
body{{font-family:Arial,sans-serif;line-height:1.45;color:#172033;max-width:1200px;margin:auto;padding:28px;background:#f6f8fb}}
h1{{margin-bottom:4px}} h2{{margin-top:34px;color:#1d4c62}} h3{{margin-top:22px}}
.subtitle,.note{{color:#4d5e70}} .warning{{border-left:5px solid #c57a00;background:#fff5df;padding:14px 16px}}
.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:18px}} .card{{background:#fff;border:1px solid #dce3ec;border-radius:10px;padding:16px}}
table{{border-collapse:collapse;width:100%;background:#fff;margin:10px 0 16px}} th,td{{border:1px solid #dce3ec;padding:8px;text-align:left;vertical-align:top}} th{{background:#eaf2f5;color:#183d4e}}
.queue{{font-size:0.88rem;overflow-wrap:anywhere}} footer{{font-size:.9rem;color:#4d5e70;margin-top:32px}}
</style></head><body>
<h1>{PROGRAM} run report</h1>
<p class="subtitle">Research and molecular-surveillance use only. This file is a review aid; the accompanying TSV is the complete result.</p>
<div class="warning"><strong>Interpretation safeguard.</strong> <em>Cryptic</em> is a combined category, not an assignment to a specific cryptic lineage. Calls with missing, approximate, or unknown alleles require MLST/QC confirmation and, when biologically consequential, phylogenomic review. This software is not validated for clinical diagnosis or treatment decisions.</div>
<section><h2>Run overview</h2><div class="grid"><div class="card"><h3>Final calls</h3>{_table(["Class", "Profiles"], class_rows)}</div><div class="card"><h3>Quality flags</h3>{_table(["Flag", "Profiles"], quality_rows)}</div><div class="card"><h3>Classical-stage support</h3>{_table(["Band", "Profiles"], support_rows)}</div></div></section>
<section><h2>Priority review queue</h2><p class="note">Ordered by input-quality risk, then Cryptic calls and lower classical-stage support. Up to 100 rows are shown.</p>{_table(queue_columns, queue_rows, "queue") if queue_rows else '<p class="note">No profiles were flagged for expert review.</p>'}</section>
{_evaluation_sections(out_df)}
<section><h2>Provenance</h2>{_table(["Field", "Value"], provenance_rows, "metrics")}</section>
<footer>Stage-1 support is not a calibrated probability of the final call. Cryptic calls with fewer than four rare loci are flagged for review. The published development benchmark has reduced Clade2 recall and very limited Clade3 support; use the model card and validation documentation when interpreting a cohort.</footer>
</body></html>"""
    Path(path).write_text(document, encoding="utf-8")
