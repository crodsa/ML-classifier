# Changelog

## 1.0 — 2026-10-22

Initial public release of CdiffCladeML (`v1.0`).

- Two-stage *Clostridioides difficile* seven-locus MLST clade classifier combining a Stage-1 one-hot-encoded Random Forest (`Clade1`–`Clade5`, 500 trees, `balanced_subsample`) with a Stage-2 allele- and profile-distance Cryptic gate (`Cryptic`).
- One-command prediction workflow (`CdiffCladeML-predict` / `cdiff-predict`) generating timestamped TSV predictions, JSON run summaries, and standalone interactive HTML reports.
- Support for PubMLST and FastMLST allele and clade notations (`adk(1)`, `adk=1`, `adk:1`, approximate novel alleles `~n`, and missing loci).
- Built-in input quality flags (`PASS`, `PARTIAL_PROFILE`, `UNKNOWN_ALLELE`, `INSUFFICIENT_LOCI`), Stage-1 support bands (`HIGH`, `MEDIUM`, `LOW`), and automated expert-review flags (`review_recommended`).
- Reproducible ST-grouped training workflow (`CdiffCladeML-train` / `cdiff-train`) with cross-fitted Cryptic gate evaluation, one-vote-per-ST and isolate-weighted metrics, and SHA-256-verified model cards (`model.joblib` and `model_model_card.json`).
