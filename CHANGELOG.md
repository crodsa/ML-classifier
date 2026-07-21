# Changelog

## 5.4.0 — 2026-07-20

- Audited the independent test split by isolate, ST and seven-locus profile.
- Added ST-weighted evaluation and corrected interpretation of historical OOF metrics.
- Fixed Cryptic pseudoreplication and fold-level distance leakage in future training.
- Added cross-fitted evaluation of the Cryptic gate.
- Replaced silent input-row dropping with strict CSV/TSV and label validation.
- Added distance-matrix, model-structure, SHA-256 and sklearn-version checks.
- Added missing/approximate/unknown-allele quality flags and review guidance.
- Bundled the validated model and its model card in the wheel while retaining
  `--model` for expert model replacement.
- Pinned the exact scikit-learn serialization version.
- Added documentation, security guidance, license, tests and release metadata.

The bundled binary model remains version 5.3.0. Its predictions and independent
test metrics are unchanged; 5.4.0 is the audited software/package release.
