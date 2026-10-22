# Contributing

Thank you for improving CdiffCladeML. This project is research and
molecular-surveillance software; changes must preserve the boundary between
screening output and clinical interpretation.

## Before opening a pull request

1. Create a focused branch and explain the biological or engineering rationale.
2. Add or update tests for every behaviour change.
3. Run `python -m ruff check .` and `python -m pytest` in a clean environment.
4. Do not commit identifiable clinical data, unreviewed assemblies, credentials,
   or untrusted joblib files.

## Changes to the model or validation evidence

Model updates require a new model card, SHA-256 digest, exact software versions,
and documented data provenance. They must use leakage-free ST-grouped
validation, report one-vote-per-ST results, and be evaluated against a newly
locked cohort. Do not describe a development benchmark as independent,
prospective, clinical, or diagnostic validation.

## Reporting issues

Include the software version, model hash, command used, a minimal de-identified
input example, and the observed versus expected behaviour. Never post patient
identifiers or raw clinical metadata in a public issue.
