# Engineering audit and distribution readiness

Audit date: 2026-07-20. Resulting package version: 5.4.0. Retained binary model
version: 5.3.0.

## Achieved state

The predictor and official model form an installable, self-contained wheel. An
external model remains replaceable through `--model`. The distribution includes
the license, user guide, model card, security controls, and serialization-
compatible dependencies.

Verified artifact:

```text
dist/cdiffclademl-5.4.0-py3-none-any.whl
SHA-256: recorded in dist/SHA256SUMS
```

The wheel contains `model_v5.joblib`; its digest matches the model card:

```text
c2a8294cd473ff3f79e9ced92bc743e13706b6a0c1a5a122dbec3746b2efd115
```

## Confirmed and corrected defects

| Area | Previous state | 5.4 correction |
|---|---|---|
| Input | The reader assigned 11 column names before counting fields, so its 10/11-column check was ineffective | Strict CSV/TSV parser, optional header, per-row width checks, and explicit errors |
| Labels | Missing or unknown classes could result in empty datasets or omitted rows | Vocabulary, completeness, ST, profile, and required-class validation |
| Cryptic OOF | Each of 104 Cryptic samples occurred five times | Every Cryptic ST occurs in exactly one fold |
| CV distances | Classical allele sets were learned from the full dataset | Allele references are reconstructed inside each fold |
| Gate tuning/evaluation | The same data and thresholds were used | Separate cross-fitted diagnostic, explicitly identified as non-nested |
| Metrics | Only isolate-weighted results | Additional one-vote-per-ST output |
| Matrices | Columns, finiteness, diagonals, symmetry, and coverage were not verified | Strong validation and rejection of contradictory distances |
| Model | Joblib file without a contract or integrity check | Bundle schema, structural validation, model card, and SHA-256 |
| Compatibility | Unbounded dependencies and a scikit-learn 1.7.1→1.9.0 warning | Exact `scikit-learn==1.7.1` dependency and version-mismatch rejection by default |
| Security | Pickle loading without a trust warning | Trust guidance, pre-load digest verification, and `SECURITY.md` |
| Query quality | No explicit signal for new or approximate alleles | Counts, `quality_flag`, and `review_recommended` |
| Distribution | No README, license, tests, or model in the wheel | Package metadata, documentation, automated tests, and a self-contained wheel |

## Verification evidence

- Ruff static analysis: passed;
- Python module compilation: passed;
- 15 unit and integration tests, including report rendering and FastMLST
  parenthesized-missing-label handling: passed;
- complete synthetic training with 30 samples/30 STs, seven matrices, five
  folds, cross-fitting, save, reload, and prediction: passed;
- structural validation of 25,124 real records and seven real matrices: passed;
- regression over 4,650 real samples: no Stage 1 or final prediction changes
  relative to the original model;
- clean wheel installation with Python 3.12 and scikit-learn 1.7.1: passed;
- bundled-model discovery, digest verification, metadata check, and six-class
  smoke test: passed.

## Residual risks

1. Public distribution requires a tagged GitHub Release containing both wheel
   and source distribution, SHA-256 checksums, and successful CI for Python
   3.10–3.13. The release workflow added with this audit enforces the test,
   lint, build, and package-metadata checks, but maintainers must confirm the
   release assets before publishing a tag.
2. The production model was deliberately not retrained so the validated binary
   artifact could be preserved and replaced later. The corrected trainer passed
   a synthetic end-to-end test and validated the real data and matrices, but
   the next production model must complete the full process and use a new
   external cohort.
3. The initial audit environment used scikit-learn 1.9.0. Under the expert
   override it produced the same 4,650 predictions. The final wheel test also
   ran with the required version 1.7.1 and no override.
4. The no-overlap benchmark was consulted during several development
   iterations. It is a development benchmark, not prospective validation. A
   major release must replace it with a final locked cohort.
5. Temporal and geographic provenance, label criteria, and a formally versioned
   PubMLST snapshot are missing. They are requirements for the next model card,
   not optional details.

## Distribution decision

Package 5.4.0 is ready for public distribution as research software once the
maintainer creates the tagged release and attaches the CI-built wheel, source
distribution, and checksums. It is not ready for a clinical claim or for
presenting the current development benchmark as prospective validation.
