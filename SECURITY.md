# Security policy

## Trusted model files only

CdiffCladeML model bundles use joblib/pickle. Loading a malicious pickle can
execute arbitrary code before structural validation is possible. Never load a
model downloaded from an untrusted source.

Official releases ship the model with a JSON model card containing its
SHA-256. The predictor checks that digest before unpickling. A missing card
produces a warning; a mismatched digest stops execution.

Replacing a model requires moving both files together:

```text
my_model.joblib
my_model_model_card.json
```

Keep the exact scikit-learn version recorded by the model card. Cross-version
loading is rejected by default because it can fail or change predictions.

Do not use `--allow-version-mismatch` in production. That option exists only
for expert migration diagnostics followed by a full regression test.
