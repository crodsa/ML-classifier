# CdiffCladeML

CdiffCladeML classifies seven-locus *Clostridioides difficile* MLST profiles as
`Clade1`, `Clade2`, `Clade3`, `Clade4`, `Clade5`, or `Cryptic`. It is a
two-stage classifier: a Random Forest first assigns a classical clade, then a
gate based on allele distances and classical-model support detects Cryptic
profiles.

Routine use requires only the `CdiffCladeML-predict` command. Training is an
advanced operation intended for model maintainers.

> **Scope:** research and molecular-surveillance software. It is not a clinical
> diagnostic device and does not replace bioinformatic, epidemiological, or
> phylogenomic review.

## Quick start: classify MLST profiles

### 1. Requirements

- Python 3.10, 3.11, 3.12, or 3.13.
- The `cdiffclademl-5.4.0-py3-none-any.whl` file from this distribution.
- A CSV or TSV file containing MLST profiles.

You do not need to download a model separately. The wheel includes the
validated model and its provenance card.

### 2. Create an environment and install

On macOS or Linux, open a terminal in the directory containing the wheel:

```bash
python3 -m venv cdiffclademl-env
. cdiffclademl-env/bin/activate
python -m pip install --upgrade pip
python -m pip install ./cdiffclademl-5.4.0-py3-none-any.whl
```

On Windows PowerShell:

```powershell
py -3.11 -m venv cdiffclademl-env
.\cdiffclademl-env\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install .\cdiffclademl-5.4.0-py3-none-any.whl
```

Verify the installation:

```bash
CdiffCladeML-predict --version
CdiffCladeML-predict --help
```

### 3. Prepare the input file

CSV and TSV files are accepted; the delimiter is detected automatically. The
following header may be included but is optional:

```text
isolate,species,st,adk,atpA,dxr,glyA,recA,sodA,tpi
```

Example `profiles.csv`:

```csv
isolate,species,st,adk,atpA,dxr,glyA,recA,sodA,tpi
sample_001,cdifficile,29,1,1,2,16,1,3,1
sample_002,cdifficile,638,1,5,9,37,1,1,78
```

Input rules:

- Each row must contain exactly 10 columns in the order shown above.
- Alleles may be written as `1` or `adk(1)`.
- An approximate allele call may be written as `~1` or `adk(~1)`.
- `-`, `NA`, `N/A`, `?`, and an empty cell represent missing data.
- The program accepts MLST allele identifiers, not FASTA files, assemblies, or
  nucleotide sequences.
- For evaluation, an eleventh column named `mlst_clade` may contain the known
  label.

### 4. Run the classifier

```bash
CdiffCladeML-predict profiles.csv -o results/predictions.tsv
```

The output directory is created automatically. To display progress details:

```bash
CdiffCladeML-predict profiles.csv \
  -o results/predictions.tsv \
  --verbose
```

Each run creates:

- `results/predictions.tsv`: one detailed row per sample;
- `results/predictions_summary.json`: class counts, quality flags, runtime
  versions, and the model hash.

If the executable is not available on `PATH`, use the equivalent module form:

```bash
python -m cdiff_clademl.cli.predict profiles.csv \
  -o results/predictions.tsv
```

### 5. Interpret the results

The main output columns are:

| Field | Interpretation |
|---|---|
| `final_prediction` | Final result: `Clade1`–`Clade5` or `Cryptic` |
| `final_reason` | Stage that produced the result |
| `quality_flag` | Profile quality-control status |
| `review_recommended` | `True` when expert review is recommended |
| `stage1_confidence` | Maximum uncalibrated Random Forest support for a classical clade |
| `stage1_support_band` | Descriptive `HIGH`, `MEDIUM`, or `LOW` classical-support band |
| `cryptic_distance_call` | Whether the Cryptic gate was activated |
| `n_missing_loci` | Number of loci without an allele call |
| `n_approx_loci` | Number of approximate allele calls |
| `n_unknown_alleles` | Alleles absent from the model's distance matrices |

Possible `quality_flag` values are:

- `PASS`: complete profile with no unknown alleles;
- `PARTIAL_PROFILE`: at least one missing or approximate call;
- `UNKNOWN_ALLELE`: at least one allele is absent from the model matrices;
- `INSUFFICIENT_LOCI`: more than one locus is missing.

`stage1_confidence` and `stage1_support_band` describe classical-model support.
They are not calibrated probabilities that the final call is correct. A
`Cryptic` call may have low classical support by design.

Manually review every row with `review_recommended=True`, particularly profiles
with unknown alleles or missing data.

## Evaluate a labelled table

For an 11-column table whose final column is `mlst_clade`, generate evaluation
metrics with:

```bash
CdiffCladeML-predict labelled_profiles.tsv \
  -o results/predictions.tsv \
  --metrics_prefix results/evaluation
```

Accepted labels are `Clade1`, `Clade2`, `Clade3`, `Clade4`, `Clade5`, and
`Cryptic`. The command writes confusion matrices, per-class reports, and
aggregate metrics. It also calculates one-vote-per-ST results, which are more
appropriate when several isolates share the same ST.

## Use an updated model

An expert user can replace the bundled model without modifying the software:

```bash
CdiffCladeML-predict profiles.csv \
  --model models/model_2027.joblib \
  -o results/predictions.tsv
```

Keep these two files together:

```text
model_2027.joblib
model_2027_model_card.json
```

The model card records provenance, software versions, and a SHA-256 digest that
the classifier verifies before loading the model. Joblib uses pickle and can
execute code while loading: use only models you created or obtained from a
trusted source. Do not use `--allow-version-mismatch` during routine operation.

## Retrain the model: expert maintainers only

Training is not part of the routine workflow. It requires a table labelled with
all six categories and a directory containing one PubMLST distance matrix per
locus:

```text
pubmlst_distances/
├── adk.tsv
├── atpA.tsv
├── dxr.tsv
├── glyA.tsv
├── recA.tsv
├── sodA.tsv
└── tpi.tsv
```

Run:

```bash
CdiffCladeML-train train.csv \
  --distance_dir pubmlst_distances/ \
  --metrics_prefix results/oof \
  --n_jobs 5 \
  -o models/model_2027.joblib \
  --verbose
```

Training requires at least five distinct ST groups per category, consistent
labels and profiles, and valid matrices covering the observed alleles. The
model and its `*_model_card.json` are created together. Before replacing the
official model, repeat ST-grouped validation and evaluate an independent locked
cohort that was not used to select the architecture or thresholds.

## Troubleshooting

| Message or situation | Recommended action |
|---|---|
| `command not found` | Reactivate the environment or use `python -m cdiff_clademl.cli.predict` |
| `Input table does not exist` | Check the path and input filename |
| `must have 10 or 11 columns` | Check the delimiter, empty cells, and column order |
| `Unknown mlst_clade` | Use exactly one of the six accepted labels |
| scikit-learn version mismatch | Create a clean environment and reinstall the wheel, which pins the compatible version |
| `UNKNOWN_ALLELE` in the output | This is not an execution failure; review the profile and consider updating the model |
| Model hash error | Do not use that model/card pair; obtain an intact copy from a trusted source |

## Bundled model and validation evidence

Software release 5.4 bundles `model_v5.joblib`, model version 5.3.0. On the
independent test set, the primary one-vote-per-ST analysis produced accuracy
0.9373, balanced accuracy 0.9416, and macro-F1 0.9527. It detected 25/25
Cryptic STs; `Clade2` recall was 43/59, and `Clade3` had only three STs, so
these classes require particular caution.

See [`docs/VALIDATION.md`](docs/VALIDATION.md) for the methodology, confusion
matrices, limitations, and future model-release criteria. See
[`docs/ENGINEERING_AUDIT.md`](docs/ENGINEERING_AUDIT.md) for the software
engineering review. Model provenance and hashes are recorded in
`model_v5_model_card.json`.

## Package development and verification

```bash
uv sync
uv run ruff check .
uv run pytest --cov=cdiff_clademl
uv build
uv run twine check dist/*
```

MIT License. Read [`SECURITY.md`](SECURITY.md) before distributing third-party
joblib models.
