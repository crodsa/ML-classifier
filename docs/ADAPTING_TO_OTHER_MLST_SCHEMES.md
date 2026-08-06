# Adapting CdiffCladeML to another MLST allele-profile classification

## What is reusable, and what is not

The reusable part of this repository is the **architecture**: allele-profile
parsing, per-locus distance features, grouped validation by sequence type (ST),
model provenance, and an HTML review report. The bundled `model_v5.joblib`, its
thresholds, its PubMLST matrices, and every performance result are specific to
seven-locus *Clostridioides difficile* clade classification. They are not valid
for another species, MLST scheme, or biological question.

Treat a new use as a new model-development study. Keep it in a fork or a
separately named package; do not overwrite the bundled model or present a
derivative model as CdiffCladeML v5 performance.

## Decide whether this architecture fits the question

The current workflow has two stages:

1. Stage 1 classifies profiles among **reference classes** (`CLASSICAL_CLADES`).
2. Stage 2 identifies a biologically divergent **outlier class** (`Cryptic`) by
   combining low Stage-1 support with distance from the reference profiles.

It is appropriate when the proposed outlier has a defensible biological meaning
(for example, a divergent lineage or a quality-controlled novel-profile group)
and the reference classes are well defined. It is not a drop-in solution for:

- a purely multiclass problem with no outlier category;
- classification from allele IDs when allele-to-allele biological distances are
  unavailable or unreliable;
- species or schemes whose labels are unstable, overlapping, or unsupported by
  representative independent STs.

For a pure multiclass problem, retain the parsing, reporting, grouped-split, and
provenance components, but remove or redesign the Stage-2 gate rather than
renaming `Cryptic`.

## Implementation checklist

### 1. Fork and name the derivative project

Create a new repository and document the organism, MLST scheme version, data
source, classification reference standard, and intended use. Change the package
name, command names, `PROGRAM`, README, citation metadata, and model-card
description so users cannot confuse it with the *C. difficile* release.

### 2. Define loci and classes in the source

Update `cdiff_clademl/config.py`:

- `LOCI`: ordered loci for the target MLST scheme;
- `ALL_CLADES`: all output labels, renamed for the target problem;
- `CLASSICAL_CLADES`: the Stage-1 reference labels. The remaining outlier label
  is currently hard-coded as `Cryptic` and must be renamed or parameterised.

The present reader assumes seven loci and therefore ten input columns
(`isolate`, `species`, `st`, plus seven loci), or eleven with the label.
Update the fixed `10`/`11` column checks, error text, header documentation,
and input examples in `cdiff_clademl/io.py`, `README.md`, and tests to use
`3 + len(LOCI)` and `4 + len(LOCI)` instead.

The terms `Cryptic`, `classical`, and the associated gate are also present in:

| Component | Required adaptation |
|---|---|
| `cdiff_clademl/cli/train.py` | Redefine the outlier label, its recall target, and the grouped cross-validation design. |
| `cdiff_clademl/cli/predict.py` | Rename final reasons, output fields, and review rules that currently mention `Cryptic`. |
| `cdiff_clademl/model.py` | Reassess the Stage-2 threshold objective and whether the distance-plus-low-support gate is biologically justified. |
| `cdiff_clademl/bundle.py` | Update bundle validation so its expected Stage-1 classes match the derivative configuration. |
| `cdiff_clademl/report.py` | Replace *C. difficile*/Cryptic-specific interpretation safeguards with target-specific wording. |
| Tests and smoke data | Replace all current clade labels, loci, profiles, and expected predictions. |

Do not merely edit displayed labels: a renamed but unchanged gate has no
established sensitivity, specificity, or biological interpretation.

### 3. Prepare the labelled allele-profile table

Use one row per isolate with a stable isolate identifier and ST. For a target
scheme with `L` loci, the training table has `3 + L + 1` columns:

```text
isolate,species,st,locus_1,...,locus_L,target_label
```

The current parser accepts exact allele calls, approximate calls prefixed by
`~`, and explicit missing values. For a derivative, document whether those
conventions remain valid. Each training ST and each complete allele profile
must map to exactly one target label; contradictory labels are a data curation
problem, not a model problem.

At minimum, the present five-fold design needs five distinct ST groups per
class. That is only a technical lower bound. Use substantially more independent
STs, especially for minority classes and the outlier group, before making a
scientific performance claim.

### 4. Build per-locus allele-distance matrices

Create one tab-separated file per locus, named `<locus>.tsv`, with this schema:

```text
allele_1    allele_2    distance
1           1           0.0
1           2           0.125
```

Distances must be finite, non-negative, symmetric, and biologically justified.
They must cover the observed training alleles. In the current implementation,
the reference set used to construct distance features consists only of Stage-1
reference classes; decide explicitly whether that is correct for the derivative
question. Allele numbers alone do not encode sequence distance, so do not
substitute numeric difference between allele identifiers for validated sequence
or allele-distance information.

### 5. Train, tune, and document the model

After adapting the source and input contract, train a new bundle with the
derivative's labelled table and distance directory:

```bash
<derivative>-train training_profiles.tsv \
  --distance_dir allele_distances/ \
  --metrics_prefix results/oof \
  -o models/derivative_model.joblib
```

Use grouped splits by ST (or a better justified genomic/epidemiological grouping)
so that related profiles do not occur in both train and validation partitions.
Tune all thresholds only within development data. Record the target-label
criteria, distance generation method, source database version/date, software
versions, hashes, and threshold-selection procedure in the derivative model
card.

### 6. Validate before release

Reserve a temporal, geographic, or otherwise independent locked cohort that is
not consulted while selecting features, model family, or thresholds. Report:

- per-class precision, recall, F1, and confusion matrix;
- one-vote-per-ST metrics as the primary analysis when isolates share STs;
- the number of unique STs and profiles per class;
- errors, uncertain calls, and outlier-gate review rates;
- performance limitations for sparse classes.

The HTML report and its star review marker are a review interface, not
validation evidence. Before public release, create a new model card and release
the model only with its matching SHA-256 sidecar.
