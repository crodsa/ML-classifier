# CdiffCladeML (`v1.0.0`)

**Rapid, sequence-free clade screening for *Clostridioides difficile* seven-locus MLST allele tables (`Clade1`–`Clade5` + `Cryptic`).**

`CdiffCladeML` takes standard **FastMLST** or **PubMLST** allele-number tables (`.csv` or `.tsv`) and classifies thousands of isolates in seconds—without needing genome assemblies or nucleotide FASTA files. Every prediction comes with an interactive `report.html`, quality flags for novel (`~n`) or missing alleles, Stage-1 support bands, and an explicit `review_recommended` flag.

> **Part of a Two-Tier *C. difficile* MLST Suite:**
> - **Tier 1 — `CdiffCladeML` (this repository):** Screens tabular allele numbers in seconds (`Clade1`–`Clade5` vs pooled `Cryptic`) and flags uncertain or novel-allele profiles.
> - **Tier 2 — [`TreeTip Placer`](https://github.com/crodsa/treetip-placer):** Takes 7-locus nucleotide FASTA sequences for isolates flagged by Tier 1 (`Cryptic`, `~n` novel alleles, or `< 0.80` support) to resolve individual cryptic clades (`C-I` to `C-VII`) and verify `Clade 1`/`Clade 2` boundary strains on fixed reference trees.

---

## Quick Start

### Requirements
- **Python** `3.10`–`3.13` (works on Linux, macOS, and Windows; no external bioinformatics binaries required).

### Option A: Install the standalone release wheel (fastest)
Download `cdiffclademl-1.0.0-py3-none-any.whl` from the [`v1.0.0` GitHub Release](https://github.com/crodsa/ML-classifier/releases/tag/v1.0.0) and install it in a clean virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate          # On Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install cdiffclademl-1.0.0-py3-none-any.whl
CdiffCladeML-predict --version    # Prints: CdiffCladeML-predict 1.0.0
```

### Option B: Install from GitHub and run the bundled demo
```bash
git clone https://github.com/crodsa/ML-classifier.git
cd ML-classifier
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .

# Run the 6-isolate demo table (1 isolate from each output class)
CdiffCladeML-predict tests/data/smoke_profiles.csv
```

This creates a timestamped folder (`smoke_profiles_cdiffclademl_<UTC_TIMESTAMP>/`) containing:
- **`report.html`** — open in any web browser for summary charts and a colour-coded isolate table.
- **`predictions.tsv`** — complete spreadsheet/pipeline-ready results (one row per isolate).
- **`predictions_summary.json`** — run metadata, class counts, software version, and SHA-256 model provenance.

---

## Preparing Your Input Table

### Direct FastMLST compatibility
If you run [**FastMLST**](https://github.com/EnzoAndree/FastMLST) on genome assemblies, you can pass its CSV/TSV output **directly** to `CdiffCladeML-predict` with no editing.

### Column order (10 or 11 columns, CSV or TSV, header optional)
`CdiffCladeML` auto-detects comma (`,`) or tab (`\t`) delimiters. Each row must contain **10 columns** (or **11 columns** if you include known `mlst_clade` labels for benchmarking):

```text
isolate,species,st,adk,atpA,dxr,glyA,recA,sodA,tpi[,mlst_clade]
```

| Columns | Example values | Notes |
| :--- | :--- | :--- |
| `1–3` (`isolate`, `species`, `st`) | `G13S1,Clostridioides difficile,29` | Retained for traceability; `st` can be `NA` or `new_alleles` (`st` is not used by the classifier) |
| `4–10` (`adk`..`tpi`) | `1,1,2,16,1,3,1` or `adk(1),atpA(~1),...` | **Used for classification.** Accepts plain numbers (`1`), FastMLST format (`adk(1)`), `adk=1`, `adk:1`, novel/approximate calls (`~1` or `adk(~1)`), and missing markers (`-`, `NA`, `?`) |
| `11` (`mlst_clade`, *optional*) | `Clade1` or `mlst_clade(1.0)` | Optional ground-truth label (`Clade1`–`Clade5` or `Cryptic`) to compute accuracy metrics |

### Minimal copy-paste example
Save as `my_profiles.csv`:

```csv
isolate,species,st,adk,atpA,dxr,glyA,recA,sodA,tpi
sample_clinical_1,Clostridioides difficile,1,1,1,1,10,1,3,5
sample_soil_novel,Clostridioides difficile,new_alleles,adk(~13),atpA(46),dxr(64),glyA(31),recA(18),sodA(32),tpi(26)
```

Run:

```bash
CdiffCladeML-predict my_profiles.csv
```

To write outputs to fixed paths for an automated pipeline (never overwrites existing files):

```bash
CdiffCladeML-predict my_profiles.csv \
  -o results/predictions.tsv \
  --report-html results/report.html
```

---

## Interpreting Your Results & Two-Tier Routing

Start by opening **`report.html`** in your browser, or inspect these five key columns in **`predictions.tsv`**:

| Column | Possible values | How to interpret & what to do next |
| :--- | :--- | :--- |
| **`final_prediction`** | `Clade1`–`Clade5`, `Cryptic` | Predicted classical *C. difficile* clade, or pooled `Cryptic` (*Clostridioides* cryptic clades `C-I`–`C-VII`). |
| **`quality_flag`** | `PASS`<br>`PARTIAL_PROFILE`<br>`UNKNOWN_ALLELE`<br>`INSUFFICIENT_LOCI` | **All 7 alleles exact & known** (`PASS`).<br>**1 missing or `>=1` novel `~n` allele** (`PARTIAL_PROFILE`).<br>**Allele ID newer than bundled 918-allele cache** (`UNKNOWN_ALLELE`).<br>**`>=2` missing loci** (`INSUFFICIENT_LOCI` — re-check assembly coverage; do not rely on clade call). |
| **`stage1_support_band`** | `HIGH` (`>= 0.80`)<br>`MEDIUM` (`0.55–0.79`)<br>`LOW` (`< 0.55`) | Random Forest support among classical Clades 1–5 (`stage1_confidence`). `MEDIUM`/`LOW` on a classical call typically indicates a `Clade 1`/`Clade 2` boundary ST. |
| **`review_recommended`** | `False` / `True` | `True` whenever `stage1_support_band == LOW`, `quality_flag != PASS`, or `Cryptic` has `< 4` non-classical loci. |
| **`final_reason`** | `classical_stage1`<br>`distance+low_classical_support` | Shows whether the isolate stayed in Stage 1 (classical RF) or triggered the Stage-2 evolutionary-distance Cryptic gate. |

### Recommended operational workflow
1. **Finalize directly at Tier 1 (~81% of routine isolates):** `final_prediction` in `Clade1`–`Clade5` **AND** `stage1_support_band == HIGH` **AND** `quality_flag == PASS`.
2. **Forward to Tier 2 ([`TreeTip Placer`](https://github.com/crodsa/treetip-placer)):** Any isolate predicted as `Cryptic` (to resolve `C-I` through `C-VII`), any classical isolate with `MEDIUM`/`LOW` support (`< 0.80`, common on `Clade 1`/`Clade 2` boundary STs), or any `PARTIAL_PROFILE` / `UNKNOWN_ALLELE` isolate.

---

## How It Works & Validated Performance

Standard allele-number classifiers use one-hot encoding (`OHE`), which maps unseen alleles in divergent cryptic strains to all-zero vectors and misclassifies them into the largest clinical clade (`Clade1`). `CdiffCladeML v1.0.0` solves this with a two-stage architecture:

1. **Stage 1 (Classical Random Forest):** 500-tree balanced Random Forest trained on `25,020` classical PubMLST isolates (`81` exact alleles across `Clade1`–`Clade5`).
2. **Stage 2 (Evolutionary-Distance Cryptic Gate):** Evaluates precomputed pairwise nucleotide distances across **918 PubMLST alleles** at all 7 loci. When an isolate is genetically distant from classical alleles/profiles **and** has low Stage-1 support (`<= 0.4455` or top-two margin `<= 0.0100`, with `<= 1` missing locus), Stage 2 overrides the call to `Cryptic`.

### Benchmark summary (`v1.0.0`)
- **Held-out unseen PubMLST STs (`test_v2.csv`: 4,650 isolates across 319 unseen STs, zero overlap with training STs):**
  - **Isolate accuracy:** `96.95%` (`4,508/4,650`; macro-F1 `0.9353`)
  - **Equal-weight ST accuracy:** `93.73%` (`299/319`; macro-F1 `0.9527`)
  - **Cryptic ST recall:** `100.00%` (`25/25`, vs `0.00%` [`0/25`] with Stage-1 RF alone)
  - **When paired with `TreeTip Placer` (routed two-tier mode):** `98.12%` ST accuracy (`313/319`) and `99.59%` isolate accuracy (`4,631/4,650`) while finalizing `80.77%` of isolates directly in Tier 1.
- **Tropical soil surveillance (48 Costa Rican *Clostridioides* genomes, ENA BioProject [`PRJEB103842`](https://www.ebi.ac.uk/ena/browser/view/PRJEB103842)):** `100.0%` (`48/48`) coarse clade concordance with whole-genome `FastANI` clusters (`7/7` classical, `41/41` cryptic) and automatic `PARTIAL_PROFILE` flagging of all `12/12` novel-allele (`~n`) environmental isolates.

See [`docs/VALIDATION.md`](docs/VALIDATION.md) and [`model_model_card.json`](model_model_card.json) for full per-clade confusion matrices, thresholds, and SHA-256 digests.

---

## CLI Reference & Advanced Usage

Both `CdiffCladeML-predict` (alias `cdiff-predict`) and `CdiffCladeML-train` (alias `cdiff-train`) are installed automatically:

```bash
# Routine prediction (creates timestamped output folder)
CdiffCladeML-predict profiles.csv

# Benchmark a labelled 11-column table and export isolate- + ST-level metrics
CdiffCladeML-predict labelled_profiles.tsv \
  -o eval/predictions.tsv \
  --metrics_prefix eval/metrics

# Retrain a custom model bundle on an updated PubMLST snapshot + 7 locus distance TSVs
CdiffCladeML-train training_profiles.tsv \
  --distance_dir distance_matrices/ \
  -o models/custom_model.joblib \
  --metrics_prefix models/oof
```

### Quick troubleshooting
- **`Every row ... must have 10 or 11 columns`:** Ensure your CSV/TSV has no extra metadata columns and follows `isolate,species,st,adk,atpA,dxr,glyA,recA,sodA,tpi[,mlst_clade]`.
- **`Model requires scikit-learn 1.7.1`:** `model.joblib` verifies exact `scikit-learn==1.7.1` compatibility and SHA-256 integrity on load. Reinstall inside a dedicated virtual environment (`python -m pip install --force-reinstall .`).
- **Additional docs:** [`docs/VALIDATION.md`](docs/VALIDATION.md) · [`docs/MODEL_PROVENANCE.md`](docs/MODEL_PROVENANCE.md) · [`docs/ADAPTING_TO_OTHER_MLST_SCHEMES.md`](docs/ADAPTING_TO_OTHER_MLST_SCHEMES.md) · [`SECURITY.md`](SECURITY.md)

---

## Citation & License

If you use `CdiffCladeML` in your research or surveillance workflow, please cite [`CITATION.cff`](CITATION.cff) and the accompanying manuscript:

> Rodríguez C. *Two-tier multilocus sequence typing clade classification and fixed-tree phylogenetic placement resolve classical and cryptic Clostridioides difficile lineages.* (2026). Companion Tier-2 repository: https://github.com/crodsa/treetip-placer · Environmental validation genomes: ENA BioProject [`PRJEB103842`](https://www.ebi.ac.uk/ena/browser/view/PRJEB103842).

Distributed under the **[MIT License](LICENSE)** (`model.joblib` SHA-256: `55c20798d584d09b854417bda67879b33247c0ab7629d2205dbb8102b9c70564`).
