# Validation and methodology summary

## Overview

The bundled `model.joblib` (`v1.0.0`) classifies seven-locus *Clostridioides difficile* MLST allele profiles into `Clade1`–`Clade5` or a combined `Cryptic` category. Validation was conducted on a strictly non-overlapping held-out PubMLST benchmark (`test_v2.csv`; 4,650 isolates across 319 unique sequence types [STs]). The recommended use is research and genomic surveillance; the software is not intended for clinical diagnostic decisions.

## Benchmark split independence

Direct verification of the training (`train_v2.csv`) and held-out test (`test_v2.csv`) datasets confirmed zero overlap across isolates, STs, and seven-locus profiles:

| Unit | Training | Held-out test | Overlap |
|---|---:|---:|---:|
| Isolates | 25,124 | 4,650 | 0 |
| Sequence types (STs) | 130 | 319 | 0 |
| Seven-locus allele profiles | 145 total (42 classical) | 319 | 0 |

Because the 4,650 held-out isolates represent 319 STs with varying sampling depth in PubMLST, the primary evaluation gives equal weight to each unique ST (`n = 319`), while the isolate-weighted evaluation (`n = 4,650`) is reported alongside it.

## Held-out benchmark performance

### One vote per ST (n = 319)

| True class | STs | Correct | Recall | 95% Wilson CI |
|---|---:|---:|---:|---:|
| Clade1 | 179 | 176 | 0.9832 | 0.9519–0.9943 |
| Clade2 | 59 | 43 | 0.7288 | 0.6040–0.8256 |
| Clade3 | 3 | 3 | 1.0000 | 0.4385–1.0000 |
| Clade4 | 37 | 37 | 1.0000 | 0.9059–1.0000 |
| Clade5 | 16 | 15 | 0.9375 | 0.7167–0.9889 |
| Cryptic | 25 | 25 | 1.0000 | 0.8668–1.0000 |

Aggregate ST-weighted metrics: accuracy `0.9373` (95% Wilson CI `0.9052–0.9591`), balanced accuracy `0.9416`, macro-F1 `0.9527`, and weighted F1 `0.9346`.

Confusion matrix by ST (true classes in rows, predictions in columns):

| | Clade1 | Clade2 | Clade3 | Clade4 | Clade5 | Cryptic |
|---|---:|---:|---:|---:|---:|---:|
| Clade1 | 176 | 3 | 0 | 0 | 0 | 0 |
| Clade2 | 16 | 43 | 0 | 0 | 0 | 0 |
| Clade3 | 0 | 0 | 3 | 0 | 0 | 0 |
| Clade4 | 0 | 0 | 0 | 37 | 0 | 0 |
| Clade5 | 0 | 0 | 0 | 0 | 15 | 1 |
| Cryptic | 0 | 0 | 0 | 0 | 0 | 25 |

### Isolate-weighted analysis (n = 4,650)

Aggregate isolate-weighted metrics: accuracy `0.9695` (`4,508/4,650`), balanced accuracy `0.9161`, macro-F1 `0.9353`, and weighted F1 `0.9663`.

## Cross-validation and distance-matrix verification

`CdiffCladeML-train` partitions training data by ST (`StratifiedGroupKFold` for classical clades and `GroupKFold` for Cryptic STs), reconstructs classical allele reference sets inside each fold to prevent cross-validation leakage, and cross-fits the Stage-2 Cryptic gate before fitting final deployment rules across all out-of-fold predictions. All seven per-locus distance matrices (`adk`, `atpA`, `dxr`, `glyA`, `recA`, `sodA`, `tpi`; 918 total alleles) are verified for zero diagonals, symmetry, finite non-negative distances, and allele coverage.

## Scope and limitations

- `Clade2` STs with novel or unshared classical alleles can fall near the `Clade1`/`Clade2` boundary; 14 of the 16 misclassified `Clade2` STs receive `MEDIUM` or `LOW` Stage-1 support and are resolved when paired with sequence-based placement in `TreeTip Placer`.
- `Clade3` has three STs in the held-out test set, resulting in a wide confidence interval.
- The `Cryptic` output category is a combined screening label and does not resolve individual cryptic clades (`C-I`–`C-VII`); use `TreeTip Placer` on concatenated seven-locus FASTA sequences for fine-scale cryptic lineage assignment.
- Stage-1 Random Forest probabilities describe support among the five classical clades and are not calibrated posterior probabilities of the final call.
