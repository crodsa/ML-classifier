# Model and distance-matrix provenance

## Distribution statement

The repository author created the training curation, the seven per-locus MLST distance matrices, and the bundled `model.joblib` (`v1.0.0`) model, and distributes the derived model bundle and embedded distance matrices under the MIT terms in `LICENSE`. The underlying curated training (`train_v2.csv`) and held-out test (`test_v2.csv`) tables, unannotated PubMLST screen (`nan.mlst`), and per-locus allele caches are publicly distributed in the companion Zenodo / Supplementary Data 1 reproducibility archive (`v1.0.0`; CC BY 4.0), alongside environmental validation genome accessions under ENA BioProject `PRJEB103842`.

Dataset hashes, cohort sizes, and model/distance-matrix SHA-256 digests are recorded in `model_model_card.json` for full cryptographic traceability.

## Distributed artifacts

| Artifact | License / Status | Verification |
| --- | --- | --- |
| Source code and documentation | MIT | `LICENSE` and package metadata |
| `model.joblib` (`v1.0.0`) | MIT | SHA-256 recorded in `model_model_card.json` |
| Seven per-locus distance matrices (embedded in `model.joblib`) | MIT | Per-locus SHA-256 digests in `model_model_card.json` |

## Scientific scope

`CdiffCladeML` (`v1.0.0`) is research and molecular-surveillance software for seven-locus *Clostridioides difficile* MLST profile screening; see `docs/VALIDATION.md` for benchmark metrics and interpretation guidelines.
