# Validation and methodology audit

Review date: 2026-07-20.

## Conclusion

The bundled model provides initial evidence of utility on unseen MLST profiles,
particularly for the combined `Cryptic` category. This evidence comes from
`test_v2.csv`, not from the historical out-of-fold (OOF) metrics, and must be
interpreted as a development benchmark because the split informed iterative
development. The recommended use is research and surveillance; the evidence
does not support clinical use or a precise performance claim for `Clade3`.

## Development-benchmark split independence

Direct review of `train_v2.csv` and `test_v2.csv` verified:

| Unit | Training | Test | Overlap |
|---|---:|---:|---:|
| Isolates | 25,124 | 4,650 | 0 |
| STs | 130 | 319 | 0 |
| Seven-locus allele profiles | 145 total (42 classical) | 319 | 0 |

The 4,650 test isolates represent 319 repeated STs. Therefore, the primary
analysis gives one vote to each ST, while the isolate-weighted analysis is
reported as secondary.

## Development-benchmark performance

### One vote per ST (n = 319)

| True class | STs | Correct | Recall | 95% Wilson CI |
|---|---:|---:|---:|---:|
| Clade1 | 179 | 176 | 0.9832 | 0.9519–0.9943 |
| Clade2 | 59 | 43 | 0.7288 | 0.6040–0.8256 |
| Clade3 | 3 | 3 | 1.0000 | 0.4385–1.0000 |
| Clade4 | 37 | 37 | 1.0000 | 0.9059–1.0000 |
| Clade5 | 16 | 15 | 0.9375 | 0.7167–0.9889 |
| Cryptic | 25 | 25 | 1.0000 | 0.8668–1.0000 |

Aggregate metrics: accuracy 0.9373 (95% Wilson CI 0.9052–0.9591), balanced
accuracy 0.9416, macro-F1 0.9527, and weighted F1 0.9346. The binomial interval
by ST does not account for residual phylogenetic dependence among related STs,
so it may remain optimistic when clonal complexes are present.

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

Accuracy was 0.9695, balanced accuracy 0.9161, macro-F1 0.9353, and weighted
F1 0.9663. These figures describe the isolate set, but their intervals must not
be calculated under the assumption of 4,650 independent observations.

## Critical finding in the historical OOF metrics

The historical `_oof_meta.tsv` file contains 25,540 rows for 25,124 samples.
The 25,020 classical samples occur once, but each of the 104 `Cryptic` samples
occurs in all five folds, producing a reported support of 520. In addition, the
Stage 2 thresholds were selected and evaluated on that same table, and the
classical allele sets used by some distance features came from the entire
training set.

Consequences:

1. the Cryptic category was pseudoreplicated;
2. distance features contained cross-validation leakage;
3. post-threshold-selection performance was optimistic.

Those metrics are retained only as a historical tuning diagnostic. Trainer 5.4
assigns each ST, including every Cryptic ST, to one fold; learns allele
references within each fold; and applies cross-fitting to the gate so that a
sample is not evaluated directly with the thresholds fitted from its own
prediction. This diagnostic is not a strict nested cross-validation of the
complete system and does not replace an external cohort. Final deployment
thresholds are fitted afterward using all OOF predictions.

## Data and distance-matrix quality

In the reviewed datasets, each ST and each allele profile maps to exactly one
label. All seven matrices have zero diagonals, finite non-negative distances,
and complete coverage of training and test alleles. The 5.4 loader now verifies
these invariants and rejects contradictory duplicate distances.

The classical training collection is highly imbalanced in ST diversity:
Clade1=18, Clade2=6, Clade3=6, Clade4=6, and Clade5=6. It contains 88 Cryptic
STs. A large isolate count for a few STs does not represent equivalent
phylogenetic diversity. Future updates should prioritize additional independent
STs from the minority classical clades, not merely more isolates from existing
STs.

## Limitations and criteria for the next model release

- The no-overlap test split was consulted repeatedly while developing several
  versions and is therefore a development benchmark, not independent final
  validation. A future major release requires a final locked temporal or
  geographic cohort that is not used for architecture or threshold selection.
- `Clade3` has only three STs in the test set. Its perfect result is not a
  stable performance estimate.
- The `Cryptic` label combines cryptic lineages and does not identify an
  individual cryptic lineage.
- Random Forest probabilities are uncalibrated. Support bands describe the
  classical stage and are not posterior probabilities of the final call.
- A new allele identifier provides no nucleotide sequence. Its fallback
  distance is heuristic and triggers an expert-review flag.
- The retained v5.3.0 model reports the number of rare loci but does not use a
  four-rare-locus criterion as a hard Cryptic-gate condition. The predictor
  flags Cryptic calls below that threshold for review; making it a hard gate
  requires retraining and a new locked-cohort evaluation.
- The next dataset must document the PubMLST version and date, clade-labelling
  criteria, geographic and temporal provenance, and deduplication procedure.

Before replacing the bundled model, require a model card with hashes and
software versions, leakage-free grouped cross-validation, one-vote-per-ST
metrics, evaluation on a locked cohort, and explicit review of `Clade2` and
`Clade5` errors.
