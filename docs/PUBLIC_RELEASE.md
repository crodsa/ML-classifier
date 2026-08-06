# Public-release requirements

## Scope

CdiffCladeML is open research software for screening seven-locus *C.
difficile* MLST profiles. It is not a clinical diagnostic device, does not
assign individual cryptic lineages, and does not replace sequence-based or
phylogenomic review.

## Data and model provenance

The source repository distributes the software, bundled model and model card;
it does not distribute the training or development-benchmark tables. The model
card records their hashes, cohort sizes and distance-matrix hashes. Any future
model release must document the PubMLST snapshot/date, clade-labelling criteria,
geographic and temporal provenance, deduplication, access conditions, and the
exact training and locked-evaluation identifiers or archive locations.

## Release checklist

Before publishing a GitHub Release:

1. Run the CI workflow successfully on all supported Python versions.
2. Build both the wheel and source distribution; run `twine check` and publish
   SHA-256 checksums.
3. Attach the wheel, source distribution and checksums to a signed/tagged
   release, then test installation in a fresh environment.
4. Verify that the model SHA-256 matches its model card and run the bundled
   six-profile regression test.
5. Confirm that documentation calls the 4,650-isolate set a development
   benchmark and clearly states the Clade2/Clade3 limitations.
6. Review all examples and issue templates for potentially identifiable data.
