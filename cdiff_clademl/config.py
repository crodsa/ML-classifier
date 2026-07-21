"""
Shared constants: loci, clade lists, and derived column name sets.
All downstream modules import from here — never hard-code these values elsewhere.
"""

LOCI = ["adk", "atpA", "dxr", "glyA", "recA", "sodA", "tpi"]

BASE_COLS = [f"{locus}_base" for locus in LOCI]
STATUS_COLS = [f"{locus}_status" for locus in LOCI]
TOKEN_COLS = [f"{locus}_token" for locus in LOCI]

ALL_CLADES = ["Clade1", "Clade2", "Clade3", "Clade4", "Clade5", "Cryptic"]
CLASSICAL_CLADES = ["Clade1", "Clade2", "Clade3", "Clade4", "Clade5"]

# Columns written to every prediction output
META_FEATURE_COLS = [
    "stage1_predicted_clade",
    "stage1_confidence",
    "stage1_margin",
    "stage1_entropy",
    "sum_min_dist_any_classical",
    "nearest_profile_dist_any_classical",
    "n_missing_loci",
    "n_approx_loci",
    "n_unknown_alleles",
    "rare_locus_count_report",
]

# Confidence band thresholds applied to stage1_confidence in predictions
CONFIDENCE_BAND_THRESHOLDS = {
    "HIGH": 0.80,  # stage1_confidence >= 0.80
    "MEDIUM": 0.55,  # 0.55 <= confidence < 0.80
    # Below 0.55 → LOW
}
