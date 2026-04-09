# GIF Framework — Genomic Intelligence Framework
# Copyright 2026 Javier Gamboa
# Licensed under the Apache License, Version 2.0
"""
P-Score (Persistence) calculation for the GIF Framework.

Scoring levels:
  Level 1  — Genetic Markers (0-75 cap)
  Level 1b — Operational Persistence (0-15 pts within the 75 cap)
  Level 2  — Temporal / Spatial (0-50 pts, if metadata available)
  Bayesian imputation when metadata is absent

Final score normalised by P_NORM_MAX (75) to [0, 100].
"""

from typing import Dict, Optional, Tuple

from gif.scoring.weights import P_NORM_MAX

# =============================================================================
# P-SCORE CONSTANTS — Level 1: Genetic Markers
# =============================================================================

P_SCORE_QACEDELTA1 = 15       # qacEDelta1 — THE MOST IMPORTANT marker
P_SCORE_QACH = 10             # qacH
P_SCORE_BCRABC = 20           # bcrABC complete operon (all 3 genes)
P_SCORE_SSI1 = 15             # SSI qualification (>= 6/10 genes across SSI-1/SSI-2/GAD)
P_SCORE_CADMIUM = 5           # Cadmium resistance (cadA variants) — AEM 2024
P_SCORE_LIX_OPERON = 5        # lix operon (DISABLED — genes don't exist in Listeria)
P_SCORE_PROPHAGES = 5         # >2 attenuated prophages
P_SCORE_CRISPR_ABSENT = 5     # CRISPR absence facilitates exogenous DNA uptake
P_SCORE_GENETIC_CAP = 75      # Paper v2.4: max theoretical 75 pts

# SSI combined threshold: >= 6/10 genes (SSI-1: 5, SSI-2: 2, GAD: 3)
SSI_THRESHOLD = 6

# =============================================================================
# P-SCORE CONSTANTS — Level 1b: Operational Persistence (0-15 pts within cap)
# Only ACCESSORY genes that discriminate between strains.
# Core genes (dnaK, groEL, sigB, kat, sod) -> 0 pts (present in 100% of strains)
# =============================================================================

P_SCORE_BAPL = 2               # bapL surface adhesion protein (variable)
P_SCORE_AGR_COMPLETE = 1       # Complete agr QS operon (>= 3 genes)
P_SCORE_EPS_CLUSTER = 2        # EPS biosynthesis cluster (>= 4 pss genes)
P_SCORE_ACID_NON_SSI = 2       # Non-SSI acid tolerance (arcA, gadD2, gadT2)
P_SCORE_ACCESSORY_EFFLUX = 3   # Accessory efflux pumps (emrC/emrE/qacC/sugE)
P_SCORE_COLD_SHOCK = 2         # Cold-shock proteins >= 2/3 cspL/cspB/cspD (v2.3)
P_SCORE_OXIDATIVE = 2          # Oxidative stress tolerance kat + >= 3 others (v2.3)
P_SCORE_IS1216 = 2             # IS1216 transposon (dissemination)
P_SCORE_LGI1 = 2              # LGI-1 genomic island (dissemination)
P_SCORE_ACTIVE_PROPHAGE = 1    # Active prophage markers (tRNA_int + phage_portal)
P_SCORE_MOTILITY = 1           # Complete flagellar system flaA+motA+motB (v2.3)

# =============================================================================
# P-SCORE CONSTANTS — Level 2: Temporal / Spatial
# =============================================================================

P_SCORE_BAYESIAN_IMPUTATION_FACTOR = 0.60  # Recalibrated v2.3 — Expert Review SS5.3

P_SCORE_TIMESPAN_GT52W = 30    # >52 weeks (>1 year)
P_SCORE_TIMESPAN_13_52W = 20   # 13-52 weeks
P_SCORE_TIMESPAN_4_12W = 10    # 4-12 weeks
P_SCORE_DETECTIONS_GTE5 = 10   # >= 5 independent detections
P_SCORE_DETECTIONS_GTE3 = 5    # >= 3 independent detections
P_SCORE_SDS_GT05 = 10          # SDS > 0.5
P_SCORE_SDS_GTE02 = 5          # SDS >= 0.2


def calculate_p_score(
    markers: dict,
    metadata: Optional[dict] = None,
) -> Tuple[float, Dict]:
    """
    Calculate P-Score (Persistence) according to the GIF specification.

    Args:
        markers: Dict with genomic marker fields.  Relevant keys:
            - qac_e_delta1_detected (bool)
            - qac_h_detected (bool)
            - bcr_abc_detected (bool)
            - ssi1_detected (bool)  — kept for legacy; combined detection below
            - ssi1_genes_count (int)
            - ssi2_genes_count (int)
            - gad_genes_count (int)
            - cadmium_resistance_detected (bool)
            - prophages_count (int)
            - crispr_detected (bool)
            - biofilm_bapL_detected (bool)
            - biofilm_agr_complete (bool)
            - biofilm_eps_cluster_count (int)
            - biocide_acid_non_ssi (list)
            - biocide_accessory_efflux (list)
            - cold_shock_detected (bool)
            - oxidative_stress_tolerant (bool)
            - mobility_IS1216_detected (bool)
            - mobility_LGI1_detected (bool)
            - mobility_active_prophage (bool)
            - flagellar_motility_complete (bool)

        metadata: Optional dict with temporal/spatial data:
            - timespan_weeks (int)
            - independent_detections (int)
            - unique_zones (int)
            - total_zones_sampled (int)

    Returns:
        (score 0-100, details dict)
    """
    details: Dict = {}

    # =========================================================================
    # LEVEL 1: Genetic Markers (cap at 75)
    # =========================================================================
    genetic_score = 0

    # qacEDelta1
    qac_e = markers.get("qac_e_delta1_detected", False)
    if qac_e:
        genetic_score += P_SCORE_QACEDELTA1
    details["qac_e_delta1_detected"] = qac_e
    details["qac_e_delta1_score"] = P_SCORE_QACEDELTA1 if qac_e else 0

    # qacH
    qac_h = markers.get("qac_h_detected", False)
    if qac_h:
        genetic_score += P_SCORE_QACH
    details["qac_h_detected"] = qac_h
    details["qac_h_score"] = P_SCORE_QACH if qac_h else 0

    # bcrABC (complete operon — all 3 genes)
    bcr = markers.get("bcr_abc_detected", False)
    if bcr:
        genetic_score += P_SCORE_BCRABC
    details["bcr_abc_detected"] = bcr
    details["bcr_abc_score"] = P_SCORE_BCRABC if bcr else 0

    # SSI Detection (combined SSI-1 + SSI-2 + GAD, >= 6/10 threshold)
    ssi1_count = markers.get("ssi1_genes_count", 0)
    ssi2_count = markers.get("ssi2_genes_count", 0)
    gad_count = markers.get("gad_genes_count", 0)
    total_ssi_genes = ssi1_count + ssi2_count + gad_count
    ssi_qualifies = total_ssi_genes >= SSI_THRESHOLD

    if ssi_qualifies:
        genetic_score += P_SCORE_SSI1

    details["ssi1_detected"] = ssi_qualifies
    details["ssi1_genes_count"] = ssi1_count
    details["ssi2_genes_count"] = ssi2_count
    details["gad_genes_count"] = gad_count
    details["total_ssi_genes"] = total_ssi_genes
    details["ssi1_score"] = P_SCORE_SSI1 if ssi_qualifies else 0

    # Cadmium resistance (cadA1-cadA4)
    cadmium = markers.get("cadmium_resistance_detected", False)
    cadmium_score = P_SCORE_CADMIUM if cadmium else 0
    if cadmium:
        genetic_score += cadmium_score
    details["cadmium_resistance_detected"] = cadmium
    details["cadmium_score"] = cadmium_score

    # Adaptation factors
    adaptation_score = 0

    # lix operon — DISABLED (genes don't exist in Listeria nomenclature)
    details["lix_operon_detected"] = False

    # Prophages (>2 attenuated prophages)
    prophages = markers.get("prophages_count", 0)
    if prophages > 2:
        adaptation_score += P_SCORE_PROPHAGES
    details["prophages_count"] = prophages

    # CRISPR absence bonus (Paper v2.4 — restored)
    # Only ~41% of L. monocytogenes carry functional CRISPR-Cas
    crispr = markers.get("crispr_detected", False)
    if not crispr:
        adaptation_score += P_SCORE_CRISPR_ABSENT
    details["crispr_detected"] = crispr
    details["crispr_absent_bonus"] = P_SCORE_CRISPR_ABSENT if not crispr else 0

    genetic_score += adaptation_score
    details["adaptation_score"] = adaptation_score

    # =========================================================================
    # LEVEL 1b: Operational Persistence (0-15 pts within existing cap)
    # =========================================================================
    operational_score = 0

    # Biofilm accessory module (0-5 pts)
    biofilm_op = 0
    if markers.get("biofilm_bapL_detected", False):
        biofilm_op += P_SCORE_BAPL
    if markers.get("biofilm_agr_complete", False):
        biofilm_op += P_SCORE_AGR_COMPLETE
    if markers.get("biofilm_eps_cluster_count", 0) >= 4:
        biofilm_op += P_SCORE_EPS_CLUSTER
    biofilm_op = min(5, biofilm_op)
    operational_score += biofilm_op
    details["operational_biofilm_score"] = biofilm_op

    # Environmental survival module (0-5 pts)
    survival_op = 0
    acid_non_ssi = markers.get("biocide_acid_non_ssi", [])
    if len(acid_non_ssi) >= 2:
        survival_op += P_SCORE_ACID_NON_SSI
    accessory_efflux = markers.get("biocide_accessory_efflux", [])
    if len(accessory_efflux) >= 1:
        survival_op += P_SCORE_ACCESSORY_EFFLUX
    if markers.get("cold_shock_detected", False):
        survival_op += P_SCORE_COLD_SHOCK
    if markers.get("oxidative_stress_tolerant", False):
        survival_op += P_SCORE_OXIDATIVE
    survival_op = min(5, survival_op)
    operational_score += survival_op
    details["operational_survival_score"] = survival_op
    details["cold_shock_detected"] = markers.get("cold_shock_detected", False)
    details["oxidative_stress_tolerant"] = markers.get("oxidative_stress_tolerant", False)

    # Dissemination potential module (0-5 pts)
    dissemination_op = 0
    if markers.get("mobility_IS1216_detected", False):
        dissemination_op += P_SCORE_IS1216
    if markers.get("mobility_LGI1_detected", False):
        dissemination_op += P_SCORE_LGI1
    if markers.get("mobility_active_prophage", False):
        dissemination_op += P_SCORE_ACTIVE_PROPHAGE
    if markers.get("flagellar_motility_complete", False):
        dissemination_op += P_SCORE_MOTILITY
    dissemination_op = min(5, dissemination_op)
    operational_score += dissemination_op
    details["operational_dissemination_score"] = dissemination_op
    details["flagellar_motility_complete"] = markers.get("flagellar_motility_complete", False)

    details["operational_persistence_score"] = operational_score
    genetic_score += operational_score

    # Cap Level 1 at 75 (Paper v2.4: max natural 75 pts)
    genetic_score = min(genetic_score, P_SCORE_GENETIC_CAP)
    details["genetic_score_raw"] = genetic_score
    details["genetic_score_capped"] = genetic_score

    # =========================================================================
    # LEVEL 2: Temporal / Spatial (if metadata available)
    # =========================================================================
    temporal_spatial_score = 0

    if metadata is not None:
        details["has_metadata"] = True

        # Metric 1: Colonisation timespan
        timespan_weeks = metadata.get("timespan_weeks", 0)
        timespan_score = 0
        if timespan_weeks > 52:
            timespan_score = P_SCORE_TIMESPAN_GT52W
        elif timespan_weeks >= 13:
            timespan_score = P_SCORE_TIMESPAN_13_52W
        elif timespan_weeks >= 4:
            timespan_score = P_SCORE_TIMESPAN_4_12W

        temporal_spatial_score += timespan_score
        details["timespan_weeks"] = timespan_weeks
        details["timespan_score"] = timespan_score

        # Metric 2: Independent detections
        detections = metadata.get("independent_detections", 0)
        detection_score = 0
        if detections >= 5:
            detection_score = P_SCORE_DETECTIONS_GTE5
        elif detections >= 3:
            detection_score = P_SCORE_DETECTIONS_GTE3

        temporal_spatial_score += detection_score
        details["independent_detections"] = detections
        details["detection_score"] = detection_score

        # Metric 3: Spatial Dispersion Score (SDS)
        unique_zones = metadata.get("unique_zones", 0)
        total_zones = metadata.get("total_zones_sampled", 1)
        sds = unique_zones / total_zones if total_zones > 0 else 0

        sds_score = 0
        if sds > 0.5:
            sds_score = P_SCORE_SDS_GT05
        elif sds >= 0.2:
            sds_score = P_SCORE_SDS_GTE02

        temporal_spatial_score += sds_score
        details["unique_zones"] = unique_zones
        details["total_zones_sampled"] = total_zones
        details["spatial_dispersion_score"] = round(sds, 2)
        details["sds_score"] = sds_score
        details["temporal_spatial_score"] = temporal_spatial_score

        # Combine levels (WITH metadata) and normalise by P_NORM_MAX
        total_score = genetic_score + temporal_spatial_score
        total_score = min(100.0, max(0.0, total_score / P_NORM_MAX * 100))

    else:
        # =================================================================
        # NO METADATA — BAYESIAN IMPUTATION (recalibrated v2.3)
        # =================================================================
        details["has_metadata"] = False
        details["imputation_method"] = "bayesian"
        details["imputation_factor"] = P_SCORE_BAYESIAN_IMPUTATION_FACTOR

        # expected_temporal = genetic * factor
        expected_temporal = (
            (genetic_score / P_SCORE_GENETIC_CAP)
            * (P_SCORE_GENETIC_CAP * P_SCORE_BAYESIAN_IMPUTATION_FACTOR)
        )

        details["expected_temporal_score"] = round(expected_temporal, 1)
        details["temporal_spatial_score"] = round(expected_temporal, 1)

        total_score = genetic_score + expected_temporal
        total_score = min(100.0, max(0.0, total_score / P_NORM_MAX * 100))

        details["imputation_note"] = (
            f"Bayesian imputation (factor={P_SCORE_BAYESIAN_IMPUTATION_FACTOR}): "
            f"genetic={genetic_score}, "
            f"expected_temporal={expected_temporal:.1f}, total={total_score:.1f}"
        )

    total_score = round(total_score, 1)
    details["total_score"] = total_score

    return total_score, details
