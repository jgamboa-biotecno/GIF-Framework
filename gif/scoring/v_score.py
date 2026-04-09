# GIF Framework — Genomic Intelligence Framework
# Copyright 2026 Javier Gamboa
# Licensed under the Apache License, Version 2.0
"""
V-Score (Virulence) calculation for the GIF Framework.

Scoring levels:
  Level 1 — CC Classification (0-40 pts)
  Level 2 — HC Proximity to clinical isolates (0-10 pts)
  Level 3 — LIPI Islands (-10 to +35 pts)
  Level 4 — InlA Status (-40 to +20 pts)

Raw range: [-45, +105] normalised linearly to [0, 100].
"""

from typing import Dict, Optional, Tuple

# =============================================================================
# V-SCORE CONSTANTS
# =============================================================================

# Level 1: CC Classification (0-40 pts) — Scientifically validated
# Source: PMC8764371 — Virulence Stratification in L. monocytogenes
# Hypervirulent: High clinical frequency, CNS/MN tropism, LIPI-3/4
#   CC1, CC2, CC4: Maury 2016, Genome Medicine 2024
#   CC87: LIPI-4+, dominant hypervirulent clone in Asia (PMC7290528)
#   CC14: Hypervirulent in Galleria model, MN infections (Frontiers Microbiol 2024)
# Intermediate: CC6 reclassified per Genome Medicine 2024 (ST6 = "average virulence")
CC_HYPERVIRULENT = ["CC1", "CC2", "CC4", "CC87", "CC14"]   # 40 points
CC_INTERMEDIATE = ["CC3", "CC5", "CC6", "CC7", "CC11"]     # 20 points
# All others: 5 points — Hypovirulent (CC9, CC121, etc.)

# Normalisation range
V_SCORE_MAX_RAW = 105   # Max: CC(40) + HC(10) + LIPI-3(15) + LIPI-4(20) + InlA(20) = 105
V_SCORE_MIN_RAW = -45   # Min: CC_hypo(5) + LIPI-1_deleted(-10) + InlA_severely_truncated(-40) = -45

# HC proximity scoring thresholds (allelic distance to nearest clinical isolate)
# Scoring identical to phiercc_service.py
HC_PROXIMITY_THRESHOLDS = [
    (7, 10),     # AD_min <= 7:   epidemic cluster — virtually identical to clinical cases
    (20, 7),     # AD_min 8-20:   related sublineage — moderate divergence, virulence preserved
    (50, 4),     # AD_min 21-50:  same CC, substantial divergence — intermediate clinical risk
    (100, 2),    # AD_min 51-100: endemic circulation — reduced potential virulence
    (150, 1),    # AD_min 101-150: same CC, highly divergent — expected low virulence
]
# AD_min > 150: 0 pts — different CCs, no recent evolutionary relationship

# LIPI-4 context-valid clonal complexes (cellobiose PTS, neuro/placental tropism)
LIPI4_VALID_CCS = ["CC4", "CC87", "CC382", "CC619"]

# Lineage-to-CC mapping for LIPI-3 context validation
# In the CLI context, lineage is derived from the markers dict.
# LIPI-3 is only valid in Lineage I.
LINEAGE_I_LABELS = ["Lineage_I", "I"]


def _score_hc_proximity(cgmlst_closest_distance: Optional[int]) -> Tuple[int, Dict]:
    """
    Score HC proximity to clinical isolates (0-10 pts).

    In CLI mode this uses a pre-computed ``cgmlst_closest_distance`` value
    rather than running a live allelic profile comparison.

    Args:
        cgmlst_closest_distance: Minimum allelic distance to a clinical
            isolate in the reference database.  ``None`` means unknown.

    Returns:
        (score, details_dict)
    """
    if cgmlst_closest_distance is None:
        return 0, {
            "method": "hc_proximity_v2",
            "success": False,
            "reason": "No cgMLST distance provided",
        }

    distance = cgmlst_closest_distance
    score = 0
    category = "distant"

    for threshold, pts in HC_PROXIMITY_THRESHOLDS:
        if distance <= threshold:
            score = pts
            category = f"AD_min<={threshold}"
            break

    return score, {
        "method": "hc_proximity_v2",
        "success": True,
        "distance": distance,
        "score": score,
        "category": category,
    }


def calculate_v_score(markers: dict) -> Tuple[float, Dict]:
    """
    Calculate V-Score (Virulence) according to the GIF specification.

    Args:
        markers: Dict with genomic marker fields.  Required keys:
            - clonal_complex (str|None)
            - lineage (str|None)
            - lipi1_intact (bool|None)
            - lipi3_complete (bool)
            - lipi3_genes_detected (int)
            - lipi4_complete (bool)
            - lipi4_genes_detected (int)
            - inla_status (str): "complete", "truncated",
              "severely_truncated", "not_found", or "unknown"
            - inla_protein_length (int|None)
            - cgmlst_closest_distance (int|None)

    Returns:
        (normalised_score 0-100, details dict)
    """
    score = 0
    details: Dict = {}

    cc = markers.get("clonal_complex")
    lineage = markers.get("lineage")

    # =========================================================================
    # LEVEL 1: CC Classification (0-40 points)
    # =========================================================================
    if cc in CC_HYPERVIRULENT:
        cc_score = 40
        cc_category = "hypervirulent"
    elif cc in CC_INTERMEDIATE:
        cc_score = 20
        cc_category = "intermediate"
    else:
        cc_score = 5
        cc_category = "hypovirulent"

    score += cc_score
    details["level1_cc_score"] = cc_score
    details["level1_cc_category"] = cc_category
    details["cc"] = cc
    details["lineage"] = lineage or ""

    # =========================================================================
    # LEVEL 2 (Metrica 1.2): HC Proximity to clinical isolates (0-10 pts)
    # =========================================================================
    hc_proximity_score, hc_proximity_details = _score_hc_proximity(
        markers.get("cgmlst_closest_distance")
    )
    score += hc_proximity_score
    details["level2_hc_proximity_score"] = hc_proximity_score
    details["level2_hc_proximity_details"] = hc_proximity_details

    # =========================================================================
    # LEVEL 3: LIPI Islands (penalisation / bonus)
    # =========================================================================

    # LIPI-1: Universal; only penalise if explicitly detected as deleted (-10 pts)
    lipi1_intact = markers.get("lipi1_intact")
    lipi1_penalty = 0
    if lipi1_intact is False:  # Explicitly False, not None
        lipi1_penalty = -10
        score += lipi1_penalty

    details["lipi1_intact"] = lipi1_intact
    details["lipi1_penalty"] = lipi1_penalty

    # LIPI-3: +15 points ONLY if Lineage I
    # (Cotter et al. 2008; Maury et al. 2016)
    lipi3_complete = markers.get("lipi3_complete", False)
    lipi3_context_valid = lineage in LINEAGE_I_LABELS
    lipi3_score = 0
    if lipi3_complete and lipi3_context_valid:
        lipi3_score = 15
        score += lipi3_score

    details["lipi3_complete"] = lipi3_complete
    details["lipi3_genes_detected"] = markers.get("lipi3_genes_detected", 0)
    details["lipi3_context_valid"] = lipi3_context_valid
    details["lipi3_score"] = lipi3_score

    # LIPI-4: +20 points if CC carries LIPI-4 (cellobiose PTS, neuro/placental tropism)
    # CC4: Maury 2016; CC87: PMC7290528; CC382/CC619: Frontiers Mol Biosci 2023
    lipi4_complete = markers.get("lipi4_complete", False)
    lipi4_context_valid = cc in LIPI4_VALID_CCS
    lipi4_score = 0
    if lipi4_complete and lipi4_context_valid:
        lipi4_score = 20
        score += lipi4_score

    details["lipi4_complete"] = lipi4_complete
    details["lipi4_genes_detected"] = markers.get("lipi4_genes_detected", 0)
    details["lipi4_context_valid"] = lipi4_context_valid
    details["lipi4_score"] = lipi4_score

    # =========================================================================
    # LEVEL 4: InlA Status (-40 to +20 points)
    # =========================================================================
    inla_status = markers.get("inla_status", "unknown")
    inla_score = 0

    if inla_status == "complete":
        inla_score = 20
    elif inla_status == "truncated":
        inla_score = -15   # Paper v2.4: partial attenuation (Nightingale et al. 2005)
    elif inla_status == "severely_truncated":
        inla_score = -40
    elif inla_status == "not_found":
        inla_score = 0     # Paper v2.4: neutral — do not penalise missing data
    # "unknown" -> 0

    score += inla_score

    details["inla_status"] = inla_status
    details["inla_protein_length"] = markers.get("inla_protein_length")
    details["inla_score"] = inla_score

    # =========================================================================
    # NORMALISATION: Raw score to [0, 100]
    # Formula: ((raw - min) / (max - min)) * 100
    # =========================================================================
    details["raw_score"] = score

    score_range = V_SCORE_MAX_RAW - V_SCORE_MIN_RAW  # 105 - (-45) = 150
    normalised_score = ((score - V_SCORE_MIN_RAW) / score_range) * 100.0
    normalised_score = min(100.0, max(0.0, normalised_score))
    normalised_score = round(normalised_score, 1)

    details["normalized_score"] = normalised_score

    return normalised_score, details
