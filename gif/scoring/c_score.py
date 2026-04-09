# GIF Framework — Genomic Intelligence Framework
# Copyright 2026 Javier Gamboa
# Licensed under the Apache License, Version 2.0
"""
C-Score (Clonality) calculation for the GIF Framework.

C-Score v3.0 structure:
  Level 1   — HC-based Clonality (0-70 pts)
  Level 2.1 — Cluster Size (0-30 pts)
  Level 2.2 — Clinical Association Index / CAI (informational in CLI mode)

The C-Score requires pre-computed cgMLST values supplied via the markers dict:
  - cgmlst_closest_distance  (int|None)
  - cgmlst_cluster_size      (int, default 0)

When no cgMLST data is available, the function falls back to PubMLST ST-count
heuristics.
"""

from typing import Dict, Optional, Tuple

# =============================================================================
# C-SCORE CONSTANTS
# =============================================================================

# Level 1: HC-based Clonality (0-70 pts)
# Scientific references:
#   ECDC multi-country threshold (<= 4 AD)
#   Institut Pasteur CT (<= 7 AD)
#   Ruppitsch et al. 2015 (<= 10 AD)
HC_CLONALITY_THRESHOLDS = [
    (4, 70),     # AD <= 4:   ECDC multi-country outbreak
    (7, 60),     # AD 5-7:    Institut Pasteur CT — confirmed outbreak cluster
    (10, 50),    # AD 8-10:   Ruppitsch outbreak definition
    (20, 35),    # AD 11-20:  related sublineage
    (50, 20),    # AD 21-50:  same CC, substantial divergence
    (100, 10),   # AD 51-100: endemic circulation
    (150, 5),    # AD 101-150: highly divergent within CC
]
# AD > 150: 0 pts

# Level 2: Cluster Size (0-30 pts)
# Measures clonal EXPANSION (different concept from Level 1 proximity)
CLUSTER_SIZE_THRESHOLDS = [
    (100, 30),   # >100 isolates in cluster at HC10
    (50, 25),    # 51-100
    (20, 20),    # 21-50
    (10, 15),    # 11-20
    (5, 10),     # 6-10
    (2, 5),      # 2-5
]
# Singleton (1): 0 pts

# Legacy PubMLST fallback constants
C_SCORE_STS_GT50 = 70     # >50 STs in PubMLST
C_SCORE_STS_20_50 = 40    # 20-50 STs
C_SCORE_STS_1_19 = 10     # 1-19 STs
C_SCORE_STS_NOVEL = 5     # Novel ST

# PubMLST fallback ST counts (approximate, from literature)
_FALLBACK_ST_COUNTS = {
    "CC1": 120, "CC2": 80, "CC3": 45, "CC4": 90, "CC5": 35,
    "CC6": 75, "CC7": 40, "CC8": 55, "CC9": 60, "CC11": 30,
    "CC14": 25, "CC37": 20, "CC121": 50, "CC155": 35,
}


def _score_hc_clonality(distance: Optional[int]) -> Tuple[int, Dict]:
    """Score Level 1 HC-based Clonality (0-70 pts) from allelic distance."""
    if distance is None:
        return 0, {
            "method": "hc_clonality_v3",
            "success": False,
            "reason": "No cgMLST distance provided",
        }

    score = 0
    category = "distant"

    for threshold, pts in HC_CLONALITY_THRESHOLDS:
        if distance <= threshold:
            score = pts
            category = f"AD<={threshold}"
            break

    return score, {
        "method": "hc_clonality_v3",
        "success": True,
        "distance": distance,
        "score": score,
        "category": category,
    }


def _score_cluster_size(cluster_size: int) -> Tuple[int, Dict]:
    """Score Level 2 Cluster Size (0-30 pts)."""
    if cluster_size <= 1:
        return 0, {
            "method": "cluster_size_v3",
            "success": True,
            "cluster_size": cluster_size,
            "score": 0,
            "category": "singleton",
        }

    score = 0
    category = "small"

    for threshold, pts in CLUSTER_SIZE_THRESHOLDS:
        if cluster_size > threshold:
            score = pts
            category = f">{threshold}"
            break

    # Handle the 2-5 range explicitly (smallest threshold is 2)
    if score == 0 and cluster_size >= 2:
        score = 5
        category = ">=2"

    return score, {
        "method": "cluster_size_v3",
        "success": True,
        "cluster_size": cluster_size,
        "score": score,
        "category": category,
    }


def _fallback_pubmlst(cc: Optional[str]) -> Tuple[int, Dict]:
    """Fallback scoring based on PubMLST ST counts when no cgMLST data."""
    if not cc:
        return 0, {
            "method": "pubmlst_fallback",
            "success": False,
            "reason": "No clonal complex provided",
        }

    st_count = _FALLBACK_ST_COUNTS.get(cc, 10)

    if st_count > 50:
        score = C_SCORE_STS_GT50
        category = ">50_STs"
    elif st_count >= 20:
        score = C_SCORE_STS_20_50
        category = "20-50_STs"
    elif st_count >= 1:
        score = C_SCORE_STS_1_19
        category = "1-19_STs"
    else:
        score = C_SCORE_STS_NOVEL
        category = "novel_ST"

    return score, {
        "method": "pubmlst_fallback",
        "success": True,
        "clonal_complex": cc,
        "st_count": st_count,
        "score": score,
        "category": category,
    }


def calculate_c_score(markers: dict) -> Tuple[float, Dict]:
    """
    Calculate C-Score (Clonality) according to the GIF specification v3.0.

    In CLI mode, pre-computed values are expected in the markers dict:
      - cgmlst_closest_distance (int|None): minimum allelic distance to
        any isolate in the reference database
      - cgmlst_cluster_size (int): number of isolates within the tightest
        HC cluster (e.g. HC10)

    When neither value is available, the function falls back to an approximate
    score based on PubMLST ST counts for the given clonal complex.

    Args:
        markers: Dict with genomic marker fields.

    Returns:
        (score 0-100, details dict)
    """
    score = 0
    details: Dict = {}

    cc = markers.get("clonal_complex")
    details["clonal_complex"] = cc

    cgmlst_distance = markers.get("cgmlst_closest_distance")
    cgmlst_cluster_size = markers.get("cgmlst_cluster_size", 0)

    has_cgmlst = cgmlst_distance is not None
    details["cgmlst_data_available"] = has_cgmlst

    if has_cgmlst:
        # =================================================================
        # Level 1: HC-based Clonality (0-70 pts)
        # =================================================================
        hc_score, hc_details = _score_hc_clonality(cgmlst_distance)
        score += hc_score
        details["level1_hc_clonality_score"] = hc_score
        details["level1_hc_clonality_details"] = hc_details

        # =================================================================
        # Level 2: Cluster Size (0-30 pts)
        # =================================================================
        cs_score, cs_details = _score_cluster_size(cgmlst_cluster_size)
        score += cs_score
        details["level2_cluster_size_score"] = cs_score
        details["level2_cluster_size_details"] = cs_details
    else:
        # =================================================================
        # Fallback: PubMLST ST-count heuristic
        # =================================================================
        fb_score, fb_details = _fallback_pubmlst(cc)
        score += fb_score
        details["fallback_pubmlst_score"] = fb_score
        details["fallback_pubmlst_details"] = fb_details

    # Clamp to [0, 100]
    score = min(100, max(0, score))
    details["total_score"] = score

    return float(score), details
