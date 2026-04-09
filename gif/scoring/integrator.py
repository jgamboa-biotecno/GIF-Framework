# GIF Framework — Genomic Intelligence Framework
# Copyright 2026 Javier Gamboa
# Licensed under the Apache License, Version 2.0
"""
GIF Score integrator and trophic strategy classification.

Combines V-Score, P-Score, C-Score, and R-Score using context-dependent
weights, then classifies the isolate's trophic strategy.
"""

from typing import Dict, Optional, Tuple

from gif.scoring.weights import (
    CONTEXT_PROFILES,
    DEFAULT_CONTEXT,
    classify_risk_level,
)
from gif.scoring.v_score import calculate_v_score
from gif.scoring.p_score import calculate_p_score
from gif.scoring.c_score import calculate_c_score
from gif.scoring.r_score import calculate_r_score

# =============================================================================
# TROPHIC STRATEGY CLASSIFICATION
# =============================================================================

TROPHIC_DESCRIPTIONS = {
    "Nosotroph": "Pathogenic niche specialist (high V, low P)",
    "Saprotroph": "Industrial niche specialist (low V, high P)",
    "Amphitroph": "Dual-niche (functional V + environmental P)",
    "Unassigned": "Intermediate profile, no defined trophic strategy",
}


def classify_trophic_strategy(
    v_score: float,
    p_score: float,
) -> Tuple[str, str]:
    """
    Classify trophic strategy based on V-Score and P-Score.

    Classification is PHENOTYPIC (by GIF scores), NOT by clonal complex.
    A CC121 with complete InlA can be amphitrophic; a CC5 with truncated
    InlA can be saprotrophic.

    Args:
        v_score: Normalised V-Score (0-100).
        p_score: Normalised P-Score (0-100).

    Returns:
        (strategy_name, human-readable description)
    """
    if v_score > 65 and p_score < 35:
        strategy = "Nosotroph"
    elif v_score < 30 and p_score > 45:
        strategy = "Saprotroph"
    elif v_score >= 35 and p_score >= 40:
        strategy = "Amphitroph"
    else:
        strategy = "Unassigned"

    return strategy, TROPHIC_DESCRIPTIONS[strategy]


# =============================================================================
# MAIN INTEGRATOR
# =============================================================================

def calculate_gif_score(
    markers: dict,
    metadata: Optional[dict] = None,
    region: Optional[str] = None,
    context: str = "industrial",
) -> dict:
    """
    Calculate the complete GIF Score according to the official specification.

    Default (industrial): GIF = V*0.30 + P*0.40 + C*0.20 + R*0.10
    Clinical:             GIF = V*0.40 + P*0.20 + C*0.30 + R*0.10

    Risk Tiers (v3.1.0):
      CRITICAL >= 85
      HIGH     >= 65
      MODERATE >= 40
      LOW      <  40

    Args:
        markers: Dict with all genomic marker fields (see individual score
            modules for required keys).
        metadata: Optional dict with temporal/spatial persistence data
            (keys: timespan_weeks, independent_detections, unique_zones,
            total_zones_sampled).
        region: Optional region/country for R-Score AMR prevalence
            adjustment.
        context: Scoring context — ``"industrial"``, ``"clinical"``, or
            ``"custom"``.

    Returns:
        Dict with complete scoring breakdown::

            {
                "gif_score": float,
                "risk_tier": str,
                "v_score": float,
                "p_score": float,
                "c_score": float,
                "r_score": float,
                "virulence_details": dict,
                "persistence_details": dict,
                "context_details": dict,
                "resistance_details": dict,
                "context_profile": str,
                "weights_applied": dict,
                "gif_score_industrial": float,
                "gif_score_clinical": float,
                "trophic_strategy": str,
                "trophic_description": str,
                "is_partial": bool,
                "failed_components": list,
                "component_errors": dict,
            }
    """
    failed_components = []
    component_errors: Dict[str, str] = {}

    # -- V-Score ---------------------------------------------------------------
    try:
        v_score, v_details = calculate_v_score(markers)
    except Exception as exc:
        v_score, v_details = 0.0, {"error": str(exc), "failed": True}
        failed_components.append("v_score")
        component_errors["v_score"] = str(exc)

    # -- P-Score ---------------------------------------------------------------
    try:
        p_score, p_details = calculate_p_score(markers, metadata)
    except Exception as exc:
        p_score, p_details = 0.0, {"error": str(exc), "failed": True}
        failed_components.append("p_score")
        component_errors["p_score"] = str(exc)

    # -- C-Score ---------------------------------------------------------------
    try:
        c_score, c_details = calculate_c_score(markers)
    except Exception as exc:
        c_score, c_details = 0.0, {"error": str(exc), "failed": True}
        failed_components.append("c_score")
        component_errors["c_score"] = str(exc)

    # -- R-Score ---------------------------------------------------------------
    try:
        r_score, r_details = calculate_r_score(markers, region)
    except Exception as exc:
        r_score, r_details = 0.0, {"error": str(exc), "failed": True}
        failed_components.append("r_score")
        component_errors["r_score"] = str(exc)

    is_partial = len(failed_components) > 0

    # -- Resolve context weights -----------------------------------------------
    active_context = context if context in CONTEXT_PROFILES else DEFAULT_CONTEXT
    weights = CONTEXT_PROFILES[active_context]["weights"]

    gif_score = round(
        v_score * weights["v"]
        + p_score * weights["p"]
        + c_score * weights["c"]
        + r_score * weights["r"],
        1,
    )

    # Dual scores (always compute both standard profiles)
    ind_w = CONTEXT_PROFILES["industrial"]["weights"]
    cli_w = CONTEXT_PROFILES["clinical"]["weights"]

    gif_score_industrial = round(
        v_score * ind_w["v"] + p_score * ind_w["p"]
        + c_score * ind_w["c"] + r_score * ind_w["r"],
        1,
    )
    gif_score_clinical = round(
        v_score * cli_w["v"] + p_score * cli_w["p"]
        + c_score * cli_w["c"] + r_score * cli_w["r"],
        1,
    )

    # Risk tier
    risk_tier = classify_risk_level(gif_score)

    # Trophic strategy
    trophic_strategy, trophic_description = classify_trophic_strategy(
        v_score, p_score
    )

    return {
        "gif_score": gif_score,
        "risk_tier": risk_tier,
        "v_score": v_score,
        "p_score": p_score,
        "c_score": c_score,
        "r_score": r_score,
        "virulence_details": v_details,
        "persistence_details": p_details,
        "context_details": c_details,
        "resistance_details": r_details,
        "context_profile": active_context,
        "weights_applied": weights,
        "gif_score_industrial": gif_score_industrial,
        "gif_score_clinical": gif_score_clinical,
        "trophic_strategy": trophic_strategy,
        "trophic_description": trophic_description,
        "is_partial": is_partial,
        "failed_components": failed_components,
        "component_errors": component_errors,
    }
