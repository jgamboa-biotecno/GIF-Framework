# GIF Framework — Genomic Intelligence Framework
# Copyright 2026 Javier Gamboa
# Licensed under the Apache License, Version 2.0
"""
Context profiles, weight configuration, and risk tier classification.
"""

from typing import Any, Dict


# =============================================================================
# CONTEXT PROFILES
# Empirically validated against published external datasets.
# AUC values reported below are computed against the corresponding
# validation_dataset; see the preprints in CITATION.cff for the full
# validation methodology and discussion.
# =============================================================================

CONTEXT_PROFILES: Dict[str, Dict[str, Any]] = {
    "industrial": {
        "weights": {"v": 0.30, "p": 0.40, "c": 0.20, "r": 0.10},
        "description": "Optimised for surveillance in food-processing facilities",
        "validated_auc": 0.933,
        "validation_dataset": "Fagerlund et al. 2022 (n=513)",
    },
    "clinical": {
        "weights": {"v": 0.40, "p": 0.20, "c": 0.30, "r": 0.10},
        "description": "Optimised for clinical risk assessment and OneHealth traceability",
        "validated_auc": None,
    },
    "custom": {
        "weights": {"v": 0.25, "p": 0.25, "c": 0.25, "r": 0.25},
        "description": "User-defined custom weights",
        "validated_auc": None,
    },
}

DEFAULT_CONTEXT = "industrial"

# P-Score normalisation denominator (v3.1.0)
# Theoretical genetic max = 75; normalising by 75 instead of 125 gives better
# discrimination (AUC 0.933 vs 0.785 with P_NORM=100)
P_NORM_MAX = 75

# =============================================================================
# RISK TIER THRESHOLDS (v3.1.0 recalibrated)
# =============================================================================

RISK_CRITICAL_THRESHOLD = 85
RISK_HIGH_THRESHOLD = 65
RISK_MODERATE_THRESHOLD = 40


def classify_risk_level(gif_score: float) -> str:
    """
    Classify risk level based on GIF score using official thresholds.

    Args:
        gif_score: Final GIF score (0-100)

    Returns:
        Risk tier string: "CRITICAL", "HIGH", "MODERATE", or "LOW"
    """
    if gif_score >= RISK_CRITICAL_THRESHOLD:
        return "CRITICAL"
    elif gif_score >= RISK_HIGH_THRESHOLD:
        return "HIGH"
    elif gif_score >= RISK_MODERATE_THRESHOLD:
        return "MODERATE"
    else:
        return "LOW"
