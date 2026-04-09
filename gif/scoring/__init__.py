# GIF Framework — Genomic Intelligence Framework
# Copyright 2026 Javier Gamboa
# Licensed under the Apache License, Version 2.0
"""
GIF Scoring Module

Exports the four component scoring functions and the integrator.
"""

from gif.scoring.weights import (
    CONTEXT_PROFILES,
    DEFAULT_CONTEXT,
    P_NORM_MAX,
    RISK_CRITICAL_THRESHOLD,
    RISK_HIGH_THRESHOLD,
    RISK_MODERATE_THRESHOLD,
    classify_risk_level,
)
from gif.scoring.v_score import calculate_v_score
from gif.scoring.p_score import calculate_p_score
from gif.scoring.c_score import calculate_c_score
from gif.scoring.r_score import calculate_r_score
from gif.scoring.integrator import calculate_gif_score, classify_trophic_strategy

__all__ = [
    "CONTEXT_PROFILES",
    "DEFAULT_CONTEXT",
    "P_NORM_MAX",
    "RISK_CRITICAL_THRESHOLD",
    "RISK_HIGH_THRESHOLD",
    "RISK_MODERATE_THRESHOLD",
    "classify_risk_level",
    "calculate_v_score",
    "calculate_p_score",
    "calculate_c_score",
    "calculate_r_score",
    "calculate_gif_score",
    "classify_trophic_strategy",
]
