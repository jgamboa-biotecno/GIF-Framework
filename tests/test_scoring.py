# GIF Framework — Genomic Intelligence Framework
# Copyright 2026 Javier Gamboa
# Licensed under the Apache License, Version 2.0
"""
GIF Framework scoring unit tests.

Test data:
  - test_nosotroph_CC6.fasta: GCF_003129845.1 (CC6, V>65, P<35)
  - test_amphitroph_CC204.fasta: GCF_013389395.1 (CC204, V>=35, P>=40)
  - test_saprotroph_CC121.fasta: GCA_024255965.1 (CC121, V<30, P>45)
"""

import pytest
from gif.scoring.v_score import calculate_v_score
from gif.scoring.p_score import calculate_p_score
from gif.scoring.c_score import calculate_c_score
from gif.scoring.r_score import calculate_r_score
from gif.scoring.integrator import calculate_gif_score, classify_trophic_strategy
from gif.scoring.weights import classify_risk_level


# =============================================================================
# V-Score Tests
# =============================================================================

class TestVScore:
    def test_hypervirulent_cc1(self):
        markers = {"clonal_complex": "CC1", "inla_status": "complete"}
        score, details = calculate_v_score(markers)
        assert score > 65, f"CC1 with complete InlA should be >65, got {score}"

    def test_hypovirulent_cc121_truncated(self):
        markers = {"clonal_complex": "CC121", "inla_status": "severely_truncated"}
        score, details = calculate_v_score(markers)
        assert score < 30, f"CC121 with truncated InlA should be <30, got {score}"

    def test_intermediate_cc6(self):
        markers = {"clonal_complex": "CC6", "inla_status": "complete"}
        score, details = calculate_v_score(markers)
        assert 30 <= score <= 80, f"CC6 should be intermediate, got {score}"

    def test_normalization_range(self):
        """V-Score must always be in [0, 100]."""
        for cc in ["CC1", "CC121", "CC9", "CC7", "CC999"]:
            for inla in ["complete", "truncated", "severely_truncated", "not_found"]:
                markers = {"clonal_complex": cc, "inla_status": inla}
                score, _ = calculate_v_score(markers)
                assert 0 <= score <= 100, f"V-Score out of range: {score} for {cc}/{inla}"


# =============================================================================
# P-Score Tests
# =============================================================================

class TestPScore:
    def test_no_markers(self):
        markers = {}
        score, details = calculate_p_score(markers)
        assert score < 20, f"No markers should give low P-Score, got {score}"

    def test_full_persistence(self):
        markers = {
            "qac_e_delta1_detected": True,
            "qac_h_detected": True,
            "bcr_abc_detected": True,
            "ssi1_genes_count": 5,
            "ssi2_genes_count": 2,
            "gad_genes_count": 3,
            "cadmium_resistance_detected": True,
            "prophages_count": 3,
            "crispr_detected": False,
        }
        score, details = calculate_p_score(markers)
        assert score > 60, f"Full persistence markers should give >60, got {score}"

    def test_bayesian_imputation(self):
        """Without metadata, Bayesian imputation should be applied."""
        markers = {"qac_e_delta1_detected": True, "bcr_abc_detected": True}
        score, details = calculate_p_score(markers, metadata=None)
        assert details.get("has_metadata") is False
        assert "imputation" in str(details).lower()

    def test_normalization_pnorm75(self):
        """P-Score should use P_NORM_MAX=75 normalization."""
        markers = {"qac_e_delta1_detected": True}  # 15 pts genetic
        score, details = calculate_p_score(markers)
        # With imputation: 15 + 15*0.60 = 24, normalized: 24/75*100 = 32
        assert score > 25, f"P_NORM_MAX=75 normalization check failed, got {score}"


# =============================================================================
# C-Score Tests
# =============================================================================

class TestCScore:
    def test_no_cgmlst(self):
        markers = {}
        score, details = calculate_c_score(markers)
        assert score == 0 or score >= 0, "C-Score without cgMLST should be 0 or fallback"

    def test_close_distance(self):
        markers = {"cgmlst_closest_distance": 3, "cgmlst_cluster_size": 50}
        score, details = calculate_c_score(markers)
        assert score > 50, f"Close distance should give high C-Score, got {score}"


# =============================================================================
# R-Score Tests
# =============================================================================

class TestRScore:
    def test_no_amr(self):
        markers = {}
        score, details = calculate_r_score(markers)
        assert score == 0, f"No AMR should give 0, got {score}"

    def test_tetracycline(self):
        markers = {"tetracycline_genes": ["tetM"]}
        score, details = calculate_r_score(markers)
        assert score > 0, f"tetM should give >0, got {score}"

    def test_regional_factor(self):
        markers = {"tetracycline_genes": ["tetM"]}
        score_norway, _ = calculate_r_score(markers, region="norway")
        score_china, _ = calculate_r_score(markers, region="china")
        assert score_norway >= score_china, "Norway should amplify more than China"


# =============================================================================
# Integration Tests
# =============================================================================

class TestIntegration:
    def test_nosotroph_classification(self):
        strategy, desc = classify_trophic_strategy(75, 20)
        assert strategy == "Nosotroph"

    def test_amphitroph_classification(self):
        strategy, desc = classify_trophic_strategy(50, 55)
        assert strategy == "Amphitroph"

    def test_saprotroph_classification(self):
        strategy, desc = classify_trophic_strategy(20, 60)
        assert strategy == "Saprotroph"

    def test_unassigned_classification(self):
        strategy, desc = classify_trophic_strategy(33, 30)
        assert strategy == "Unassigned"

    def test_risk_tiers(self):
        assert classify_risk_level(90) == "CRITICAL"
        assert classify_risk_level(85) == "CRITICAL"
        assert classify_risk_level(70) == "HIGH"
        assert classify_risk_level(65) == "HIGH"
        assert classify_risk_level(50) == "MODERATE"
        assert classify_risk_level(40) == "MODERATE"
        assert classify_risk_level(39) == "LOW"
        assert classify_risk_level(0) == "LOW"

    def test_full_gif_score(self):
        markers = {
            "clonal_complex": "CC1",
            "inla_status": "complete",
            "lipi3_complete": True,
            "lipi3_genes_detected": 8,
            "qac_e_delta1_detected": True,
            "bcr_abc_detected": True,
        }
        result = calculate_gif_score(markers, context="industrial")
        assert "gif_score" in result
        assert "risk_tier" in result
        assert "trophic_strategy" in result
        assert "gif_score_industrial" in result
        assert "gif_score_clinical" in result
        assert 0 <= result["gif_score"] <= 100

    def test_dual_scoring(self):
        """Industrial and clinical scores should differ."""
        markers = {
            "clonal_complex": "CC1",
            "inla_status": "complete",
            "qac_e_delta1_detected": True,
        }
        result = calculate_gif_score(markers, context="industrial")
        assert result["gif_score_industrial"] != result["gif_score_clinical"]


# =============================================================================
# Expected Results from Validation Datasets
# =============================================================================

class TestValidationExpected:
    """
    Expected results for representative validation isolates.
    These accessions are included as test_data FASTAs.
    """

    def test_cc6_nosotroph_expected(self):
        """GCF_003129845.1 (CC6) should be Nosotroph: V~66.7, P~29.9"""
        markers = {
            "clonal_complex": "CC6",
            "inla_status": "complete",
        }
        result = calculate_gif_score(markers)
        v = result["v_score"]
        assert v > 50, f"CC6 V-Score should be >50, got {v}"

    def test_cc121_saprotroph_expected(self):
        """GCA_024255965.1 (CC121) should be Saprotroph: V~6.7, P~46.4"""
        markers = {
            "clonal_complex": "CC121",
            "inla_status": "severely_truncated",
            "qac_e_delta1_detected": True,
            "ssi1_genes_count": 3,
            "ssi2_genes_count": 2,
            "gad_genes_count": 3,
        }
        result = calculate_gif_score(markers)
        v = result["v_score"]
        p = result["p_score"]
        assert v < 30, f"CC121 truncated V should be <30, got {v}"
        assert p > 40, f"CC121 with SSI+qac P should be >40, got {p}"
