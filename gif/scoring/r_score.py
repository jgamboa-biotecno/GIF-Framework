# GIF Framework — Genomic Intelligence Framework
# Copyright 2026 Javier Gamboa
# Licensed under the Apache License, Version 2.0
"""
R-Score (Resistance) calculation for the GIF Framework.

Scoring levels:
  Level 1  — Genetic Determinants (cap at 80 pts)
  Level 1b — Biocide Resistance (0-15 pts within the 80 cap)
  Level 2  — MDR Plasmid bonus (+20 pts)

Regional amplification factor applied at the end.
Final score clamped to [0, 100].
"""

from typing import Any, Dict, Optional, Tuple

# =============================================================================
# R-SCORE CONSTANTS — Level 1: Genetic Determinants
# =============================================================================

R_SCORE_BETA_LACTAM_EACH = 10       # +10 per beta-lactam gene detected
R_SCORE_TETRACYCLINE_EACH = 10      # +10 per tet gene
R_SCORE_AMINOGLYCOSIDE = 10         # +10 binary (any detected)
R_SCORE_FLUOROQUINOLONE_EACH = 15   # +15 per gyrA/parC QRDR mutation
R_SCORE_MACROLIDE = 10              # ermB, mphB (Lancet Reg Health EU 2024)
R_SCORE_PHENICOL = 10               # fexA — Paper v2.4
R_SCORE_LINCOSAMIDE = 10            # lnuG — Paper v2.4
R_SCORE_TRIMETHOPRIM = 10           # dfrD — Paper v2.4
R_SCORE_GENETIC_CAP = 80

# Level 1b: Biocide Resistance (within the 80 cap)
R_SCORE_ACCESSORY_QAC_EFFLUX = 5    # qacC/emrC/emrE/sugE (NOT qacH/bcrABC)
R_SCORE_BC_EFFLUX = 5               # mdrL + lde co-detection
R_SCORE_RESISTANCE_PLASMID = 5      # repA_pLM80/cadA_plasmid/tetS markers

# Level 2: MDR Plasmid
R_SCORE_MULTIDRUG_PLASMID = 20      # >= 40 kb, >= 3 AMR genes

# =============================================================================
# REGIONAL AMR PREVALENCE FACTORS
# Formula: F_regional = min(MAX_FACTOR, max(1.0, 5% / P_AMR))
# Cap at 3.0x to avoid disproportionate geographic distortion
# =============================================================================

REGIONAL_AMR_MAX_FACTOR = 3.0

AMR_REGIONAL_FACTORS: Dict[str, Dict[str, Any]] = {
    # Europe — Low AMR prevalence
    "norway":         {"prevalence": 0.008, "factor": 3.0,  "source": "NORM 2023"},
    "finland":        {"prevalence": 0.012, "factor": 3.0,  "source": "Finres 2023"},
    "spain":          {"prevalence": 0.025, "factor": 2.0,  "source": "CNM 2022"},
    "france":         {"prevalence": 0.018, "factor": 2.78, "source": "CNR-Listeria 2023"},
    "germany":        {"prevalence": 0.032, "factor": 1.56, "source": "RKI 2022"},
    "italy":          {"prevalence": 0.040, "factor": 1.25, "source": "ISS 2022"},
    "uk":             {"prevalence": 0.038, "factor": 1.32, "source": "PHE 2023"},
    "united_kingdom": {"prevalence": 0.038, "factor": 1.32, "source": "PHE 2023"},

    # North America
    "usa":            {"prevalence": 0.042, "factor": 1.19, "source": "NARMS 2022"},
    "united_states":  {"prevalence": 0.042, "factor": 1.19, "source": "NARMS 2022"},
    "canada":         {"prevalence": 0.035, "factor": 1.43, "source": "CIPARS 2023"},

    # Asia — Higher AMR prevalence (no amplification)
    "china":          {"prevalence": 0.15,  "factor": 1.0,  "source": "Yan et al. 2025"},
    "india":          {"prevalence": 0.08,  "factor": 1.0,  "source": "Gupta et al. 2023"},
    "japan":          {"prevalence": 0.055, "factor": 1.0,  "source": "JANIS 2023"},

    # Regional defaults
    "europe":         {"prevalence": 0.030, "factor": 1.67, "source": "EUCAST regional"},
    "north_america":  {"prevalence": 0.040, "factor": 1.25, "source": "Regional average"},
    "asia":           {"prevalence": 0.08,  "factor": 1.0,  "source": "Regional average"},
    "south_america":  {"prevalence": 0.06,  "factor": 1.0,  "source": "WHO GLASS"},
    "africa":         {"prevalence": 0.07,  "factor": 1.0,  "source": "WHO GLASS"},
    "oceania":        {"prevalence": 0.035, "factor": 1.43, "source": "Regional estimate"},

    # Fallback
    "unknown":        {"prevalence": 0.05,  "factor": 1.0,  "source": "Conservative default"},
}

# Aliases for country-name normalisation
_REGION_ALIASES: Dict[str, str] = {
    "us": "usa",
    "united_states_of_america": "usa",
    "estados_unidos": "usa",
    "great_britain": "united_kingdom",
    "reino_unido": "united_kingdom",
    "españa": "spain",
    "deutschland": "germany",
    "alemania": "germany",
    "italia": "italy",
    "francia": "france",
    "noruega": "norway",
    "finlandia": "finland",
    "china_mainland": "china",
    "peoples_republic_of_china": "china",
}

# Country-to-continent mapping for fallback
_REGION_TO_CONTINENT: Dict[str, str] = {
    "sweden": "europe", "denmark": "europe", "netherlands": "europe",
    "belgium": "europe", "austria": "europe", "switzerland": "europe",
    "portugal": "europe", "greece": "europe", "czech_republic": "europe",
    "poland": "europe", "hungary": "europe", "croatia": "europe",
    "slovenia": "europe", "slovakia": "europe", "romania": "europe",
    "bulgaria": "europe", "lithuania": "europe", "latvia": "europe",
    "estonia": "europe", "ireland": "europe", "luxembourg": "europe",
    "malta": "europe", "cyprus": "europe",
    "mexico": "north_america", "guatemala": "north_america",
    "costa_rica": "north_america", "panama": "north_america",
    "brazil": "south_america", "argentina": "south_america",
    "chile": "south_america", "colombia": "south_america",
    "peru": "south_america", "venezuela": "south_america",
    "ecuador": "south_america", "uruguay": "south_america",
    "bolivia": "south_america", "paraguay": "south_america",
    "south_korea": "asia", "korea": "asia", "taiwan": "asia",
    "thailand": "asia", "vietnam": "asia", "indonesia": "asia",
    "malaysia": "asia", "singapore": "asia", "philippines": "asia",
    "bangladesh": "asia", "pakistan": "asia", "sri_lanka": "asia",
    "australia": "oceania", "new_zealand": "oceania",
    "south_africa": "africa", "egypt": "africa", "nigeria": "africa",
    "kenya": "africa", "morocco": "africa", "tunisia": "africa",
    "algeria": "africa", "ethiopia": "africa", "ghana": "africa",
    "tanzania": "africa",
}


def get_regional_amr_factor(region: Optional[str] = None) -> Dict[str, Any]:
    """
    Get regional AMR amplification factor for R-Score calculation.

    Args:
        region: Region code or country name (case-insensitive).

    Returns:
        Dict with ``factor``, ``prevalence``, and ``source``.
    """
    if not region:
        return AMR_REGIONAL_FACTORS["unknown"]

    norm = region.lower().strip().replace(" ", "_").replace("-", "_")

    # Direct lookup
    if norm in AMR_REGIONAL_FACTORS:
        return AMR_REGIONAL_FACTORS[norm]

    # Alias lookup
    if norm in _REGION_ALIASES:
        return AMR_REGIONAL_FACTORS[_REGION_ALIASES[norm]]

    # Continent fallback
    continent = _REGION_TO_CONTINENT.get(norm)
    if continent:
        info = AMR_REGIONAL_FACTORS[continent].copy()
        info["source"] = f"{info['source']} (country {region} mapped to {continent})"
        return info

    # Ultimate fallback
    fallback = AMR_REGIONAL_FACTORS["unknown"].copy()
    fallback["source"] = f"Unknown region '{region}' - using conservative default"
    return fallback


def calculate_r_score(
    markers: dict,
    region: Optional[str] = None,
) -> Tuple[float, Dict]:
    """
    Calculate R-Score (Resistance) according to the GIF specification.

    NOTE: qacH and bcrABC are counted in P-Score, NOT here.
    NOTE: fosX is INTRINSIC to L. monocytogenes — QC marker only, NOT scored.

    Args:
        markers: Dict with genomic marker fields.  Relevant keys:
            - beta_lactam_genes (list[str])
            - tetracycline_genes (list[str])
            - aminoglycoside_genes (list[str])
            - fluoroquinolone_mutations (list[str])
            - macrolide_genes (list[str])
            - phenicol_genes (list[str])
            - lincosamide_genes (list[str])
            - trimethoprim_genes (list[str])
            - biocide_accessory_efflux (list[str])
            - biocide_mdrL_detected (bool)  — optional, for BC efflux
            - biocide_lde_detected (bool)   — optional, for BC efflux
            - mobility_resistance_plasmid (bool) — optional
            - has_multidrug_plasmid (bool)
            - plasmid_amr_genes_count (int)

        region: Optional region/country for AMR prevalence adjustment.

    Returns:
        (score 0-100, details dict)
    """
    score = 0
    details: Dict = {}

    # =========================================================================
    # LEVEL 1: Genetic Determinants (cap at 80)
    # =========================================================================

    # Beta-lactamases: +10 each
    bl_genes = markers.get("beta_lactam_genes", [])
    bl_count = len(bl_genes)
    bl_score = bl_count * R_SCORE_BETA_LACTAM_EACH
    score += bl_score
    details["beta_lactam_genes"] = bl_genes
    details["beta_lactam_count"] = bl_count
    details["beta_lactam_score"] = bl_score

    # Tetracyclines: +10 each
    tet_genes = markers.get("tetracycline_genes", [])
    tet_count = len(tet_genes)
    tet_score = tet_count * R_SCORE_TETRACYCLINE_EACH
    score += tet_score
    details["tetracycline_genes"] = tet_genes
    details["tetracycline_count"] = tet_count
    details["tetracycline_score"] = tet_score

    # Aminoglycosides: +10 binary (any detected)
    amino_genes = markers.get("aminoglycoside_genes", [])
    amino_detected = len(amino_genes) > 0
    amino_score = R_SCORE_AMINOGLYCOSIDE if amino_detected else 0
    score += amino_score
    details["aminoglycoside_detected"] = amino_detected
    details["aminoglycoside_genes"] = amino_genes
    details["aminoglycoside_score"] = amino_score

    # Fluoroquinolones: +15 per QRDR mutation
    fq_mutations = markers.get("fluoroquinolone_mutations", [])
    fq_count = len(fq_mutations)
    fq_score = fq_count * R_SCORE_FLUOROQUINOLONE_EACH
    score += fq_score
    details["fluoroquinolone_mutations"] = fq_mutations
    details["fluoroquinolone_count"] = fq_count
    details["fluoroquinolone_score"] = fq_score

    # Emerging acquired AMR (Lancet Reg Health EU 2024)
    mac_genes = markers.get("macrolide_genes", [])
    mac_score = R_SCORE_MACROLIDE if mac_genes else 0
    score += mac_score
    details["macrolide_genes"] = mac_genes
    details["macrolide_score"] = mac_score

    phen_genes = markers.get("phenicol_genes", [])
    phen_score = R_SCORE_PHENICOL if phen_genes else 0
    score += phen_score
    details["phenicol_genes"] = phen_genes
    details["phenicol_score"] = phen_score

    linco_genes = markers.get("lincosamide_genes", [])
    linco_score = R_SCORE_LINCOSAMIDE if linco_genes else 0
    score += linco_score
    details["lincosamide_genes"] = linco_genes
    details["lincosamide_score"] = linco_score

    trim_genes = markers.get("trimethoprim_genes", [])
    trim_score = R_SCORE_TRIMETHOPRIM if trim_genes else 0
    score += trim_score
    details["trimethoprim_genes"] = trim_genes
    details["trimethoprim_score"] = trim_score

    # =========================================================================
    # LEVEL 1b: Biocide Resistance (0-15 pts within existing cap of 80)
    # =========================================================================
    biocide_r_score = 0

    # Acquired QAC efflux (qacC/emrC/emrE/sugE — NOT qacH/bcrABC)
    acc_efflux = markers.get("biocide_accessory_efflux", [])
    if len(acc_efflux) >= 1:
        biocide_r_score += R_SCORE_ACCESSORY_QAC_EFFLUX
        details["accessory_qac_efflux_genes"] = acc_efflux

    # BC efflux pumps (mdrL + lde co-detection)
    mdrL = markers.get("biocide_mdrL_detected", False)
    lde = markers.get("biocide_lde_detected", False)
    if mdrL and lde:
        biocide_r_score += R_SCORE_BC_EFFLUX

    # Resistance plasmid markers
    res_plasmid = markers.get("mobility_resistance_plasmid", False)
    if res_plasmid:
        biocide_r_score += R_SCORE_RESISTANCE_PLASMID

    details["biocide_resistance_score"] = biocide_r_score
    details["biocide_mdrL"] = mdrL
    details["biocide_lde"] = lde
    details["resistance_plasmid"] = res_plasmid
    score += biocide_r_score

    # Cap Level 1 at 80
    genetic_score = min(score, R_SCORE_GENETIC_CAP)
    details["genetic_score_raw"] = score
    details["genetic_score_capped"] = genetic_score
    score = genetic_score

    # =========================================================================
    # LEVEL 2: Multi-drug Resistance Plasmid (0-20 pts)
    # Criteria: plasmid >= 40 kb + >= 3 AMR genes + >= 3 different categories
    # =========================================================================
    plasmid_bonus = 0
    if markers.get("has_multidrug_plasmid", False):
        plasmid_bonus = R_SCORE_MULTIDRUG_PLASMID
        score += plasmid_bonus

    details["has_multidrug_plasmid"] = markers.get("has_multidrug_plasmid", False)
    details["plasmid_amr_genes_count"] = markers.get("plasmid_amr_genes_count", 0)
    details["plasmid_score"] = plasmid_bonus

    raw_score = score
    details["raw_score"] = raw_score

    # =========================================================================
    # REGIONAL AMPLIFICATION FACTOR
    # R_final = min(100, R_raw * F_regional)
    # =========================================================================
    regional_info = get_regional_amr_factor(region)
    regional_factor = min(regional_info["factor"], REGIONAL_AMR_MAX_FACTOR)
    regional_prevalence = regional_info["prevalence"]
    regional_source = regional_info["source"]

    if raw_score > 0:
        final_score = min(100.0, raw_score * regional_factor)
    else:
        final_score = 0.0

    details["regional_factor"] = regional_factor
    details["regional_prevalence"] = regional_prevalence
    details["regional_source"] = regional_source
    details["region_used"] = region or "unknown"
    details["amplified_score"] = final_score
    details["regional_adjustment_applied"] = (raw_score > 0 and regional_factor > 1.0)

    final_score = min(100.0, max(0.0, final_score))
    details["total_score"] = final_score

    return final_score, details
