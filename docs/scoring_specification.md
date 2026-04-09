# GIF Scoring Specification v3.1

## Overview

GIF Score = weighted sum of four components, normalized to 0-100.

**Industrial calibration (default):** GIF = V * 0.30 + P * 0.40 + C * 0.20 + R * 0.10
**Clinical calibration:** GIF = V * 0.40 + P * 0.20 + C * 0.30 + R * 0.10

---

## V-Score (Virulence) - 0 to 100

### Level 1: CC Phylogenetic Classification (0-40 pts)

| Category | CCs | Score |
|----------|-----|-------|
| Hypervirulent | CC1, CC2, CC4, CC87, CC14 | 40 |
| Intermediate | CC3, CC5, CC6, CC7, CC11 | 20 |
| Hypovirulent | All others | 5 |

### Level 2: LIPI Islands

| Marker | Score | Condition |
|--------|-------|-----------|
| LIPI-1 intact | 0 | All 6 genes present |
| LIPI-1 deleted | -10 | Missing genes |
| LIPI-3 complete | +15 | >= 6/8 genes |
| LIPI-4 complete | +20 | >= 4/6 genes AND CC4 or CC87 |

### Level 3: InlA Status

| Status | Score | Description |
|--------|-------|-------------|
| Complete | +20 | Full-length functional InlA |
| Truncated | -15 | Premature stop codon (partial function) |
| Severely truncated | -40 | Early truncation (loss of function) |
| Not found | 0 | No data (neutral) |

### HC Proximity (0-10 pts)

Based on minimum allelic distance to clinical isolates in reference database.

### Normalization

Raw score range: [-45, +105]
Normalized: ((raw - (-45)) / (105 - (-45))) * 100 = ((raw + 45) / 150) * 100

---

## P-Score (Persistence) - 0 to 100

### Level 1: Genetic Markers (0-75 pts, capped)

| Marker | Score | Description |
|--------|-------|-------------|
| qacEdelta1 | +15 | Quaternary ammonium resistance |
| qacH | +10 | Quaternary ammonium resistance |
| bcrABC | +20 | Benzalkonium chloride resistance operon |
| SSI (>=6/10 genes) | +15 | Stress Survival Islets |
| Cadmium resistance | +5 | cadA variants |
| Prophages > 2 | +5 | Attenuated prophages |
| CRISPR absent | +5 | Facilitates exogenous DNA acquisition |

### Level 1b: Operational Persistence (0-15 pts within cap)

- Biofilm module (0-5): bapL, agr, EPS cluster
- Survival module (0-5): acid tolerance, efflux, cold shock, oxidative stress
- Dissemination module (0-5): IS1216, LGI-1, active prophage, motility

### Level 2: Temporal/Spatial (0-50 pts, if metadata available)

| Factor | Condition | Score |
|--------|-----------|-------|
| Timespan | > 52 weeks | +30 |
| Timespan | 13-52 weeks | +20 |
| Timespan | 4-12 weeks | +10 |
| Detections | >= 5 | +10 |
| Detections | >= 3 | +5 |
| SDS | > 0.5 | +10 |
| SDS | >= 0.2 | +5 |

### Bayesian Imputation (no metadata)

When temporal metadata is unavailable:
- expected_temporal = genetic_score * 0.60
- total = genetic + expected_temporal

### Normalization

P_NORM_MAX = 75. Final score = min(100, total / 75 * 100)

---

## C-Score (Clonality) - 0 to 100

### Level 1: HC Clonality (0-70 pts)

Based on minimum allelic distance in cgMLST HC clustering.

### Level 2: Cluster Size (0-30 pts)

Based on number of isolates in HC10/HC20/HC50/HC100 clusters.

### Clinical Association Index (CAI)

Ratio of clinical to total isolates in the cluster.

---

## R-Score (Resistance) - 0 to 100

### Per-gene Scoring

| Category | Score per gene |
|----------|---------------|
| Beta-lactam | +10 |
| Tetracycline | +10 |
| Aminoglycoside | +10 |
| Fluoroquinolone (QRDR) | +15 |
| Macrolide | +10 |
| Phenicol | +10 |
| Lincosamide | +10 |
| Trimethoprim | +10 |

Genetic cap: 80 pts

### MDR Plasmid Bonus

+20 pts if plasmid >= 40kb with >= 3 AMR genes

### Regional Factor

Country-specific amplification factor based on AMR prevalence.
Formula: F = max(1.0, 5% / P_AMR), capped at 3.0x

---

## Risk Tiers

| Tier | Score Range |
|------|------------|
| CRITICAL | >= 85 |
| HIGH | 65 - 84 |
| MODERATE | 40 - 64 |
| LOW | < 40 |

---

## Trophic Classification

| Strategy | V-Score | P-Score | Description |
|----------|---------|---------|-------------|
| Nosotroph | > 65 | < 35 | Pathogenic niche specialist |
| Amphitroph | >= 35 | >= 40 | Dual-niche capability |
| Saprotroph | < 30 | > 45 | Industrial niche specialist |
| Unassigned | other | other | Intermediate profile |
