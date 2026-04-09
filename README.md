# GIF — Genomic Intelligence Framework

A multi-component genomic risk assessment tool for *Listeria monocytogenes*. GIF integrates virulence (V), persistence (P), clonality (C) and antimicrobial resistance (R) into a single 0–100 score and classifies each isolate by trophic strategy (nosotroph, amphitroph, saprotroph).

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)

## Status

Reference implementation accompanying the methodological and clinical-validation manuscripts listed in the [Citation](#citation) section. Both manuscripts are currently available as bioRxiv preprints.

## Installation

### From source

```bash
git clone https://github.com/jgamboa-biotecno/GIF-Framework.git
cd GIF-Framework
pip install -e .
```

### Docker

```bash
docker build -t gif-framework .
docker run --rm -v $(pwd):/data gif-framework score /data/genome.fasta --output /data/results/
```

## External dependencies

GIF requires the following bioinformatics tools to be available in `PATH`:

| Tool | Purpose | Install |
|------|---------|---------|
| BLAST+ | Sequence alignment | `conda install -c bioconda blast` |
| chewBBACA | cgMLST typing | `conda install -c bioconda chewbbaca` |
| mlst | 7-gene MLST | `conda install -c bioconda mlst` |
| ABRicate | AMR / virulence gene detection | `conda install -c bioconda abricate` |

## Usage

```bash
# Single genome
gif score genome.fasta

# Directory of FASTAs
gif score /data/genomes/ --output /data/results/

# From NCBI accessions
gif score --accession GCA_000026945.2 --output results/

# From a BioProject
gif score --bioproject PRJNA422580 --output results/

# Calibration context (default: industrial)
gif score genome.fasta --context clinical

# Output format
gif score genome.fasta --format tsv      # default
gif score genome.fasta --format json
gif score genome.fasta --format report   # markdown
```

Run `gif --help` for the complete option list.

## Scoring components

GIF computes four orthogonal components and integrates them into a weighted score on a 0–100 scale.

| Component | Industrial weight | Clinical weight | Description |
|-----------|------------------:|----------------:|-------------|
| **V** Virulence | 0.30 | 0.40 | Clonal complex classification, LIPI islands, InlA truncation status |
| **P** Persistence | 0.40 | 0.20 | Genetic persistence markers (SSI-1, SSI-2, BcrABC, qac) and temporal recurrence |
| **C** Clonality | 0.20 | 0.30 | cgMLST allelic distance to reference clusters and cluster expansion |
| **R** Resistance | 0.10 | 0.10 | AMR determinants and regional prevalence weighting |

## Trophic strategies

| Strategy | V-Score | P-Score | Interpretation |
|----------|---------|---------|----------------|
| Nosotroph | > 65 | < 35 | Pathogenic-niche specialist |
| Amphitroph | ≥ 35 | ≥ 40 | Dual-niche: virulence and persistence |
| Saprotroph | < 30 | > 45 | Industrial-niche specialist |
| Unassigned | other | other | Intermediate profile |

The methodological basis for the trophic classification is described in the first manuscript listed in the [Citation](#citation) section.

## Citation

If you use GIF in your work, please cite the following preprints:

- Gamboa J. (2026). *Amphitrophic Listeria monocytogenes: multi-dimensional genomic profiling reveals a third ecological strategy that challenges the virulence-persistence trade-off paradigm.* bioRxiv. https://doi.org/10.64898/2026.03.23.713700

- Gamboa J. (2026). *Amphitrophic Listeria monocytogenes causes one-third of invasive listeriosis yet remains undetected by clonal complex-based risk classification.* bioRxiv. https://doi.org/10.64898/2026.03.28.715028

Additional preprints describing pathogen-specific applications and the regulatory context:

- Gamboa J. (2026). *Beyond Outbreak Detection: Mandatory Genomic Surveillance in the EU as an Opportunity to Quantify Metal-Mediated Co-Selection of Antimicrobial Resistance from a One Health Perspective.* Preprints.org. https://doi.org/10.20944/preprints202604.0573.v1

- Gamboa J. (2026). *From Serogroup to Genome: The Unfinished Transition in EU Food Safety Criteria for Shiga Toxin-Producing Escherichia coli.* Preprints.org. https://doi.org/10.20944/preprints202604.0508.v1

A machine-readable citation file is provided in [CITATION.cff](CITATION.cff).

## License

Apache License 2.0. See [LICENSE](LICENSE) for the full text.
