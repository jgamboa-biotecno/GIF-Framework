# GIF Framework - Quick Start Guide

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

### Conda environment for external tools

```bash
conda create -n gif python=3.10
conda activate gif
conda install -c bioconda blast chewbbaca mlst abricate
pip install -e .
```

## External tools

GIF requires the following bioinformatics tools in `PATH`:

```bash
conda install -c bioconda blast chewbbaca mlst abricate
```

## First run

```bash
# Score a single genome
gif score my_genome.fasta

# Score with output directory
gif score my_genome.fasta --output results/

# Score all FASTAs in a directory
gif score genomes/ --output results/

# Check version
gif --version
```

## Output formats

```bash
gif score genome.fasta --format tsv      # Tab-separated table (default)
gif score genome.fasta --format json     # Structured JSON
gif score genome.fasta --format report   # Markdown per-isolate report
gif score genome.fasta --format all      # All formats
```

## Scoring contexts

```bash
# Industrial (default) - optimized for food processing surveillance
gif score genome.fasta --context industrial

# Clinical - optimized for patient risk and outbreak tracing
gif score genome.fasta --context clinical
```

## NCBI download

```bash
# From SRA accessions
gif score --accession SRR12345678

# From BioProject
gif score --bioproject PRJNA689484 --output results/

# From assembly accessions
gif score --accession GCF_000196035.1
```
