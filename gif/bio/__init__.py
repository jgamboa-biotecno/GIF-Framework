# GIF Framework — Genomic Intelligence Framework
# Copyright 2026 Javier Gamboa
# Licensed under the Apache License, Version 2.0

"""
gif.bio — Bioinformatics tool wrappers for GIF-CLI.

Provides wrappers around external tools for genome assembly QC, MLST/cgMLST
typing, virulence detection, AMR profiling, persistence marker detection,
and NCBI data fetching.
"""

from gif.bio.assembly import assess_assembly, validate_species
from gif.bio.mlst import run_mlst, run_cgmlst, get_lineage_from_cc
from gif.bio.virulence import detect_virulence
from gif.bio.resistance import detect_amr
from gif.bio.biocides import detect_persistence_markers
from gif.bio.fetch import fetch_assembly, fetch_sra, fetch_bioproject, resolve_accession

__all__ = [
    "assess_assembly",
    "validate_species",
    "run_mlst",
    "run_cgmlst",
    "get_lineage_from_cc",
    "detect_virulence",
    "detect_amr",
    "detect_persistence_markers",
    "fetch_assembly",
    "fetch_sra",
    "fetch_bioproject",
    "resolve_accession",
]
