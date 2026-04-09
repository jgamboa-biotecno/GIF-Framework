# GIF Framework — Genomic Intelligence Framework
# Copyright 2026 Javier Gamboa
# Licensed under the Apache License, Version 2.0
"""
gif.utils — Utility functions for GIF-CLI.

Provides input validation, FASTA discovery, and metadata parsing.
"""

from gif.utils.validators import validate_fasta, find_fastas, parse_metadata_tsv

__all__ = [
    "validate_fasta",
    "find_fastas",
    "parse_metadata_tsv",
]
