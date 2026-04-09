# GIF Framework — Genomic Intelligence Framework
# Copyright 2026 Javier Gamboa
# Licensed under the Apache License, Version 2.0
"""
Input validation and file discovery utilities for GIF-CLI.

Provides FASTA validation, recursive FASTA file discovery, and metadata TSV
parsing for P-Score Level 2 temporal scoring.
"""

from __future__ import annotations

import csv
import gzip
import logging
import os
from pathlib import Path
from typing import Any, Dict, List

logger = logging.getLogger("gif.utils.validators")

# Recognised FASTA file extensions (case-insensitive)
FASTA_EXTENSIONS = {".fasta", ".fa", ".fna", ".fas", ".fasta.gz", ".fa.gz", ".fna.gz"}

# Expected metadata TSV columns for P-Score Level 2
METADATA_COLUMNS = {
    "sample_id",
    "timespan_weeks",
    "independent_detections",
    "unique_zones",
    "total_zones",
}


def validate_fasta(path: str) -> bool:
    """
    Check whether a file appears to be a valid FASTA.

    Performs lightweight validation: checks the file exists, has a recognised
    extension, is non-empty, and the first non-blank line starts with ``>``.

    Args:
        path: File path to validate.

    Returns:
        ``True`` if the file looks like a valid FASTA, ``False`` otherwise.
    """
    p = Path(path)

    if not p.is_file():
        logger.debug("Not a file: %s", path)
        return False

    # Check extension
    suffixes = "".join(p.suffixes).lower()
    if not any(suffixes.endswith(ext) for ext in FASTA_EXTENSIONS):
        logger.debug("Unrecognised extension: %s", path)
        return False

    # Check file is non-empty and starts with '>'
    try:
        if suffixes.endswith(".gz"):
            with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
                first_line = _read_first_nonblank(fh)
        else:
            with open(path, "r", encoding="utf-8", errors="replace") as fh:
                first_line = _read_first_nonblank(fh)

        if first_line is None:
            logger.debug("Empty file: %s", path)
            return False

        if not first_line.startswith(">"):
            logger.debug("First line does not start with '>': %s", path)
            return False

    except Exception as exc:
        logger.debug("Error reading %s: %s", path, exc)
        return False

    return True


def find_fastas(path: str) -> List[str]:
    """
    Recursively find all FASTA files under a directory.

    Files are returned sorted alphabetically by name. Only files passing
    extension checks are included (content is not validated for performance
    on large directories).

    Args:
        path: Directory path to search.

    Returns:
        List of absolute paths to FASTA files found.
    """
    results: List[str] = []
    root = Path(path)

    if not root.is_dir():
        logger.warning("Not a directory: %s", path)
        return results

    for entry in sorted(root.rglob("*")):
        if entry.is_file():
            suffixes = "".join(entry.suffixes).lower()
            if any(suffixes.endswith(ext) for ext in FASTA_EXTENSIONS):
                results.append(str(entry.resolve()))

    return results


def parse_metadata_tsv(path: str) -> Dict[str, Dict[str, Any]]:
    """
    Parse a metadata TSV file for P-Score Level 2 temporal scoring.

    Expected columns (tab-separated, header row required):
      - ``sample_id`` — must match FASTA filename stems
      - ``timespan_weeks`` — duration of persistence observation
      - ``independent_detections`` — number of independent isolation events
      - ``unique_zones`` — number of distinct sampling zones positive
      - ``total_zones`` — total sampling zones in the facility

    Additional columns are silently ignored. Missing numeric values default
    to 0.

    Args:
        path: Path to the metadata TSV file.

    Returns:
        Dict mapping ``sample_id`` to a metadata dict with numeric fields.
    """
    metadata: Dict[str, Dict[str, Any]] = {}

    with open(path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        if reader.fieldnames is None:
            logger.warning("Empty metadata file: %s", path)
            return metadata

        # Warn about missing expected columns
        present = set(reader.fieldnames)
        missing = METADATA_COLUMNS - present
        if missing:
            logger.warning(
                "Metadata file %s is missing columns: %s", path, ", ".join(sorted(missing))
            )

        if "sample_id" not in present:
            logger.error("Metadata file %s has no 'sample_id' column", path)
            return metadata

        for row in reader:
            sample_id = row.get("sample_id", "").strip()
            if not sample_id:
                continue

            meta: Dict[str, Any] = {"sample_id": sample_id}

            for col in ["timespan_weeks", "independent_detections", "unique_zones", "total_zones"]:
                raw = row.get(col, "").strip()
                try:
                    meta[col] = float(raw) if raw else 0
                except ValueError:
                    logger.warning(
                        "Non-numeric value for %s in sample %s: %s",
                        col, sample_id, raw,
                    )
                    meta[col] = 0

            metadata[sample_id] = meta

    logger.info("Loaded metadata for %d samples from %s", len(metadata), path)
    return metadata


def _read_first_nonblank(fh) -> str | None:
    """Read and return the first non-blank line from a file handle."""
    for line in fh:
        stripped = line.strip()
        if stripped:
            return stripped
    return None
