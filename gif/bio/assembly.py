# GIF Framework — Genomic Intelligence Framework
# Copyright 2026 Javier Gamboa
# Licensed under the Apache License, Version 2.0

"""
Assembly QC module.

Uses BioPython to parse FASTA assemblies and compute standard quality metrics
(total length, N50, GC content, contig counts). No external tool dependencies.
"""

from __future__ import annotations

import logging
from pathlib import Path

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

logger = logging.getLogger(__name__)


def assess_assembly(fasta_path: str) -> dict:
    """Parse a FASTA assembly and compute quality metrics.

    Parameters
    ----------
    fasta_path : str
        Path to the genome assembly in FASTA format.

    Returns
    -------
    dict
        Assembly statistics with keys:
        - total_length (int): Sum of all contig lengths in bp.
        - num_contigs (int): Number of contigs/scaffolds.
        - n50 (int): N50 value in bp.
        - gc_content (float): GC fraction (0.0-1.0).
        - largest_contig (int): Length of the longest contig in bp.
    """
    path = Path(fasta_path)
    if not path.exists():
        raise FileNotFoundError(f"Assembly file not found: {fasta_path}")

    contig_lengths: list[int] = []
    total_gc_bases = 0
    total_bases = 0

    for record in SeqIO.parse(str(path), "fasta"):
        seq_len = len(record.seq)
        contig_lengths.append(seq_len)
        total_bases += seq_len
        total_gc_bases += record.seq.count("G") + record.seq.count("C")
        total_gc_bases += record.seq.count("g") + record.seq.count("c")

    if not contig_lengths:
        raise ValueError(f"No contigs found in assembly: {fasta_path}")

    # Sort descending for N50 calculation
    contig_lengths.sort(reverse=True)

    total_length = sum(contig_lengths)
    gc_content = total_gc_bases / total_bases if total_bases > 0 else 0.0

    # Calculate N50
    n50 = _compute_n50(contig_lengths, total_length)

    result = {
        "total_length": total_length,
        "num_contigs": len(contig_lengths),
        "n50": n50,
        "gc_content": round(gc_content, 4),
        "largest_contig": contig_lengths[0],
    }

    logger.info(
        "Assembly QC: %d bp across %d contigs, N50=%d, GC=%.2f%%",
        total_length,
        len(contig_lengths),
        n50,
        gc_content * 100,
    )

    return result


def validate_species(fasta_path: str) -> dict:
    """Basic species validation (placeholder).

    A full implementation would use Kraken2 or Mash for species identification.
    This placeholder returns *Listeria monocytogenes* and performs a simple
    genome-size sanity check.

    Parameters
    ----------
    fasta_path : str
        Path to the genome assembly in FASTA format.

    Returns
    -------
    dict
        Species validation result with keys:
        - species (str): Detected species name.
        - confidence (str): Confidence level ("high", "medium", "low").
        - genome_size_check (bool): Whether genome size is within expected
          range for *L. monocytogenes* (2.7-3.2 Mb).
    """
    stats = assess_assembly(fasta_path)
    total_length = stats["total_length"]

    # L. monocytogenes expected genome size: ~2.9 Mb (range 2.7-3.2 Mb)
    LM_MIN = 2_700_000
    LM_MAX = 3_200_000
    size_ok = LM_MIN <= total_length <= LM_MAX

    if size_ok and 0.36 <= stats["gc_content"] <= 0.40:
        confidence = "high"
    elif size_ok:
        confidence = "medium"
    else:
        confidence = "low"

    return {
        "species": "Listeria monocytogenes",
        "confidence": confidence,
        "genome_size_check": size_ok,
    }


# ---- private helpers -------------------------------------------------------

def _compute_n50(sorted_lengths: list[int], total_length: int) -> int:
    """Compute N50 from a descending-sorted list of contig lengths."""
    cumulative = 0
    half = total_length / 2
    for length in sorted_lengths:
        cumulative += length
        if cumulative >= half:
            return length
    return sorted_lengths[-1]
