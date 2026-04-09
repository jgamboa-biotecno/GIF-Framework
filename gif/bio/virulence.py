# GIF Framework — Genomic Intelligence Framework
# Copyright 2026 Javier Gamboa
# Licensed under the Apache License, Version 2.0

"""
Virulence marker detection module.

Detects Listeria pathogenicity islands (LIPI-1, LIPI-3, LIPI-4) and
Internalin A (InlA) status using BLASTn or ABRicate against curated
virulence databases.
"""

from __future__ import annotations

import logging
import shutil
import subprocess

from gif.utils.wsl import has_tool, run_tool
from pathlib import Path

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Gene sets for each pathogenicity island
# ---------------------------------------------------------------------------
LIPI1_GENES = {"prfA", "plcA", "hly", "mpl", "actA", "plcB"}

LIPI3_GENES = {"llsA", "llsG", "llsH", "llsX", "llsB", "llsY", "llsD", "llsP"}

LIPI4_GENES = {
    "LM9005581_RS14095",
    "LM9005581_RS14100",
    "LM9005581_RS14105",
    "LM9005581_RS14110",
    "LM9005581_RS14115",
    "LM9005581_RS14120",
}

# CCs where LIPI-4 is biologically validated
LIPI4_VALID_CCS = {"CC4", "CC87"}

# InlA reference protein length (full-length)
INLA_FULL_LENGTH = 800  # amino acids


def detect_virulence(fasta_path: str, db_dir: str, cc: str | None = None) -> dict:
    """Detect virulence markers in a Listeria genome assembly.

    Uses ABRicate (preferred) or BLASTn against a curated virulence database
    to identify LIPI-1, LIPI-3, LIPI-4 genes and InlA status.

    Parameters
    ----------
    fasta_path : str
        Path to the genome assembly in FASTA format.
    db_dir : str
        Path to the virulence reference database directory containing
        FASTA files for each marker category.
    cc : str or None
        Clonal complex of the isolate (e.g. "CC4"). Used for LIPI-4
        context validation — LIPI-4 is only scored as complete when
        detected in CC4 or CC87 backgrounds.

    Returns
    -------
    dict
        Virulence markers with keys:
        - inla_status (str): "full_length", "truncated", or "absent".
        - inla_protein_length (int | None): Estimated protein length.
        - lipi1_intact (bool): True if all 6 LIPI-1 genes detected.
        - lipi1_deleted_genes (list[str]): Missing LIPI-1 genes.
        - lipi3_complete (bool): True if all 8 LIPI-3 genes detected.
        - lipi3_genes_detected (list[str]): LIPI-3 genes found.
        - lipi4_complete (bool): True if all 6 LIPI-4 genes detected
          AND the isolate belongs to CC4/CC87.
        - lipi4_genes_detected (list[str]): LIPI-4 genes found.
    """
    fasta = Path(fasta_path)
    if not fasta.exists():
        raise FileNotFoundError(f"Assembly file not found: {fasta_path}")

    db = Path(db_dir)
    if not db.exists():
        raise FileNotFoundError(f"Virulence database not found: {db_dir}")

    # Use BLASTn against our reference FASTAs (preferred for local DBs)
    # ABRicate requires its own DB format — our refs are plain FASTAs
    if has_tool("blastn"):
        hits = _run_blastn(fasta_path, db_dir)
    elif has_tool("abricate"):
        hits = _run_abricate(fasta_path, db_dir)
    else:
        raise RuntimeError(
            "Neither 'blastn' nor 'abricate' found on PATH or in WSL. "
            "Install via: conda install -c bioconda blast abricate"
        )

    detected_genes_raw = {h["gene"].lower() for h in hits}

    # Normalize gene names: "actA_LIPI1" -> "acta", "prfA_reference" -> "prfa"
    detected_genes = set()
    for g in detected_genes_raw:
        # Take the part before the first underscore as the gene name
        base = g.split("_")[0]
        detected_genes.add(base)
        detected_genes.add(g)  # also keep original for exact matches

    # --- LIPI-1 ---
    lipi1_found = {g for g in LIPI1_GENES if g.lower() in detected_genes}
    lipi1_deleted = sorted(LIPI1_GENES - lipi1_found)
    lipi1_intact = len(lipi1_deleted) == 0

    # --- LIPI-3 ---
    lipi3_found = {g for g in LIPI3_GENES if g.lower() in detected_genes}
    lipi3_complete = len(lipi3_found) == len(LIPI3_GENES)

    # --- LIPI-4 (context-validated) ---
    lipi4_found_genes = set()
    for g in LIPI4_GENES:
        # Match by suffix since gene names may vary
        for hit_gene in detected_genes:
            if g.lower() in hit_gene or hit_gene in g.lower():
                lipi4_found_genes.add(g)
                break
    lipi4_all_present = len(lipi4_found_genes) == len(LIPI4_GENES)
    lipi4_context_valid = cc in LIPI4_VALID_CCS if cc else False
    lipi4_complete = lipi4_all_present and lipi4_context_valid

    # --- InlA ---
    inla_status, inla_length = _assess_inla(hits)

    result = {
        "inla_status": inla_status,
        "inla_protein_length": inla_length,
        "lipi1_intact": lipi1_intact,
        "lipi1_deleted_genes": lipi1_deleted,
        "lipi3_complete": lipi3_complete,
        "lipi3_genes_detected": sorted(lipi3_found),
        "lipi4_complete": lipi4_complete,
        "lipi4_genes_detected": sorted(lipi4_found_genes),
    }

    logger.info(
        "Virulence: LIPI-1=%s, LIPI-3=%s (%d/8), LIPI-4=%s (%d/6), InlA=%s",
        "intact" if lipi1_intact else f"deleted({len(lipi1_deleted)})",
        "complete" if lipi3_complete else "partial",
        len(lipi3_found),
        "complete" if lipi4_complete else "partial/absent",
        len(lipi4_found_genes),
        inla_status,
    )

    return result


# ---- private helpers -------------------------------------------------------

def _run_abricate(fasta_path: str, db_dir: str) -> list[dict]:
    """Run ABRicate and return parsed hit records."""
    cmd = [
        "abricate",
        "--db", db_dir,
        "--minid", "80",
        "--mincov", "60",
        str(fasta_path),
    ]
    logger.info("Running ABRicate: %s", " ".join(cmd))

    stdout, stderr, returncode = run_tool(cmd, timeout=600)
    if returncode != 0:
        logger.error("ABRicate stderr: %s", stderr.strip())
        raise RuntimeError(
            f"ABRicate failed (exit {returncode}): {stderr.strip()}"
        )

    return _parse_abricate_output(stdout)


def _run_blastn(fasta_path: str, db_dir: str) -> list[dict]:
    """Run BLASTn against each reference FASTA in db_dir and return hits."""
    db_path = Path(db_dir)
    all_hits: list[dict] = []

    # Look for FASTA reference files in the database directory
    ref_files = list(db_path.glob("*.fasta")) + list(db_path.glob("*.fa"))
    if not ref_files:
        raise FileNotFoundError(
            f"No reference FASTA files found in {db_dir}"
        )

    for ref in ref_files:
        cmd = [
            "blastn",
            "-query", str(fasta_path),
            "-subject", str(ref),
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen "
                       "qstart qend sstart send evalue bitscore slen",
            "-evalue", "1e-10",
            "-max_target_seqs", "50",
        ]

        stdout, stderr, returncode = run_tool(cmd, timeout=300)
        if returncode != 0:
            logger.warning("BLASTn failed for %s: %s", ref.name, stderr.strip())
            continue

        for line in stdout.strip().split("\n"):
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 13:
                continue
            pident = float(fields[2])
            aln_len = int(fields[3])
            slen = int(fields[12])
            coverage = (aln_len / slen * 100) if slen > 0 else 0

            if pident >= 80 and coverage >= 60:
                gene_name = fields[1].split("|")[0] if "|" in fields[1] else fields[1]
                all_hits.append({
                    "gene": gene_name,
                    "identity": pident,
                    "coverage": round(coverage, 2),
                    "length": aln_len,
                    "subject_length": slen,
                })

    return all_hits


def _parse_abricate_output(stdout: str) -> list[dict]:
    """Parse ABRicate tab-delimited output into hit dicts."""
    hits: list[dict] = []
    for line in stdout.strip().split("\n"):
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 10:
            continue
        hits.append({
            "gene": fields[5],         # GENE column
            "identity": float(fields[9]),  # %IDENTITY
            "coverage": float(fields[8]),  # %COVERAGE
            "length": int(fields[7]) - int(fields[6]) + 1,
            "subject_length": 0,
        })
    return hits


def _assess_inla(hits: list[dict]) -> tuple[str, int | None]:
    """Determine InlA status from BLAST/ABRicate hits.

    Returns (status, estimated_protein_length).
    """
    inla_hits = [h for h in hits if "inla" in h["gene"].lower()]

    if not inla_hits:
        return "absent", None

    # Take the best hit by identity
    best = max(inla_hits, key=lambda h: h["identity"])
    aln_len_aa = best["length"] // 3  # approximate protein length

    if aln_len_aa >= INLA_FULL_LENGTH * 0.95:
        return "full_length", aln_len_aa
    else:
        # Shorter alignment suggests premature stop codon (PMSC)
        return "truncated", aln_len_aa
