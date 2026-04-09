# GIF Framework — Genomic Intelligence Framework
# Copyright 2026 Javier Gamboa
# Licensed under the Apache License, Version 2.0

"""
Antimicrobial resistance (AMR) detection module.

Wraps AMRFinderPlus or ABRicate for AMR gene detection. Categorises
detected genes by drug class and identifies MDR plasmid signatures
(contigs >= 40 kb carrying >= 3 AMR genes from >= 2 drug categories).
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

from gif.utils.wsl import has_tool, run_tool

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Gene -> drug category mapping (canonical AMR genes in Listeria)
# ---------------------------------------------------------------------------
GENE_CATEGORY: dict[str, str] = {
    # Beta-lactam
    "bla": "beta_lactam", "blaZ": "beta_lactam", "penA": "beta_lactam",
    # Tetracycline
    "tetA": "tetracycline", "tetB": "tetracycline", "tetC": "tetracycline",
    "tetK": "tetracycline", "tetL": "tetracycline", "tetM": "tetracycline",
    "tetO": "tetracycline", "tetS": "tetracycline", "tetW": "tetracycline",
    # Aminoglycoside
    "aadA": "aminoglycoside", "aadB": "aminoglycoside",
    "aadE": "aminoglycoside", "aac(6')-Ie-aph(2'')-Ia": "aminoglycoside",
    "aacA4": "aminoglycoside", "aph(3')-III": "aminoglycoside",
    "ant(6)-Ia": "aminoglycoside", "str": "aminoglycoside",
    # Fluoroquinolone
    "qnrA": "fluoroquinolone", "qnrB": "fluoroquinolone",
    "qnrS": "fluoroquinolone", "gyrA": "fluoroquinolone",
    # Macrolide
    "ermA": "macrolide", "ermB": "macrolide", "ermC": "macrolide",
    "mefA": "macrolide", "msrA": "macrolide", "ereA": "macrolide",
    "ereB": "macrolide",
    # Phenicol
    "cat": "phenicol", "catA": "phenicol", "catB": "phenicol",
    "floR": "phenicol", "fexA": "phenicol", "cfr": "phenicol",
    # Lincosamide
    "lnuA": "lincosamide", "lnuB": "lincosamide", "lnuG": "lincosamide",
    "vgaA": "lincosamide",
    # Trimethoprim
    "dfrA": "trimethoprim", "dfrD": "trimethoprim", "dfrG": "trimethoprim",
    "dfrK": "trimethoprim",
}

# All recognised drug categories
DRUG_CATEGORIES = [
    "beta_lactam", "tetracycline", "aminoglycoside", "fluoroquinolone",
    "macrolide", "phenicol", "lincosamide", "trimethoprim",
]

# MDR plasmid detection thresholds
MDR_CONTIG_MIN_LENGTH = 40_000   # 40 kb
MDR_MIN_GENES = 3
MDR_MIN_CATEGORIES = 2


def detect_amr(fasta_path: str) -> dict:
    """Detect antimicrobial resistance genes in a genome assembly.

    Uses AMRFinderPlus (preferred) or ABRicate with the NCBI/ResFinder
    database. Categorises genes by drug class and checks for MDR plasmid
    signatures.

    Parameters
    ----------
    fasta_path : str
        Path to the genome assembly in FASTA format.

    Returns
    -------
    dict
        AMR detection results with keys:
        - genes_detected (list[dict]): Each dict has gene, category,
          identity, coverage, contig.
        - categories_detected (list[str]): Drug categories with hits.
        - by_category (dict[str, list[str]]): Gene names grouped by
          drug category.
        - total_genes (int): Total AMR genes found.
        - mdr_plasmid_detected (bool): True if MDR plasmid criteria met.
        - mdr_plasmid_details (dict | None): Details of the MDR contig
          if detected.
    """
    path = Path(fasta_path)
    if not path.exists():
        raise FileNotFoundError(f"Assembly file not found: {fasta_path}")

    # Try AMRFinderPlus first, then ABRicate
    if has_tool("amrfinder"):
        hits = _run_amrfinder(fasta_path)
    elif has_tool("abricate"):
        hits = _run_abricate_amr(fasta_path)
    else:
        raise RuntimeError(
            "Neither 'amrfinder' nor 'abricate' found on PATH or in WSL. "
            "Install via: conda install -c bioconda ncbi-amrfinderplus abricate"
        )

    # Categorise genes
    genes_detected: list[dict] = []
    by_category: dict[str, list[str]] = defaultdict(list)
    contig_genes: dict[str, list[dict]] = defaultdict(list)

    for hit in hits:
        gene = hit["gene"]
        contig = hit.get("contig", "unknown")
        category = _categorise_gene(gene)

        entry = {
            "gene": gene,
            "category": category,
            "identity": hit.get("identity", 0),
            "coverage": hit.get("coverage", 0),
            "contig": contig,
        }
        genes_detected.append(entry)
        by_category[category].append(gene)
        contig_genes[contig].append(entry)

    categories_detected = sorted(
        cat for cat in DRUG_CATEGORIES if cat in by_category
    )

    # MDR plasmid detection
    mdr_detected, mdr_details = _detect_mdr_plasmid(contig_genes, fasta_path)

    result = {
        "genes_detected": genes_detected,
        "categories_detected": categories_detected,
        "by_category": dict(by_category),
        "total_genes": len(genes_detected),
        "mdr_plasmid_detected": mdr_detected,
        "mdr_plasmid_details": mdr_details,
    }

    logger.info(
        "AMR: %d genes in %d categories, MDR plasmid=%s",
        len(genes_detected),
        len(categories_detected),
        mdr_detected,
    )

    return result


# ---- private helpers -------------------------------------------------------

def _run_amrfinder(fasta_path: str) -> list[dict]:
    """Run AMRFinderPlus and return parsed hits."""
    cmd = [
        "amrfinder",
        "--nucleotide", str(fasta_path),
        "--plus",
    ]
    logger.info("Running AMRFinderPlus: %s", " ".join(cmd))

    stdout, stderr, returncode = run_tool(cmd, timeout=600)
    if returncode != 0:
        logger.error("amrfinder stderr: %s", stderr.strip())
        raise RuntimeError(
            f"AMRFinderPlus failed (exit {returncode}): {stderr.strip()}"
        )

    hits: list[dict] = []
    lines = stdout.strip().split("\n")
    if len(lines) < 2:
        return hits

    # Header is first line; data starts at line 2
    header = lines[0].split("\t")
    for line in lines[1:]:
        if not line.strip():
            continue
        fields = line.split("\t")
        record = dict(zip(header, fields))

        hits.append({
            "gene": record.get("Gene symbol", record.get("Name", "")),
            "identity": float(record.get("% Identity to reference sequence", 0)),
            "coverage": float(record.get("% Coverage of reference sequence", 0)),
            "contig": record.get("Contig id", "unknown"),
        })

    return hits


def _run_abricate_amr(fasta_path: str) -> list[dict]:
    """Run ABRicate with ResFinder database and return parsed hits."""
    cmd = [
        "abricate",
        "--db", "resfinder",
        "--minid", "80",
        "--mincov", "60",
        str(fasta_path),
    ]
    logger.info("Running ABRicate (resfinder): %s", " ".join(cmd))

    stdout, stderr, returncode = run_tool(cmd, timeout=600)
    if returncode != 0:
        logger.error("ABRicate stderr: %s", stderr.strip())
        raise RuntimeError(
            f"ABRicate failed (exit {returncode}): {stderr.strip()}"
        )

    hits: list[dict] = []
    for line in stdout.strip().split("\n"):
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 10:
            continue
        try:
            identity = float(fields[9])
        except (ValueError, IndexError):
            identity = 0.0
        try:
            cov_str = fields[8].split("/")[0] if "/" in fields[8] else fields[8]
            coverage = float(cov_str)
        except (ValueError, IndexError):
            coverage = 0.0
        hits.append({
            "gene": fields[5],
            "identity": identity,
            "coverage": coverage,
            "contig": fields[1],
        })

    return hits


def _categorise_gene(gene: str) -> str:
    """Map a gene name to its drug category.

    Performs exact match first, then prefix matching for gene families.
    """
    if gene in GENE_CATEGORY:
        return GENE_CATEGORY[gene]

    # Try prefix matching (e.g. "tetM_1" -> "tetM" -> "tetracycline")
    gene_base = gene.split("_")[0].split("-")[0]
    if gene_base in GENE_CATEGORY:
        return GENE_CATEGORY[gene_base]

    # Try 3-character prefix for gene families (tet, erm, aad, etc.)
    for known_gene, category in GENE_CATEGORY.items():
        if gene.lower().startswith(known_gene.lower()):
            return category

    return "other"


def _detect_mdr_plasmid(
    contig_genes: dict[str, list[dict]],
    fasta_path: str,
) -> tuple[bool, dict | None]:
    """Check if any contig meets MDR plasmid criteria.

    MDR plasmid = contig >= 40 kb with >= 3 AMR genes from >= 2 categories.
    """
    from Bio import SeqIO

    # Build contig length lookup
    contig_lengths: dict[str, int] = {}
    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            contig_lengths[record.id] = len(record.seq)
    except Exception as exc:
        logger.warning("Could not parse contigs for MDR check: %s", exc)
        return False, None

    for contig_id, genes in contig_genes.items():
        length = contig_lengths.get(contig_id, 0)
        if length < MDR_CONTIG_MIN_LENGTH:
            continue

        categories = {g["category"] for g in genes}
        categories.discard("other")

        if len(genes) >= MDR_MIN_GENES and len(categories) >= MDR_MIN_CATEGORIES:
            details = {
                "contig": contig_id,
                "contig_length": length,
                "amr_genes": [g["gene"] for g in genes],
                "categories": sorted(categories),
                "gene_count": len(genes),
                "category_count": len(categories),
            }
            logger.info(
                "MDR plasmid detected: contig %s (%d bp), %d genes, %d categories",
                contig_id, length, len(genes), len(categories),
            )
            return True, details

    return False, None
