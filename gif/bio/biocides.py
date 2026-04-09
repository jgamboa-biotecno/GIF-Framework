# GIF Framework — Genomic Intelligence Framework
# Copyright 2026 Javier Gamboa
# Licensed under the Apache License, Version 2.0

"""
Biocide resistance and persistence marker detection.

Detects quaternary ammonium compound (QAC) efflux pumps, stress survival
islets (SSI-1, SSI-2), acid tolerance (GAD system), heavy metal resistance
(cadmium), CRISPR arrays, and prophage regions in Listeria genomes.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path

from gif.utils.wsl import has_tool, run_tool

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Marker gene sets
# ---------------------------------------------------------------------------

# QAC resistance
QAC_GENES = {"qacH", "qacEdelta1"}

# Benzalkonium chloride resistance operon
BCR_GENES = {"bcrA", "bcrB", "bcrC"}

# Stress Survival Islet 1 (SSI-1): salt/acid/bile tolerance
SSI1_GENES = {"lmo0444", "lmo0445", "lmo0446", "lmo0447", "lmo0448"}

# Stress Survival Islet 2 (SSI-2): alkaline/oxidative stress (Lineage II)
SSI2_GENES = {"lin0464", "lin0465"}

# Glutamate decarboxylase acid resistance system
GAD_GENES = {"gadD1", "gadD2", "gadD3", "gadE", "gadT1", "gadT2"}

# Cadmium resistance determinants
CADMIUM_GENES = {"cadA1", "cadA2", "cadA3", "cadA4", "cadA5", "cadC"}

# All persistence gene names (flattened for lookup)
ALL_PERSISTENCE_GENES = (
    QAC_GENES | BCR_GENES | SSI1_GENES | SSI2_GENES
    | GAD_GENES | CADMIUM_GENES
)


def detect_persistence_markers(fasta_path: str, db_dir: str) -> dict:
    """Detect biocide resistance and environmental persistence markers.

    Uses BLASTn against a curated persistence marker database. Optionally
    detects CRISPR arrays (via minced) and prophage regions (via PHASTER
    or simple contig-length heuristics).

    Parameters
    ----------
    fasta_path : str
        Path to the genome assembly in FASTA format.
    db_dir : str
        Path to the persistence marker reference database directory.

    Returns
    -------
    dict
        Persistence markers with keys:
        - qac_genes (list[str]): QAC efflux pump genes detected.
        - bcr_complete (bool): True if all 3 bcrABC genes present.
        - bcr_genes (list[str]): Detected bcrABC components.
        - ssi1_complete (bool): True if all 5 SSI-1 genes present.
        - ssi1_genes (list[str]): Detected SSI-1 genes.
        - ssi2_complete (bool): True if both SSI-2 genes present.
        - ssi2_genes (list[str]): Detected SSI-2 genes.
        - gad_genes (list[str]): Detected GAD system genes.
        - gad_system_present (bool): True if core GAD genes detected.
        - cadmium_genes (list[str]): Detected cadmium resistance genes.
        - cadmium_resistance (bool): True if any cadA variant found.
        - crispr_arrays (int): Number of CRISPR arrays detected.
        - prophage_count (int): Number of prophage regions detected.
        - total_persistence_genes (int): Total persistence genes found.
    """
    fasta = Path(fasta_path)
    if not fasta.exists():
        raise FileNotFoundError(f"Assembly file not found: {fasta_path}")

    db = Path(db_dir)
    if not db.exists():
        raise FileNotFoundError(f"Persistence database not found: {db_dir}")

    # Detect persistence genes via BLAST
    detected_genes = _blast_persistence_markers(fasta_path, db_dir)

    # Normalize BLAST hit names: "qacH_reference" -> "qach", "lmo0444_reference" -> "lmo0444"
    detected_lower = set()
    for g in detected_genes:
        gl = g.lower()
        detected_lower.add(gl)
        # Strip common suffixes from reference FASTA names
        for suffix in ("_reference", "_ref", "_full"):
            if gl.endswith(suffix):
                detected_lower.add(gl[: -len(suffix)])
        # Also keep just the first part before underscore
        base = gl.split("_")[0]
        detected_lower.add(base)

    qac_found = sorted(g for g in QAC_GENES if g.lower() in detected_lower)
    bcr_found = sorted(g for g in BCR_GENES if g.lower() in detected_lower)
    ssi1_found = sorted(g for g in SSI1_GENES if g.lower() in detected_lower)
    ssi2_found = sorted(g for g in SSI2_GENES if g.lower() in detected_lower)
    gad_found = sorted(g for g in GAD_GENES if g.lower() in detected_lower)
    cadmium_found = sorted(g for g in CADMIUM_GENES if g.lower() in detected_lower)

    # CRISPR detection
    crispr_count = _detect_crispr(fasta_path)

    # Prophage detection
    prophage_count = _detect_prophages(fasta_path)

    total = (
        len(qac_found) + len(bcr_found) + len(ssi1_found) + len(ssi2_found)
        + len(gad_found) + len(cadmium_found)
    )

    result = {
        "qac_genes": qac_found,
        "bcr_complete": len(bcr_found) == len(BCR_GENES),
        "bcr_genes": bcr_found,
        "ssi1_complete": len(ssi1_found) == len(SSI1_GENES),
        "ssi1_genes": ssi1_found,
        "ssi2_complete": len(ssi2_found) == len(SSI2_GENES),
        "ssi2_genes": ssi2_found,
        "gad_genes": gad_found,
        "gad_system_present": any(
            g.lower() in detected_lower for g in ("gadD1", "gadD2", "gadD3")
        ),
        "cadmium_genes": cadmium_found,
        "cadmium_resistance": len(cadmium_found) > 0,
        "crispr_arrays": crispr_count,
        "prophage_count": prophage_count,
        "total_persistence_genes": total,
    }

    logger.info(
        "Persistence: QAC=%d, bcrABC=%s, SSI-1=%s, SSI-2=%s, "
        "GAD=%d, Cd=%d, CRISPR=%d, prophage=%d",
        len(qac_found),
        "complete" if result["bcr_complete"] else f"partial({len(bcr_found)})",
        "complete" if result["ssi1_complete"] else f"partial({len(ssi1_found)})",
        "complete" if result["ssi2_complete"] else f"partial({len(ssi2_found)})",
        len(gad_found),
        len(cadmium_found),
        crispr_count,
        prophage_count,
    )

    return result


# ---- private helpers -------------------------------------------------------

def _blast_persistence_markers(fasta_path: str, db_dir: str) -> set[str]:
    """Run BLASTn against the persistence marker database."""
    if not has_tool("blastn"):
        raise RuntimeError(
            "Tool 'blastn' not found on PATH or in WSL. "
            "Install via: conda install -c bioconda blast"
        )

    db_path = Path(db_dir)
    ref_files = list(db_path.glob("*.fasta")) + list(db_path.glob("*.fa"))
    if not ref_files:
        raise FileNotFoundError(
            f"No reference FASTA files found in {db_dir}"
        )

    detected: set[str] = set()

    for ref in ref_files:
        cmd = [
            "blastn",
            "-query", str(fasta_path),
            "-subject", str(ref),
            "-outfmt", "6 qseqid sseqid pident length slen",
            "-evalue", "1e-10",
            "-max_target_seqs", "100",
        ]

        stdout, stderr, returncode = run_tool(cmd, timeout=300)
        if returncode != 0:
            logger.warning("BLASTn failed for %s: %s", ref.name, stderr.strip())
            continue

        for line in stdout.strip().split("\n"):
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                continue

            pident = float(fields[2])
            aln_len = int(fields[3])
            slen = int(fields[4])
            coverage = (aln_len / slen * 100) if slen > 0 else 0

            if pident >= 80 and coverage >= 60:
                gene_name = fields[1].split("|")[0] if "|" in fields[1] else fields[1]
                detected.add(gene_name)

    return detected


def _detect_crispr(fasta_path: str) -> int:
    """Detect CRISPR arrays using MinCED if available.

    Returns the number of CRISPR arrays found, or 0 if MinCED is not
    installed.
    """
    if not has_tool("minced"):
        logger.debug("minced not found; skipping CRISPR detection")
        return 0

    cmd = ["minced", str(fasta_path)]
    try:
        stdout, stderr, returncode = run_tool(cmd, timeout=300)
    except subprocess.TimeoutExpired:
        logger.warning("minced timed out")
        return 0

    if returncode != 0:
        logger.warning("minced failed: %s", stderr.strip())
        return 0

    # Count lines matching "CRISPR" in the summary output
    count = 0
    for line in stdout.split("\n"):
        if line.strip().startswith("CRISPR"):
            count += 1

    return count


def _detect_prophages(fasta_path: str) -> int:
    """Estimate prophage count using a simple contig-length heuristic.

    A proper implementation would use PHASTER or Phigaro. This heuristic
    counts contigs between 20-60 kb with atypical GC content as potential
    prophage regions.

    Returns the estimated number of prophage regions.
    """
    from Bio import SeqIO

    prophage_candidates = 0

    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            seq_len = len(record.seq)
            # Typical Listeria prophages: 20-60 kb
            if 20_000 <= seq_len <= 60_000:
                seq = str(record.seq).upper()
                gc = (seq.count("G") + seq.count("C")) / seq_len if seq_len > 0 else 0
                # Listeria genome GC ~ 38%; prophages often have different GC
                if gc < 0.34 or gc > 0.42:
                    prophage_candidates += 1
    except Exception as exc:
        logger.warning("Prophage heuristic failed: %s", exc)

    return prophage_candidates
