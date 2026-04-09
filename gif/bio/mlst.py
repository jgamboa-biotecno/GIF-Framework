# GIF Framework — Genomic Intelligence Framework
# Copyright 2026 Javier Gamboa
# Licensed under the Apache License, Version 2.0

"""
MLST and cgMLST typing wrappers.

Wraps the Seemann ``mlst`` tool for 7-gene MLST and ``chewBBACA`` for
cgMLST allele calling. Includes ST-to-CC and CC-to-lineage lookup tables.
"""

from __future__ import annotations

import csv
import logging
import shutil
import subprocess
from pathlib import Path

from gif.utils.wsl import require_tool, run_tool

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# ST -> CC lookup (most clinically relevant Listeria monocytogenes CCs)
# ---------------------------------------------------------------------------
ST_TO_CC: dict[int, str] = {
    1: "CC1",
    2: "CC2",
    3: "CC3",
    4: "CC4",
    5: "CC5",
    6: "CC6",
    7: "CC7",
    8: "CC8",
    9: "CC9",
    11: "CC11",
    14: "CC14",
    87: "CC87",
    121: "CC121",
    155: "CC155",
    204: "CC204",
    224: "CC224",
}

# ---------------------------------------------------------------------------
# CC -> Lineage mapping (Ragon et al. / Moura et al.)
# ---------------------------------------------------------------------------
_CC_LINEAGE: dict[str, str] = {
    # Lineage I — serotypes 1/2b, 4b; hypervirulent CCs
    "CC1": "Lineage_I",
    "CC2": "Lineage_I",
    "CC4": "Lineage_I",
    "CC6": "Lineage_I",
    "CC14": "Lineage_I",
    "CC87": "Lineage_I",
    # Lineage II — serotypes 1/2a, 1/2c; food-adapted CCs
    "CC3": "Lineage_II",
    "CC5": "Lineage_II",
    "CC7": "Lineage_II",
    "CC8": "Lineage_II",
    "CC9": "Lineage_II",
    "CC11": "Lineage_II",
    "CC121": "Lineage_II",
    "CC155": "Lineage_II",
    "CC204": "Lineage_II",
    "CC224": "Lineage_II",
    # Lineage III — rare, animal-associated
    "CC73": "Lineage_III",
}


def run_mlst(fasta_path: str) -> dict:
    """Run 7-gene MLST using the Seemann ``mlst`` tool.

    Parameters
    ----------
    fasta_path : str
        Path to the genome assembly in FASTA format.

    Returns
    -------
    dict
        MLST result with keys:
        - scheme (str): MLST scheme used (e.g. "lmonocytogenes").
        - st (int | None): Sequence type number.
        - cc (str | None): Clonal complex derived from ST.
        - alleles (dict[str, str]): Per-locus allele assignments.
    """
    require_tool("mlst", "conda install -c bioconda mlst")

    path = Path(fasta_path)
    if not path.exists():
        raise FileNotFoundError(f"Assembly file not found: {fasta_path}")

    cmd = ["mlst", "--scheme", "listeria_2", str(path)]
    logger.info("Running MLST: %s", " ".join(cmd))

    stdout, stderr, returncode = run_tool(cmd, timeout=300)

    if returncode != 0:
        logger.error("mlst stderr: %s", stderr.strip())
        raise RuntimeError(f"mlst failed (exit {returncode}): {stderr.strip()}")

    return _parse_mlst_output(stdout.strip())


def run_cgmlst(fasta_path: str, schema_dir: str) -> dict:
    """Run cgMLST allele calling using chewBBACA.

    Parameters
    ----------
    fasta_path : str
        Path to the genome assembly in FASTA format.
    schema_dir : str
        Path to the chewBBACA cgMLST schema directory.

    Returns
    -------
    dict
        cgMLST result with keys:
        - profile (dict[str, str]): Locus name -> allele ID mapping.
        - loci_called (int): Number of loci with valid allele calls.
        - loci_total (int): Total number of loci in the schema.
        - completion_rate (float): Fraction of loci successfully called.
    """
    _require_tool("chewBBACA.py", "pip install chewbbaca")

    fasta = Path(fasta_path)
    schema = Path(schema_dir)
    if not fasta.exists():
        raise FileNotFoundError(f"Assembly file not found: {fasta_path}")
    if not schema.exists():
        raise FileNotFoundError(f"cgMLST schema not found: {schema_dir}")

    # chewBBACA requires a directory of genomes or a list file
    work_dir = fasta.parent / "chewbbaca_work"
    work_dir.mkdir(exist_ok=True)
    output_dir = work_dir / "results"

    cmd = [
        "chewBBACA.py", "AlleleCall",
        "-i", str(fasta.parent),
        "-g", str(schema),
        "-o", str(output_dir),
        "--cpu", "4",
    ]
    logger.info("Running chewBBACA: %s", " ".join(cmd))

    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

    if proc.returncode != 0:
        logger.error("chewBBACA stderr: %s", proc.stderr.strip())
        raise RuntimeError(
            f"chewBBACA failed (exit {proc.returncode}): {proc.stderr.strip()}"
        )

    return _parse_chewbbaca_results(output_dir)


def get_lineage_from_cc(cc: str) -> str:
    """Map a clonal complex to its phylogenetic lineage.

    Parameters
    ----------
    cc : str
        Clonal complex identifier (e.g. "CC1", "CC121").

    Returns
    -------
    str
        Lineage designation ("Lineage_I", "Lineage_II", "Lineage_III",
        or "Unknown").
    """
    return _CC_LINEAGE.get(cc, "Unknown")


# ---- private helpers -------------------------------------------------------

def _require_tool(name: str, install_hint: str) -> None:
    """Raise RuntimeError if an external tool is not on PATH."""
    if shutil.which(name) is None:
        raise RuntimeError(
            f"Tool '{name}' not found on PATH. Install via: {install_hint}"
        )


def _parse_mlst_output(line: str) -> dict:
    """Parse a single line of ``mlst`` tab-delimited output.

    Format: <filename>\\t<scheme>\\t<ST>\\t<allele1>\\t...\\t<alleleN>
    """
    parts = line.split("\t")
    if len(parts) < 3:
        raise ValueError(f"Unexpected mlst output format: {line!r}")

    scheme = parts[1]

    # ST may be "-" if novel
    raw_st = parts[2]
    try:
        st = int(raw_st)
    except (ValueError, TypeError):
        st = None

    # Alleles are columns 3+; each is "locus(allele_id)"
    allele_names = [
        "abcZ", "bglA", "cat", "dapE", "dat", "ldh", "lhkA",
    ]
    alleles: dict[str, str] = {}
    for i, allele_val in enumerate(parts[3:]):
        locus = allele_names[i] if i < len(allele_names) else f"locus_{i+1}"
        alleles[locus] = allele_val

    cc = ST_TO_CC.get(st) if st is not None else None

    return {
        "scheme": scheme,
        "st": st,
        "cc": cc,
        "alleles": alleles,
    }


def _parse_chewbbaca_results(output_dir: Path) -> dict:
    """Parse the chewBBACA AlleleCall results_alleles.tsv file."""
    # chewBBACA outputs into a timestamped subdirectory
    results_dirs = sorted(output_dir.glob("results_*"))
    if not results_dirs:
        raise FileNotFoundError(
            f"No chewBBACA results found in {output_dir}"
        )

    latest = results_dirs[-1]
    tsv_file = latest / "results_alleles.tsv"
    if not tsv_file.exists():
        raise FileNotFoundError(f"results_alleles.tsv not found in {latest}")

    profile: dict[str, str] = {}
    loci_called = 0
    loci_total = 0

    with open(tsv_file, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            # First row is our single sample; columns are locus names
            for locus, allele in row.items():
                if locus in ("FILE", "Sample"):
                    continue
                loci_total += 1
                # Valid alleles are integers; missing/novel are prefixed
                # (e.g., "INF-", "LNF", "PLOT", "NIPH", "ASM", "ALM")
                allele = allele.strip()
                profile[locus] = allele
                try:
                    int(allele)
                    loci_called += 1
                except ValueError:
                    # INF- prefix means inferred — still a valid call
                    if allele.startswith("INF-"):
                        loci_called += 1
            break  # Only one sample row expected

    completion_rate = loci_called / loci_total if loci_total > 0 else 0.0

    logger.info(
        "cgMLST: %d/%d loci called (%.1f%%)",
        loci_called,
        loci_total,
        completion_rate * 100,
    )

    return {
        "profile": profile,
        "loci_called": loci_called,
        "loci_total": loci_total,
        "completion_rate": round(completion_rate, 4),
    }
