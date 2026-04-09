# GIF Framework — Genomic Intelligence Framework
# Copyright 2026 Javier Gamboa
# Licensed under the Apache License, Version 2.0

"""
NCBI genome fetching module.

Downloads genome assemblies from NCBI by GCA/GCF accession, SRA runs,
or BioProject identifiers. Handles accession type resolution, SRA
download + assembly, and batch BioProject retrieval.
"""

from __future__ import annotations

import gzip
import logging
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

import requests

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# NCBI Entrez API base URLs
# ---------------------------------------------------------------------------
ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ESEARCH_URL = f"{ENTREZ_BASE}/esearch.fcgi"
EFETCH_URL = f"{ENTREZ_BASE}/efetch.fcgi"
ELINK_URL = f"{ENTREZ_BASE}/elink.fcgi"
ESUMMARY_URL = f"{ENTREZ_BASE}/esummary.fcgi"

# NCBI Datasets API for assembly download
DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"

# Accession patterns
_RE_ASSEMBLY = re.compile(r"^GC[AF]_\d{9}\.\d+$")
_RE_SRA = re.compile(r"^[SED]RR\d+$")
_RE_BIOSAMPLE = re.compile(r"^SAM[NED]\d+$")
_RE_BIOPROJECT = re.compile(r"^PRJ[A-Z]{2}\d+$")


def resolve_accession(accession: str) -> tuple[str, str]:
    """Detect accession type and resolve to a usable identifier.

    Parameters
    ----------
    accession : str
        Any NCBI accession (GCA/GCF, SRR/ERR, SAMN, PRJNA, etc.).

    Returns
    -------
    tuple[str, str]
        (accession_type, resolved_accession) where accession_type is one
        of "assembly", "sra", "biosample", "bioproject".

    Raises
    ------
    ValueError
        If the accession format is not recognised.
    """
    accession = accession.strip()

    if _RE_ASSEMBLY.match(accession):
        return "assembly", accession

    if _RE_SRA.match(accession):
        return "sra", accession

    if _RE_BIOSAMPLE.match(accession):
        # Resolve BioSample -> assembly accession via Entrez
        assembly_acc = _biosample_to_assembly(accession)
        if assembly_acc:
            return "assembly", assembly_acc
        # Fallback: try to find SRA runs
        sra_acc = _biosample_to_sra(accession)
        if sra_acc:
            return "sra", sra_acc
        raise ValueError(
            f"Could not resolve BioSample {accession} to an assembly or SRA run"
        )

    if _RE_BIOPROJECT.match(accession):
        return "bioproject", accession

    raise ValueError(
        f"Unrecognised accession format: {accession!r}. "
        "Expected GCA_/GCF_ (assembly), SRR/ERR (SRA), "
        "SAMN (BioSample), or PRJNA (BioProject)."
    )


def fetch_assembly(accession: str, output_dir: str) -> str:
    """Download a genome assembly from NCBI by GCA/GCF accession.

    Downloads the assembly FASTA via the NCBI Datasets API, decompresses
    it, and returns the local file path.

    Parameters
    ----------
    accession : str
        Assembly accession (e.g. "GCF_000196035.1").
    output_dir : str
        Directory to save the downloaded assembly.

    Returns
    -------
    str
        Path to the downloaded (decompressed) FASTA file.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    output_fasta = out / f"{accession}.fasta"
    if output_fasta.exists():
        logger.info("Assembly already exists: %s", output_fasta)
        return str(output_fasta)

    # Try NCBI Datasets CLI first
    if shutil.which("datasets"):
        return _fetch_via_datasets_cli(accession, out, output_fasta)

    # Fallback: direct FTP download
    return _fetch_via_ftp(accession, out, output_fasta)


def fetch_sra(accession: str, output_dir: str) -> str:
    """Download an SRA run, extract FASTQ, and assemble with SPAdes.

    Parameters
    ----------
    accession : str
        SRA run accession (e.g. "SRR1234567").
    output_dir : str
        Directory for intermediate and output files.

    Returns
    -------
    str
        Path to the assembled genome FASTA.
    """
    _require_tool("fasterq-dump", "conda install -c bioconda sra-tools")
    _require_tool("spades.py", "conda install -c bioconda spades")

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    assembly_path = out / f"{accession}_assembly" / "scaffolds.fasta"
    if assembly_path.exists():
        logger.info("Assembly already exists: %s", assembly_path)
        return str(assembly_path)

    # Step 1: Download FASTQ
    logger.info("Downloading SRA run %s...", accession)
    fastq_dir = out / f"{accession}_fastq"
    fastq_dir.mkdir(exist_ok=True)

    cmd_fasterq = [
        "fasterq-dump",
        "--split-files",
        "--outdir", str(fastq_dir),
        accession,
    ]
    proc = subprocess.run(cmd_fasterq, capture_output=True, text=True, timeout=3600)
    if proc.returncode != 0:
        raise RuntimeError(
            f"fasterq-dump failed for {accession}: {proc.stderr.strip()}"
        )

    # Detect paired vs single reads
    r1 = fastq_dir / f"{accession}_1.fastq"
    r2 = fastq_dir / f"{accession}_2.fastq"
    single = fastq_dir / f"{accession}.fastq"

    # Step 2: Assemble with SPAdes
    assembly_dir = out / f"{accession}_assembly"
    if r1.exists() and r2.exists():
        cmd_spades = [
            "spades.py",
            "--careful",
            "-1", str(r1),
            "-2", str(r2),
            "-o", str(assembly_dir),
            "-t", "4",
            "-m", "8",
        ]
    elif single.exists():
        cmd_spades = [
            "spades.py",
            "--careful",
            "-s", str(single),
            "-o", str(assembly_dir),
            "-t", "4",
            "-m", "8",
        ]
    else:
        raise FileNotFoundError(
            f"No FASTQ files found after fasterq-dump for {accession}"
        )

    logger.info("Assembling %s with SPAdes...", accession)
    proc = subprocess.run(cmd_spades, capture_output=True, text=True, timeout=7200)
    if proc.returncode != 0:
        raise RuntimeError(
            f"SPAdes assembly failed for {accession}: {proc.stderr.strip()}"
        )

    if not assembly_path.exists():
        # Try contigs.fasta as fallback
        contigs = assembly_dir / "contigs.fasta"
        if contigs.exists():
            return str(contigs)
        raise FileNotFoundError(
            f"Assembly output not found for {accession} in {assembly_dir}"
        )

    return str(assembly_path)


def fetch_bioproject(bioproject: str, output_dir: str) -> list[str]:
    """Download all assemblies/runs from an NCBI BioProject.

    Parameters
    ----------
    bioproject : str
        BioProject accession (e.g. "PRJNA12345").
    output_dir : str
        Directory to save downloaded assemblies.

    Returns
    -------
    list[str]
        Paths to all downloaded/assembled genomes.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Get all SRA runs linked to this BioProject
    run_accessions = _get_bioproject_runs(bioproject)

    if not run_accessions:
        raise ValueError(f"No SRA runs found for BioProject {bioproject}")

    paths: list[str] = []
    for i, run_acc in enumerate(run_accessions, 1):
        logger.info(
            "Fetching %s (%d/%d from %s)...",
            run_acc, i, len(run_accessions), bioproject,
        )
        try:
            acc_type, resolved = resolve_accession(run_acc)
            if acc_type == "assembly":
                path = fetch_assembly(resolved, str(out))
            else:
                path = fetch_sra(resolved, str(out))
            paths.append(path)
        except Exception as exc:
            logger.error("Failed to fetch %s: %s", run_acc, exc)

    return paths


# ---- private helpers -------------------------------------------------------

def _require_tool(name: str, install_hint: str) -> None:
    """Raise RuntimeError if an external tool is not on PATH."""
    if shutil.which(name) is None:
        raise RuntimeError(
            f"Tool '{name}' not found on PATH. Install via: {install_hint}"
        )


def _fetch_via_datasets_cli(
    accession: str, out: Path, output_fasta: Path
) -> str:
    """Download assembly using NCBI datasets CLI."""
    zip_path = out / f"{accession}.zip"
    cmd = [
        "datasets", "download", "genome", "accession", accession,
        "--include", "genome",
        "--filename", str(zip_path),
    ]
    logger.info("Downloading via datasets CLI: %s", accession)

    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if proc.returncode != 0:
        raise RuntimeError(
            f"datasets download failed: {proc.stderr.strip()}"
        )

    # Extract the FASTA from the zip
    import zipfile
    extract_dir = out / f"{accession}_extract"
    with zipfile.ZipFile(zip_path, "r") as zf:
        fasta_files = [f for f in zf.namelist() if f.endswith((".fna", ".fasta"))]
        if not fasta_files:
            raise FileNotFoundError(
                f"No FASTA file found in datasets download for {accession}"
            )
        zf.extract(fasta_files[0], extract_dir)
        extracted = extract_dir / fasta_files[0]

    shutil.move(str(extracted), str(output_fasta))

    # Cleanup
    zip_path.unlink(missing_ok=True)
    shutil.rmtree(extract_dir, ignore_errors=True)

    logger.info("Assembly saved: %s", output_fasta)
    return str(output_fasta)


def _fetch_via_ftp(accession: str, out: Path, output_fasta: Path) -> str:
    """Download assembly via NCBI FTP."""
    # Build FTP path from accession
    # GCF_000196035.1 -> GCF/000/196/035/GCF_000196035.1
    prefix = accession[:3]
    digits = accession.split("_")[1].split(".")[0]
    d1, d2, d3 = digits[:3], digits[3:6], digits[6:9]

    base_url = (
        f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{d1}/{d2}/{d3}"
    )

    # List the directory to find the assembly
    logger.info("Fetching assembly listing from %s", base_url)

    # Try direct download with common naming pattern
    asm_name = accession
    fasta_url = (
        f"{base_url}/{asm_name}/{asm_name}_genomic.fna.gz"
    )

    try:
        resp = requests.get(fasta_url, stream=True, timeout=120)
        resp.raise_for_status()
    except requests.HTTPError:
        # Try to find the correct subdirectory name
        try:
            dir_resp = requests.get(base_url, timeout=60)
            dir_resp.raise_for_status()
            # Parse HTML directory listing for the assembly name
            matches = re.findall(
                rf'href="({re.escape(accession)}[^"]*)"',
                dir_resp.text,
            )
            if matches:
                asm_dir = matches[0].rstrip("/")
                fasta_url = f"{base_url}/{asm_dir}/{asm_dir}_genomic.fna.gz"
                resp = requests.get(fasta_url, stream=True, timeout=120)
                resp.raise_for_status()
            else:
                raise FileNotFoundError(
                    f"Assembly {accession} not found on NCBI FTP"
                )
        except Exception as exc:
            raise FileNotFoundError(
                f"Could not download assembly {accession}: {exc}"
            ) from exc

    # Save and decompress
    gz_path = out / f"{accession}_genomic.fna.gz"
    with open(gz_path, "wb") as fh:
        for chunk in resp.iter_content(chunk_size=8192):
            fh.write(chunk)

    with gzip.open(gz_path, "rt") as gz_in:
        with open(output_fasta, "w") as fasta_out:
            fasta_out.write(gz_in.read())

    gz_path.unlink(missing_ok=True)

    logger.info("Assembly saved: %s", output_fasta)
    return str(output_fasta)


def _biosample_to_assembly(biosample: str) -> str | None:
    """Resolve a BioSample accession to an assembly accession via Entrez."""
    try:
        # Search for assembly linked to this BioSample
        params: dict[str, Any] = {
            "db": "assembly",
            "term": biosample,
            "retmax": 1,
            "retmode": "json",
        }
        resp = requests.get(ESEARCH_URL, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()

        id_list = data.get("esearchresult", {}).get("idlist", [])
        if not id_list:
            return None

        # Get assembly accession from summary
        sum_params: dict[str, Any] = {
            "db": "assembly",
            "id": id_list[0],
            "retmode": "json",
        }
        sum_resp = requests.get(ESUMMARY_URL, params=sum_params, timeout=30)
        sum_resp.raise_for_status()
        summary = sum_resp.json()

        result = summary.get("result", {})
        uid = id_list[0]
        doc = result.get(uid, {})
        assembly_acc = doc.get("assemblyaccession")

        return assembly_acc if assembly_acc else None

    except Exception as exc:
        logger.warning("BioSample->Assembly resolution failed: %s", exc)
        return None


def _biosample_to_sra(biosample: str) -> str | None:
    """Resolve a BioSample accession to an SRA run accession."""
    try:
        params: dict[str, Any] = {
            "db": "sra",
            "term": biosample,
            "retmax": 1,
            "retmode": "json",
        }
        resp = requests.get(ESEARCH_URL, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()

        id_list = data.get("esearchresult", {}).get("idlist", [])
        if not id_list:
            return None

        # Fetch SRA record to get the run accession
        fetch_params: dict[str, Any] = {
            "db": "sra",
            "id": id_list[0],
            "rettype": "runinfo",
            "retmode": "text",
        }
        fetch_resp = requests.get(EFETCH_URL, params=fetch_params, timeout=30)
        fetch_resp.raise_for_status()

        # Parse CSV-like runinfo output
        lines = fetch_resp.text.strip().split("\n")
        if len(lines) >= 2:
            # First line is header, second line has the run accession
            run_acc = lines[1].split(",")[0]
            if _RE_SRA.match(run_acc):
                return run_acc

        return None

    except Exception as exc:
        logger.warning("BioSample->SRA resolution failed: %s", exc)
        return None


def _get_bioproject_runs(bioproject: str) -> list[str]:
    """Get all SRA run accessions linked to a BioProject."""
    try:
        params: dict[str, Any] = {
            "db": "sra",
            "term": f"{bioproject}[BioProject]",
            "retmax": 1000,
            "retmode": "json",
        }
        resp = requests.get(ESEARCH_URL, params=params, timeout=60)
        resp.raise_for_status()
        data = resp.json()

        id_list = data.get("esearchresult", {}).get("idlist", [])
        if not id_list:
            return []

        # Fetch run accessions in batch
        run_accessions: list[str] = []
        batch_size = 100
        for i in range(0, len(id_list), batch_size):
            batch = id_list[i : i + batch_size]
            fetch_params: dict[str, Any] = {
                "db": "sra",
                "id": ",".join(batch),
                "rettype": "runinfo",
                "retmode": "text",
            }
            fetch_resp = requests.get(EFETCH_URL, params=fetch_params, timeout=60)
            fetch_resp.raise_for_status()

            lines = fetch_resp.text.strip().split("\n")
            for line in lines[1:]:  # Skip header
                if not line.strip():
                    continue
                run_acc = line.split(",")[0]
                if _RE_SRA.match(run_acc):
                    run_accessions.append(run_acc)

        return run_accessions

    except Exception as exc:
        logger.error("BioProject run retrieval failed: %s", exc)
        return []
