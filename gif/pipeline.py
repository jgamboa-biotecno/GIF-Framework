# GIF Framework — Genomic Intelligence Framework
# Copyright 2026 Javier Gamboa
# Licensed under the Apache License, Version 2.0
"""
GIF-CLI Pipeline Orchestrator.

Coordinates the full GIF scoring pipeline for one or many FASTA assemblies:
assembly QC, MLST typing, virulence/persistence/AMR detection, scoring, and
output generation.
"""

from __future__ import annotations

import logging
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

from rich.console import Console
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
)

from gif.bio.assembly import assess_assembly
from gif.bio.mlst import run_mlst, get_lineage_from_cc
from gif.bio.virulence import detect_virulence
from gif.bio.biocides import detect_persistence_markers
from gif.bio.resistance import detect_amr
from gif.scoring.integrator import calculate_gif_score, classify_trophic_strategy
from gif.output.tsv_writer import write_tsv
from gif.output.json_writer import write_json
from gif.output.summary import print_summary, print_batch_summary
from gif.output.report_md import write_report

logger = logging.getLogger("gif.pipeline")
console = Console(stderr=True)

GIF_VERSION = "1.0.0"

# Reference data directory (packaged with the CLI)
_REF_DIR = str(Path(__file__).parent / "db" / "reference_data")


# ═══════════════════════════════════════════════════════════════════════════
# Single FASTA processing
# ═══════════════════════════════════════════════════════════════════════════

def process_fasta(
    fasta_path: str,
    metadata: Optional[Dict[str, Any]] = None,
    context: str = "industrial",
    region: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Run the full GIF scoring pipeline on a single FASTA assembly.

    Args:
        fasta_path: Path to the FASTA assembly file.
        metadata: Optional metadata dict with temporal/epidemiological fields
                  for P-Score Level 2 (timespan_weeks, independent_detections,
                  unique_zones, total_zones).
        context: Scoring context profile (``"industrial"`` or ``"clinical"``).
        region: Optional region/country string for R-Score regional adjustment.

    Returns:
        Complete result dict with assembly QC, typing, component scores,
        GIF scores, trophic classification, and provenance metadata.
    """
    sample_id = Path(fasta_path).stem
    filename = Path(fasta_path).name

    result: Dict[str, Any] = {
        "sample_id": sample_id,
        "filename": filename,
        "gif_version": GIF_VERSION,
        "error": None,
    }

    try:
        # ── Step 1: Assembly QC ──────────────────────────────────────────
        logger.info("Assessing assembly: %s", fasta_path)
        assembly_result = assess_assembly(fasta_path)
        result["assembly"] = assembly_result

        # ── Step 2: MLST typing ─────────────────────────────────────────
        logger.info("Running MLST: %s", sample_id)
        try:
            mlst_result = run_mlst(fasta_path)
        except RuntimeError as e:
            logger.warning("MLST skipped (%s): %s", sample_id, e)
            mlst_result = {}
        cc = mlst_result.get("cc", "")
        st = mlst_result.get("st", "")
        lineage = get_lineage_from_cc(cc)  # returns str
        # Derive serogroup from lineage
        serogroup_map = {"Lineage_I": "IVb", "Lineage_II": "IIa", "Lineage_III": "IVa"}
        typing_result = {
            "species": mlst_result.get("species", "Listeria monocytogenes"),
            "lineage": lineage if lineage != "Unknown" else "",
            "serogroup": serogroup_map.get(lineage, ""),
            "st": st,
            "cc": cc,
        }
        result["typing"] = typing_result

        # ── Step 3: Virulence detection ─────────────────────────────────
        logger.info("Detecting virulence markers: %s", sample_id)
        try:
            vir_db = str(Path(_REF_DIR) / "virulence")
            virulence_result = detect_virulence(fasta_path, vir_db, cc=cc)
        except RuntimeError as e:
            logger.warning("Virulence detection skipped (%s): %s", sample_id, e)
            virulence_result = {}

        # ── Step 4: Persistence marker detection ────────────────────────
        logger.info("Detecting persistence markers: %s", sample_id)
        try:
            per_db = str(Path(_REF_DIR) / "persistence")
            persistence_result = detect_persistence_markers(fasta_path, per_db)
        except RuntimeError as e:
            logger.warning("Persistence detection skipped (%s): %s", sample_id, e)
            persistence_result = {}

        # ── Step 5: AMR detection ───────────────────────────────────────
        logger.info("Detecting AMR genes: %s", sample_id)
        try:
            amr_result = detect_amr(fasta_path)
        except RuntimeError as e:
            logger.warning("AMR detection skipped (%s): %s", sample_id, e)
            amr_result = {}

        # ── Step 6: Build markers dict ──────────────────────────────────
        markers = _build_markers(
            virulence_result, persistence_result, amr_result, mlst_result,
        )

        # ── Step 7: Calculate GIF score ─────────────────────────────────
        logger.info("Calculating GIF scores: %s", sample_id)
        gif_result = calculate_gif_score(
            markers=markers,
            metadata=metadata,
            context=context,
            region=region,
        )

        # Merge gif_result into result
        result["scores"] = {
            "v_score": gif_result.get("v_score", 0),
            "p_score": gif_result.get("p_score", 0),
            "c_score": gif_result.get("c_score", 0),
            "r_score": gif_result.get("r_score", 0),
            "gif_score": gif_result.get("gif_score", 0),
            "gif_score_industrial": gif_result.get("gif_score_industrial", 0),
            "gif_score_clinical": gif_result.get("gif_score_clinical", 0),
            "risk_tier": gif_result.get("risk_tier", "LOW"),
        }
        result["context"] = {
            "profile": gif_result.get("context_profile", context),
            "weights_applied": gif_result.get("weights_applied", {}),
        }
        result["details"] = {
            "virulence": gif_result.get("virulence_details", {}),
            "persistence": gif_result.get("persistence_details", {}),
            "clonality": gif_result.get("context_details", {}),
            "resistance": gif_result.get("resistance_details", {}),
        }
        result["trophic"] = {
            "strategy": gif_result.get("trophic_strategy", "Unassigned"),
            "description": gif_result.get("trophic_description", ""),
        }
        result["is_partial"] = gif_result.get("is_partial", False)
        result["failed_components"] = gif_result.get("failed_components", [])

    except Exception as exc:
        logger.error("Pipeline error for %s: %s", sample_id, exc, exc_info=True)
        result["error"] = str(exc)

    return result


# ═══════════════════════════════════════════════════════════════════════════
# Batch processing
# ═══════════════════════════════════════════════════════════════════════════

def process_batch(
    input_paths: List[str],
    output_dir: str,
    metadata_file: Optional[str] = None,
    context: str = "industrial",
    fmt: str = "tsv",
    region: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """
    Process multiple FASTA assemblies through the GIF pipeline.

    Args:
        input_paths: List of FASTA file paths.
        output_dir: Directory for output files.
        metadata_file: Optional path to a metadata TSV file.
        context: Scoring context (``"industrial"`` or ``"clinical"``).
        fmt: Output format: ``"tsv"``, ``"json"``, ``"report"``, or ``"all"``.
        region: Optional region for R-Score adjustment.

    Returns:
        List of result dicts (one per isolate).
    """
    from gif.utils.validators import parse_metadata_tsv

    os.makedirs(output_dir, exist_ok=True)

    # Load metadata if provided
    metadata_map: Dict[str, Dict[str, Any]] = {}
    if metadata_file:
        metadata_map = parse_metadata_tsv(metadata_file)

    n = len(input_paths)

    results: List[Dict[str, Any]] = []

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("Scoring genomes", total=n)

        for fasta_path in input_paths:
            sample_id = Path(fasta_path).stem
            sample_metadata = metadata_map.get(sample_id)

            progress.update(task, description=f"Scoring {sample_id}")

            result = process_fasta(
                fasta_path=fasta_path,
                metadata=sample_metadata,
                context=context,
                region=region,
            )

            if result.get("error"):
                console.print(
                    f"  [red]ERROR[/red] {sample_id}: {result['error']}"
                )
            else:
                # Print individual summary to terminal
                print_summary(result)

            results.append(result)
            progress.advance(task)

    # ── Write outputs ────────────────────────────────────────────────────
    successful = [r for r in results if not r.get("error")]

    if successful:
        formats = ["tsv", "json", "report"] if fmt == "all" else [fmt]

        if "tsv" in formats:
            tsv_path = os.path.join(output_dir, "gif_results.tsv")
            write_tsv(successful, tsv_path)
            console.print(f"  [green]TSV[/green]    {tsv_path}")

        if "json" in formats:
            json_path = os.path.join(output_dir, "gif_results.json")
            write_json(successful, json_path)
            console.print(f"  [green]JSON[/green]   {json_path}")

        if "report" in formats:
            report_dir = os.path.join(output_dir, "reports")
            os.makedirs(report_dir, exist_ok=True)
            for r in successful:
                sid = r.get("sample_id", "unknown")
                report_path = os.path.join(report_dir, f"{sid}_report.md")
                write_report(r, report_path)
            console.print(f"  [green]Reports[/green] {report_dir}/")

    # Batch summary
    print_batch_summary(results)

    return results


# ═══════════════════════════════════════════════════════════════════════════
# Internal helpers
# ═══════════════════════════════════════════════════════════════════════════

def _build_markers(
    virulence: Dict[str, Any],
    persistence: Dict[str, Any],
    amr: Dict[str, Any],
    mlst: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Merge detection results into a unified markers dict consumed by the
    scoring integrator.

    Keys must match what v_score, p_score, c_score, r_score expect.
    """
    markers: Dict[str, Any] = {}

    # ── MLST / typing ──────────────────────────────────────────────────
    markers["clonal_complex"] = mlst.get("cc", "")
    markers["sequence_type"] = mlst.get("st", "")

    # ── Virulence (v_score expects) ────────────────────────────────────
    markers["inla_status"] = virulence.get("inla_status", "unknown")
    markers["inla_protein_length"] = virulence.get("inla_protein_length")
    markers["lipi1_intact"] = virulence.get("lipi1_intact", False)
    markers["lipi1_deleted_genes"] = virulence.get("lipi1_deleted_genes", [])
    markers["lipi3_complete"] = virulence.get("lipi3_complete", False)
    markers["lipi3_genes_detected"] = len(virulence.get("lipi3_genes_detected", []))
    markers["lipi4_complete"] = virulence.get("lipi4_complete", False)
    markers["lipi4_genes_detected"] = len(virulence.get("lipi4_genes_detected", []))

    # ── Persistence (p_score expects) ──────────────────────────────────
    qac_genes = persistence.get("qac_genes", [])
    markers["qac_e_delta1_detected"] = any(
        "qacedelta1" in g.lower() or "qacedelta1" in g.lower() for g in qac_genes
    )
    markers["qac_h_detected"] = any("qach" in g.lower() for g in qac_genes)
    markers["bcr_abc_detected"] = persistence.get("bcr_complete", False)
    markers["bcr_genes_count"] = len(persistence.get("bcr_genes", []))

    ssi1_genes = persistence.get("ssi1_genes", [])
    ssi2_genes = persistence.get("ssi2_genes", [])
    gad_genes = persistence.get("gad_genes", [])
    markers["ssi1_genes_count"] = len(ssi1_genes)
    markers["ssi2_genes_count"] = len(ssi2_genes)
    markers["gad_genes_count"] = len(gad_genes)
    markers["ssi1_detected"] = persistence.get("ssi1_complete", False)
    markers["ssi2_detected"] = persistence.get("ssi2_complete", False)

    markers["cadmium_resistance_detected"] = persistence.get("cadmium_resistance", False)
    markers["cadmium_genes"] = persistence.get("cadmium_genes", [])

    crispr_arrays = persistence.get("crispr_arrays", 0)
    markers["crispr_detected"] = crispr_arrays > 0
    markers["prophages_count"] = persistence.get("prophage_count", 0)

    # ── AMR (r_score expects) ──────────────────────────────────────────
    by_cat = amr.get("by_category", {})
    markers["beta_lactam_genes"] = by_cat.get("beta-lactam", [])
    markers["tetracycline_genes"] = by_cat.get("tetracycline", [])
    markers["aminoglycoside_genes"] = by_cat.get("aminoglycoside", [])
    markers["aminoglycoside_detected"] = len(markers["aminoglycoside_genes"]) > 0
    markers["fluoroquinolone_mutations"] = by_cat.get("fluoroquinolone", [])
    markers["macrolide_genes"] = by_cat.get("macrolide", [])
    markers["phenicol_genes"] = by_cat.get("phenicol", [])
    markers["lincosamide_genes"] = by_cat.get("lincosamide", [])
    markers["trimethoprim_genes"] = by_cat.get("trimethoprim", [])
    markers["has_multidrug_plasmid"] = amr.get("mdr_plasmid", False)
    markers["plasmid_amr_genes_count"] = amr.get("mdr_plasmid_genes", 0)

    return markers
