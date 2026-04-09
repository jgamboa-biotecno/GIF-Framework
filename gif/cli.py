# GIF Framework — Genomic Intelligence Framework
# Copyright 2026 Javier Gamboa
# Licensed under the Apache License, Version 2.0
"""
GIF-CLI — Command-line interface for the Genomic Intelligence Framework.

Entry point for scoring *Listeria monocytogenes* genomes from FASTA files,
NCBI accessions, or BioProjects.
"""

from __future__ import annotations

import logging
import os
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import click
from rich.console import Console

console = Console(stderr=True)
GIF_VERSION = "1.0.0"


# ═══════════════════════════════════════════════════════════════════════════
# CLI Group
# ═══════════════════════════════════════════════════════════════════════════

@click.group()
@click.version_option(
    GIF_VERSION,
    prog_name="GIF Framework",
    message="%(prog)s v%(version)s",
)
def main():
    """GIF Framework \u2014 Genomic Intelligence Framework for Listeria monocytogenes risk assessment."""
    pass


# ═══════════════════════════════════════════════════════════════════════════
# score command
# ═══════════════════════════════════════════════════════════════════════════

@main.command()
@click.argument("inputs", nargs=-1, type=click.Path())
@click.option(
    "--output", "-o",
    type=click.Path(),
    default=".",
    help="Output directory for results.",
)
@click.option(
    "--accession",
    multiple=True,
    help="NCBI accession(s) to download and score.",
)
@click.option(
    "--accession-list",
    type=click.Path(exists=True),
    help="File containing NCBI accessions (one per line).",
)
@click.option(
    "--bioproject",
    help="NCBI BioProject accession to fetch all assemblies from.",
)
@click.option(
    "--context",
    type=click.Choice(["industrial", "clinical"]),
    default="industrial",
    help="Scoring context profile (default: industrial).",
)
@click.option(
    "--format", "fmt",
    type=click.Choice(["tsv", "json", "report", "all"]),
    default="tsv",
    help="Output format (default: tsv).",
)
@click.option(
    "--region",
    help="Region or country for R-Score regional AMR adjustment.",
)
@click.option(
    "--metadata",
    type=click.Path(exists=True),
    help="Metadata TSV for P-Score Level 2 temporal scoring.",
)
@click.option(
    "--tmp-dir",
    default="/tmp/gif-cli",
    help="Temporary directory for downloaded files.",
)
def score(
    inputs: Tuple[str, ...],
    output: str,
    accession: Tuple[str, ...],
    accession_list: Optional[str],
    bioproject: Optional[str],
    context: str,
    fmt: str,
    region: Optional[str],
    metadata: Optional[str],
    tmp_dir: str,
) -> None:
    """Calculate GIF scores for Listeria monocytogenes genomes.

    Accepts one or more FASTA files, directories, NCBI accessions, or a
    BioProject. Results are written to the output directory in the requested
    format.

    \b
    Examples:
      gif score genome.fasta
      gif score *.fasta -o results/ --format all
      gif score --accession GCA_000026945.2 --context clinical
      gif score --bioproject PRJNA422580 -o prjna_results/
      gif score genomes/ --metadata facility_history.tsv --region Spain
    """
    from gif.utils.validators import validate_fasta, find_fastas
    from gif.pipeline import process_batch, process_fasta
    from gif.output.summary import print_summary

    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        handlers=[logging.StreamHandler(sys.stderr)],
    )

    fasta_paths: List[str] = []

    # ── Collect FASTA paths from positional arguments ────────────────────
    for inp in inputs:
        p = Path(inp)
        if p.is_file():
            if validate_fasta(str(p)):
                fasta_paths.append(str(p.resolve()))
            else:
                console.print(f"[yellow]Skipping invalid FASTA: {inp}[/yellow]")
        elif p.is_dir():
            found = find_fastas(str(p))
            if found:
                fasta_paths.extend(found)
            else:
                console.print(f"[yellow]No FASTA files found in: {inp}[/yellow]")
        else:
            console.print(f"[red]Path not found: {inp}[/red]")

    # ── Resolve NCBI accessions ──────────────────────────────────────────
    all_accessions: List[str] = list(accession)

    if accession_list:
        with open(accession_list, "r") as fh:
            for line in fh:
                acc = line.strip()
                if acc and not acc.startswith("#"):
                    all_accessions.append(acc)

    if all_accessions:
        from gif.bio.fetch import resolve_accession, fetch_assembly
        os.makedirs(tmp_dir, exist_ok=True)
        for acc in all_accessions:
            try:
                console.print(f"  Resolving accession: {acc}")
                resolved = resolve_accession(acc)
                fasta_path = fetch_assembly(resolved, tmp_dir)
                if fasta_path and validate_fasta(fasta_path):
                    fasta_paths.append(fasta_path)
                else:
                    console.print(f"[red]Failed to fetch: {acc}[/red]")
            except Exception as exc:
                console.print(f"[red]Error resolving {acc}: {exc}[/red]")

    # ── Resolve BioProject ───────────────────────────────────────────────
    if bioproject:
        from gif.bio.fetch import fetch_bioproject

        try:
            os.makedirs(tmp_dir, exist_ok=True)
            bp_paths = fetch_bioproject(bioproject, tmp_dir)
            fasta_paths.extend(bp_paths)
        except Exception as exc:
            console.print(f"[red]BioProject fetch error: {exc}[/red]")

    # ── Validate we have inputs ──────────────────────────────────────────
    if not fasta_paths:
        console.print(
            "[red]No valid FASTA inputs found. Provide files, directories, "
            "accessions, or a BioProject.[/red]"
        )
        raise SystemExit(1)

    # Deduplicate while preserving order
    seen: set[str] = set()
    unique_paths: List[str] = []
    for p in fasta_paths:
        norm = os.path.normpath(p)
        if norm not in seen:
            seen.add(norm)
            unique_paths.append(p)

    console.print(f"\n  [bold]{len(unique_paths)} genome(s)[/bold] to score "
                  f"(context: {context})\n")

    # ── Single file: simple path ─────────────────────────────────────────
    if len(unique_paths) == 1 and fmt != "all":
        result = process_fasta(
            fasta_path=unique_paths[0],
            metadata=_load_single_metadata(metadata, unique_paths[0]),
            context=context,
            region=region,
        )

        if result.get("error"):
            console.print(f"[red]Error: {result['error']}[/red]")
            raise SystemExit(1)

        print_summary(result)

        # Write output
        output_dir = os.path.abspath(output)
        os.makedirs(output_dir, exist_ok=True)

        if fmt == "tsv":
            from gif.output.tsv_writer import write_tsv
            out_path = os.path.join(output_dir, "gif_results.tsv")
            write_tsv([result], out_path)
            console.print(f"  [green]Results saved to {out_path}[/green]")
        elif fmt == "json":
            from gif.output.json_writer import write_json
            out_path = os.path.join(output_dir, "gif_results.json")
            write_json([result], out_path)
            console.print(f"  [green]Results saved to {out_path}[/green]")
        elif fmt == "report":
            from gif.output.report_md import write_report
            sid = result.get("sample_id", "unknown")
            out_path = os.path.join(output_dir, f"{sid}_report.md")
            write_report(result, out_path)
            console.print(f"  [green]Report saved to {out_path}[/green]")

        return

    # ── Batch path ───────────────────────────────────────────────────────
    output_dir = os.path.abspath(output)
    process_batch(
        input_paths=unique_paths,
        output_dir=output_dir,
        metadata_file=metadata,
        context=context,
        fmt=fmt,
        region=region,
    )


# ═══════════════════════════════════════════════════════════════════════════
# info command
# ═══════════════════════════════════════════════════════════════════════════

@main.command()
@click.argument("query", required=False, default=None)
@click.option("--trophic", help="Show info about a trophic strategy.")
@click.option("--version", "show_version", is_flag=True, help="Show GIF version info.")
@click.option("--spec", is_flag=True, help="Show scoring specification summary.")
def info(
    query: Optional[str],
    trophic: Optional[str],
    show_version: bool,
    spec: bool,
) -> None:
    """Information about CC profiles, trophic strategies, and scoring specification.

    \b
    Examples:
      gif info CC1
      gif info --trophic Nosotroph
      gif info --spec
      gif info --version
    """
    from gif.scoring.weights import (
        CONTEXT_PROFILES,
        RISK_CRITICAL_THRESHOLD,
        RISK_HIGH_THRESHOLD,
        RISK_MODERATE_THRESHOLD,
    )

    if show_version:
        console.print(f"[bold]GIF Framework[/bold] v{GIF_VERSION}")
        return

    if spec:
        console.print("\n[bold underline]GIF Scoring Specification[/bold underline]\n")
        console.print("  Components: V-Score (Virulence), P-Score (Persistence), "
                      "C-Score (Clonality), R-Score (AMR)")
        console.print()
        for name, profile in CONTEXT_PROFILES.items():
            w = profile["weights"]
            desc = profile["description"]
            console.print(f"  [bold]{name}[/bold]: V={w['v']:.0%} P={w['p']:.0%} "
                          f"C={w['c']:.0%} R={w['r']:.0%}")
            console.print(f"    {desc}")
            if profile.get("validated_auc"):
                console.print(f"    Validated AUC: {profile['validated_auc']}")
            console.print()

        console.print("  [bold]Risk Tiers:[/bold]")
        console.print(f"    CRITICAL  >= {RISK_CRITICAL_THRESHOLD}")
        console.print(f"    HIGH      >= {RISK_HIGH_THRESHOLD}")
        console.print(f"    MODERATE  >= {RISK_MODERATE_THRESHOLD}")
        console.print(f"    LOW       <  {RISK_MODERATE_THRESHOLD}")
        console.print()
        return

    if trophic:
        _show_trophic_info(trophic)
        return

    if query:
        _show_cc_info(query)
        return

    # No arguments: show general help
    console.print(
        "\n[bold]GIF Framework v1.0[/bold] \u2014 Genomic Intelligence Framework\n"
    )
    console.print("  Use [bold]gif info CC1[/bold] for CC profile information")
    console.print("  Use [bold]gif info --trophic Nosotroph[/bold] for trophic strategy details")
    console.print("  Use [bold]gif info --spec[/bold] for scoring specification")
    console.print("  Use [bold]gif info --version[/bold] for version information")
    console.print()


def _show_trophic_info(strategy: str) -> None:
    """Display information about a trophic strategy."""
    strategies = {
        "nosotroph": {
            "name": "Nosotroph",
            "emoji": "\U0001f534",
            "description": (
                "Disease-adapted lineages with high clinical potential. "
                "Characterised by hypervirulent CCs (CC1, CC2, CC4, CC6) with "
                "intact InlA, full LIPI complement, and high V-Scores. "
                "Predominant in clinical listeriosis cases."
            ),
            "typical_ccs": "CC1, CC2, CC4, CC6",
            "typical_v_score": ">70",
        },
        "amphitroph": {
            "name": "Amphitroph",
            "emoji": "\U0001f536",
            "description": (
                "Dual-capability lineages that persist in food-processing "
                "environments AND retain clinical virulence. Characterised by "
                "moderate-to-high V-Score combined with persistence markers "
                "(SSI-1, qacH, bcrABC). Most concerning for food safety as they "
                "combine facility colonisation with outbreak potential."
            ),
            "typical_ccs": "CC9, CC121, CC8, CC3",
            "typical_v_score": "40-70 (with high P-Score)",
        },
        "saprotroph": {
            "name": "Saprotroph",
            "emoji": "\U0001f7e2",
            "description": (
                "Environmental lineages with limited clinical significance. "
                "Low virulence scores, often with truncated InlA. Typically "
                "isolated from natural environments rather than clinical cases."
            ),
            "typical_ccs": "Various low-frequency CCs",
            "typical_v_score": "<40",
        },
    }

    key = strategy.lower()
    if key not in strategies:
        console.print(f"[yellow]Unknown trophic strategy: {strategy}[/yellow]")
        console.print("Valid strategies: Nosotroph, Amphitroph, Saprotroph")
        return

    info = strategies[key]
    console.print(f"\n  {info['emoji']} [bold]{info['name']}[/bold]\n")
    console.print(f"  {info['description']}\n")
    console.print(f"  Typical CCs: {info['typical_ccs']}")
    console.print(f"  Typical V-Score: {info['typical_v_score']}")
    console.print()


def _show_cc_info(query: str) -> None:
    """Display information about a clonal complex."""
    # Normalise: "CC1" -> "1", "cc121" -> "121"
    cc_num = query.upper().replace("CC", "").strip()

    # Known CC profiles (representative subset)
    cc_profiles = {
        "1": ("Lineage I", "4b", "Nosotroph", "Hypervirulent clinical lineage"),
        "2": ("Lineage I", "4b", "Nosotroph", "Major outbreak lineage (e.g., 2011 cantaloupe)"),
        "3": ("Lineage I", "1/2b", "Amphitroph", "Clinical + environmental dual capability"),
        "4": ("Lineage I", "4b", "Nosotroph", "Hypervirulent, common in Southern Europe"),
        "5": ("Lineage I", "1/2b", "Amphitroph", "Moderate virulence with persistence"),
        "6": ("Lineage I", "4b", "Nosotroph", "Hypervirulent clinical lineage"),
        "7": ("Lineage II", "1/2a", "Amphitroph", "Environmental with clinical reports"),
        "8": ("Lineage II", "1/2a", "Amphitroph", "Food processing persistent lineage"),
        "9": ("Lineage II", "1/2c", "Amphitroph", "Major food-processing lineage, SSI-1+"),
        "11": ("Lineage II", "1/2a", "Amphitroph", "Persistent in dairy processing"),
        "14": ("Lineage II", "1/2a", "Saprotroph", "Environmental lineage"),
        "37": ("Lineage II", "1/2a", "Amphitroph", "Dairy and RTE lineage"),
        "87": ("Lineage II", "1/2a", "Amphitroph", "Food processing persistent lineage"),
        "101": ("Lineage II", "1/2a", "Saprotroph", "Environmental, rarely clinical"),
        "121": ("Lineage II", "1/2a", "Amphitroph", "Most persistent food-industry lineage, SSI-1+"),
        "155": ("Lineage II", "1/2a", "Saprotroph", "Environmental with truncated InlA"),
        "204": ("Lineage II", "1/2a", "Amphitroph", "Persistent in meat processing"),
        "321": ("Lineage II", "1/2a", "Saprotroph", "Rarely isolated environmental lineage"),
    }

    if cc_num in cc_profiles:
        lineage, serogroup, trophic, desc = cc_profiles[cc_num]
        console.print(f"\n  [bold]CC{cc_num}[/bold]")
        console.print(f"  Lineage:    {lineage}")
        console.print(f"  Serogroup:  {serogroup}")
        console.print(f"  Trophic:    {trophic}")
        console.print(f"  {desc}")
        console.print()
    else:
        console.print(f"\n  [yellow]CC{cc_num} not found in built-in profiles.[/yellow]")
        console.print()


# ═══════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════

def _load_single_metadata(
    metadata_path: Optional[str],
    fasta_path: str,
) -> Optional[dict]:
    """Load metadata for a single sample from a metadata TSV, if provided."""
    if not metadata_path:
        return None
    from gif.utils.validators import parse_metadata_tsv
    meta_map = parse_metadata_tsv(metadata_path)
    sample_id = Path(fasta_path).stem
    return meta_map.get(sample_id)


# ═══════════════════════════════════════════════════════════════════════════
# Entry point
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    main()
