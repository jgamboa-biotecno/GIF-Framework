# GIF Framework — Genomic Intelligence Framework
# Copyright 2026 Javier Gamboa
# Licensed under the Apache License, Version 2.0
"""
WSL-aware tool execution for GIF-CLI.

On Windows, bioinformatics tools (blastn, mlst, chewBBACA, abricate, etc.)
run via WSL. On Linux/macOS they run natively. This module provides a
transparent wrapper that handles both cases.
"""

import os
import shutil
import subprocess
import logging
from typing import Tuple, Optional

logger = logging.getLogger("gif.utils.wsl")

IS_WINDOWS = os.name == "nt"
WSL_DISTRO = "Ubuntu-22.04"
# Try miniconda3 bioinfo env first, fallback to miniforge3
CONDA_ACTIVATE = (
    "source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null && "
    "conda activate bioinfo 2>/dev/null || "
    "source ~/miniforge3/bin/activate 2>/dev/null || true"
)


def windows_to_wsl_path(windows_path: str) -> str:
    """Convert a Windows path to a WSL-compatible /mnt/ path."""
    path = os.path.abspath(windows_path).replace("\\", "/")
    if len(path) >= 2 and path[1] == ":":
        drive = path[0].lower()
        path = f"/mnt/{drive}{path[2:]}"
    return path


def has_tool(name: str) -> bool:
    """Check if a tool is available (natively or via WSL)."""
    if shutil.which(name) is not None:
        return True
    if IS_WINDOWS:
        try:
            result = subprocess.run(
                ["wsl", "-d", WSL_DISTRO, "--", "bash", "-c",
                 f"{CONDA_ACTIVATE} && which {name}"],
                capture_output=True, text=True, timeout=10,
            )
            return result.returncode == 0
        except Exception:
            return False
    return False


def require_tool(name: str, install_hint: str) -> None:
    """Raise RuntimeError if a tool is not available anywhere."""
    if not has_tool(name):
        raise RuntimeError(
            f"Tool '{name}' not found on PATH"
            + (" or in WSL" if IS_WINDOWS else "")
            + f". Install via: {install_hint}"
        )


def run_tool(
    cmd: list[str],
    timeout: int = 600,
    convert_paths: bool = True,
) -> Tuple[str, str, int]:
    """
    Run an external bioinformatics tool, transparently via WSL on Windows.

    Parameters
    ----------
    cmd : list[str]
        Command and arguments (e.g. ["blastn", "-query", "input.fasta", ...]).
    timeout : int
        Timeout in seconds.
    convert_paths : bool
        If True and on Windows, convert any path-like arguments to WSL paths.

    Returns
    -------
    Tuple of (stdout, stderr, returncode).
    """
    tool_name = cmd[0]

    if shutil.which(tool_name) is not None:
        # Tool available natively — run directly
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return result.stdout, result.stderr, result.returncode

    if IS_WINDOWS:
        # Run via WSL
        if convert_paths:
            converted_cmd = []
            for arg in cmd:
                # Heuristic: if arg looks like a Windows path, convert it
                if (os.path.sep in arg or "/" in arg) and (
                    os.path.exists(arg) or ":" in arg
                ):
                    converted_cmd.append(windows_to_wsl_path(arg))
                else:
                    converted_cmd.append(arg)
            cmd = converted_cmd

        import shlex
        shell_cmd = " ".join(shlex.quote(a) for a in cmd)
        full_command = f"{CONDA_ACTIVATE} && {shell_cmd}"

        logger.debug("WSL command: %s", full_command)
        result = subprocess.run(
            ["wsl", "-d", WSL_DISTRO, "--", "bash", "-c", full_command],
            capture_output=True, text=True, timeout=timeout,
        )
        return result.stdout, result.stderr, result.returncode

    raise RuntimeError(
        f"Tool '{tool_name}' not found. Install via conda or use Docker."
    )
