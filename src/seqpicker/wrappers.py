"""Wrappers for external tools used in sequence database reduction."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional, Union

from seqpicker.utils import run_command


def run_cdhit(
    input_fasta: Union[str, Path],
    output_fasta: Union[str, Path],
    additional_args: Optional[str] = None,
) -> None:
    """Run CD-HIT to obtain representative sequences.

    Args:
        input_fasta: Path to input FASTA file
        output_fasta: Path to output file
        additional_args: Additional arguments to CD-HIT
    """
    cmd = f"cd-hit -i {input_fasta} -o {output_fasta}"
    if additional_args:
        cmd += f" {additional_args}"
    run_command(cmd, suppress_output=True)


def run_mafft(
    input_fasta: Union[str, Path],
    output_file: Union[str, Path],
    processes: int = -1,
    parallel: bool = True,
    additional_args: Optional[str] = None,
) -> None:
    """Run MAFFT for multiple sequence alignment.

    Args:
        input_fasta: Path to input FASTA file
        output_file: Path to output file
        processes: Number of processes to use (-1 for all available)
        parallel: Whether to use parallel processing
        additional_args: Additional arguments to MAFFT
    """
    if parallel:
        thread_str = f"--thread {processes}"
    else:
        thread_str = ""

    cmd = f"mafft {thread_str}"
    if additional_args:
        cmd += f" {additional_args}"
    cmd += f" {input_fasta} > {output_file}"

    run_command(cmd, suppress_output=False)
