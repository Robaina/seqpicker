"""Core functionality for reducing protein sequence datasets."""

from __future__ import annotations

import logging
import shutil
from pathlib import Path
from typing import Optional, Union

from seqpicker import repset
from seqpicker.utils import TemporaryDirectoryPath, TemporaryFilePath, run_command
from seqpicker.wrappers import run_cdhit, run_mafft


# Update the core.py code to use the new repset module location
def reduce_database_redundancy(
    input_fasta: Union[str, Path],
    output_fasta: Optional[Union[str, Path]] = None,
    cdhit: bool = True,
    maxsize: Optional[int] = None,
    cdhit_args: Optional[str] = None,
    mixture_weight: float = 0.5,
) -> Path:
    """Reduce redundancy of protein sequence database.

    Uses CD-HIT and/or facility location optimization to select representative sequences.

    Args:
        input_fasta: Path to input FASTA file
        output_fasta: Path to output reduced FASTA file
        cdhit: Whether to use CD-HIT for initial clustering
        maxsize: Maximum number of sequences in final database
        cdhit_args: Additional arguments to CD-HIT
        mixture_weight: Weight between facility location and redundancy objectives

    Returns:
        Path to output FASTA file containing representative sequences
    """
    input_fasta = Path(input_fasta).resolve()
    if output_fasta is None:
        output_fasta = input_fasta.with_stem(f"{input_fasta.stem}_reduced")
    else:
        output_fasta = Path(output_fasta).resolve()

    # Rest of the implementation remains the same, but uses the new repset module path
    ...
