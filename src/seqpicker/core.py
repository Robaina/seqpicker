"""Core functionality for reducing protein sequence datasets."""

from __future__ import annotations

import logging
import shutil
from pathlib import Path
from typing import Optional, Union

from seqpicker import repset
from seqpicker.utils import TemporaryFilePath, run_command
from seqpicker.wrappers import run_cdhit, run_mafft

logger = logging.getLogger(__name__)

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

    if not cdhit and maxsize is None:
        logger.info("No reduction algorithm selected, copying input to output")
        shutil.copy(input_fasta, output_fasta)
        return output_fasta

    with TemporaryFilePath() as tempaln, \
         TemporaryFilePath() as tempfasta, \
         TemporaryFilePath() as tempfasta2, \
         TemporaryFilePath() as tempident:

        current_fasta = tempfasta
        shutil.copy(input_fasta, current_fasta)

        if cdhit:
            logger.info("Running CD-HIT for initial redundancy reduction")
            run_cdhit(
                input_fasta=current_fasta,
                output_fasta=tempfasta2,
                additional_args=cdhit_args,
            )
            Path(f"{tempfasta2}.clstr").unlink(missing_ok=True)
            shutil.move(tempfasta2, current_fasta)

        if maxsize is not None:
            logger.info("Running RepSet to select representative sequences")
            # Generate alignment
            run_mafft(
                input_fasta=current_fasta,
                output_file=tempaln,
                processes=-1,
                parallel=True,
                additional_args="--retree 1 --maxiterate 0 --quiet",
            )

            # Get pairwise identities
            run_command(
                f"esl-alipid {tempaln} > {tempident}",
                suppress_output=False
            )

            # Build sequence database using original repset method
            logger.info("Reading identity matrix...")
            db = repset.get_pident(tempident)
            
            # Create mixture objective with original functions
            objective = repset.MixtureObjective(
                [repset.summaxacross, repset.sumsumwithin],
                [mixture_weight, 1.0 - mixture_weight]
            )
            
            # Setup similarity functions
            sim = [repset.fraciden, repset.fraciden]
            
            # Select representatives using original selection function
            logger.info("Selecting representative sequences...")
            representatives = repset.accelerated_greedy_selection(
                db=db,
                objective=objective,
                sim=sim,
                repset_size=maxsize
            )

            # Filter FASTA to just representatives
            filter_fasta_by_ids(
                input_fasta=current_fasta,
                record_ids=representatives,
                output_fasta=tempfasta2
            )
            shutil.move(tempfasta2, current_fasta)

        shutil.move(current_fasta, output_fasta)
        logger.info(f"Successfully wrote representative sequences to: {output_fasta}")

def filter_fasta_by_ids(
    input_fasta: Union[str, Path],
    record_ids: list[str],
    output_fasta: Optional[Union[str, Path]] = None,
) -> Path:
    """Filter records in FASTA file matching provided IDs.

    Args:
        input_fasta: Path to input FASTA file
        record_ids: List of record IDs to filter FASTA file by
        output_fasta: Optional path to output file. If None, appends '_filtered' to input name

    Returns:
        Path to the filtered output FASTA file
    """
    input_fasta = Path(input_fasta).resolve()
    if output_fasta is None:
        output_fasta = input_fasta.with_stem(f"{input_fasta.stem}_filtered")
    else:
        output_fasta = Path(output_fasta).resolve()

    # Remove any existing index file
    index_file = Path(str(input_fasta) + ".fxi")
    index_file.unlink(missing_ok=True)

    # Write IDs to temporary file for seqkit
    record_ids = set(record_ids)
    with TemporaryFilePath() as tmp_ids:
        tmp_ids.write_text("\n".join(record_ids))
        run_command(
            f"seqkit grep -i -f {tmp_ids} {input_fasta} -o {output_fasta}",
            suppress_output=True,
        )
    
    return output_fasta
