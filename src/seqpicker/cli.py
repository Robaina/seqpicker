"""Command line interface for seqpicker package."""

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional

from seqpicker.core import reduce_database_redundancy

logger = logging.getLogger(__name__)


def setup_logging(verbose: bool = False) -> None:
    """Adjust logging level for verbose mode."""
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)


def parse_args(args: Optional[list[str]] = None) -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Select representative protein sequences from large datasets",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "input", type=Path, help="Input FASTA file containing protein sequences"
    )

    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output FASTA file for representative sequences",
    )

    parser.add_argument(
        "--maxsize", type=int, help="Maximum number of sequences in final dataset"
    )

    # Method selection
    method_group = parser.add_argument_group("Method Selection")
    method_exclusive = method_group.add_mutually_exclusive_group()

    method_exclusive.add_argument(
        "--cdhit-only",
        action="store_true",
        help="Use only CD-HIT for redundancy reduction",
    )

    method_exclusive.add_argument(
        "--repset-only", action="store_true", help="Use only RepSet selection algorithm"
    )

    # CD-HIT options
    cdhit_group = parser.add_argument_group("CD-HIT Options")
    cdhit_group.add_argument(
        "--similarity",
        type=float,
        default=0.9,
        help="CD-HIT sequence similarity threshold",
    )

    cdhit_group.add_argument(
        "--cdhit-args", type=str, help="Additional arguments to pass to CD-HIT"
    )

    # RepSet options
    repset_group = parser.add_argument_group("RepSet Options")
    repset_group.add_argument(
        "--mixture-weight",
        type=float,
        default=0.5,
        help="Weight between facility location (1.0) and redundancy (0.0)",
    )

    # Other options
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose output"
    )

    parsed_args = parser.parse_args(args)

    # Validate arguments
    if parsed_args.similarity < 0 or parsed_args.similarity > 1:
        parser.error("--similarity must be between 0 and 1")

    if parsed_args.mixture_weight < 0 or parsed_args.mixture_weight > 1:
        parser.error("--mixture-weight must be between 0 and 1")

    return parsed_args


def main(args: Optional[list[str]] = None) -> int:
    """Main entry point for the seqpicker command."""
    parsed_args = parse_args(args)
    setup_logging(parsed_args.verbose)

    try:
        # Determine which methods to use
        use_cdhit = not parsed_args.repset_only
        use_repset = not parsed_args.cdhit_only

        if not use_cdhit and not use_repset:
            logger.error("Must use at least one reduction method")
            return 1

        # Build CD-HIT args if needed
        cdhit_args = parsed_args.cdhit_args or ""
        if use_cdhit and "-c" not in cdhit_args:
            cdhit_args = f"-c {parsed_args.similarity} {cdhit_args}"

        # Run reduction
        output_file = reduce_database_redundancy(
            input_fasta=parsed_args.input,
            output_fasta=parsed_args.output,
            cdhit=use_cdhit,
            maxsize=parsed_args.maxsize if use_repset else None,
            cdhit_args=cdhit_args if use_cdhit else None,
            mixture_weight=parsed_args.mixture_weight if use_repset else 0.5,
        )

        logger.info(f"Successfully wrote representative sequences to: {output_file}")
        return 0

    except Exception as e:
        logger.error(f"Error: {str(e)}")
        if parsed_args.verbose:
            logger.exception("Detailed error trace:")
        return 1


if __name__ == "__main__":
    sys.exit(main())
