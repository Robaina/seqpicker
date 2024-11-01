import logging
import sys

from .core import reduce_database_redundancy, filter_fasta_by_ids


# Configure logging when the package is imported
def _setup_logging():
    """Setup logging configuration."""
    # Get the root logger
    root_logger = logging.getLogger()

    # Remove any existing handlers to avoid duplicate messages
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    # Create and add a stream handler that writes to stdout
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter("%(asctime)s | %(levelname)s: %(message)s"))
    root_logger.addHandler(handler)
    root_logger.setLevel(logging.INFO)


# Run setup when package is imported
_setup_logging()
