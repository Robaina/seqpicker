"""Tests for utility functions."""

import os
from pathlib import Path

import pytest

from seqpicker.utils import TemporaryFilePath, TemporaryDirectoryPath, run_command


def test_temporary_file_path(temp_dir: Path):
    """Test TemporaryFilePath context manager."""
    # Test basic functionality
    with TemporaryFilePath() as temp_file:
        assert isinstance(temp_file, Path)
