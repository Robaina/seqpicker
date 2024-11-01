"""Tests for utility functions."""

from pathlib import Path

from seqpicker.utils import TemporaryFilePath


def test_temporary_file_path():
    """Test TemporaryFilePath context manager."""
    # Test basic functionality
    with TemporaryFilePath() as temp_file:
        assert isinstance(temp_file, Path)
        temp_file.write_text("Test content")
        assert temp_file.exists()
    # Verify that the file is deleted after exiting the context
    assert not temp_file.exists()
