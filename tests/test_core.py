"""Tests for core sequence reduction functionality."""

from pathlib import Path

import pytest

from seqpicker.core import reduce_database_redundancy, filter_fasta_by_ids
from seqpicker.utils import TemporaryFilePath


def test_filter_fasta_by_ids(sample_fasta: Path, temp_dir: Path):
    """Test filtering FASTA file by sequence IDs."""
    output_file = temp_dir / "filtered.fasta"
    record_ids = ["seq1", "seq4", "seq6"]

    result = filter_fasta_by_ids(
        input_fasta=sample_fasta, record_ids=record_ids, output_fasta=output_file
    )

    assert result == output_file
    assert result.exists()

    # Check content
    content = result.read_text()
    assert all(f">{seq_id}" in content for seq_id in record_ids)
    assert len(content.split(">")) - 1 == len(record_ids)


def test_reduce_database_redundancy_cdhit_only(sample_fasta: Path, temp_dir: Path):
    """Test database reduction using only CD-HIT."""
    output_file = temp_dir / "reduced.fasta"

    result = reduce_database_redundancy(
        input_fasta=sample_fasta,
        output_fasta=output_file,
        cdhit=True,
        maxsize=None,
        cdhit_args="-c 0.9",  # 90% identity threshold
    )

    assert result == output_file
    assert result.exists()

    # Should reduce similar sequences
    content = result.read_text()
    seq_count = len(content.split(">")) - 1
    assert seq_count < 10  # Original has 10 sequences


def test_reduce_database_redundancy_repset_only(
    sample_fasta: Path, sample_identity_matrix: Path, temp_dir: Path
):
    """Test database reduction using only RepSet."""
    output_file = temp_dir / "reduced.fasta"

    result = reduce_database_redundancy(
        input_fasta=sample_fasta,
        output_fasta=output_file,
        cdhit=False,
        maxsize=5,
        mixture_weight=0.5,
    )

    assert result == output_file
    assert result.exists()

    # Should respect maxsize
    content = result.read_text()
    seq_count = len(content.split(">")) - 1
    assert seq_count <= 5


def test_reduce_database_redundancy_combined(
    sample_fasta: Path, sample_identity_matrix: Path, temp_dir: Path
):
    """Test database reduction using both CD-HIT and RepSet."""
    output_file = temp_dir / "reduced.fasta"

    result = reduce_database_redundancy(
        input_fasta=sample_fasta,
        output_fasta=output_file,
        cdhit=True,
        maxsize=3,
        cdhit_args="-c 0.9",
        mixture_weight=0.5,
    )

    assert result == output_file
    assert result.exists()

    # Check final sequence count
    content = result.read_text()
    seq_count = len(content.split(">")) - 1
    assert seq_count <= 3


def test_reduce_database_redundancy_error_handling(temp_dir: Path):
    """Test error handling in database reduction."""
    # Test with non-existent file
    with pytest.raises(FileNotFoundError):
        reduce_database_redundancy(
            input_fasta="nonexistent.fasta", output_fasta=temp_dir / "out.fasta"
        )

    # Test with invalid mixture weight
    with pytest.raises(ValueError):
        reduce_database_redundancy(
            input_fasta=sample_fasta,
            output_fasta=temp_dir / "out.fasta",
            mixture_weight=2.0,
        )


def test_reduce_database_redundancy_preserves_sequences(
    sample_fasta: Path, temp_dir: Path
):
    """Test that sequence content is preserved during reduction."""
    output_file = temp_dir / "reduced.fasta"

    result = reduce_database_redundancy(
        input_fasta=sample_fasta, output_fasta=output_file, cdhit=True, maxsize=5
    )

    # Check that sequences in output are unchanged from input
    original_content = sample_fasta.read_text()
    reduced_content = result.read_text()

    for record in reduced_content.split(">")[1:]:  # Skip empty first split
        header, sequence = record.split("\n", 1)
        seq_id = header.split()[0]

        # Find this sequence in original
        original_seq = None
        for orig_record in original_content.split(">")[1:]:
            orig_header, orig_sequence = orig_record.split("\n", 1)
            if orig_header.split()[0] == seq_id:
                original_seq = orig_sequence.strip()
                break

        assert original_seq == sequence.strip()
