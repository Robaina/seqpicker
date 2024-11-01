"""Test fixtures for seqpicker."""

from pathlib import Path

import pytest


@pytest.fixture
def test_data_dir() -> Path:
    """Return path to test data directory."""
    return Path(__file__).parent / "data"


@pytest.fixture
def sample_fasta(test_data_dir: Path) -> Path:
    """Create a sample FASTA file for testing.

    Contains 10 protein sequences with varying similarity.
    """
    fasta_content = """>seq1
MSLLPTPTVLPTLAPPLTFQPTLAEVLPKPVTVLTTLPVHLMKRVDQPVAPTLTPALAPK
>seq2
MSLLPTPTVLPTLAPPLTFQPTLAEVLPKPVTVLTTLPVHLMKRVDQPVAPTLTPALAPK
>seq3
MPLLPTPTVLPTLAPPLTFQPTLAEVLPKPVTVLTTLPVHLMKRVDQPVAPTLTPALAPK
>seq4
MKSINRTILLSLLSCFVLSQVIFQGQNLGFKQSSPLAFMFNKQPQNVIFSASFTTKTKSP
>seq5
MKSINRTILLSLLSCFVLSQVIFQGENLGFKQSSPLAFMFNKQPQNVIFSASFTTKTKSP
>seq6
MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNSSSGGKNGQGEPARVRCSHLL
>seq7
MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNSSSGGKNGQGEPARVRCSHLL
>seq8
MALLPTPTVLPTLAPPLTFQPTLAEVLPKPVTVLTTLPVHLMKRVDQPVAPTLTPALAPK
>seq9
MALLPTPTVLPTLAPALTFQPTLAEVLPKPVTVLTTLPVHLMKRVDQPVAPTLTPALAPK
>seq10
MKSINRTILLSLLSCFVLSQVVFQGQNLGFKQSSPLAFMFNKQPQNVIFSASFTTKTKSP"""

    fasta_file = test_data_dir / "sample.fasta"
    fasta_file.parent.mkdir(parents=True, exist_ok=True)
    fasta_file.write_text(fasta_content)
    return fasta_file


@pytest.fixture
def sample_identity_matrix(test_data_dir: Path) -> Path:
    """Create a sample pairwise identity matrix."""
    matrix_content = """# p1              p2              %id     nid   denomid  %match  nmatch  denommatch
seq1             seq2             100.0   62    62       100.0   62      62
seq1             seq3             98.4    61    62       98.4    61      62
seq1             seq8             96.8    60    62       96.8    60      62
seq2             seq3             98.4    61    62       98.4    61      62
seq3             seq8             95.2    59    62       95.2    59      62
seq4             seq5             98.4    61    62       98.4    61      62
seq4             seq10            96.8    60    62       96.8    60      62
seq6             seq7             100.0   62    62       100.0   62      62
seq8             seq9             98.4    61    62       98.4    61      62"""

    matrix_file = test_data_dir / "identities.txt"
    matrix_file.parent.mkdir(parents=True, exist_ok=True)
    matrix_file.write_text(matrix_content)
    return matrix_file


@pytest.fixture
def temp_dir(tmp_path: Path) -> Path:
    """Create a temporary directory for test outputs."""
    return tmp_path
