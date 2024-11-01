"""Tests for representative set selection algorithm."""

from pathlib import Path

import pytest
import numpy as np

from seqpicker.repset import (
    get_pident,
    fraciden,
    summaxacross,
    sumsumwithin,
    MixtureObjective,
    accelerated_greedy_selection,
)


def test_parse_identities(sample_identity_matrix: Path):
    """Test parsing of pairwise identity matrix."""
    db = get_pident(sample_identity_matrix)

    # Check basic structure
    assert isinstance(db, dict)
    assert all(isinstance(v, dict) for v in db.values())
    assert all({"neighbors", "in_neighbors"} <= set(v.keys()) for v in db.values())

    # Check specific relationships
    assert "seq2" in db["seq1"]["neighbors"]
    assert db["seq1"]["neighbors"]["seq2"]["pct_identical"] == 100.0
    assert "seq1" in db["seq2"]["in_neighbors"]
    assert db["seq2"]["in_neighbors"]["seq1"]["pct_identical"] == 100.0


def test_similarity_functions():
    """Test similarity calculation functions."""
    # Test fraciden function
    assert fraciden(-100, 75.0) == 0.75
    assert fraciden(-100, 100.0) == 1.0
    assert fraciden(-100, 0.0) == 0.0


def test_facility_location_objective(sample_identity_matrix: Path):
    """Test facility location objective function."""
    db = get_pident(sample_identity_matrix)
    obj = summaxacross

    # Test evaluation
    eval_result = obj["eval"](db, ["seq1", "seq4"], fraciden)
    assert isinstance(eval_result, float)
    assert eval_result > 0

    # Test difference calculation
    data = obj["base_data"](db, fraciden)
    diff = obj["diff"](db, "seq1", fraciden, data)
    assert isinstance(diff, float)


def test_redundancy_objective(sample_identity_matrix: Path):
    """Test redundancy minimization objective function."""
    db = get_pident(sample_identity_matrix)
    obj = sumsumwithin

    # Test evaluation
    eval_result = obj["eval"](db, ["seq1", "seq2"], fraciden)
    assert isinstance(eval_result, float)
    assert eval_result <= 0  # Should be negative due to penalties

    # Test difference calculation
    data = obj["base_data"](db, fraciden)
    diff = obj["diff"](db, "seq1", fraciden, data)
    assert isinstance(diff, float)


def test_mixture_objective(sample_identity_matrix: Path):
    """Test mixture of objectives."""
    facility_loc = summaxacross
    redundancy = sumsumwithin

    # Test valid initialization
    obj = MixtureObjective(objectives=[facility_loc, redundancy], weights=[0.7, 0.3])

    # Test invalid weights
    with pytest.raises(ValueError):
        MixtureObjective(objectives=[facility_loc, redundancy], weights=[0.7, 0.7])

    # Test evaluation
    db = get_pident(sample_identity_matrix)
    sims = [fraciden, fraciden]
    eval_result = obj.eval(db, ["seq1", "seq4"], sims)
    assert isinstance(eval_result, float)


def test_select_representatives(sample_identity_matrix: Path):
    """Test representative sequence selection."""
    db = get_pident(sample_identity_matrix)

    # Test with different max_size values
    for max_size in [3, 5]:
        representatives = accelerated_greedy_selection(
            db=db,
            objective=MixtureObjective(
                objectives=[summaxacross, sumsumwithin], weights=[0.5, 0.5]
            ),
            sim=[fraciden, fraciden],
            repset_size=max_size,
        )

        assert isinstance(representatives, list)
        assert len(representatives) <= max_size
        assert len(set(representatives)) == len(representatives)  # No duplicates
        assert all(r in db for r in representatives)  # All valid sequences

    # Test with different mixture weights
    for weight in [0.0, 0.5, 1.0]:
        obj = MixtureObjective(
            objectives=[summaxacross, sumsumwithin], weights=[weight, 1.0 - weight]
        )
        representatives = accelerated_greedy_selection(
            db=db,
            objective=obj,
            sim=[fraciden, fraciden],
            repset_size=3,
        )
        assert isinstance(representatives, list)
        assert len(representatives) <= 3


def test_select_representatives_error_handling(sample_identity_matrix: Path):
    """Test error handling in representative selection."""
    db = get_pident(sample_identity_matrix)

    # Test invalid mixture weight
    with pytest.raises(ValueError):
        MixtureObjective(objectives=[summaxacross, sumsumwithin], weights=[1.5, -0.5])

    # Test invalid max_size
    with pytest.raises(AssertionError):
        accelerated_greedy_selection(
            db=db,
            objective=MixtureObjective(
                objectives=[summaxacross, sumsumwithin], weights=[0.5, 0.5]
            ),
            sim=[fraciden, fraciden],
            repset_size=0,
        )

    # Test empty database
    with pytest.raises(ValueError):
        accelerated_greedy_selection(
            db={},
            objective=MixtureObjective(
                objectives=[summaxacross, sumsumwithin], weights=[0.5, 0.5]
            ),
            sim=[fraciden, fraciden],
            repset_size=3,
        )
