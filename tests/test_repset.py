"""Tests for representative set selection algorithm."""

from pathlib import Path

import pytest
import numpy as np

from seqpicker.repset import (
    parse_identities,
    similarity_from_db,
    fraction_identity,
    create_facility_location_objective,
    create_redundancy_objective,
    MixtureObjective,
    select_representatives,
)


def test_parse_identities(sample_identity_matrix: Path):
    """Test parsing of pairwise identity matrix."""
    db = parse_identities(sample_identity_matrix)

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
    # Test fraction identity
    assert fraction_identity(-100, 75.0) == 0.75
    assert fraction_identity(-100, 100.0) == 1.0
    assert fraction_identity(-100, 0.0) == 0.0


def test_facility_location_objective(sample_identity_matrix: Path):
    """Test facility location objective function."""
    db = parse_identities(sample_identity_matrix)
    obj = create_facility_location_objective()

    # Test evaluation
    eval_result = obj.eval_func(db, ["seq1", "seq4"], fraction_identity)
    assert isinstance(eval_result, float)
    assert eval_result > 0

    # Test difference calculation
    data = obj.base_data_func(db, fraction_identity)
    diff = obj.diff_func(db, "seq1", fraction_identity, data)
    assert isinstance(diff, float)


def test_redundancy_objective(sample_identity_matrix: Path):
    """Test redundancy minimization objective function."""
    db = parse_identities(sample_identity_matrix)
    obj = create_redundancy_objective()

    # Test evaluation
    eval_result = obj.eval_func(db, ["seq1", "seq2"], fraction_identity)
    assert isinstance(eval_result, float)
    assert eval_result <= 0  # Should be negative due to penalties

    # Test difference calculation
    data = obj.base_data_func(db, fraction_identity)
    diff = obj.diff_func(db, "seq1", fraction_identity, data)
    assert isinstance(diff, float)


def test_mixture_objective(sample_identity_matrix: Path):
    """Test mixture of objectives."""
    facility_loc = create_facility_location_objective()
    redundancy = create_redundancy_objective()

    # Test valid initialization
    obj = MixtureObjective(objectives=[facility_loc, redundancy], weights=[0.7, 0.3])

    # Test invalid weights
    with pytest.raises(ValueError):
        MixtureObjective(objectives=[facility_loc, redundancy], weights=[0.7, 0.7])

    # Test evaluation
    db = parse_identities(sample_identity_matrix)
    sims = [fraction_identity, fraction_identity]
    eval_result = obj.eval(db, ["seq1", "seq4"], sims)
    assert isinstance(eval_result, float)


def test_select_representatives(sample_identity_matrix: Path):
    """Test representative sequence selection."""
    db = parse_identities(sample_identity_matrix)

    # Test with different max_size values
    for max_size in [3, 5]:
        representatives = select_representatives(
            db=db, mixture_weight=0.5, max_size=max_size
        )

        assert isinstance(representatives, list)
        assert len(representatives) <= max_size
        assert len(set(representatives)) == len(representatives)  # No duplicates
        assert all(r in db for r in representatives)  # All valid sequences

    # Test with different mixture weights
    for weight in [0.0, 0.5, 1.0]:
        representatives = select_representatives(
            db=db, mixture_weight=weight, max_size=3
        )
        assert isinstance(representatives, list)
        assert len(representatives) <= 3


def test_select_representatives_error_handling(sample_identity_matrix: Path):
    """Test error handling in representative selection."""
    db = parse_identities(sample_identity_matrix)

    # Test invalid mixture weight
    with pytest.raises(ValueError):
        select_representatives(db=db, mixture_weight=1.5)

    # Test invalid max_size
    with pytest.raises(ValueError):
        select_representatives(db=db, max_size=0)

    # Test empty database
    with pytest.raises(ValueError):
        select_representatives({})
