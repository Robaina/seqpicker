"""Sequence reduction algorithm using facility location and redundancy minimization.

This module implements the representative set selection algorithm based on:
https://doi.org/10.1002/prot.25461

The algorithm uses facility location (summaxacross) and redundancy minimization 
(sumsumwithin) objectives to select a representative subset of sequences.
"""

from __future__ import annotations

import heapq
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import (
    Callable,
    Dict,
    List,
    Optional,
    TypeVar,
    Protocol,
    TypedDict,
    Union,
    cast,
)

import pandas as pd

# Type definitions
T = TypeVar("T")
SimScore = float
SeqId = str

logger = logging.getLogger(__name__)


class NeighborInfo(TypedDict):
    """Type definition for neighbor relationship information."""

    log10_e: float
    pct_identical: float


class NeighborDict(TypedDict):
    """Type definition for neighbor dictionaries."""

    neighbors: Dict[SeqId, NeighborInfo]
    in_neighbors: Dict[SeqId, NeighborInfo]


# Complete database type
Database = Dict[SeqId, NeighborDict]


class SimilarityFunction(Protocol):
    """Protocol defining similarity function interface."""

    def __call__(self, log10_e: float, pct_identical: float) -> SimScore: ...


@dataclass
class ObjectiveFunction:
    """Base class for submodular objective functions.

    Attributes:
        name: Identifier for the objective function
        eval_func: Function to evaluate objective on a set of sequences
        diff_func: Function to compute marginal gain of adding a sequence
        update_func: Function to update data structures after adding a sequence
        base_data_func: Function to initialize data structures
        full_data_func: Function to compute complete data representation
    """

    name: str
    eval_func: Callable[[Database, List[SeqId], SimilarityFunction], float]
    diff_func: Callable[[Database, SeqId, SimilarityFunction, T], float]
    update_func: Callable[[Database, SeqId, SimilarityFunction, T], T]
    base_data_func: Callable[[Database, SimilarityFunction], T]
    full_data_func: Callable[[Database, SimilarityFunction], T]


class MixtureObjective:
    """Mixture of multiple objective functions with weights.

    This class combines multiple objective functions with corresponding weights
    to create a single composite objective function.
    """

    def __init__(
        self,
        objectives: List[ObjectiveFunction],
        weights: List[float],
    ) -> None:
        """Initialize mixture objective.

        Args:
            objectives: List of objective functions to combine
            weights: List of weights for each objective

        Raises:
            ValueError: If lengths of objectives and weights don't match
        """
        if len(objectives) != len(weights):
            raise ValueError("Number of objectives must match number of weights")

        self.objectives = objectives
        self.weights = weights
        self.name = self._generate_name()

    def _generate_name(self) -> str:
        """Generate descriptive name for the mixture objective."""
        return "mix-" + "-".join(
            f"{obj.name}({weight})"
            for obj, weight in zip(self.objectives, self.weights)
        )

    def evaluate(
        self,
        db: Database,
        seq_ids: List[SeqId],
        sim_funcs: List[SimilarityFunction],
    ) -> float:
        """Evaluate the mixture objective.

        Args:
            db: Sequence database
            seq_ids: List of sequence IDs to evaluate
            sim_funcs: List of similarity functions

        Returns:
            Combined objective value
        """
        return sum(
            weight * obj.eval_func(db, seq_ids, sim)
            for obj, sim, weight in zip(self.objectives, sim_funcs, self.weights)
        )

    def compute_marginal_gain(
        self,
        db: Database,
        seq_id: SeqId,
        sim_funcs: List[SimilarityFunction],
        objective_data: List[T],
    ) -> float:
        """Compute marginal gain from adding a sequence.

        Args:
            db: Sequence database
            seq_id: ID of sequence to evaluate
            sim_funcs: List of similarity functions
            objective_data: Current objective function data

        Returns:
            Combined marginal gain value
        """
        return sum(
            weight * obj.diff_func(db, seq_id, sim, data)
            for obj, sim, weight, data in zip(
                self.objectives, sim_funcs, self.weights, objective_data
            )
        )

    def update_data(
        self,
        db: Database,
        seq_id: SeqId,
        sim_funcs: List[SimilarityFunction],
        objective_data: List[T],
    ) -> List[T]:
        """Update objective data after adding a sequence.

        Args:
            db: Sequence database
            seq_id: ID of sequence that was added
            sim_funcs: List of similarity functions
            objective_data: Current objective function data

        Returns:
            Updated objective data
        """
        return [
            obj.update_func(db, seq_id, sim, data)
            for obj, sim, data in zip(self.objectives, sim_funcs, objective_data)
        ]


def compute_similarity(log10_e: float, pct_identical: float) -> SimScore:
    """Compute similarity score from percent identity.

    Args:
        log10_e: Log10 of e-value (unused)
        pct_identical: Percent identity between sequences

    Returns:
        Normalized similarity score between 0 and 1
    """
    return float(pct_identical) / 100


def parse_pairwise_identities(pi_file: Union[str, Path]) -> Database:
    """Parse pairwise sequence identities from esl-alipid output.

    Args:
        pi_file: Path to esl-alipid output file

    Returns:
        Database dictionary containing pairwise identities

    Raises:
        FileNotFoundError: If input file doesn't exist
        pd.errors.EmptyDataError: If input file is empty
    """
    logger.info("Reading pairwise identity data")
    pi_path = Path(pi_file)

    if not pi_path.exists():
        raise FileNotFoundError(f"Identity file not found: {pi_file}")

    # Read and validate input data
    try:
        df = pd.read_csv(
            pi_path,
            sep=r"\s+",
            skiprows=1,
            header=None,
            names=[
                "seqname1",
                "seqname2",
                "%id",
                "nid",
                "denomid",
                "%match",
                "nmatch",
                "denommatch",
            ],
        )
    except pd.errors.EmptyDataError as e:
        raise pd.errors.EmptyDataError("Identity file is empty") from e

    # Build database structure
    db: Database = {}

    for _, row in df.iterrows():
        seq_id1 = row.seqname1.split("/")[0]
        seq_id2 = row.seqname2.split("/")[0]
        pident = row["%id"]

        # Initialize database entries
        for seq_id in (seq_id1, seq_id2):
            if seq_id not in db:
                db[seq_id] = {"neighbors": {}, "in_neighbors": {}}

        # Store bidirectional relationships
        neighbor_info: NeighborInfo = {
            "log10_e": -100,  # Placeholder e-value
            "pct_identical": pident,
        }

        db[seq_id2]["neighbors"][seq_id1] = neighbor_info.copy()
        db[seq_id1]["in_neighbors"][seq_id2] = neighbor_info.copy()

    return db


def select_representatives(
    db: Database,
    objective: MixtureObjective,
    sim_funcs: List[SimilarityFunction],
    max_size: Optional[int] = None,
    approx_ratio: float = 1.0,
    initial_priority: float = float("-inf"),
    gain_denominator_offset: float = 0.01,
    min_relative_gain: Optional[float] = None,
) -> List[SeqId]:
    """Select representative sequences using accelerated greedy selection.

    This function implements an accelerated greedy algorithm to select a subset
    of representative sequences that maximize the given objective function.

    Args:
        db: Database containing sequence relationships
        objective: Mixture objective function
        sim_funcs: List of similarity functions
        max_size: Maximum number of representatives to select
        approx_ratio: Approximation ratio for acceleration (1.0 = exact)
        initial_priority: Initial priority value for candidates in queue
        gain_denominator_offset: Small constant to prevent division by zero when
            computing relative gains (default: 0.01)
        min_relative_gain: Minimum relative gain ratio required for selection.
            If None, uses (approx_ratio - 1.0)

    Returns:
        List of selected representative sequence IDs

    Raises:
        ValueError: If database is empty or sim_funcs don't match objective
    """
    if not db:
        raise ValueError("Empty sequence database")

    if len(sim_funcs) != len(objective.objectives):
        raise ValueError(
            "Number of similarity functions must match number of objectives"
        )

    max_representatives = float("inf") if max_size is None else max_size
    min_relative_gain = (
        min_relative_gain if min_relative_gain is not None else (approx_ratio - 1.0)
    )

    logger.info(f"Selecting up to {max_representatives} representatives")

    # Initialize selection algorithm
    representatives: List[SeqId] = []
    priority_queue = [(initial_priority, seq_id) for seq_id in db]
    heapq.heapify(priority_queue)

    # Initialize objective tracking
    objective_data = [
        obj.base_data_func(db, sim) for obj, sim in zip(objective.objectives, sim_funcs)
    ]

    evaluations_since_last_selection = 0

    # Main selection loop
    while len(representatives) < max_representatives and len(priority_queue) > 1:
        possible_gain, seq_id = heapq.heappop(priority_queue)
        actual_gain = objective.compute_marginal_gain(
            db, seq_id, sim_funcs, objective_data
        )

        next_possible_gain = -priority_queue[0][0] if priority_queue else float("-inf")
        evaluations_since_last_selection += 1

        # Check if we should select this sequence
        relative_gain_ratio = (actual_gain - next_possible_gain) / (
            abs(actual_gain) + gain_denominator_offset
        )

        if relative_gain_ratio >= min_relative_gain:
            representatives.append(seq_id)
            objective_data = objective.update_data(
                db, seq_id, sim_funcs, objective_data
            )
            evaluations_since_last_selection = 0
        else:
            heapq.heappush(priority_queue, (-actual_gain, seq_id))

    # Add final sequence if needed and possible
    if len(priority_queue) == 1 and len(representatives) < max_representatives:
        representatives.append(priority_queue[0][1])

    return representatives


def write_representatives(
    representatives: List[SeqId],
    output_file: Union[str, Path],
) -> None:
    """Write selected representative sequence IDs to file.

    Args:
        representatives: List of representative sequence IDs
        output_file: Path to output file

    Raises:
        OSError: If unable to create output directory or write file
    """
    output_path = Path(output_file)

    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text("\n".join(representatives))
    except OSError as e:
        raise OSError(f"Failed to write representatives to {output_file}") from e
