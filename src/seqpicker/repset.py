#!/bin/env python

"""
1. Compatible with Python3
2. Reads percent identity from els-alipid output file instead of calling psi-blast
3. Minimized version
"""

import argparse
import copy
import heapq
import logging
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, required=True, help="Output directory")
    parser.add_argument(
        "--seqs", type=Path, required=True, help="Input sequences in FASTA format"
    )
    parser.add_argument(
        "--pi", type=Path, required=True, help="Input text file with PI from esl-alipid"
    )
    parser.add_argument(
        "--mixture",
        type=float,
        default=0.5,
        help=(
            "Mixture parameter determining the relative weight of facility-location "
            "relative to sum-redundancy. Default=0.5"
        ),
    )
    parser.add_argument(
        "--size", type=float, default=float("inf"), help="Repset size. Default=inf"
    )
    args = parser.parse_args()
    workdir: Path = args.outdir

    if not 0.0 <= args.mixture <= 1.0:
        raise ValueError("Mixture parameter must be between 0.0 and 1.0")

    if not workdir.exists():
        workdir.mkdir(parents=True)

    # Logging setup
    logging.basicConfig(format="%(asctime)s %(levelname)s:%(message)s")
    logger = logging.getLogger("log")
    logger.setLevel(logging.DEBUG)

    # Create file handler which logs even debug messages
    fh = logging.FileHandler(workdir / "stdout.txt")
    fh.setLevel(logging.DEBUG)

    # Create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)

    # Create formatter and add it to the handlers
    formatter = logging.Formatter("%(asctime)s %(levelname)s - %(message)s")
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)

    # Add the handlers to logger
    logger.addHandler(ch)
    logger.addHandler(fh)
else:
    logger = logging.getLogger("log")


def get_pident(pi_file: Path) -> dict:
    """
    Parse esl-alipid output file.

    Args:
        pi_file (Path): Path to the esl-alipid output file.

    Returns:
        dict: A database containing sequence identities.
    """
    logger.info("Building dataframe")
    df = pd.read_csv(pi_file, sep=r"\s+", skiprows=1, header=None)
    df.columns = [
        "seqname1",
        "seqname2",
        "%id",
        "nid",
        "denomid",
        "%match",
        "nmatch",
        "denommatch",
    ]
    logger.info("Dataframe built")

    db = {}
    for _, row in df.iterrows():
        seq_id1 = row.seqname1.split("/")[0]
        seq_id2 = row.seqname2.split("/")[0]
        log10_e = -100.0
        pident = row["%id"]

        if seq_id1 not in db:
            db[seq_id1] = {"neighbors": {}, "in_neighbors": {}}
        if seq_id2 not in db:
            db[seq_id2] = {"neighbors": {}, "in_neighbors": {}}

        db[seq_id2]["neighbors"][seq_id1] = {
            "log10_e": log10_e,
            "pct_identical": pident,
        }
        db[seq_id1]["in_neighbors"][seq_id2] = {
            "log10_e": log10_e,
            "pct_identical": pident,
        }

    del df
    return db


###############################################################
###############################################################
# Submodular optimization functions
###############################################################
###############################################################

#############################
# Similarity functions
#############################


def sim_from_db(db, sim, seq_id1, seq_id2):
    d = db[seq_id1]["neighbors"][seq_id2]
    return sim_from_neighbor(sim, d)


def sim_from_neighbor(sim, d):
    return sim(d["log10_e"], d["pct_identical"])


def fraciden(log10_e: float, pct_identical: float) -> float:
    """
    Compute fractional identity, ignoring log10_e.

    Args:
        log10_e (float): Not used.
        pct_identical (float): Percentage identity.

    Returns:
        float: Fractional identity as a decimal.
    """
    return pct_identical / 100.0


def rankpropsim(log10_e, pct_identical):
    return np.exp(-np.power(10, log10_e) / 100.0)


def rankpropsim_loge(log10_e, pct_identical):
    return np.exp(-log10_e / 100.0)


def logloge(log10_e, pct_identical):
    if (-log10_e) <= 0.1:
        return 0.0
    elif (-log10_e) >= 1000:
        return 3.0
    else:
        return np.log10(-log10_e)


def oneprankpropsim(log10_e, pct_identical):
    return 1.0 + 1e-3 * rankpropsim_loge(log10_e, pct_identical)


def prodevaliden(log10_e, pct_identical):
    return fraciden(log10_e, pct_identical) * logloge(log10_e, pct_identical) / 3


def one(log10_e, pct_identical):
    return 1.0


def p90(log10_e, pct_identical):
    return float(pct_identical >= 0.9)


#############################
# Objective functions
# -----------------
# An objective is a dictionary
# {"eval": db, seq_ids, sim -> float, # The value of an objective function
#  "diff": db, seq_ids, sim, data -> float, # The difference in value after adding seq_id
#  "negdiff": db, seq_ids, sim, data -> float, # The difference in value after removing seq_id
#  "update": db, seq_ids, sim, data -> data, # Update the data structure to add seq_id as a representative
#  "negupdate": db, seq_ids, sim, data -> data, # Update the data structure to remove seq_id as a representative
#  "base_data": db, sim -> data, # The data structure corresponding to no representatives chosen
#  "full_data": db, sim -> data, # The data structure corresponding to all representatives chosen
#  "name": name}
# db: Database
# sim: Similiarity function
# data: function-specific data structure which may be modified over the course of an optimization algorithm
#############################

######################
# summaxacross
# AKA facility location
######################


def summaxacross_eval(db, seq_ids, sim):
    max_sim = {seq_id: 0 for seq_id in db}
    for chosen_seq_id in seq_ids:
        for neighbor_seq_id, d in db[chosen_seq_id]["in_neighbors"].items():
            if neighbor_seq_id in max_sim:
                sim_val = sim(d["log10_e"], d["pct_identical"])
                if sim_val > max_sim[neighbor_seq_id]:
                    max_sim[neighbor_seq_id] = sim_val
                # max_sim[neighbor_seq_id] = max(max_sim[neighbor_seq_id], sim(d["log10_e"], d["pct_identical"]))
            else:
                pass
                # raise Exception("Found node with neighbor not in set")

    return sum(max_sim.values())


# summaxacross data:
# Who each example is represented by
# "examples": {seq_id: (representive, val) }
# Who each represntative represents
# "representatives" {seq_id: {example: val}}


def summaxacross_base_data(db, sim):
    return {"examples": {seq_id: (None, 0) for seq_id in db}, "representatives": {}}


def summaxacross_full_data(db, sim):
    return {
        "examples": {
            seq_id: (seq_id, sim_from_db(db, sim, seq_id, seq_id)) for seq_id in db
        },
        "representatives": {
            seq_id: {seq_id: sim_from_db(db, sim, seq_id, seq_id)} for seq_id in db
        },
    }


def summaxacross_diff(db, seq_id, sim, data):
    diff = 0
    for neighbor_seq_id, d in db[seq_id]["in_neighbors"].items():
        if neighbor_seq_id in data["examples"]:
            sim_val = sim(d["log10_e"], d["pct_identical"])
            if sim_val > data["examples"][neighbor_seq_id][1]:
                diff += sim_val - data["examples"][neighbor_seq_id][1]
        else:
            pass
            # raise Exception("Found node with neighbor not in set")
    return diff


def summaxacross_update(db, seq_id, sim, data):
    data = copy.deepcopy(data)
    data["representatives"][seq_id] = {}
    for neighbor_seq_id, d in db[seq_id]["in_neighbors"].items():
        if neighbor_seq_id in data["examples"]:
            sim_val = sim(d["log10_e"], d["pct_identical"])
            if sim_val > data["examples"][neighbor_seq_id][1]:
                data["examples"][neighbor_seq_id] = (seq_id, sim_val)
                data["representatives"][seq_id][neighbor_seq_id] = sim_val
        else:
            pass
            # raise Exception("Found node with neighbor not in set")
    return data


# O(D^2)
def summaxacross_negdiff(db, seq_id, sim, data):
    diff = 0
    # For each neighbor_seq_id that was previously represented by seq_id
    new_representatives = set(data["representatives"].keys()) - set([seq_id])
    for neighbor_seq_id, d in data["representatives"][seq_id].items():
        # Find the next-best representative for neighbor_seq_id
        candidate_ids = (
            set(db[neighbor_seq_id]["neighbors"].keys()) & new_representatives
        )
        if len(candidate_ids) == 0:
            diff += -d
        else:
            best_id = max(
                candidate_ids, key=lambda x: sim_from_db(db, sim, neighbor_seq_id, x)
            )
            diff += sim_from_db(db, sim, neighbor_seq_id, best_id) - d
    return diff


# O(D^2)
def summaxacross_negupdate(db, seq_id, sim, data):
    data = copy.deepcopy(data)
    new_representatives = set(data["representatives"].keys()) - set([seq_id])
    for neighbor_seq_id, d in data["representatives"][seq_id].items():
        # Find the next-best representative for neighbor_seq_id
        candidate_ids = (
            set(db[neighbor_seq_id]["neighbors"].keys()) & new_representatives
        )
        if len(candidate_ids) == 0:
            data["examples"][neighbor_seq_id] = (None, 0)
        else:
            best_id = max(
                candidate_ids, key=lambda x: sim_from_db(db, sim, neighbor_seq_id, x)
            )
            data["examples"][neighbor_seq_id] = (
                best_id,
                sim_from_db(db, sim, neighbor_seq_id, best_id),
            )
            data["representatives"][best_id][neighbor_seq_id] = sim_from_db(
                db, sim, neighbor_seq_id, best_id
            )
    del data["representatives"][seq_id]
    return data


summaxacross = {
    "eval": summaxacross_eval,
    "diff": summaxacross_diff,
    "negdiff": summaxacross_negdiff,
    "update": summaxacross_update,
    "negupdate": summaxacross_negupdate,
    "base_data": summaxacross_base_data,
    "full_data": summaxacross_full_data,
    "name": "summaxacross",
}

######################
# summaxwithin
# AKA negfacloc
######################


def summaxwithin_eval(db, seq_ids, sim):
    max_sim = {seq_id: 0 for seq_id in seq_ids}
    for chosen_seq_id in seq_ids:
        for neighbor_seq_id, d in db[chosen_seq_id]["in_neighbors"].items():
            if neighbor_seq_id == chosen_seq_id:
                continue
            if neighbor_seq_id in max_sim:
                sim_val = sim(d["log10_e"], d["pct_identical"])
                if sim_val > max_sim[neighbor_seq_id]:
                    max_sim[neighbor_seq_id] = sim_val
            else:
                pass
    return -sum(max_sim.values())


# summaxwithin data:
# Who each example is represented by
# "examples": {seq_id: (representive, val) }
# Who each represntative represents
# "representatives" {seq_id: {example: val}}


def summaxwithin_base_data(db, sim):
    return {"examples": {seq_id: (None, 0) for seq_id in db}, "representatives": {}}


def summaxwithin_full_data(db, sim):
    data = {}
    data["examples"] = {}
    data["representatives"] = {seq_id: {} for seq_id in db}
    for seq_id in db:
        neighbors = {
            neighbor_seq_id: d
            for neighbor_seq_id, d in db[seq_id]["neighbors"].items()
            if neighbor_seq_id != seq_id
        }
        if len(neighbors) == 0:
            data["examples"][seq_id] = (None, 0)
        else:
            d = max(neighbors.items(), key=lambda d: sim_from_neighbor(sim, d[1]))
            data["examples"][seq_id] = (d[0], sim_from_neighbor(sim, d[1]))
            data["representatives"][d[0]][seq_id] = sim_from_neighbor(sim, d[1])
    return data


def summaxwithin_diff(db, seq_id, sim, data):
    diff = 0
    # Difference introduced in other representatives
    for neighbor_seq_id, d in db[seq_id]["in_neighbors"].items():
        if neighbor_seq_id == seq_id:
            continue
        if neighbor_seq_id in data["representatives"]:
            sim_val = sim(d["log10_e"], d["pct_identical"])
            if sim_val > data["examples"][neighbor_seq_id][1]:
                # adding a penalty of sim_val, removing old penalty
                diff -= sim_val - data["examples"][neighbor_seq_id][1]
    # Difference from adding this representative
    neighbors = {
        neighbor_seq_id: d
        for neighbor_seq_id, d in db[seq_id]["neighbors"].items()
        if neighbor_seq_id != seq_id
    }
    if len(neighbors) == 0:
        diff -= 0
    else:
        d = max(neighbors.items(), key=lambda d: sim_from_neighbor(sim, d[1]))
        diff -= sim_from_neighbor(sim, d[1])
    return diff


def summaxwithin_update(db, seq_id, sim, data):
    data = copy.deepcopy(data)
    # Find best repr for seq_id
    candidate_ids = (
        set(db[seq_id]["neighbors"].keys()) & set(data["representatives"].keys())
    ) - set([seq_id])
    if len(candidate_ids) == 0:
        data["examples"][seq_id] = (None, 0)
    else:
        best_id = max(candidate_ids, key=lambda x: sim_from_db(db, sim, seq_id, x))
        data["examples"][seq_id] = (best_id, sim_from_db(db, sim, seq_id, best_id))
        data["representatives"][best_id][seq_id] = sim_from_db(db, sim, seq_id, best_id)
    # Find ids represented by seq_id
    data["representatives"][seq_id] = {}
    for neighbor_seq_id, d in db[seq_id]["in_neighbors"].items():
        if neighbor_seq_id in data["examples"]:
            if neighbor_seq_id == seq_id:
                continue
            sim_val = sim(d["log10_e"], d["pct_identical"])
            if sim_val > data["examples"][neighbor_seq_id][1]:
                data["examples"][neighbor_seq_id] = (seq_id, sim_val)
                data["representatives"][seq_id][neighbor_seq_id] = sim_val
        else:
            pass
            # raise Exception("Found node with neighbor not in set")
    return data


# O(D^2)
def summaxwithin_negdiff(db, seq_id, sim, data):
    diff = 0
    # Difference introduced in other representatives
    # For each neighbor_seq_id that was previously represented by seq_id
    new_representatives = set(data["representatives"].keys()) - set([seq_id])
    for neighbor_seq_id, d in data["representatives"][seq_id].items():
        # Find the next-best representative for neighbor_seq_id
        candidate_ids = (
            set(db[neighbor_seq_id]["neighbors"].keys()) & new_representatives
        )
        if len(candidate_ids) == 0:
            diff += d  # removing a penalty of -d
        else:
            best_id = max(
                candidate_ids, key=lambda x: sim_from_db(db, sim, neighbor_seq_id, x)
            )
            # removing a penalty of d, adding a new penalty of -sim(neighbor, best)
            diff += d - sim_from_db(db, sim, neighbor_seq_id, best_id)
    # Difference from adding this representative
    diff += data["examples"][seq_id][1]  # removing a penalty of -sim
    return diff


# O(D^2)
def summaxwithin_negupdate(db, seq_id, sim, data):
    data = copy.deepcopy(data)
    del data["examples"][seq_id]
    new_representatives = set(data["representatives"].keys()) - set([seq_id])
    for neighbor_seq_id, d in data["representatives"][seq_id].items():
        # Find the next-best representative for neighbor_seq_id
        candidate_ids = (
            set(db[neighbor_seq_id]["neighbors"].keys()) & new_representatives
        )
        if len(candidate_ids) == 0:
            data["examples"][neighbor_seq_id] = (None, 0)
        else:
            best_id = max(
                candidate_ids, key=lambda x: sim_from_db(db, sim, neighbor_seq_id, x)
            )
            data["examples"][neighbor_seq_id] = (
                best_id,
                sim_from_db(db, sim, neighbor_seq_id, best_id),
            )
            data["representatives"][best_id][neighbor_seq_id] = sim_from_db(
                db, sim, neighbor_seq_id, best_id
            )
    del data["representatives"][seq_id]
    return data


summaxwithin = {
    "eval": summaxwithin_eval,
    "diff": summaxwithin_diff,
    "negdiff": summaxwithin_negdiff,
    "update": summaxwithin_update,
    "negupdate": summaxwithin_negupdate,
    "base_data": summaxwithin_base_data,
    "full_data": summaxwithin_full_data,
    "name": "summaxwithin",
}


######################
# sumsumwithin
######################


def bisim(db, sim, seq_id1, seq_id2):
    ret = 0
    if seq_id2 in db[seq_id1]["neighbors"]:
        d = db[seq_id1]["neighbors"][seq_id2]
        ret += sim(d["log10_e"], d["pct_identical"])
    if seq_id1 in db[seq_id2]["neighbors"]:
        d = db[seq_id2]["neighbors"][seq_id1]
        ret += sim(d["log10_e"], d["pct_identical"])
    return ret


def sumsumwithin_eval(db, seq_ids, sim):
    seq_ids = set(seq_ids)
    s = 0
    for chosen_id in seq_ids:
        for neighbor, d in db[chosen_id]["neighbors"].items():
            if chosen_id == neighbor:
                continue
            if neighbor in seq_ids:
                s += -sim(d["log10_e"], d["pct_identical"])
    return s


def sumsumwithin_base_data(db, sim):
    return set()


def sumsumwithin_full_data(db, sim):
    return set(db.keys())


def sumsumwithin_diff(db, seq_id, sim, data):
    diff = 0
    data = data | set([seq_id])
    for neighbor, d in db[seq_id]["neighbors"].items():
        if seq_id == neighbor:
            continue
        if neighbor not in data:
            continue
        diff += -sim_from_neighbor(sim, d)
        # neighbor_bisim = bisim(db, sim, seq_id, neighbor)
        # diff += -neighbor_bisim
    for neighbor, d in db[seq_id]["in_neighbors"].items():
        if seq_id == neighbor:
            continue
        if neighbor not in data:
            continue
        diff += -sim_from_neighbor(sim, d)
    return diff


def sumsumwithin_update(db, seq_id, sim, data):
    data.add(seq_id)
    return data


def sumsumwithin_negdiff(db, seq_id, sim, data):
    diff = 0
    # data = data - set([seq_id])
    for neighbor, d in db[seq_id]["neighbors"].items():
        if seq_id == neighbor:
            continue
        if neighbor not in data:
            continue
        # neighbor_bisim = bisim(db, sim, seq_id, neighbor)
        # diff -= -neighbor_bisim
        diff += sim_from_neighbor(sim, d)  # removing a penalty
    for neighbor, d in db[seq_id]["in_neighbors"].items():
        if seq_id == neighbor:
            continue
        if neighbor not in data:
            continue
        diff += sim_from_neighbor(sim, d)  # removing a penalty
    return diff


def sumsumwithin_negupdate(db, seq_id, sim, data):
    data.remove(seq_id)
    return data


sumsumwithin = {
    "eval": sumsumwithin_eval,
    "diff": sumsumwithin_diff,
    "negdiff": sumsumwithin_negdiff,
    "update": sumsumwithin_update,
    "negupdate": sumsumwithin_negupdate,
    "base_data": sumsumwithin_base_data,
    "full_data": sumsumwithin_full_data,
    "name": "sumsumwithin",
}

######################
# sumsumacross
######################


def sumsumacross_eval(db, seq_ids, sim):
    seq_ids = set(seq_ids)
    s = 0
    for chosen_id in seq_ids:
        for neighbor, d in db[chosen_id]["neighbors"].items():
            s += -sim(d["log10_e"], d["pct_identical"])
    return s


def sumsumacross_base_data(db, sim):
    return None


def sumsumacross_full_data(db, sim):
    return None


def sumsumacross_diff(db, seq_id, sim, data):
    diff = 0
    for neighbor, d in db[seq_id]["neighbors"].items():
        if seq_id == neighbor:
            continue
        diff += -sim(d["log10_e"], d["pct_identical"])
    return diff


def sumsumacross_negdiff(db, seq_id, sim, data):
    diff = 0
    for neighbor, d in db[seq_id]["neighbors"].items():
        if seq_id == neighbor:
            continue
        diff -= -sim(d["log10_e"], d["pct_identical"])
    return diff


def sumsumacross_update(db, seq_id, sim, data):
    # raise Exception("Not used")
    return None


def sumsumacross_negupdate(db, seq_id, sim, data):
    # raise Exception("Not used")
    return None


sumsumacross = {
    "eval": sumsumacross_eval,
    "diff": sumsumacross_diff,
    "negdiff": sumsumacross_negdiff,
    "update": sumsumacross_update,
    "negupdate": sumsumacross_negupdate,
    "base_data": sumsumacross_base_data,
    "full_data": sumsumacross_full_data,
    "name": "sumsumacross",
}

######################
# MixtureObjective
# ------------------------
# Create a mixture objective with:
# MixtureObjective([summaxacross, sumsumwithin], [0.1, 1.2])
# Must be used with a sim of the form
# [sim1, sim2]
# (same number of sims as objectives)
######################


class MixtureObjective(object):
    def __init__(self, objectives, weights):
        self.objectives = objectives
        self.weights = weights
        self.name = "mix-" + "-".join(
            [
                "{0}({1})".format(objective["name"], self.weights[i])
                for i, objective in enumerate(self.objectives)
            ]
        )

    def __getitem__(self, key):
        return self.__getattribute__(key)

    def __contains__(self, item):
        all_contain = True
        for i, objective in enumerate(self.objectives):
            all_contain = all_contain and (item in objective)
        return all_contain

    def eval(self, db, seq_ids, sims):
        s = 0
        for i, objective in enumerate(self.objectives):
            s += self.weights[i] * objective["eval"](db, seq_ids, sims[i])
        return s

    def diff(self, db, seq_id, sims, datas):
        s = 0
        for i, objective in enumerate(self.objectives):
            s += self.weights[i] * objective["diff"](db, seq_id, sims[i], datas[i])
        return s

    def negdiff(self, db, seq_id, sims, datas):
        s = 0
        for i, objective in enumerate(self.objectives):
            s += self.weights[i] * objective["negdiff"](db, seq_id, sims[i], datas[i])
        return s

    def update(self, db, seq_id, sims, datas):
        new_datas = []
        for i, objective in enumerate(self.objectives):
            new_datas.append(objective["update"](db, seq_id, sims[i], datas[i]))
        return new_datas

    def negupdate(self, db, seq_id, sims, datas):
        new_datas = []
        for i, objective in enumerate(self.objectives):
            new_datas.append(objective["negupdate"](db, seq_id, sims[i], datas[i]))
        return new_datas

    def base_data(self, db, sims):
        datas = []
        for i, objective in enumerate(self.objectives):
            datas.append(objective["base_data"](db, sims[i]))
        return datas

    def full_data(self, db, sims):
        datas = []
        for i, objective in enumerate(self.objectives):
            datas.append(objective["full_data"](db, sims[i]))
        return datas


#############################
# Optimization algorithms
# -------------------------
# Each returns either a specific
# subset or an order.
#############################


# returns an order
def accelerated_greedy_selection(
    db: dict,
    objective,
    sim,
    max_evals: float = float("inf"),
    diff_approx_ratio: float = 1.0,
    repset_size: float = float("inf"),
    target_obj_val: float = float("inf"),
) -> list:
    """
    Perform accelerated greedy selection algorithm.

    Args:
        db (dict): The database.
        objective: The objective function.
        sim: Similarity function.
        max_evals (float, optional): Maximum number of evaluations. Defaults to infinity.
        diff_approx_ratio (float, optional): Difference approximation ratio. Defaults to 1.0.
        repset_size (float, optional): Representative set size. Defaults to infinity.
        target_obj_val (float, optional): Target objective value. Defaults to infinity.

    Returns:
        list: The representative set order.
    """
    logger.info(f"Repset size: {repset_size}")

    assert diff_approx_ratio <= 1.0
    repset = []
    pq = [(-float("inf"), seq_id) for seq_id in db]
    objective_data = objective["base_data"](db, sim)
    cur_objective = 0
    num_evals = 0
    while (
        (len(repset) < repset_size)
        and (len(pq) > 1)
        and (cur_objective < target_obj_val)
    ):
        possible_diff, seq_id = heapq.heappop(pq)
        diff = objective["diff"](db, seq_id, sim, objective_data)
        next_diff = -pq[0][0]
        num_evals += 1
        if (num_evals >= max_evals) or (
            ((diff - next_diff) / (abs(diff) + 0.01)) >= (diff_approx_ratio - 1.0)
        ):
            repset.append(seq_id)
            objective_data = objective["update"](db, seq_id, sim, objective_data)
            cur_objective += diff
            num_evals = 0
        else:
            heapq.heappush(pq, (-diff, seq_id))
    if len(pq) == 1:
        repset.append(pq[0][1])
    return repset


def write_results_table(repset_order: list, output_file: Path) -> None:
    """
    Write repset results to file.

    Args:
        repset_order (list): List of sequence IDs in the representative set.
        output_file (Path): Path to the output file.
    """
    with open(output_file, "w") as f:
        for seq_id in repset_order:
            f.write(f"{seq_id}\n")


###############################################################
###############################################################
# Run optimization and output
###############################################################
###############################################################

if __name__ == "__main__":
    # db = run_psiblast(workdir, args.seqs)
    logger.info("Reading PI database...")
    db = get_pident(
        pi_file=args.pi
    )  # , seqs=args.seqs) # Make db with PI from esl-alipid
    logger.info("Finished building database...")
    objective = MixtureObjective(
        [summaxacross, sumsumwithin], [args.mixture, 1.0 - args.mixture]
    )
    logger.info(
        "Starting mixture of summaxacross and sumsumwithin with weight %s...",
        args.mixture,
    )
    sim = [fraciden, fraciden]
    repset_order = accelerated_greedy_selection(
        db, objective, sim, repset_size=args.size
    )  # Call main algorithm

    write_results_table(repset_order, workdir / "repset.txt")
