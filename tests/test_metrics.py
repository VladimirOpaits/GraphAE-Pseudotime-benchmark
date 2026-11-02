import os

import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

from bench.data import DATA, clone_pairs, load_larry
from bench.metrics import sister_order

PAIRS = (np.array([0, 2, 4]), np.array([1, 3, 5]), np.array([0, 1, 2]))
BRANCH = np.array(["A"] * 6, dtype=object)


def acc(pt, min_pairs=1):
    res = sister_order(np.asarray(pt, dtype=float), BRANCH, PAIRS, n_boot=50, min_pairs=min_pairs)
    return res.set_index("branch")


def test_later_sister_must_have_larger_pseudotime():
    assert acc(np.arange(6)).loc["A", "acc"] == 1.0
    assert acc(-np.arange(6)).loc["A", "acc"] == 0.0
    assert acc(np.zeros(6)).loc["A", "acc"] == 0.5


def test_small_branch_excluded_from_mean():
    res = acc(np.arange(6), min_pairs=1000)
    assert not res.loc["A", "used"]
    assert "mean" not in res.index


def test_clone_pairs_are_cross_day_and_early_first():
    day = np.array([2, 4, 4, 6, 2])
    rows, cols = [0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 1, 3]
    clones = csr_matrix((np.ones(6), (cols, rows)), shape=(5, 2))
    i, j, _ = clone_pairs(day, clones)
    assert len(i) == 5
    assert np.all(day[i] < day[j])


@pytest.mark.skipif(not os.path.exists(DATA + "stateFate_inVitro_metadata.txt.gz"), reason="no LARRY data")
def test_larry_reference_values():
    meta, clones = load_larry()
    pairs = clone_pairs(meta["Time point"].to_numpy(), clones)
    ps = pd.read_csv(DATA + "stateFate_inVitro_neutrophil_pseudotime.txt.gz", sep="\t")
    pt = np.full(len(meta), np.nan)
    pt[ps["Cell index"].to_numpy()] = ps["pseudotime"].to_numpy()
    branch = np.where(np.isnan(pt), None, "Neutrophil")
    authors = sister_order(pt, branch, pairs).set_index("branch").loc["mean"]
    rnd = sister_order(np.random.default_rng(0).random(len(meta)), branch, pairs).set_index("branch").loc["mean"]
    assert 0.92 <= authors["acc"] <= 0.96
    assert authors["violation"] < 0.02
    assert 0.47 <= rnd["acc"] <= 0.53
    assert 0.15 <= rnd["violation"] <= 0.19
