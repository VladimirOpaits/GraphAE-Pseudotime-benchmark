import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr

from bench.frontends import make_frontend, save_frontend
from bench.graph import add_graph
from bench.trajectory import export_csv, run_dpt


@pytest.fixture
def curve(tmp_path):
    rng = np.random.default_rng(0)
    n = 300
    t = np.sort(rng.uniform(0, 1, n))
    x = np.column_stack([t, np.sin(6 * t), rng.normal(0, 0.01, (n, 3))]) + rng.normal(0, 0.02, (n, 5))
    obs = pd.DataFrame({"label": np.where(t < 0.2, "Start", "Other"), "time": np.where(t < 0.2, 2.0, 4.0)},
                       index=[f"c{i:04d}" for i in range(n)])
    data = make_frontend(x, obs)
    add_graph(data, 15)
    path = str(tmp_path / "curve.h5ad")
    save_frontend(data, path)
    return path, data, t


def test_pseudotime_aligned_ordered_and_rooted(curve, tmp_path):
    path, data, t = curve
    pt, root = run_dpt(path, str(tmp_path / "a"), "Start", 2.0, neighbors_key="k15", n_dcs=2)
    assert list(pt.index) == list(data.obs_names)
    assert data.obs.loc[root, "label"] == "Start"
    late = t > 0.2
    assert spearmanr(pt.to_numpy()[late], t[late])[0] > 0.95


def test_stored_graph_equals_graph_built_in_script(curve, tmp_path):
    path, _, _ = curve
    pt1, r1 = run_dpt(path, str(tmp_path / "a"), "Start", 2.0, seed=3, neighbors_key="k15")
    pt2, r2 = run_dpt(path, str(tmp_path / "b"), "Start", 2.0, seed=3, n_neighbors=15)
    assert r1 == r2
    assert np.allclose(pt1.to_numpy(), pt2.to_numpy())


def test_export_csv_keeps_cell_order(curve, tmp_path):
    _, data, _ = curve
    paths = export_csv(data, str(tmp_path / "csv"))
    coords = pd.read_csv(paths["coords"], index_col=0)
    assert list(coords.index) == list(data.obs_names)
    assert np.allclose(coords.to_numpy(), data.obsm["X_latent"])
