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
    root = data.obs_names[0]
    pt = run_dpt(path, str(tmp_path / "a"), [root], neighbors_key="k15", n_dcs=2)
    assert list(pt.index) == list(data.obs_names)
    assert list(pt.columns) == [root]
    assert pt.loc[root, root] == 0
    late = t > 0.2
    assert spearmanr(pt[root].to_numpy()[late], t[late])[0] > 0.95


def test_stored_graph_equals_graph_built_in_script(curve, tmp_path):
    path, data, _ = curve
    roots = list(data.obs_names[[0, 5]])
    pt1 = run_dpt(path, str(tmp_path / "a"), roots, neighbors_key="k15")
    pt2 = run_dpt(path, str(tmp_path / "b"), roots, n_neighbors=15)
    assert np.allclose(pt1.to_numpy(), pt2.to_numpy())


def test_many_roots_equal_single_runs(curve, tmp_path):
    path, data, _ = curve
    roots = list(data.obs_names[[0, 7, 30]])
    many = run_dpt(path, str(tmp_path / "many"), roots, neighbors_key="k15")
    for r in roots:
        one = run_dpt(path, str(tmp_path / r), [r], neighbors_key="k15")
        assert np.allclose(many[r].to_numpy(), one[r].to_numpy())


def test_export_csv_keeps_cell_order(curve, tmp_path):
    _, data, _ = curve
    paths = export_csv(data, str(tmp_path / "csv"))
    coords = pd.read_csv(paths["coords"], index_col=0)
    assert list(coords.index) == list(data.obs_names)
    assert np.allclose(coords.to_numpy(), data.obsm["X_latent"])
