import numpy as np
import pandas as pd
import scipy.sparse as sp

from bench.frontends import make_frontend, pca, smooth
from bench.graph import get_graph


def frontend(n=200, d=10):
    x = np.random.default_rng(0).normal(size=(n, d))
    obs = pd.DataFrame({"label": "A", "time": 2.0}, index=[f"c{i}" for i in range(n)])
    return make_frontend(x, obs)


def test_pca_is_deterministic():
    rng = np.random.default_rng(0)
    X = sp.csr_matrix(rng.poisson(rng.gamma(0.5, 2.0, 500), (300, 500)).astype(np.float32))
    a = pca(X.copy(), n_hvg=200, n_pcs=10)
    b = pca(X.copy(), n_hvg=200, n_pcs=10)
    assert np.array_equal(a, b)


def test_graph_is_deterministic():
    _, c1 = get_graph(frontend(), 15)
    _, c2 = get_graph(frontend(), 15)
    assert (c1 != c2).nnz == 0


def test_smooth_averages_neighbours():
    data = frontend()
    x = data.obsm["X_latent"]
    _, conn = get_graph(data, 10)
    assert np.allclose(smooth(np.ones_like(x), conn, t=3), 1)
    assert smooth(x, conn, t=3).var() < x.var()
