import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from bench.data import DATA

LATENT = "X_latent"
ROOT = "runs/frontends"


def load_counts(cells=None):
    X = sp.load_npz(DATA + "cache/normed_counts_f32.npz").tocsr() # формат компрессированных sparse матриц
    return X if cells is None else X[cells]


def larry_obs(meta, cells=None):
    idx = np.arange(len(meta)) if cells is None else np.asarray(cells)
    m = meta.iloc[idx]
    return pd.DataFrame(
        {"label": m["Cell type annotation"].to_numpy(), "time": m["Time point"].to_numpy().astype(float)},
        index=pd.Index(idx.astype(str), name="cell_id"))


def pca(X, n_hvg=2000, n_pcs=50):
    data = sc.AnnData(X)
    sc.pp.log1p(data)
    sc.pp.highly_variable_genes(data, n_top_genes=n_hvg)
    data = data[:, data.var["highly_variable"]].copy()
    sc.pp.scale(data, max_value=10)
    sc.pp.pca(data, n_comps=n_pcs, svd_solver="covariance_eigh")
    return data.obsm["X_pca"].astype(np.float32)


def diffusion_operator(conn, self_loops=True):
    w = conn + sp.identity(conn.shape[0], format="csr") if self_loops else conn
    deg = np.asarray(w.sum(1)).ravel()
    return (sp.diags(1 / deg) @ w).tocsr()


def smooth(latent, conn, t=3):
    p = diffusion_operator(conn)
    x = np.asarray(latent, dtype=np.float32)
    for _ in range(t):
        x = p @ x
    return x.astype(np.float32)


def make_frontend(latent, obs, **params):
    data = ad.AnnData(obs=obs.copy())
    data.obsm[LATENT] = np.asarray(latent, dtype=np.float32)
    data.uns["frontend"] = params
    return data


def frontend_path(name, root=ROOT):
    return os.path.join(root, f"{name}.h5ad")


def save_frontend(data, path):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    data.write_h5ad(path)


def load_frontend(path):
    return ad.read_h5ad(path)
