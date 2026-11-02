import hashlib
import os

import anndata as ad
import numpy as np
import scanpy as sc
import scipy.sparse as sp


def latent_hash(latent):
    return hashlib.sha1(np.ascontiguousarray(latent, dtype=np.float32).tobytes()).hexdigest()


def build_graph(latent, k, seed):
    data = ad.AnnData(np.zeros((len(latent), 1), dtype=np.float32))
    data.obsm["X_latent"] = np.asarray(latent, dtype=np.float32)
    sc.pp.neighbors(data, n_neighbors=k, use_rep="X_latent", random_state=seed)
    return data.obsp["distances"].tocsr(), data.obsp["connectivities"].tocsr()


def save_graph(path, dist, conn, cell_ids, k, seed, h):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    np.savez_compressed(
        path,
        dist_data=dist.data, dist_indices=dist.indices, dist_indptr=dist.indptr,
        conn_data=conn.data, conn_indices=conn.indices, conn_indptr=conn.indptr,
        cell_ids=np.asarray(cell_ids, dtype=str), k=k, seed=seed, latent_hash=h,
    )


def load_graph(path):
    z = np.load(path)
    n = len(z["cell_ids"])
    dist = sp.csr_matrix((z["dist_data"], z["dist_indices"], z["dist_indptr"]), shape=(n, n))
    conn = sp.csr_matrix((z["conn_data"], z["conn_indices"], z["conn_indptr"]), shape=(n, n))
    meta = {"cell_ids": z["cell_ids"], "k": int(z["k"]), "seed": int(z["seed"]), "latent_hash": str(z["latent_hash"])}
    return dist, conn, meta


def knn_graph(latent, cell_ids, k, seed, path):
    h = latent_hash(latent)
    if os.path.exists(path):
        dist, conn, meta = load_graph(path)
        same = (meta["latent_hash"] == h and meta["k"] == k and meta["seed"] == seed
                and np.array_equal(meta["cell_ids"], np.asarray(cell_ids, dtype=str)))
        if same:
            return dist, conn
        print(f"stale graph, rebuilding: {path}")
    dist, conn = build_graph(latent, k, seed)
    save_graph(path, dist, conn, cell_ids, k, seed, h)
    return dist, conn


def graph_path(root, name, k, seed):
    return os.path.join(root, name, f"knn_k{k}_s{seed}.npz")
