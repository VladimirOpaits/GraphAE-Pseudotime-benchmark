import argparse
import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


def read_table(path):
    return pd.read_csv(path, index_col=0, dtype={0: str})


def pick_col(df, col):
    col = col or df.columns[0]
    if col not in df.columns:
        raise SystemExit(f"Column '{col}' not found")
    return df[col]


def load_graph(path, cell_ids):
    z = np.load(path)
    if not np.array_equal(z["cell_ids"], np.asarray(cell_ids, dtype=str)):
        raise SystemExit("Graph cell ids do not match coords")
    n = len(cell_ids)
    dist = sp.csr_matrix((z["dist_data"], z["dist_indices"], z["dist_indptr"]), shape=(n, n))
    conn = sp.csr_matrix((z["conn_data"], z["conn_indices"], z["conn_indptr"]), shape=(n, n))
    return dist, conn, int(z["k"]), int(z["seed"])


def set_graph(data, dist, conn, k, seed):
    data.obsp["distances"] = dist
    data.obsp["connectivities"] = conn
    data.uns["neighbors"] = {
        "connectivities_key": "connectivities",
        "distances_key": "distances",
        "params": {"n_neighbors": k, "method": "umap", "random_state": seed, "metric": "euclidean", "use_rep": "X_latent"},
    }


def pick_root(latent, pool, seed, top_frac, min_top):
    center = latent[pool].mean(0)
    dist = ((latent[pool] - center) ** 2).sum(1)
    k = max(min_top, int(top_frac * len(pool)))
    near = pool[np.argsort(dist)[:k]]
    return int(np.random.default_rng(seed).choice(near))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--coords", required=True)
    p.add_argument("--labels", required=True) # тип клетки
    p.add_argument("--time", required=True)   # день каждой клетки
    p.add_argument("--start", required=True)  # метка типа для выбора корня
    p.add_argument("--start-time", type=float, required=True)  # почти всегда берем 2.0
    p.add_argument("--outdir", required=True)
    p.add_argument("--graph")
    p.add_argument("--label-col")             # названия в файлах
    p.add_argument("--time-col")              # названия в файлах
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--root-seed", type=int)
    p.add_argument("--n-dcs", type=int, default=10)  # количество компонент, слишком много - шум
    p.add_argument("--n-neighbors", type=int, default=15)  # размер kNN графа
    p.add_argument("--top-frac", type=float, default=0.05)
    p.add_argument("--min-top", type=int, default=20)
    p.add_argument("--prefix", default="dpt")
    a = p.parse_args()

    coords = read_table(a.coords)
    labels = pick_col(read_table(a.labels), a.label_col).reindex(coords.index) # 
    time = pick_col(read_table(a.time), a.time_col).reindex(coords.index)
    latent = coords.to_numpy(dtype=np.float32)

    pool = np.flatnonzero((labels.to_numpy() == a.start) & (time.to_numpy() == a.start_time))
    if not len(pool):
        raise SystemExit(f"No cells with label '{a.start}' at time {a.start_time}")
    root_seed = a.seed if a.root_seed is None else a.root_seed
    root = pick_root(latent, pool, root_seed, a.top_frac, a.min_top)

    data = ad.AnnData(np.zeros((len(coords), 1), dtype=np.float32), obs=pd.DataFrame(index=coords.index))
    data.obsm["X_latent"] = latent
    if a.graph:
        set_graph(data, *load_graph(a.graph, coords.index))
    else:
        sc.pp.neighbors(data, n_neighbors=a.n_neighbors, use_rep="X_latent", random_state=a.seed)
    sc.tl.diffmap(data, n_comps=a.n_dcs + 5)
    data.uns["iroot"] = root
    sc.tl.dpt(data, n_dcs=a.n_dcs)

    pt = data.obs["dpt_pseudotime"].to_numpy().copy()
    pt[~np.isfinite(pt)] = np.nan
    os.makedirs(a.outdir, exist_ok=True)
    out = pd.DataFrame({"cell_id": coords.index, "label": labels.to_numpy(), "dpt": pt})
    out.to_csv(os.path.join(a.outdir, f"{a.prefix}_pseudotime.csv"), index=False)
    with open(os.path.join(a.outdir, f"{a.prefix}_root.txt"), "w") as f:
        f.write(f"{coords.index[root]}\n")


if __name__ == "__main__":
    main()
