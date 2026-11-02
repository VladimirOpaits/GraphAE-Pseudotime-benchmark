import argparse
import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


def pick_root(latent, pool, seed, top_frac, min_top):
    center = latent[pool].mean(0)                       # стартовое среднее
    dist = ((latent[pool] - center) ** 2).sum(1)        # все дистанции в типе
    k = max(min_top, int(top_frac * len(pool)))         # набираем из них самые центральные
    near = pool[np.argsort(dist)[:k]]                   # сортировка
    return int(np.random.default_rng(seed).choice(near))  # берем конкретным seed


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--adata", required=True)
    p.add_argument("--start", required=True)  # метка типа для выбора корня
    p.add_argument("--start-time", type=float, required=True)  # почти всегда берем 2.0
    p.add_argument("--outdir", required=True)
    p.add_argument("--latent-key", default="X_latent")
    p.add_argument("--neighbors-key")
    p.add_argument("--label-col", default="label")  # названия в файлах
    p.add_argument("--time-col", default="time")    # названия в файлах
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--n-dcs", type=int, default=10)  # количество компонент, слишком много - шум
    p.add_argument("--n-neighbors", type=int, default=15)  # размер kNN графа
    p.add_argument("--top-frac", type=float, default=0.05)
    p.add_argument("--min-top", type=int, default=20)
    p.add_argument("--prefix", default="dpt")
    a = p.parse_args()

    data = ad.read_h5ad(a.adata)
    for col in [a.label_col, a.time_col]:
        if col not in data.obs:
            raise SystemExit(f"Column '{col}' not in obs")
    if a.latent_key not in data.obsm:
        raise SystemExit(f"'{a.latent_key}' not in obsm")
    latent = np.asarray(data.obsm[a.latent_key], dtype=np.float32)
    labels = data.obs[a.label_col].to_numpy()
    time = data.obs[a.time_col].to_numpy()

    pool = np.flatnonzero((labels == a.start) & (time == a.start_time))
    if not len(pool):
        raise SystemExit(f"No cells with label '{a.start}' at time {a.start_time}")
    root = pick_root(latent, pool, a.seed, a.top_frac, a.min_top)

    key = a.neighbors_key
    if key is None:
        key = "neighbors"
        sc.pp.neighbors(data, n_neighbors=a.n_neighbors, use_rep=a.latent_key, transformer="sklearn")
    elif key not in data.uns:
        raise SystemExit(f"Graph '{key}' not in uns")
    sc.tl.diffmap(data, n_comps=a.n_dcs + 5, neighbors_key=key)
    data.uns["iroot"] = root
    sc.tl.dpt(data, n_dcs=a.n_dcs, neighbors_key=key)

    pt = data.obs["dpt_pseudotime"].to_numpy().copy()
    pt[~np.isfinite(pt)] = np.nan
    os.makedirs(a.outdir, exist_ok=True)
    out = pd.DataFrame({"cell_id": data.obs_names, "label": labels, "dpt": pt})
    out.to_csv(os.path.join(a.outdir, f"{a.prefix}_pseudotime.csv"), index=False)
    with open(os.path.join(a.outdir, f"{a.prefix}_root.txt"), "w") as f:
        f.write(f"{data.obs_names[root]}\n")


if __name__ == "__main__":
    main()
