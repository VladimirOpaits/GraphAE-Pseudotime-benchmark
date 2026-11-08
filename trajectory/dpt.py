import argparse
import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--adata", required=True)
    p.add_argument("--roots", nargs="+", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--latent-key", default="X_latent")
    p.add_argument("--neighbors-key")
    p.add_argument("--label-col", default="label")  # названия в файлах
    p.add_argument("--n-dcs", type=int, default=10)  # количество компонент, слишком много - шум
    p.add_argument("--n-neighbors", type=int, default=15)  # размер kNN графа
    p.add_argument("--prefix", default="dpt")
    a = p.parse_args()

    data = ad.read_h5ad(a.adata)
    if a.latent_key not in data.obsm:
        raise SystemExit(f"'{a.latent_key}' not in obsm")
    missing = [r for r in a.roots if r not in data.obs_names]
    if missing:
        raise SystemExit(f"Roots not in obs_names: {', '.join(missing)}")

    key = a.neighbors_key
    if key is None:
        key = "neighbors"
        sc.pp.neighbors(data, n_neighbors=a.n_neighbors, use_rep=a.latent_key, transformer="sklearn")
    elif key not in data.uns:
        raise SystemExit(f"Graph '{key}' not in uns")
    sc.tl.diffmap(data, n_comps=a.n_dcs + 5, neighbors_key=key)

    out = pd.DataFrame(index=pd.Index(data.obs_names, name="cell_id"))
    if a.label_col in data.obs:
        out["label"] = data.obs[a.label_col].to_numpy()
    for r in a.roots:
        data.uns["iroot"] = int(data.obs_names.get_loc(r))
        sc.tl.dpt(data, n_dcs=a.n_dcs, neighbors_key=key)
        pt = data.obs["dpt_pseudotime"].to_numpy().copy()
        pt[~np.isfinite(pt)] = np.nan
        out[r] = pt

    os.makedirs(a.outdir, exist_ok=True)
    out.to_csv(os.path.join(a.outdir, f"{a.prefix}_pseudotime.csv"))


if __name__ == "__main__":
    main()
