import argparse
import os

import anndata as ad
import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--adata", required=True)
    p.add_argument("--roots", nargs="+", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--latent-key", default="X_latent")
    p.add_argument("--label-col", default="label")
    p.add_argument("--prefix", default="euclid")
    a = p.parse_args()

    data = ad.read_h5ad(a.adata)
    if a.latent_key not in data.obsm:
        raise SystemExit(f"'{a.latent_key}' not in obsm")
    missing = [r for r in a.roots if r not in data.obs_names]
    if missing:
        raise SystemExit(f"Roots not in obs_names: {', '.join(missing)}")

    x = np.asarray(data.obsm[a.latent_key], dtype=np.float64)
    out = pd.DataFrame(index=pd.Index(data.obs_names, name="cell_id"))
    if a.label_col in data.obs:
        out["label"] = data.obs[a.label_col].to_numpy()
    for r in a.roots:
        out[r] = np.linalg.norm(x - x[data.obs_names.get_loc(r)], axis=1)

    os.makedirs(a.outdir, exist_ok=True)
    out.to_csv(os.path.join(a.outdir, f"{a.prefix}_pseudotime.csv"))


if __name__ == "__main__":
    main()
