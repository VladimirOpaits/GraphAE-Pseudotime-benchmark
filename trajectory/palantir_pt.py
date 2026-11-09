import argparse
import os

import anndata as ad
import numpy as np
import pandas as pd
import palantir


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--adata", required=True)
    p.add_argument("--roots", nargs="+", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--latent-key", default="X_latent")
    p.add_argument("--label-col", default="label")
    p.add_argument("--n-components", type=int, default=10)
    p.add_argument("--n-neighbors", type=int, default=30)
    p.add_argument("--num-waypoints", type=int, default=1200)
    p.add_argument("--seed", type=int, default=20)
    p.add_argument("--n-jobs", type=int, default=-1)
    p.add_argument("--prefix", default="palantir")
    a = p.parse_args()

    data = ad.read_h5ad(a.adata)
    if a.latent_key not in data.obsm:
        raise SystemExit(f"'{a.latent_key}' not in obsm")
    missing = [r for r in a.roots if r not in data.obs_names]
    if missing:
        raise SystemExit(f"Roots not in obs_names: {', '.join(missing)}")

    data.obsm[a.latent_key] = np.asarray(data.obsm[a.latent_key], dtype=np.float64)
    palantir.utils.run_diffusion_maps(data, n_components=a.n_components, knn=a.n_neighbors,
                                      pca_key=a.latent_key, seed=a.seed)
    palantir.utils.determine_multiscale_space(data)

    out = pd.DataFrame(index=pd.Index(data.obs_names, name="cell_id"))
    if a.label_col in data.obs:
        out["label"] = data.obs[a.label_col].to_numpy()
    os.makedirs(a.outdir, exist_ok=True)
    for r in a.roots:
        palantir.core.run_palantir(data, r, knn=a.n_neighbors, num_waypoints=a.num_waypoints,
                                   use_early_cell_as_start=True, n_jobs=a.n_jobs, seed=a.seed,
                                   save_as_df=True)
        out[r] = data.obs["palantir_pseudotime"].to_numpy()
        fate = pd.DataFrame(data.obsm["palantir_fate_probabilities"], index=out.index)
        fate.to_csv(os.path.join(a.outdir, f"{a.prefix}_fate_{r}.csv"))

    out.to_csv(os.path.join(a.outdir, f"{a.prefix}_pseudotime.csv"))


if __name__ == "__main__":
    main()
