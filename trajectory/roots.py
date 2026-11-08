import argparse
import anndata as ad
import numpy as np

def pick_roots(latent, pool, n_roots, pool_size, seed):
    center = latent[pool].mean(0)
    dist = ((latent[pool] - center) ** 2).sum(1)
    near = pool[np.argsort(dist)[:pool_size]]
    return np.random.default_rng(seed).choice(near, size=n_roots, replace=False)



p = argparse.ArgumentParser()
p.add_argument("--adata", required=True)
p.add_argument("--out", required=True)
p.add_argument("--latent-key", default="X_latent")
p.add_argument("--label-col", default="label")
p.add_argument("--time-col", default="time")
p.add_argument("--start", default="Undifferentiated")
p.add_argument("--start-time", type=float, default=2.0)
p.add_argument("--n-roots", type=int, default=5)
p.add_argument("--seed", type=int, default=0)
p.add_argument("--pool-size", type = int, default = 100)
a = p.parse_args()

data = ad.read_h5ad(a.adata)
for col in [a.label_col, a.time_col]:
    if col not in data.obs: raise SystemExit(f"Column '{col}' not in obs")
if a.latent_key not in data.obsm: raise SystemExit(f"'{a.latent_key}' not in obsm")

latent = np.asarray(data.obsm[a.latent_key], dtype=np.float32)
labels = data.obs[a.label_col].to_numpy()
time = data.obs[a.time_col].to_numpy()
pool = np.flatnonzero((labels == a.start) & (time == a.start_time))
if len(pool) < a.n_roots: raise SystemExit(f"Only {len(pool)} cells with label '{a.start}' at time {a.start_time}")

roots = pick_roots(latent, pool, a.n_roots, a.seed, a.top_frac, a.min_top)
with open(a.out, "w") as f: f.write("\n".join(data.obs_names[roots]) + "\n")