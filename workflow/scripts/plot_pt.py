import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--pt", required=True)
p.add_argument("--meta", required=True)
p.add_argument("--out", required=True)
p.add_argument("--max-cells", type=int, default=30000)
p.add_argument("--seed", type=int, default=0)
a = p.parse_args()

pt = pd.read_csv(a.pt, index_col="cell_id")
meta = pd.read_csv(a.meta, sep="\t", usecols=["SPRING-x", "SPRING-y", "Cell type annotation"])
m = meta.iloc[pt.index.to_numpy()]
roots = [c for c in pt.columns if c != "label"]
rng = np.random.default_rng(a.seed)
show = np.sort(rng.choice(len(pt), min(a.max_cells, len(pt)), replace=False))
x, y = m["SPRING-x"].to_numpy()[show], m["SPRING-y"].to_numpy()[show]

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
labels = m["Cell type annotation"].to_numpy()[show]
for lab in sorted(set(labels)):
    k = labels == lab
    axes[0].scatter(x[k], y[k], s=1, alpha=0.5, label=lab, rasterized=True)
axes[0].legend(markerscale=6, fontsize=7, frameon=False)
axes[0].set_title("labels")

r = pd.Series(pt[roots[0]].to_numpy()).rank(pct=True).to_numpy()[show]
o = np.argsort(np.nan_to_num(r, nan=-1))
sc_ = axes[1].scatter(x[o], y[o], c=r[o], s=1, cmap="plasma", rasterized=True)
fig.colorbar(sc_, ax=axes[1], label="rank")
for rt in roots:
    axes[1].scatter(meta["SPRING-x"].iloc[int(rt)], meta["SPRING-y"].iloc[int(rt)], s=60, c="cyan",
                    edgecolors="black", marker="*")
axes[1].set_title(f"pseudotime rank, root {roots[0]} (stars: all roots)")
for ax in axes:
    ax.set_xticks([])
    ax.set_yticks([])
fig.suptitle(a.pt)
fig.tight_layout()
fig.savefig(a.out, dpi=150)
