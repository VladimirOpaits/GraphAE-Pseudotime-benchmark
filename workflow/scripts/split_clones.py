import argparse

import numpy as np
import pandas as pd
import scipy.io as sio

p = argparse.ArgumentParser()
p.add_argument("--meta", required=True)
p.add_argument("--clones", required=True)
p.add_argument("--out", required=True)
p.add_argument("--dev-frac", type=float, default=0.3)
p.add_argument("--seed", type=int, default=0)
a = p.parse_args()

labels = pd.read_csv(a.meta, sep="\t")["Cell type annotation"].to_numpy()
cl = sio.mmread(a.clones).tocsc()

rows = []
for c in range(cl.shape[1]):
    lab = labels[cl.indices[cl.indptr[c]:cl.indptr[c + 1]]]
    mature = lab[lab != "Undifferentiated"]
    fate = pd.Series(mature).mode()[0] if len(mature) else "Undifferentiated"
    rows.append((c, len(lab), fate))
df = pd.DataFrame(rows, columns=["clone", "size", "fate"])
df["size_bin"] = pd.cut(df["size"], [0, 2, 4, 8, 16, np.inf], labels=["1-2", "3-4", "5-8", "9-16", "17+"])

rng = np.random.default_rng(a.seed)
df["split"] = "test"
for _, g in df.groupby(["fate", "size_bin"], observed=True):
    dev = rng.choice(g.index.to_numpy(), size=int(round(a.dev_frac * len(g))), replace=False)
    df.loc[dev, "split"] = "dev"
df.to_csv(a.out, index=False)
