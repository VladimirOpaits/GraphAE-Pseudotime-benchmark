import argparse

import numpy as np
import pandas as pd

from bench.data import clone_pairs, load_branches, load_larry
from bench.metrics import sister_order

p = argparse.ArgumentParser()
p.add_argument("--pt", required=True)
p.add_argument("--split", required=True)
p.add_argument("--part", choices=["dev", "test"], required=True)
p.add_argument("--out", required=True)
p.add_argument("--n-boot", type=int, default=1000)
a = p.parse_args()

meta, clones = load_larry()
n = len(meta)
branch = load_branches(n)
i, j, c = clone_pairs(meta["Time point"].to_numpy(), clones)
split = pd.read_csv(a.split)
keep = np.isin(c, split.loc[split["split"] == a.part, "clone"].to_numpy())
pairs = (i[keep], j[keep], c[keep])

pt = pd.read_csv(a.pt, index_col="cell_id")
ids = pt.index.to_numpy().astype(int)
rows = []
for r in [col for col in pt.columns if col != "label"]:
    x = np.full(n, np.nan)
    x[ids] = pt[r].to_numpy()
    res = sister_order(x, branch, pairs, n_boot=a.n_boot)
    res.insert(0, "root", r)
    res.insert(1, "n_nan", int(np.isnan(x).sum()))
    rows.append(res)
pd.concat(rows).to_csv(a.out, index=False)
