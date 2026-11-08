import time
import numpy as np, pandas as pd
from bench.data import load_larry
from bench.frontends import frontend_path, load_frontend, make_frontend, smooth
from bench.graph import get_graph

meta, _ = load_larry()
base = load_frontend(frontend_path("pca"))
lib = meta["Library"].to_numpy()
pop = np.where(meta["Starting population"].str.endswith("+"), "LSK", "LK")
strat = pd.Series(list(zip(meta["Time point"], meta["Well"], pop, meta["Cell type annotation"]))).astype(str).to_numpy()
s_codes, s_idx = np.unique(strat, return_inverse=True)
l_codes, l_idx = np.unique(lib, return_inverse=True)
n_s = np.bincount(s_idx)
n_sl = pd.Series(1, index=pd.MultiIndex.from_arrays([s_idx, l_idx])).groupby(level=[0, 1]).size()
libs_per_s = n_sl.groupby(level=0).size().reindex(range(len(s_codes)), fill_value=0).to_numpy()
ok_s = (n_s >= 300) & (libs_per_s >= 2)
exp_cell = n_sl.reindex(pd.MultiIndex.from_arrays([s_idx, l_idx])).to_numpy() / n_s[s_idx]
day_type = pd.Series(list(zip(meta["Time point"], meta["Cell type annotation"]))).astype(str).to_numpy()

def mixing(dist):
    d = dist.tocsr()
    rows = np.repeat(np.arange(d.shape[0]), np.diff(d.indptr)); cols = d.indices
    in_s = s_idx[rows] == s_idx[cols]
    same = in_s & (l_idx[rows] == l_idx[cols])
    n_in = np.bincount(rows[in_s], minlength=d.shape[0]); n_same = np.bincount(rows[same], minlength=d.shape[0])
    has = n_in > 0
    obs_cell = np.where(has, n_same / np.maximum(n_in, 1), np.nan)
    df = pd.DataFrame({"s": s_idx, "obs": obs_cell, "exp": exp_cell, "pop": pop})[has & ok_s[s_idx]]
    g = df.groupby(["pop", "s"]).agg(obs=("obs", "mean"), exp=("exp", "mean"))
    r = (g.obs / g.exp).groupby(level=0).median()
    return r.get("LSK"), r.get("LK")

def eta2(X, groups, mask):
    out = []
    for p in ["LSK", "LK"]:
        m = mask & (pop == p)
        x = X[m]; g = groups[m]; s = s_idx[m]
        df = pd.DataFrame(x); df["s"] = s; df["g"] = g
        tot = df.groupby("s").transform(lambda c: (c - c.mean()) ** 2).drop(columns="g").sum()
        mu_sg = df.groupby(["s", "g"]).transform("mean").drop(columns=[], errors="ignore")
        mu_s = df.drop(columns="g").groupby("s").transform("mean")
        between = ((mu_sg[list(range(X.shape[1]))] - mu_s[list(range(X.shape[1]))]) ** 2).sum()
        out.append((between / tot).to_numpy())
    return out

def eta2_type(X):
    out = []
    for p in ["LSK", "LK"]:
        m = pop == p
        df = pd.DataFrame(X[m]); df["d"] = meta["Time point"].to_numpy()[m]; df["t"] = day_type[m]
        mu_d = df.drop(columns="t").groupby("d").transform("mean")
        mu_t = df.drop(columns="d").groupby("t").transform("mean")
        cols = list(range(X.shape[1]))
        tot = ((df[cols] - mu_d[cols]) ** 2).sum()
        between = ((mu_t[cols] - mu_d[cols]) ** 2).sum()
        out.append((between / tot).to_numpy())
    return out

X0 = np.asarray(base.obsm["X_latent"])
_, conn0 = get_graph(base, 15)
mask = ok_s[s_idx]
lat = {"pca": X0}
for t in [1, 3, 10]:
    lat[f"smooth_t{t}"] = smooth(X0, conn0, t=t)

summary, comp = [], {}
for name, X in lat.items():
    t0 = time.time()
    data = base if name == "pca" else make_frontend(X, base.obs)
    dist, _ = get_graph(data, 15)
    r_lsk, r_lk = mixing(dist)
    e_lsk, e_lk = eta2(X, l_idx, mask)
    ty_lsk, ty_lk = eta2_type(X)
    w = X.var(0) / X.var(0).sum()
    summary.append((name, r_lsk, r_lk, (w * e_lsk).sum(), (w * e_lk).sum(), (w * ty_lk).sum()))
    comp[name] = (e_lk, ty_lk, w)
    print(name, "done", round(time.time() - t0), "s")

print("\nkNN same-library ratio (median over strata) and variance-weighted eta2 within strata")
print(pd.DataFrame(summary, columns=["latent", "ratio_LSK", "ratio_LK", "eta2_lib_LSK", "eta2_lib_LK", "eta2_type_LK"]).round(3).to_string(index=False))

e_lk, ty_lk, w = comp["pca"]
order = np.argsort(-e_lk)
print("\nPCA, LK: components with most library variance")
print(pd.DataFrame({"pc": order[:10] + 1, "eta2_lib": e_lk[order[:10]], "eta2_type": ty_lk[order[:10]], "var_share": w[order[:10]]}).round(3).to_string(index=False))
print("\nPCA, LK: eta2_lib by component block")
for a, b in [(0, 5), (5, 10), (10, 20), (20, 30), (30, 50)]:
    print(f"PC{a+1}-{b}: lib {e_lk[a:b].mean():.3f}  type {ty_lk[a:b].mean():.3f}  var share {w[a:b].sum():.3f}")
