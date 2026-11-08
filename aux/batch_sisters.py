import numpy as np, pandas as pd
from bench.data import load_larry
from bench.frontends import frontend_path, load_frontend

meta, clones = load_larry()
data = load_frontend(frontend_path("pca"))
X = np.asarray(data.obsm["X_latent"])
lib = meta["Library"].to_numpy(); day = meta["Time point"].to_numpy()
well = meta["Well"].to_numpy(); pop = meta["Starting population"].to_numpy()
typ = meta["Cell type annotation"].to_numpy()
scheme = np.where(pd.Series(lib).str.startswith(("LK_", "LSK_")), "named", "bare")
rng = np.random.default_rng(0)
key = pd.Series(list(zip(day, well, pop, typ))).astype(str).to_numpy()
key_lib = pd.Series(list(zip(day, well, pop, typ, lib))).astype(str).to_numpy()
by_key = pd.Series(np.arange(len(meta))).groupby(key).apply(np.array).to_dict()
by_key_lib = pd.Series(np.arange(len(meta))).groupby(key_lib).apply(np.array).to_dict()
cell_clone = np.full(len(meta), -1); cl = clones.tocsc()
for c in range(cl.shape[1]):
    cell_clone[cl.indices[cl.indptr[c]:cl.indptr[c + 1]]] = c
rows = []
for c in range(cl.shape[1]):
    cells = cl.indices[cl.indptr[c]:cl.indptr[c + 1]]
    for k in np.unique(key[cells]):
        g = cells[key[cells] == k]
        if len(g) < 2: continue
        a, b = np.triu_indices(len(g), 1)
        if len(a) > 50:
            s = rng.choice(len(a), 50, replace=False); a, b = a[s], b[s]
        for i, j in zip(g[a], g[b]):
            same = lib[i] == lib[j]
            pool = by_key_lib[key_lib[j]] if same else by_key[key[j]]
            if not same: pool = pool[lib[pool] != lib[i]]
            pool = pool[cell_clone[pool] != c]
            if len(pool) == 0: continue
            rows.append((("LSK" if pop[i].endswith("+") else "LK"), day[i], same, np.linalg.norm(X[i] - X[j]), np.linalg.norm(X[i] - X[rng.choice(pool)])))
df = pd.DataFrame(rows, columns=["pop", "day", "same_lib", "sister", "random"])
print("TEST 1 by population")
print(df.groupby(["pop", "day", "same_lib"]).agg(n=("sister", "size"), sister=("sister", "median"), random=("random", "median")).round(2).to_string())

D = data.obsp["k15_distances"].tocsr()
print("TEST 2 within day x well x pop x type (types with >=300 cells in stratum)")
for p in np.unique(pop):
    strat = pd.Series(list(zip(day, well, pop, typ))).astype(str).to_numpy()
    out = []
    for s in np.unique(strat[pop == p]):
        m = np.flatnonzero(strat == s)
        if len(m) < 300 or len(np.unique(lib[m])) < 2: continue
        libs, cnt = np.unique(lib[m], return_counts=True); frac = dict(zip(libs, cnt / cnt.sum()))
        sch, sc = np.unique(scheme[m], return_counts=True); sfrac = dict(zip(sch, sc / sc.sum()))
        ol, el, os_, es = [], [], [], []
        for i in m:
            nb = D.indices[D.indptr[i]:D.indptr[i + 1]]; nb = nb[strat[nb] == s]
            if len(nb) == 0: continue
            ol.append((lib[nb] == lib[i]).mean()); el.append(frac[lib[i]])
            if len(sch) > 1: os_.append((scheme[nb] == scheme[i]).mean()); es.append(sfrac[scheme[i]])
        out.append((s, len(m), np.mean(ol) / np.mean(el), np.mean(os_) / np.mean(es) if os_ else np.nan))
    o = pd.DataFrame(out, columns=["stratum", "n", "lib_ratio", "scheme_ratio"])
    print(p, f"strata {len(o)}  lib ratio median {o.lib_ratio.median():.2f} [{o.lib_ratio.min():.2f}-{o.lib_ratio.max():.2f}]  scheme ratio median {o.scheme_ratio.median():.2f}")
    print(o.sort_values("lib_ratio", ascending=False).head(6).round(2).to_string(index=False))
