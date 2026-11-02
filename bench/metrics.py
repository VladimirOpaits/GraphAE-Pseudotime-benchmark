import numpy as np
import pandas as pd

from bench.data import DATA, clone_pairs, load_larry


def _pct_rank_by_branch(pt, branch):
    out = np.full(len(pt), np.nan)
    br = np.asarray(branch, dtype=object)
    ok = pd.notna(br) & ~np.isnan(pt)
    for b in pd.unique(br[ok]):
        m = ok & (br == b)
        out[m] = pd.Series(pt[m]).rank(pct=True).to_numpy()
    return out


def sister_order(pt, branch, pairs, n_boot=1000, min_pairs=200, seed=0):
    i, j, c = pairs
    bi, bj = np.asarray(branch, dtype=object)[i], np.asarray(branch, dtype=object)[j]
    ok = (bi == bj) & pd.notna(bi) & ~np.isnan(pt[i]) & ~np.isnan(pt[j])
    i, j, c, br = i[ok], j[ok], c[ok], bi[ok]
    r = _pct_rank_by_branch(pt, branch)
    d = r[j] - r[i]
    acc = np.where(pt[j] > pt[i], 1.0, np.where(pt[j] == pt[i], 0.5, 0.0))
    df = pd.DataFrame({"clone": c, "branch": br, "acc": acc,
                       "violation": np.maximum(0.0, -d), "gap": d})
    stats = ["acc", "violation", "gap"]
    g = df.groupby(["branch", "clone"]).agg(
        n=("acc", "size"), **{k: (k, "sum") for k in stats}).reset_index()
    rng = np.random.default_rng(seed)
    rows, boots = [], {}
    for b, gb in g.groupby("branch"):
        n = gb["n"].to_numpy()
        idx = rng.integers(0, len(gb), (n_boot, len(gb)))
        row = dict(branch=b, n_pairs=int(n.sum()), n_clones=len(gb),
                   used=n.sum() >= min_pairs)
        boots[b] = {}
        for k in stats:
            sk = gb[k].to_numpy()
            bs = sk[idx].sum(1) / n[idx].sum(1)
            boots[b][k] = bs
            row[k] = sk.sum() / n.sum()
            row[k + "_lo"], row[k + "_hi"] = np.percentile(bs, [2.5, 97.5])
        rows.append(row)
    res = pd.DataFrame(rows)
    used = res[res["used"]]
    if len(used):
        row = dict(branch="mean", n_pairs=int(used["n_pairs"].sum()),
                   n_clones=int(used["n_clones"].sum()), used=True)
        for k in stats:
            mb = np.mean([boots[b][k] for b in used["branch"]], axis=0)
            row[k] = used[k].mean()
            row[k + "_lo"], row[k + "_hi"] = np.percentile(mb, [2.5, 97.5])
        res.loc[len(res)] = row
    cols = ["branch", "n_pairs", "n_clones", "used"] + [
        x for k in stats for x in (k, k + "_lo", k + "_hi")]
    return res[cols]


def report(res):
    cols = ["branch", "n_pairs", "acc", "violation", "gap"]
    out = res[cols].copy()
    for k in ["acc", "violation", "gap"]:
        out[k] = [f"{v:.3f} [{lo:.3f}, {hi:.3f}]" for v, lo, hi in
                  zip(res[k], res[k + "_lo"], res[k + "_hi"])]
    return out.to_string(index=False)


if __name__ == "__main__":
    meta, clones = load_larry()
    pairs = clone_pairs(meta["Time point"].to_numpy(), clones)
    ps = pd.read_csv(DATA + "stateFate_inVitro_neutrophil_pseudotime.txt.gz", sep="\t")
    pt = np.full(len(meta), np.nan)
    pt[ps["Cell index"].to_numpy()] = ps["pseudotime"].to_numpy()
    branch = np.where(np.isnan(pt), None, "Neutrophil")
    rnd = np.random.default_rng(0).random(len(meta))
    print("authors' pseudotime\n" + report(sister_order(pt, branch, pairs)))
    print("random pseudotime\n" + report(sister_order(rnd, branch, pairs)))
