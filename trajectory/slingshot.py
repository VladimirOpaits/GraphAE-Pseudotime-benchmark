import argparse
import os
import subprocess
import tempfile

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
    p.add_argument("--label-col", default="label")
    p.add_argument("--n-neighbors", type=int, default=15)
    p.add_argument("--resolution", type=float, default=1.0)
    p.add_argument("--cluster-by", choices=["leiden", "leiden_start", "label"], default="leiden")
    p.add_argument("--time-col", default="time")
    p.add_argument("--start", default="Undifferentiated")
    p.add_argument("--start-time", type=float, default=2.0)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--prefix", default="slingshot")
    p.add_argument("--plot", action="store_true")
    p.add_argument("--meta", default="rt_data/stateFate_inVitro_metadata.txt.gz")
    a = p.parse_args()

    data = ad.read_h5ad(a.adata)
    if a.latent_key not in data.obsm:
        raise SystemExit(f"'{a.latent_key}' not in obsm")
    missing = [r for r in a.roots if r not in data.obs_names]
    if missing:
        raise SystemExit(f"Roots not in obs_names: {', '.join(missing)}")

    if a.cluster_by == "label":
        cl = data.obs[a.label_col].astype(str).to_numpy()
    else:
        sc.pp.neighbors(data, n_neighbors=a.n_neighbors, use_rep=a.latent_key, transformer="sklearn")
        sc.tl.leiden(data, resolution=a.resolution, random_state=a.seed, flavor="igraph",
                     n_iterations=2, directed=False)
        cl = data.obs["leiden"].astype(str).to_numpy()
        if a.cluster_by == "leiden_start":
            m = (data.obs[a.label_col].to_numpy() == a.start) & (data.obs[a.time_col].to_numpy() == a.start_time)
            cl = np.where(m, "start", cl)
    starts = {r: cl[data.obs_names.get_loc(r)] for r in a.roots}
    uniq = sorted(set(starts.values()))

    x = np.ascontiguousarray(data.obsm[a.latent_key], dtype=np.float64)
    script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "slingshot.R")
    os.makedirs(a.outdir, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=a.outdir) as tmp:
        x.tofile(os.path.join(tmp, "x.bin"))
        with open(os.path.join(tmp, "cl.txt"), "w") as f:
            f.write("\n".join(cl) + "\n")
        cmd = ["Rscript", script, "--x", os.path.join(tmp, "x.bin"), "--n", str(x.shape[0]),
               "--d", str(x.shape[1]), "--clusters", os.path.join(tmp, "cl.txt"),
               "--starts", ",".join(uniq), "--outdir", tmp, "--seed", str(a.seed), "--prefix", a.prefix]
        if a.plot:
            meta = pd.read_csv(a.meta, sep="\t", usecols=["SPRING-x", "SPRING-y"])
            xy = np.ascontiguousarray(meta.to_numpy()[data.obs_names.astype(int)], dtype=np.float64)
            xy.tofile(os.path.join(tmp, "xy.bin"))
            with open(os.path.join(tmp, "lab.txt"), "w") as f:
                f.write("\n".join(data.obs[a.label_col].astype(str)) + "\n")
            cmd += ["--plot-xy", os.path.join(tmp, "xy.bin"), "--plot-labels", os.path.join(tmp, "lab.txt"),
                    "--plot-dir", a.outdir]
        subprocess.run(cmd, check=True)
        pts = {s: np.fromfile(os.path.join(tmp, f"pt_{s}.bin")) for s in uniq}
        lins = {s: open(os.path.join(tmp, f"lineages_{s}.txt")).read().strip().split("\n") for s in uniq}

    out = pd.DataFrame(index=pd.Index(data.obs_names, name="cell_id"))
    if a.label_col in data.obs:
        out["label"] = data.obs[a.label_col].to_numpy()
    for r in a.roots:
        out[r] = pts[starts[r]]
    out.to_csv(os.path.join(a.outdir, f"{a.prefix}_pseudotime.csv"))
    pd.DataFrame({"root": a.roots, "start_cluster": [starts[r] for r in a.roots],
                  "n_lineages": [len(lins[starts[r]]) for r in a.roots],
                  "lineages": [" | ".join(lins[starts[r]]) for r in a.roots],
                  "n_clusters": len(set(cl))}).to_csv(os.path.join(a.outdir, f"{a.prefix}_info.csv"), index=False)


if __name__ == "__main__":
    main()
