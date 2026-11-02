import os
import subprocess
import sys

import numpy as np
import pandas as pd

from bench.frontends import LATENT


def export_csv(data, workdir, label_col="label", time_col="time"):
    os.makedirs(workdir, exist_ok=True)
    paths = {k: os.path.join(workdir, f"{k}.csv") for k in ["coords", "labels", "time"]}
    ids = data.obs_names.rename("cell_id")
    pd.DataFrame(np.asarray(data.obsm[LATENT]), index=ids).to_csv(paths["coords"])
    pd.DataFrame({"label": data.obs[label_col].to_numpy()}, index=ids).to_csv(paths["labels"])
    pd.DataFrame({"time": data.obs[time_col].to_numpy()}, index=ids).to_csv(paths["time"])
    return paths


def run_dpt(adata_path, outdir, start, start_time, seed=0, neighbors_key=None, **kw):
    cmd = [sys.executable, "trajectory/dpt.py", "--adata", adata_path,
           "--start", start, "--start-time", str(start_time),
           "--outdir", outdir, "--seed", str(seed)]
    if neighbors_key is not None:
        cmd += ["--neighbors-key", neighbors_key]
    for k, v in kw.items():
        cmd += [f"--{k.replace('_', '-')}", str(v)]
    subprocess.run(cmd, check=True)
    out = pd.read_csv(os.path.join(outdir, "dpt_pseudotime.csv"), index_col="cell_id", dtype={"cell_id": str})
    with open(os.path.join(outdir, "dpt_root.txt")) as f:
        root = f.read().strip()
    return out["dpt"], root
