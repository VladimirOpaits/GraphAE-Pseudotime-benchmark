import os
import subprocess
import sys

import pandas as pd


def write_inputs(workdir, latent, cell_ids, labels, time):
    os.makedirs(workdir, exist_ok=True)
    paths = {k: os.path.join(workdir, f"{k}.csv") for k in ["coords", "labels", "time"]}
    pd.DataFrame(latent, index=cell_ids).rename_axis("cell_id").to_csv(paths["coords"])
    pd.DataFrame({"label": labels}, index=cell_ids).rename_axis("cell_id").to_csv(paths["labels"])
    pd.DataFrame({"time": time}, index=cell_ids).rename_axis("cell_id").to_csv(paths["time"])
    return paths


def run_dpt(workdir, latent, cell_ids, labels, time, start, start_time, seed=0, root_seed=None, graph=None, **kw):
    paths = write_inputs(workdir, latent, cell_ids, labels, time)
    cmd = [sys.executable, "trajectory/dpt.py",
           "--coords", paths["coords"], "--labels", paths["labels"], "--time", paths["time"],
           "--start", start, "--start-time", str(start_time),
           "--outdir", workdir, "--seed", str(seed)]
    if root_seed is not None:
        cmd += ["--root-seed", str(root_seed)]
    if graph is not None:
        cmd += ["--graph", graph]
    for k, v in kw.items():
        cmd += [f"--{k.replace('_', '-')}", str(v)]
    subprocess.run(cmd, check=True)
    out = pd.read_csv(os.path.join(workdir, "dpt_pseudotime.csv"), index_col="cell_id", dtype={"cell_id": str})
    with open(os.path.join(workdir, "dpt_root.txt")) as f:
        root = f.read().strip()
    return out["dpt"].reindex(cell_ids).to_numpy(), root
