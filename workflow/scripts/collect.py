import argparse
import re

import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--scores", nargs="+", required=True)
p.add_argument("--out", required=True)
a = p.parse_args()

rows = []
for path in a.scores:
    m = re.search(r"score/(?P<part>[^/]+)/(?P<fe>[^/]+)/(?P<backend>[^/]+)/(?:(?P<rest>.*)/)?[^/]+_score\.csv$", path)
    rest = m["rest"].split("/") if m["rest"] else []
    d = pd.read_csv(path)
    d.insert(0, "part", m["part"])
    d.insert(1, "frontend", m["fe"])
    d.insert(2, "backend", m["backend"])
    d.insert(3, "setting", "/".join(rest) if rest else "-")
    rows.append(d)
long = pd.concat(rows)
long.to_csv(a.out.replace(".csv", "_long.csv"), index=False)

key = ["part", "frontend", "backend", "setting", "branch"]
summary = long.groupby(key).agg(
    n_roots=("root", "size"), acc=("acc", "mean"), acc_sd=("acc", "std"),
    acc_lo=("acc_lo", "mean"), acc_hi=("acc_hi", "mean"),
    violation=("violation", "mean"), violation_sd=("violation", "std"),
    n_pairs=("n_pairs", "first"), n_nan=("n_nan", "max")).reset_index()
summary.to_csv(a.out, index=False)
print(summary[summary["branch"] == "mean"].round(3).to_string(index=False))
