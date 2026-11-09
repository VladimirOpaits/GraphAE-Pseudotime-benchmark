import argparse

import numpy as np

from bench.frontends import LATENT, load_frontend, make_frontend, save_frontend, smooth
from bench.graph import get_graph

p = argparse.ArgumentParser()
p.add_argument("--pca", required=True)
p.add_argument("--out", required=True)
p.add_argument("--k", type=int, required=True)
p.add_argument("--t", type=int, required=True)
a = p.parse_args()

data = load_frontend(a.pca)
_, conn = get_graph(data, a.k)
x = smooth(np.asarray(data.obsm[LATENT]), conn, t=a.t)
save_frontend(make_frontend(x, data.obs, method="pca_smooth", t=a.t, graph_k=a.k), a.out)
