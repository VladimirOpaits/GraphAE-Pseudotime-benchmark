import numpy as np
import pandas as pd
import scipy.io as sio

DATA = "rt_data/"


def load_larry():
    meta = pd.read_csv(DATA + "stateFate_inVitro_metadata.txt.gz", sep="\t")
    clones = sio.mmread(DATA + "stateFate_inVitro_clone_matrix.mtx.gz").tocsr()
    return meta, clones


def clone_pairs(day, clones):
    clones = clones.tocsc()
    ii, jj, cc = [], [], []
    for c in range(clones.shape[1]):
        cells = clones.indices[clones.indptr[c]:clones.indptr[c + 1]]
        if len(cells) < 2:
            continue
        d = day[cells]
        a, b = np.triu_indices(len(cells), 1)
        keep = d[a] != d[b]
        a, b = cells[a[keep]], cells[b[keep]]
        swap = day[a] > day[b]
        a, b = np.where(swap, b, a), np.where(swap, a, b)
        ii.append(a)
        jj.append(b)
        cc.append(np.full(len(a), c))
    return np.concatenate(ii), np.concatenate(jj), np.concatenate(cc)
