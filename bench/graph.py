import scanpy as sc

from bench.frontends import LATENT


def graph_key(k):
    return f"k{k}"


def add_graph(data, k):
    key = graph_key(k)
    if key not in data.uns:
        sc.pp.neighbors(data, n_neighbors=k, use_rep=LATENT, key_added=key, transformer="sklearn")
    return key


def get_graph(data, k):
    key = add_graph(data, k)
    return data.obsp[f"{key}_distances"], data.obsp[f"{key}_connectivities"]
