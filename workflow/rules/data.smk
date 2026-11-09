rule split_clones:
    input:
        meta="rt_data/stateFate_inVitro_metadata.txt.gz",
        clones="rt_data/stateFate_inVitro_clone_matrix.mtx.gz"
    output:
        "runs/fixed/clone_split.csv"
    params:
        dev_frac=config["split"]["dev_frac"],
        seed=config["split"]["seed"]
    shell:
        "python workflow/scripts/split_clones.py --meta {input.meta} --clones {input.clones} "
        "--out {output} --dev-frac {params.dev_frac} --seed {params.seed}"


rule roots:
    input:
        "runs/frontends/pca.h5ad"
    output:
        "runs/fixed/roots.txt"
    params:
        n_roots=config["roots"]["n_roots"],
        pool_size=config["roots"]["pool_size"],
        seed=config["roots"]["seed"]
    shell:
        "python trajectory/roots.py --adata {input} --out {output} "
        "--n-roots {params.n_roots} --pool-size {params.pool_size} --seed {params.seed}"


rule data_all:
    input:
        "runs/fixed/clone_split.csv",
        "runs/fixed/roots.txt"
