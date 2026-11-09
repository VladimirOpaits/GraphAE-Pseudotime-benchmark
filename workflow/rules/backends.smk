wildcard_constraints:
    run=r"[^/]+",
    fe=r"[^/]+",
    k=r"\d+",
    n=r"\d+",
    res=r"[0-9.]+"


rule euclid:
    input:
        adata="runs/frontends/{fe}.h5ad",
        roots="runs/fixed/roots.txt"
    output:
        "runs/{run}/pt/{fe}/euclid/euclid_pseudotime.csv"
    params:
        n_roots=config["roots"]["n_used"]
    resources:
        mem_gb=2
    shell:
        "python trajectory/euclid.py --adata {input.adata} --roots $(head -n {params.n_roots} {input.roots}) "
        "--outdir $(dirname {output})"


rule dpt:
    input:
        adata="runs/frontends/{fe}.h5ad",
        roots="runs/fixed/roots.txt"
    output:
        "runs/{run}/pt/{fe}/dpt/k{k}/dpt_pseudotime.csv"
    params:
        n_roots=config["roots"]["n_used"],
        n_dcs=config["dpt"]["n_dcs"]
    resources:
        mem_gb=3
    shell:
        "python trajectory/dpt.py --adata {input.adata} --roots $(head -n {params.n_roots} {input.roots}) "
        "--n-neighbors {wildcards.k} --n-dcs {params.n_dcs} --outdir $(dirname {output})"


rule palantir:
    input:
        adata="runs/frontends/{fe}.h5ad",
        roots="runs/fixed/roots.txt"
    output:
        "runs/{run}/pt/{fe}/palantir/k{k}/palantir_pseudotime.csv"
    params:
        n_roots=config["roots"]["n_used"],
        n_components=config["palantir"]["n_components"],
        waypoints=config["palantir"]["num_waypoints"],
        seed=config["palantir"]["seed"]
    threads: 8
    resources:
        mem_gb=5
    shell:
        "python trajectory/palantir_pt.py --adata {input.adata} --roots $(head -n {params.n_roots} {input.roots}) "
        "--n-neighbors {wildcards.k} --n-components {params.n_components} "
        "--num-waypoints {params.waypoints} --seed {params.seed} --n-jobs {threads} "
        "--outdir $(dirname {output})"


rule slingshot:
    input:
        adata="runs/frontends/{fe}.h5ad",
        roots="runs/fixed/roots.txt"
    output:
        "runs/{run}/pt/{fe}/slingshot/k{k}/res{res}/slingshot_pseudotime.csv"
    params:
        n_roots=config["roots"]["n_used"],
        cluster_by=config["slingshot"]["cluster_by"],
        seed=config["slingshot"]["seed"]
    resources:
        mem_gb=4
    shell:
        "python trajectory/slingshot.py --adata {input.adata} --roots $(head -n {params.n_roots} {input.roots}) "
        "--n-neighbors {wildcards.k} --resolution {wildcards.res} --cluster-by {params.cluster_by} "
        "--seed {params.seed} --outdir $(dirname {output})"
