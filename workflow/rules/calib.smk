def calib_scores():
    c, run = config["calib"], config["run"]
    base = f"runs/{run}/score/dev/{c['frontend']}"
    out = [f"{base}/euclid/euclid_score.csv"]
    out += [f"{base}/dpt/k{c['k']}/dpt_score.csv"]
    out += [f"{base}/palantir/k{c['k']}/palantir_score.csv"]
    out += [f"{base}/slingshot/k{c['k']}/res{r}/slingshot_score.csv" for r in c["slingshot_resolution"]]
    return out


rule calib_summary:
    input:
        calib_scores()
    output:
        f"runs/{config['run']}/calib_summary.csv"
    shell:
        "python workflow/scripts/collect.py --scores {input} --out {output}"


rule calib_all:
    input:
        rules.calib_summary.output


rule calib_plots:
    input:
        [p.replace("/score/dev/", "/plots/").replace("_score.csv", "_pt.png") for p in calib_scores()]


def kcross_scores():
    run, ks = config["run"], config["kcross"]["k"]
    fe = "pca_smooth_t{t}_k{kf}".replace("{t}", str(config["kcross"]["t"]))
    return [f"runs/{run}/score/dev/{fe.replace('{kf}', str(kf))}/dpt/k{kb}/dpt_score.csv"
            for kf in ks for kb in ks]


rule kcross_summary:
    input:
        kcross_scores()
    output:
        f"runs/{config['run']}/kcross_summary.csv"
    shell:
        "python workflow/scripts/collect.py --scores {input} --out {output}"


rule palantir_calib:
    input:
        adata="runs/frontends/{fe}.h5ad",
        roots="runs/fixed/roots.txt"
    output:
        "runs/{run}/pt/{fe}/palantir/k{k}/nc{n}/palantir_pseudotime.csv"
    params:
        n_roots=config["roots"]["n_used"],
        waypoints=config["palantir"]["num_waypoints"],
        seed=config["palantir"]["seed"]
    threads: 8
    resources:
        mem_gb=5
    shell:
        "python trajectory/palantir_pt.py --adata {input.adata} --roots $(head -n {params.n_roots} {input.roots}) "
        "--n-neighbors {wildcards.k} --n-components {wildcards.n} "
        "--num-waypoints {params.waypoints} --seed {params.seed} --n-jobs {threads} "
        "--outdir $(dirname {output})"
