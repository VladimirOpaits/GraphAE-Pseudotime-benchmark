rule score:
    input:
        pt="runs/{run}/pt/{path}/{backend}_pseudotime.csv",
        split="runs/fixed/clone_split.csv"
    output:
        "runs/{run}/score/{part}/{path}/{backend}_score.csv"
    wildcard_constraints:
        part="dev|test",
        backend="euclid|dpt|palantir|slingshot"
    params:
        n_boot=config["score"]["n_boot"]
    resources:
        mem_gb=3
    shell:
        "PYTHONPATH=. python workflow/scripts/score.py --pt {input.pt} --split {input.split} "
        "--part {wildcards.part} --n-boot {params.n_boot} --out {output}"


rule plot_pt:
    input:
        pt="runs/{run}/pt/{path}/{backend}_pseudotime.csv",
        meta="rt_data/stateFate_inVitro_metadata.txt.gz"
    output:
        "runs/{run}/plots/{path}/{backend}_pt.png"
    wildcard_constraints:
        backend="euclid|dpt|palantir|slingshot"
    resources:
        mem_gb=2
    shell:
        "python workflow/scripts/plot_pt.py --pt {input.pt} --meta {input.meta} --out {output}"
