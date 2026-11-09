rule smooth_pca:
    input:
        "runs/frontends/pca.h5ad"
    output:
        "runs/frontends/pca_smooth_t{t}_k{kf}.h5ad"
    wildcard_constraints:
        t=r"\d+",
        kf=r"\d+"
    resources:
        mem_gb=3
    shell:
        "PYTHONPATH=. python workflow/scripts/smooth_pca.py --pca {input} --out {output} "
        "--k {wildcards.kf} --t {wildcards.t}"
