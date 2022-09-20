# GTDBTK taxonomy
rule gtdbtk:
    input:
        rules.dRep.output.final
    output:
        directory(os.path.join(RESULTS_DIR, "Bins/gtdbtk_final"))
    log:
        os.path.join(RESULTS_DIR, "logs/gtdbtk.log")
    conda:
        os.path.join(ENV_DIR, "gtdbtk.yaml")
    params:
        config["gtdbtk"]["path"]
    threads:
        config["gtdbtk"]["threads"]
    message:
        "Running GTDB on MAGs"
    shell:
        "(date && export GTDBTK_DATA_PATH={params} && gtdbtk classify_wf --cpus {threads} -x fa --genome_dir {input} --out_dir {output} && date) &> >(tee {log})"

# Checking bin quality
rule checkm_final:
    input:
        rules.dRep.output.final
    output:
        directory(os.path.join(RESULTS_DIR, "Bins/checkm_final"))
    conda:
        os.path.join(ENV_DIR, "checkm.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/checkm_final/checkm.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/checkm_final/checkm.err.log")
    threads:
        config["checkM"]["threads"]
    message:
        "Running Final Checkm on dereplicated output"
    shell:
        """
        (date && checkm lineage_wf -t {threads} -f {output}.tsv --tab_table -x fa {input} {output} && date) 2> {log.err} > {log.out}
        """
