import os
import glob
import shutil
# import pandas 

configfile:"config/config.yaml"

READS_DIR=config['reads_dir']
ASSEMBLY_DIR=config['assembly_dir']
RESULTS_DIR=config['results_dir']
ENV_DIR=config['env_dir']
DB_DIR=config['db_dir']


SAMPLES = [line.strip() for line in open("config/sample_list.txt").readlines()]
SAMPLES_2=SAMPLES


rule all:
    input:
        expand(os.path.join(RESULTS_DIR,"Assembly/{sample}_filter.fasta.sa"), sample = SAMPLES),
        expand(os.path.join(RESULTS_DIR, "Bam/{sample}/{sample}_{sample_2}.bam"),sample = SAMPLES, sample_2 = SAMPLES_2),
        expand(os.path.join(RESULTS_DIR, "Bins/{sample}/Metabat"), sample = SAMPLES),
        expand(os.path.join(RESULTS_DIR, "Bins/{sample}/Concoct/concoct_clustering_merged.csv"),sample = SAMPLES),
        expand(os.path.join(RESULTS_DIR, "Bins/{sample}/Metabinner/metabinner_result.tsv"), sample = SAMPLES),
        expand(os.path.join(RESULTS_DIR, "Bins/{sample}/DasTool/concoct_das.tsv"), sample = SAMPLES),
        expand(os.path.join(RESULTS_DIR, "Bins/{sample}/DasTool/das_DASTool_summary.tsv"), sample = SAMPLES),
        expand(os.path.join(RESULTS_DIR, "Bins/dRep/dereplicated_genomes")),
        expand(os.path.join(RESULTS_DIR, "Bins/checkm_final")),
        expand(os.path.join(RESULTS_DIR, "Bins/gtdbtk_final"))


include: "rules/coverage.smk"
include: "rules/bins.smk"
include: "rules/dereplicate.smk"
include: "rules/taxqual.smk"
