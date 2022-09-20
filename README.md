# Uganda MAGs

## Pipeline description

This pipeline starts with different individual assemblies (`fasta` files) and their respective reads (`fastq.gz` files). It will map all reads against all the assemblies using [BWA](https://github.com/lh3/bwa).
Then per sample, [Metabat2](https://bitbucket.org/berkeleylab/metabat/src), [Concoct](https://github.com/BinPro/CONCOCT) and [Metabinner](https://github.com/ziyewang/MetaBinner) are ran to create genomic bins. These bins are then consolidated using [DasTool](https://github.com/cmks/DAS_Tool). After that, the bins of all samples are dereplicated with [dRep](https://github.com/MrOlm/drep) to form MAGs. [CheckM](https://github.com/Ecogenomics/CheckM) and [GtdbTk](https://github.com/Ecogenomics/GTDBTk) are then used to estimate the quality and taxonomy of the MAGs

## Setup

* Place your preprocessed/trim reads (e.g. `sample_R2.fastq.gz` and `sample_R2.fastq.gz` files) in a `reads` folder
* Place the individual assemblies (e.g. `sample.fasta`) into an `assembly` folder
* Modify the `config.yaml` file to change the different paths and eventually the different options
* Modify the `sample_list.txt` file to include your samples

### Without Slurm

`snakemake -rp --cores 28 --use-conda --rerun-incomplete`

### With Slurm

This part was mainly taken from [@susheelbhanu](https://github.com/susheelbhanu/) [nomis_pipeline](https://github.com/susheelbhanu/nomis_pipeline)

* Modify the `slurm.yaml` file by checking `partition`, `qos` and `account` that heavily depends on your system
* Modify the `sbatch.sh` file by checking `#SBATCH -p`, `#SBATCH --qos=` and `#SBATCH -A` options that heavily depends on your system

`sbatch config/sbatch.sh`
