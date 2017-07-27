#!/bin/bash
#SBATCH --job-name masterSnake
#SBATCH --output job_masterSnake_%j.out
#SBATCH --error job_masterSnake_%j.err
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --qos long

module load snakemake/3.5.4
module load fastqc/0.11.5
module load STAR/2.5.0a
module load subread/1.4.6-p3
module load MultiQC/0.7

# The snakemake file with pipeline rules
SM=pre_bam.sm

# Cluster specification for each job
CLUS_CONF=cluster.conf

# This configuration file can be create with config_generator.py
SM_JSON=snakemake_tars_conf_example.json

snakemake -j 12 --cluster-config $CLUS_CONF --cluster "sbatch -t {cluster.time}" -s $SM --configfile $SM_JSON
