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
module load samtools/1.6
# deeptools dependencies
module load bedtools
module load BigWig_tools
module load bx-python
module load deepTools/2.4.2

# The snakemake file with pipeline rules
SM=~/Programming/pyNextGen/snake_make_recipes/pre_bam.sm 

# Cluster specification for each job
CLUS_CONF=~/Programming/pyNextGen/snake_make_recipes/cluster.conf

# This configuration file can be create with config_generator.py
SM_JSON=snakemake.json

# Logs directory (HAS TO FINISH WITH "/")
#LOGS_DIR=~/logs/snakemake/$(date +%Y-%m-%d-%H-%M-%S)/;
# Another Good possible locations for logs could be just the working directory
LOGS_DIR=logs_cluster/X$(date +%Y-%m-%d-%H-%M-%S)/;
mkdir $LOGS_DIR

snakemake -j 12 --printshellcmds --cluster-config $CLUS_CONF --cluster "sbatch --qos {cluster.qos} --ntasks {cluster.tasks} --cpus-per-task {cluster.cpus} --job-name {cluster.jobName} --output $LOGS_DIR{cluster.output} --error $LOGS_DIR{cluster.error}" -s $SM --configfile $SM_JSON
