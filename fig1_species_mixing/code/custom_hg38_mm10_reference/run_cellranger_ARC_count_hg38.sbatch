#!/bin/bash 
#
# bash script to run cellranger-arc count on sherlock
# Assumes you are in the run directory when running
# Usage: sbatch run_cellranger_ARC_count.sbatch <sample_id> <fastq_path> <ref_genome_path>
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling  
#################
#set a job name  
#SBATCH --job-name=hg38_count
#################  
#a file for job output, you can check job progress
#SBATCH --output=hg38_count.out
#################
# a file for errors from the job
#SBATCH --error=hg38_count.err
#################
#time you think you need; default is one hour
#in minutes in this case, hh:mm:ss
#SBATCH --time=48:00:00
#################
#quality of service; think of it as job priority
#SBATCH --partition=wjg,biochem,sfgf,owners
#SBATCH --qos=normal
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#################
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#################

srun cellranger-arc count --id=multiome_counts_hg38 --reference=/oak/stanford/groups/wjg/amyfchen/cellranger/cellranger_multiome_software/genomes/refdata-cellranger-arc-GRCh38-2020-A --libraries multi_samplesheet.csv --localcores=20
