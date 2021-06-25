#!/bin/bash 
#
#SBATCH --job-name=multiome_count
#SBATCH --output=04_data/custom_multiome_reference/logs/hg38_mm10_count-%j.out
#SBATCH --error=04_data/custom_multiome_reference/logs/hg38_mm10_count-%j.err
#SBATCH --time=48:00:00

# We want a new biochem node so we can have 2TB of local storage
#SBATCH --constraint=IB:HDR
#SBATCH --partition=biochem

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128GB

set -euxo pipefail

PROJ_DIR=/oak/stanford/groups/wjg/bparks/TF_ADTs_Amy/
SAMPLE_NAME=Dec2020_multiome_hg38_mm10

cd $L_SCRATCH_JOB

$PROJ_DIR/03_code/cellranger-arc-1.0.1/cellranger-arc count \
    --id="${SAMPLE_NAME}" \
    --reference="${PROJ_DIR}"/04_data/custom_multiome_reference/hg38_mm10 \
    --libraries="${PROJ_DIR}"/03_code/custom_multiome_reference/multi_samplesheet.csv \
    --localcores=16 \
    --localmem=100 

cp -r "${SAMPLE_NAME}/outs/" "${PROJ_DIR}/04_data/${SAMPLE_NAME}"