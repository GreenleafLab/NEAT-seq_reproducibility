#!/bin/bash 
#
#SBATCH --job-name=multiome_count
#SBATCH --output=revisions/outputs/logs/GM12878_count-%j.out
#SBATCH --error=revisions/outputs/logs/GM12878_count-%j.err
#SBATCH --time=48:00:00

# We want a new biochem node so we can have 2TB of local storage
#SBATCH --constraint=IB:HDR
#SBATCH --partition=biochem

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128GB



# This is the code used for the GM12878 multiome data processing, but it takes
# many hours to run and so the code is just for illustrative purposes, not set up to
# run right out of the box.
# The only important output from this is the filtered_feature_bc_matrix which is
# provided in the code repository as input data.

set -euxo pipefail

PROJ_DIR=/oak/stanford/groups/wjg/bparks/NEAT-seq_reproducibility/

CELLRANGER_PATH=/oak/stanford/groups/wjg/bparks/multiome/03_code/cellranger_arc_v1/installation/cellranger-arc-1.0.1/cellranger-arc
REF_DIR=/oak/stanford/groups/wjg/bparks/multiome/03_code/cellranger_arc_v1/installation/refdata-cellranger-arc-GRCh38-2020-A

SAMPLE_NAME=GM12878

cd $L_SCRATCH_JOB

$CELLRANGER_PATH count \
    --id="${SAMPLE_NAME}" \
    --reference=$REF_DIR \
    --libraries="${PROJ_DIR}"/revisions/code/GM12878_cellranger/sample_sheet.csv \
    --localcores=16 \
    --localmem=100 

cp -r "${SAMPLE_NAME}/outs/" "${PROJ_DIR}/revisions/outputs/${SAMPLE_NAME}"