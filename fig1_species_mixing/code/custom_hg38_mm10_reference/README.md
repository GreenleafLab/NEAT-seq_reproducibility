This directory contains the Snakemake pipeline used to construct the combined
hg38 & mm10 reference. 

The code is not possible to run from scratch as it contains paths which are not
present in this repository, however it should demonstrate the parameters used
as well URLs for the original data sources

## Files
- `cellranger_count.sbatch`: `cellranger-arc count` command and parameters
- `hg38_mm10.config`: config file for `cellranger-arc mkref`
- `multi_samplesheet.csv`: Samplesheet for `cellranger-arc count`
- `Snakefile`: Snakemake pipeline to generate chimeric reference genome based on the published 10x genomes `refdata-cellranger-arc-GRCh38-2020-A` and `refdata-cellranger-arc-mm10-2020-A`