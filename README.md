# NEAT-seq reproducibility
Code to reproduce results from the NEAT-seq manuscript. (Link will be added upon publication)

## Getting started
The exact R package versions used for running analysis code are saved in `renv.lock`. These can be 
installed into a reproducible R environment using the `renv` package's `renv::restore()` function.

Input data from GEO should be downloaded before running any analysis scripts. This can be done
with the script `code_utils/download_data.py` and will download about 1.5GB of data

### Code layout
Each main text figure and the corresponding supplemental figure is in its own directory. Details on the code + data files are provided in the folders for each figure. All scripts should be run from the project root, rather than inside any of the figure folders.

### Running the code
Most of the code files run in <5 minutes, but can consume substantial RAM. It is recommended to run on a machine with >10GB of available memory.

**Library Installation**
- R: `install.packages("renv"); renv::restore()`
    - Additional packages that are required for certian scripts:
      `Seurat`, `immunogenomics/presto`, `limma`, `seriation`, `org.Hs.eg.db`, `clusterProfiler` (may require `rvcheck@1.0.8` for clusterProfiler install to work)
- python: `pip install pysam snakemake`
    - Some ArchR analysis further requires `macs2`, though it may be easiest
      to `pip install macs3` then alias `macs2` to point to `macs3`
- other: `conda install samtools`


**Data download**
```shell
python code_utils/download_data.py
```

**Figure 1**
```shell
Rscript fig1_species_mixing/code/barnyard_analysis.R
```
**Figure 2**
```shell
Rscript fig2_CD4_Tcells/code/hematopoiesis_projection.R
Rscript fig2_CD4_Tcells/code/CD4_HTO_singlet_ADT_counts.R
Rscript fig2_CD4_Tcells/code/CD4_ArchR_plots.R # This is the slow step
Rscript fig2_CD4_Tcells/code/ADT_normalization.R
Rscript fig2_CD4_Tcells/code/ArchR_CD4cells_add25xADT.R
Rscript fig2_CD4_Tcells/code/ADT_vs_RNA_correlations.R
Rscript fig2_CD4_Tcells/code/Seurat_markers.R
Rscript fig2_CD4_Tcells/code/Seurat_GATA3_differential_analysis.R
```
**Figure3**
```shell
Rscript fig3_correlation_analysis/code/correlation_analysis.R
Rscript fig3_correlation_analysis/code/correlation_plots.R

python fig3_correlation_analysis/code/extract_snp_information.py
Rscript fig3_correlation_analysis/code/snp_plots.R

Rscript fig3_correlation_analysis/code/trackplots.R
Rscript fig3_correlation_analysis/code/heatmaps.R
```

**Revisions**
```shell
snakemake -c5 -s Supplementary_figures/code/K562_bulk_ATAC/K562_download.snake
Rscript Supplementary_figures/code/bulk_correlation.R

Rscript Supplementary_figures/code/ADT_normalization_tests.R
Rscript Supplementary_figures/code/ArchR_CD4cells_add5xADT_250.R
Rscript Supplementary_figures/code/peak_gene_GO.R
```

## To-do checklist upon SRA publication
[ ] Add instructions for accessing ATAC-seq bam files from SRA 
