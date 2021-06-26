# NEAT-seq reproducibility
Code to reproduce results from the NEAT-seq manuscript. (Link will be added upon publication)

## Getting started
The exact R package versions used for running analysis code are saved in `renv.lock`. These can be 
installed into a reproducible R environment using the `renv` package's `renv::restore()` function.

### Code layout
Each main text figure and the corresponding supplemental figure is in its own directory. Details on the code + data files are provided in the folders for each figure. All scripts should be run from the project root, rather than inside any of the figure folders.

### Running the code
Most of the code files run in <5 minutes, but can consume substantial RAM. It is recommended to run on a machine with >10GB of available memory.

**Library Installation**
- R: `install.packages("renv"); renv::restore()`
- python: `pip install pysam`

**Figure 1**
```shell
Rscript fig1_species_mixing/barnyard_analysis.R
```
**Figure 2**
```shell
Rscript fig2_CD4_Tcells/code/hematopoiesis_projection.R
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


## To-do checklist upon GEO publication

[ ] Write script/function to download files from GEO prior to loading  
[ ] Update code everywhere that loads GEO data  
[ ] Add instructions for accessing ATAC-seq bam files from SRA 
