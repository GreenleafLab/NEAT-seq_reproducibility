# Figure 2 CD4 T cells
Code to reproduce the CD4 T cell analysis

## Data files
- ArchR project files at different filtering stages:
  1. `ArchR_HTOsinglets.rds` - Cells passing QC, and HTO doublets removed
  2. `ArchR_HTOsinglets_CD4only.rds` - HTO doublets removed, non-CD4 T cells removed
  3. `ArchR_HTOsinglets_CD4only_25XADT.rds` - HTO doublets removed, non-CD4 T cells removed, only cells with more concentrated ADT staining titration
- `hematopoiesies_peaks_hg38.bed`: Peaks from healthy hematopoiesis analysis in Granja, Klemm, McGinnis et al. 2019, lifted over from hg19 to hg38

## Code files
- `hematopoiesis_projection.R`: Identify non-CD4 cell types using published hematopoiesis scATAC-seq data
    - Inputs: 
        - ATAC fragment files from GEO
        - Published dimensionality reduction from Granga, Klemm, McGinnis et al. 2019
    - Outputs:
        - `hematopoiesis_projection.tsv.gz`. Label transfers from Granga, Klemm, McGinnis et al. 2019 into the NEAT-seq dataset, including UMAP coordinates of the nearest neighbor in the healthy hematopoiesis UMAP.