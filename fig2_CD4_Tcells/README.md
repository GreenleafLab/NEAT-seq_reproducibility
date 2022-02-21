# Figure 2 CD4 T cells
Code to reproduce the CD4 T cell analysis

## Data files
- ArchR project files at different filtering stages:
  1. `ArchR_HTOsinglets.rds` - Cells passing QC, and HTO doublets removed
  2. `ArchR_HTOsinglets_CD4only.rds` - HTO doublets removed, non-CD4 T cells removed
  3. `ArchR_HTOsinglets_CD4only_25XADT.rds` - HTO doublets removed, non-CD4 T cells removed, only cells with more concentrated ADT staining titration
- `hematopoiesies_peaks_hg38.bed`: Peaks from healthy hematopoiesis analysis in Granja, Klemm, McGinnis et al. 2019, lifted over from hg19 to hg38
- MPAL_alignment.tsv - Label transfers from Granga, Klemm, McGinnis et al. 2019 into the NEAT-seq dataset.
- helios_pwms.Rds - Helios PWMs from Transfac database

## Code files
- `hematopoiesis_projection.R`: Identify non-CD4 cell types using published hematopoiesis scATAC-seq data
    - Inputs: 
        - ATAC fragment files from GEO
        - Published dimensionality reduction from Granga, Klemm, McGinnis et al. 2019
    - Outputs:
        - `hematopoiesis_projection.tsv.gz`. Label transfers from Granga, Klemm, McGinnis et al. 2019 into the NEAT-seq dataset, including UMAP coordinates of the nearest neighbor in the healthy hematopoiesis UMAP.

- `CD4_HTO_singlet_ADT_counts.R`: Identifies singlets based on HTO and NPC hashing (not filtered for cells passing ATAC/RNA QC)
	- Inputs:
		- HTO and ADT counts files from GEO
	- Outputs:
		- TF_ADT_counts_singlets_from_NPCandHTO.csv

- `CD4_ArchR_plots.R`: Creates ArchR project of cells passing QC and with HTO doublets removed and plots all UMAP plots except ADT UMAP plots
	- Inputs:
		- ATAC fragment files from GEO
		- TF_ADT_counts_singlets_from_NPCandHTO.csv (output of CD4_HTO_singlet_ADT_counts.R)
		- MPAL_alignment.tsv from fig2_CD4_Tcells/data folder
		- helios_pwms.Rds from fig2_CD4_Tcells/data folder
		- CD4_RNA_counts.rds from GEO
	- Outputs:
		- directory with ArchR project (ArchR_HTOsinglets_CD4only)
		- PDFs of UMAP plots

- `ADT_normalization.R`: Creates tables of normalized ADT counts for the two antibody panel concentrations used (not filtered for cells passing QC)
	- Inputs:
		- TF_ADT_counts_singlets_from_NPCandHTO.csv (output of CD4_HTO_singlet_ADT_counts.R)
	- Outputs:
		- ADT_5x_NPCnorm.csv for "Concentration 1" in manuscript
		- ADT_25x_NPCnorm.csv for "Concentration 2" in manuscript

- `ArchR_CD4cells_add25xADT.R`: Creates ArchR project of CD4 singlets stained with antibody panel Concentration 2 and plots ADT UMAP and density plots
	- Inputs:
		- ArchR_HTOsinglets_CD4only ArchR project
		- ADT_25x_NPCnorm.csv from ADT_normalization.R
	- Outputs:
		- directory with ArchR project (ArchR_HTOsinglets_CD4only_ADT)
		- CD4cells_RNA_25XADT.csv with normalized ADT data for CD4 singlets
		- PDFs of UMAP and density plots

- `ADT_vs_RNA_correlations.R`: Plots scatterplots of ADT vs RNA data across CD4 singlets
	- Inputs:
		- CD4cells_RNA_25XADT.csv from ArchR_CD4cells_add25xADT.R
	- Outputs:
		- PDF of scatterplots

- `Seurat_markers.R`: Creates Seurat object of CD4 singlets and plots CD4 T cell subset marker expression
	- Inputs:
		- filtered feature barcode matrix files from GEO
		- CD4cells.csv from GEO
	- Outputs:
		- Tmem_allCD4cells.rds Seurat object
		- PDF of RNA UMAP colored by ATAC-defined clusters
		- PDf of marker expression in ATAC clusters

- `Seurat_GATA3_differential_analysis.R`: Performs differential RNA analysis on GATA3 high RNA and high vs low ADT cells
	- Inputs:
		- Tmem_allCD4cells.rds Seurat object from output of Seurat_markers.R
		- CD4cells_RNA_25XADT.csv from ArchR_CD4cells_add25xADT.R
	- Outputs:
		- Tmem_CD4cells_withADT.rds containing Seurat object with all cells stained with antibody panel Concentration 2
		- PDF of scatterplot of GATA3 ADT vs RNA counts with high and low cutoffs indicated
		- PDF of differential expression analysis volcano plot
