# Supplementary Figures
Code to reproduce supplementary figures not generated in code from main figures


## Code files
- ADT_normalization_tests.R: Comparison of different normalization methods for ADT data
	- Inputs:
		- TF_ADT_counts_singlets_from_NPCandHTO.csv (output of CD4_HTO_singlet_ADT_counts.R)
		- CD4cells_RNA_25XADT.csv from ArchR_CD4cells_add25xADT.R
		- ArchR_HTOsinglets_CD4only_ADT ArchR project from ArchR_CD4cells_add25xADT.R
	- Outputs:
		- PDF of UMAP plots of ADT levels after each normalization strategy 

- ArchR_CD4cells_add5xADT_250.R:  Creates ArchR project of CD4 singlets stained with antibody panel Concentration 1 and plots ADT density plots
	- Inputs:
		- ArchR_HTOsinglets_CD4only ArchR project
		- ADT_5x_NPCnorm.csv from ADT_normalization.R
	- Outputs:
		- CD4cells_RNA_5XADT.csv with normalized ADT data for CD4 singlets
		- PDFs of density plots

- peak_gene_GO.R: GO enrichment analysis for candidate TF-driven peak-gene linkages
	- Inputs:
		- putative_targets.csv.gz from GEO
		- CD4_RNA_counts.rds from GEO
	- Outputs:
		- PDF of enriched biological process terms from clusterProfiler