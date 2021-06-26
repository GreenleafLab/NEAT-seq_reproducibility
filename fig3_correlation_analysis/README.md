# Figure 3 Correlation analysis
Code to reproduce the results from the correlation analysis and browser track plots in Figures 3b-d and S3a-c.

## Data files
- `genes.gtf.gz`: Gene coordinates from the `cellranger-arc` reference `refdata-cellranger-arc-GRCh38-2020-A`
- `Helios_PFM_1.tsv`: ChIP-seq derived motif from TransFac, accession M09745. Derived from ENCODE dataset ENCSR680UQE

Note that the figure 3 code also relies on some of the data files from the figure 2 data folder.

## Code files
- `correlation_analysis.R`: Re-generate correlation tables from GEO
    - Inputs: RNA, ADT, and ATAC counts matrices from GEO. Filtered ArchR projects from figure 2 data
    - Outputs: 
        - `corr_peak.csv.gz`: Peak-TF correlations
        - `corr_rna.csv.gz`: RNA-TF correlations
        - `corr_gene_peak.csv.gz`: Peak-RNA correlations
        - `top_links.csv.gz`: Putative direct targets with Peak-RNA-TF correlations
- `correlation_plots.R`: Re-generate figures S3a+b
    - Inputs: Correlation tables from GEO
    - Outputs: 
        - `gene_adt_corr_plots.pdf`: RNA-TF correlation plots
        - `peak_adt_corr_plots.pdf`: Peak-TF correlation plots
- `correlation_utils.R`: Functions for Spearman correlation statistics on sparse matrices, pseudobulk calculations, MAGIC imputation, and loading of PWM and GTF files.
- `extract_snp_information.py`:
    - Inputs: Raw ATAC-seq BAM files from GEO/SRA upload
    - Outputs:
        - `lane1_fragments.tsv` and `lane2_fragments.tsv`: Fragment files for all reads overlapping SNP rs62088464, annotated by which allele was sequenced.
- `heatmaps.R`: Recreate heatmaps from figure 3b
    - Inputs: RNA, Peak, and ADT matrices from GEO, along with putative direct targets table from GEO
    - Outputs: 
        - `heatmap_plot.pdf`: Heatmaps from figure 3b
- `rcpp_utils.cc`: C++ code required for certain functionality in `correlation_utils.R`
- `snp_plots.R`: Recreate allele-specific accessibility barplot from figure 3d
    - Inputs: `lane1_fragments.tsv` and `lane2_fragments.tsv`
    - Outputs:
        - `snp_plots.pdf`: Barplot from figure 3d
- `trackplot_utils.R`: Lightweight, modular library for constructing aligned trackplots using the `patchwork` framework.
- `trackplots.R`: Recreate trackplots from figures 3c+d, and S3c
    - Inputs: RNA, Peak, and ADT matrices from GEO, along with putative direct targets table from GEO
    - Outputs: 
        - `trackplots.pdf`: Plots from figures 3c+d and S3c
