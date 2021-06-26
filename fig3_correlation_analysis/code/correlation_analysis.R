# Runs in 2 min, 40s taking a maximum of 10GB of memory
suppressPackageStartupMessages({
  library(ArchR)
  
   library(chromVARmotifs)
   library(TFBSTools)
   library(BSgenome.Hsapiens.UCSC.hg38)
   library(motifmatchr)
  
   library(tidyverse)
})

input_data_path_matrices <- "/oak/stanford/groups/wjg/bparks/TF_ADTs_Amy/GEO_upload/CD4_Tcell/"
input_data_path <- "/oak/stanford/groups/wjg/amyfchen/GEO_submission_June2021"
output_path <- "fig3_correlation_analysis/outputs"

source("fig3_correlation_analysis/code/correlation_utils.R")
Rcpp::sourceCpp("fig3_correlation_analysis/code/rcpp_utils.cc")

tf_ensg_ids <- c(
  "RORgT" = "ENSG00000143365",
  "Tbet" = "ENSG00000073861",
  "FOXP3" = "ENSG00000049768",
  "GATA3" = "ENSG00000107485",
  "Helios" = "ENSG00000030419"
)

proj <- readRDS("fig2_CD4_Tcells/data/ArchR_HTOsinglets_CD4only_25XADT.rds")
proj_all <- readRDS("fig2_CD4_Tcells/data/ArchR_HTOsinglets_CD4only.rds")

rna_counts <- readRDS(file.path(input_data_path_matrices, "CD4_RNA_counts.rds"))
peak_counts <- readRDS(file.path(input_data_path_matrices, "CD4_Peak_matrix.rds"))

rna <- rna_counts %>%
  t() %>% {log10(10000 * . / rowMeans(.) + 1)} %>% as("dgCMatrix") %>%
  .[,colSums(.) > 0]
peaks <- peak_counts / proj_all$ReadsInTSS

adt <- bind_rows(
    lane1 = read_csv(file.path(input_data_path, "CD4_lane1/ADT_counts_lane1.csv")),
    lane2 = read_csv(file.path(input_data_path, "CD4_lane2/ADT_counts_lane2.csv")),
    .id = "lane"
) %>%
  mutate(cell_id = str_c(lane, "#", cell, "-1")) %>%
  mutate(across(c("FOXP3", "Tbet", "GATA3", "Helios", "RORgT"), ~ log10(250*.x/pmax(NPC1,NPC2) + 1))) %>%
  select(cell_id, c("FOXP3", "Tbet", "GATA3", "Helios", "RORgT")) %>%
  column_to_rownames("cell_id") %>%
  as.matrix() %>%
  .[proj$cellNames,]

peak_motif_matches <-  list(
  "RORgT"=human_pwms_v2$ENSG00000143365_LINE3014_RORC_D,
  "Tbet"=human_pwms_v2$ENSG00000073861_LINE3513_TBX21_D_N2,
  "FOXP3"=human_pwms_v2$ENSG00000049768_LINE1933_FOXP3_D_N1,
  "GATA3"=human_pwms_v2$ENSG00000107485_LINE2094_GATA3_D_N6,
  "Helios" = load_pwm("fig3_correlation_analysis/data/Helios_PFM_1.tsv", "Helios")
) %>% do.call(PWMatrixList, .) %>%
  motifmatchr::matchMotifs(getPeakSet(proj), genome=BSgenome.Hsapiens.UCSC.hg38, out="matches") %>%
  assay() %>% 
  set_rownames(as.character(as.character(getPeakSet(proj)))) %>%
  .[colnames(peaks),] %>%
  as.matrix() %>% as_tibble(rownames = "peak_id") %>%
  pivot_longer(!peak_id, names_to="adt", values_to="has_motif")
  
gene_names <- read_tsv(
    file.path(input_data_path, "CD4_lane1/filtered_feature_bc_matrix/lane1_features.tsv.gz"),
    col_names=c("gene_id", "gene_name", "feature_type", "chr", "start", "end"),
    col_types="ccccii"
  ) %>%
  filter(feature_type == "Gene Expression") %>%
  pull(gene_name, name=gene_id)

gene_coords <- fread("fig3_correlation_analysis/data/genes.gtf.gz",
                           skip="chr", col.names = c("seqname", "source", "feature", "start", "end",
                                                     "score", "strand", "frame", "attributes")) %>%
  as_tibble() %>%
  filter(feature == "gene") %>%
  mutate(gene_id = str_match(attributes, 'gene_id "([^"]*)"')[,2]) %>%
  filter(gene_id %in% colnames(rna)) %>%
  column_to_rownames("gene_id") %>%
  makeGRangesFromDataFrame() %>%
  .[colnames(rna)]


expressed_genes <- colSums(rna[rownames(adt),] > 0) > 10 | colnames(rna) %in% tf_ensg_ids
accessible_peaks <- colSums(peaks[rownames(adt),] > 0) > 10


set.seed(125134)
sample_order_all <- sample(seq_along(proj_all$cellNames))
pseudobulk_mat <- knn_pseudobulk(getReducedDims(proj_all), starter_cells = sample_order_all, k=100, max_overlap = 0.8) %>%
  .[1:500,]


corr_rna <- correlation_test(
    rna[rownames(adt), expressed_genes], adt, 
    method="spearman", use_permutations=FALSE
  ) %>%
  dplyr::rename(gene_id=feature_x, adt=feature_y) %>%
  group_by(adt) %>%
  mutate(
    p.ttest.adj = p.adjust(p.ttest, method = "BH"),
    avg_expression = colMeans(rna)[gene_id],
    gene_name = gene_names[gene_id]
  ) %>%
  ungroup() %>%
  select(gene_name, everything())

corr_peak <- correlation_test(
    peaks[rownames(adt), accessible_peaks], adt, 
    method="spearman", use_permutations=FALSE
  ) %>%
  dplyr::rename(peak_id=feature_x, adt=feature_y) %>%
  group_by(adt) %>%
  mutate(
    p.ttest.adj = p.adjust(p.ttest, method = "BH"),
    avg_accessibility = colMeans(peaks)[peak_id]
  ) %>%
  ungroup() %>%
  inner_join(peak_motif_matches, by=c("peak_id", "adt")) %>%
  select(peak_id, everything())

# Peak-gene correlations
corr_gene_peak <- findOverlaps(
  promoters(gene_coords, upstream=100000, downstream=100000),
  as(colnames(peaks), "GRanges")
) %>% as_tibble() %>%
  dplyr::rename(gene_idx = queryHits, peak_idx = subjectHits) %>%
  mutate(
    gene_id = colnames(rna)[gene_idx],
    peak_id = colnames(peaks)[peak_idx],
    distance = distance(
      resize(gene_coords, 1, fix="start")[gene_idx], 
      resize(as(colnames(peaks), "GRanges"), 1, fix="center")[peak_idx]
    ),
    spearman = sparse.cor.partial(
      as.matrix(pseudobulk_mat %*% rna), as.matrix(pseudobulk_mat %*% peaks), 
      gene_idx, peak_idx, method="spearman"
    ),
    p.ttest = cor.pval(spearman, nrow(pseudobulk_mat)),
    p.ttest.adj = p.adjust(p.ttest, method="BH")
  )

top_peaks <- corr_peak %>%
  filter(has_motif) %>%
  group_by(adt) %>%
  filter(spearman > quantile(spearman, 0.8))
top_genes <- corr_rna %>%
  group_by(adt) %>%
  filter(spearman > quantile(spearman, 0.8))

top_links <- corr_gene_peak %>%
  select(!c(gene_idx, peak_idx)) %>%
  filter(p.ttest.adj < 0.05, spearman > 0) %>%
  inner_join(top_genes, 
             by=c("gene_id"), suffix=c(".peak_gene", "")) %>%
  inner_join(top_peaks, 
             by=c("peak_id", "adt"), suffix=c(".gene_adt", ".peak_adt")) %>%
  as.data.table()

dir.create(output_path, showWarnings=FALSE)
fwrite(corr_peak, file.path(output_path, "corr_peak.csv.gz"), showProgress = FALSE, compress="gzip")
fwrite(corr_rna, file.path(output_path, "corr_rna.csv.gz"), showProgress = FALSE, compress="gzip")
fwrite(corr_gene_peak, file.path(output_path, "corr_gene_peak.csv.gz"), showProgress = FALSE, compress="gzip")
fwrite(top_links, file.path(output_path, "top_links.csv.gz"), showProgress = FALSE, compress="gzip")

