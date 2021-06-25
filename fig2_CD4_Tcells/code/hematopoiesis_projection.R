# Takes about 3 minutes and 8.5GB of memory to run
suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
  library(RcppHNSW)
  library(cowplot)
  library(rtracklayer)
})
addArchRThreads(2)

source("code_utils/projection_utils.R")
source("code_utils/fragmentsS4.R")

input_data_dir <- "/oak/stanford/groups/wjg/amyfchen/GEO_submission_June2021"
mapping_output <- "fig2_CD4_Tcells/outputs/hematopoiesis_projection.tsv"
dir.create(dirname(mapping_output), showWarnings = FALSE)

#Adapted & simplified from https://github.com/GreenleafLab/MPAL-Single-Cell-2019
projectLSI <- function(mat, LSI) {
  if(LSI$binarize) {
    message(paste0("Binarizing matrix..."))
    mat@x[mat@x > 0] <- 1
  }
  colSm <- Matrix::colSums(mat)
  
  message("Computing Term Frequency IDF...")
  freqs <- t(t(mat)/colSm)
  idf   <- as(log(1 + ncol(mat) / LSI$rowSm), "sparseVector")
  tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs
  
  matSVD <- as.matrix(t(tfidf) %*% LSI$svd$u)
  colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))
  return(matSVD)
}




proj <- readRDS("fig2_CD4_Tcells/data/ArchR_HTOsinglets.rds")

fragment_files <- c(
  "lane1"=file.path(input_data_dir, "CD4_lane1/lane1_atac_fragments.tsv.gz"),
  "lane2"=file.path(input_data_dir, "CD4_lane2/lane2_atac_fragments.tsv.gz")
)

# Takes 1.5 minutes to load fragments + insertions; 
# Total memory usage for the process afterwards is ~7GB
fragments <- map(fragment_files, ~fread(.x, col.names = c("chr", "start", "end", "cell_id", "duplicates"))) %>%
  rbindlist(idcol="Sample") %>%
  mutate(cell_id = str_c(Sample, "#", cell_id)) %>%
  filter(cell_id %in% proj$cellNames)

insertions <- Insertions(
  makeGRangesFromDataFrame(fragments),
  fragments$cell_id,
  cell_levels=proj$cellNames
)
rm(fragments)

# Downloaded from: https://jeffgranja.s3.amazonaws.com/MPAL-10x/Supplementary_Data/Healthy-Data/scATAC-Healthy-Hematopoiesis-191120.rds
# As linked from https://github.com/GreenleafLab/MPAL-Single-Cell-2019
# Takes about 1.5 minutes to download + parse for me
mpal <- readRDS(url("https://jeffgranja.s3.amazonaws.com/MPAL-10x/Supplementary_Data/Healthy-Data/scATAC-Healthy-Hematopoiesis-191120.rds"))


cluster_colors_table <- metadata(mpal)$colorMap$Clusters %>%
  tibble(
    Clusters = names(.),
    color = .
  ) %>%
  group_by(Clusters) %>% slice_head(n=1) %>% ungroup() %>% #Handle Cluster25 listed twice
  inner_join(as_tibble(colData(mpal)), by="Clusters") %>%
  select(Clusters, BioClassification, color) %>%
  distinct()

cluster_colors <- select(cluster_colors_table, BioClassification, color) %>%
  pull(color, name=BioClassification)

ref_peaks <- str_split_fixed(mpal@metadata$variablePeaks, "_", 3) %>%
  set_colnames(c("chr", "start", "end")) %>%
  as_tibble() %>%
  makeGRangesFromDataFrame()

if (FALSE) {
  # I've included the liftover coordinates as a data file, but
  # this is the code for how the liftover was generated
  as_tibble(ref_peaks) %>%
    select(seqnames, start, end) %>%
    mutate(id=seq_len(nrow(.))) %>%
    write_tsv("04_data/TM_030321/hematopoiesis_peaks_hg19.bed", col_names=FALSE)
  
  system("liftOver 04_data/TM_030321/hematopoiesis_peaks_hg19.bed 01_raw_data/hg19ToHg38.over.chain.gz 04_data/TM_030321/hematopoiesis_peaks_hg38.bed 04_data/TM_030321/hematopoiesis_peaks_hg19_unmapped.bed -bedPlus=3")
}

ref_peaks_hg38 <- read_tsv("fig2_CD4_Tcells/data/hematopoiesis_peaks_hg38.bed", 
                           col_names = c("chr", "start", "end", "id"),
                           col_types="ciii") %>%
  filter(str_detect(chr, "chr[0-9XY]+$")) %>%
  full_join(tibble(id=seq_along(ref_peaks)), by="id") %>%
  mutate(
    chr=replace_na(chr, "chr1"),
    start=replace_na(start, 0),
    end=replace_na(end, 500)
  ) %>%
  arrange(id) %>%
  makeGRangesFromDataFrame()

peak_mat <- peakMatrix(insertions, ref_peaks_hg38)
lsi_coords <- projectLSI(peak_mat, mpal@metadata$LSI)

knn <- get_knn(mpal@metadata$matSVD, lsi_coords, k=10, ef=200)

knn_data <- tibble(
  cell_id = rownames(lsi_coords),
  UMAP1 = mpal$UMAP1[knn$idx[,1]],
  UMAP2 = mpal$UMAP2[knn$idx[,2]],
  cell_type = transfer_labels(knn, mpal$BioClassification),
)

select(knn_data, cell_id, UMAP1, UMAP2, cell_type) %>% 
  write_tsv(mapping_output)



