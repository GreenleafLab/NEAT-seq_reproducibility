# Takes 1 minute to run, uses 8GB RAM
suppressPackageStartupMessages({
  library(ArchR)
  
  library(tidyverse)
  library(dtplyr)
  library(ggrepel)
  library(glue)
  
  library(Matrix)
  library(patchwork)
})

source("fig3_correlation_analysis/code/correlation_utils.R")
Rcpp::sourceCpp("fig3_correlation_analysis/code/rcpp_utils.cc")


input_data_path <- "geo_download"

output_path <- "fig3_correlation_analysis/outputs/heatmap_plot.pdf"

genes_to_highlight <- list(
  RORgT = c("CCSAP", "CCL20", "CMTM6", "GLB1", "USP46", "SEPTIN11", "ARHGAP10", "ELOVL4", "SASH1", "AL139393.3", "MAP3K4", "CCR6", "EEPD1", "LRP12", "VIM", "TSPAN15", "VCL", "ARL3", "ADAM12", "MCAM", "ELK3", "FRY", "WDFY2", "TMED8", "EHD4", "MPRIP", "NR1D1", "HLF", "AC010754.1", "TTC39C-AS1", "EPG5", "KDSR", "IL4I1", "B4GALT5", "YWHAH"),
  GATA3 = c("SLC25A33", "NCAPH", "CCR4", "CCR8", "GPR146", "GNAQ", "STAM", "IATPR", "SORL1", "ANKRD13A", "FRY", "NDFIP2", "IL4R", "MORC4", "TSEN54"),
  Tbet = c("GFI1", "ARHGAP26")
)

proj <- readRDS("fig2_CD4_Tcells/data/ArchR_HTOsinglets_CD4only_25XADT.rds")
proj_all <- readRDS("fig2_CD4_Tcells/data/ArchR_HTOsinglets_CD4only.rds")

gene_names <- read_tsv(
    file.path(input_data_path, "GSM5396333_lane1_features.tsv.gz"),
    col_names=c("gene_id", "gene_name", "feature_type", "chr", "start", "end"),
    col_types="ccccii"
  ) %>%
  filter(feature_type == "Gene Expression") %>%
  pull(gene_name, name=gene_id)

rna_counts <- readRDS(file.path(input_data_path, "GSM5396333_CD4_RNA_counts.rds"))
peak_counts <- readRDS(file.path(input_data_path, "GSM5396336_CD4_Peak_matrix.rds"))

rna <- rna_counts %>%
  t() %>% {log10(10000 * . / rowMeans(.) + 1)} %>% as("dgCMatrix") %>%
  .[,colSums(.) > 0]
peaks <- peak_counts / proj_all$ReadsInTSS

adt <- bind_rows(
    lane1 = read_csv(file.path(input_data_path, "GSM5396330_ADT_counts_lane1.csv.gz")),
    lane2 = read_csv(file.path(input_data_path, "GSM5396334_ADT_counts_lane2.csv.gz")),
    .id = "lane"
) %>%
  mutate(cell_id = str_c(lane, "#", cell, "-1")) %>%
  mutate(across(c("FOXP3", "Tbet", "GATA3", "Helios", "RORgT"), ~ log10(250*.x/pmax(NPC1,NPC2) + 1))) %>%
  select(cell_id, c("FOXP3", "Tbet", "GATA3", "Helios", "RORgT")) %>%
  column_to_rownames("cell_id") %>%
  as.matrix() %>%
  .[proj$cellNames,]


set.seed(125134)
sample_order <- sample(seq_along(proj$cellNames))

pseudobulk_mat <- knn_pseudobulk(getReducedDims(proj), starter_cells = sample_order, k=100, max_overlap = 0.3) %>%
    .[1:500,]

clamp <- function(x, lo, hi) {
  pmin(pmax(x, lo), hi)
}

#' Convert matrix to a 3-column data frame holding elements of (row, column, value) 
matrix_to_data_frame <- function(matrix) {
  # Output colnames as integers if they are not already given in the matrix
  mat <- matrix
  colnames(mat) <- seq(ncol(mat))
  rownames(mat) <- seq(nrow(mat))
  mat %>% as.data.frame %>% tibble::rownames_to_column("row") %>%
    pivot_longer(-row, names_to="column", values_to="value") %>%
    mutate(row = as.integer(row), column = as.integer(column))
}



top_links <- read_csv(file.path(input_data_path, "GSM5396330_putative_targets.csv.gz"))

heatmaps <- list()
shared_plot <- NULL
for (tf in c("RORgT", "GATA3", "Tbet")) {
    
  links <- top_links %>% as_tibble() %>% filter(adt == !!tf)
  
  
  
  adt_vec <- as.matrix(pseudobulk_mat %*% adt[,tf])
  adt_order <- order(adt_vec)
  adt_norm <- scale(adt_vec)
  
  rna_mat <- as.matrix(pseudobulk_mat %*% rna[rownames(adt),links$gene_id])
  rna_mat_norm <- scale(rna_mat)
  
  peak_mat <- as.matrix(pseudobulk_mat %*% peaks[rownames(adt), links$peak_id])
  peak_mat_norm <- scale(peak_mat)
  
  
  
  link_order <- seriation::seriate(rbind(rna_mat_norm,peak_mat_norm), margin=2) %>%
    seriation::get_order()
  
  adt_frame <- matrix_to_data_frame(t(adt_norm[adt_order,,drop=FALSE])) %>%
    mutate(
      value = clamp(value, -2, 2)
    )
  
  rna_frame <- matrix_to_data_frame(t(rna_mat_norm[adt_order,link_order])) %>%
    mutate(
      value = clamp(value, -2, 2)
    )
  
  peak_frame <- matrix_to_data_frame(t(peak_mat_norm[adt_order,link_order])) %>%
    mutate(
      value = clamp(value, -2, 2)
    )
  
  adt_plot <- ggplot(adt_frame, aes(column, -row)) +
    geom_raster(aes(fill=value)) + 
    scale_x_discrete(breaks=NULL) +
    scale_y_continuous(breaks=NULL) +
    scale_fill_gradientn(colors=viridisLite::magma(50)) +
    theme_void() +
    theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")) +
    labs(fill="ADT")
  
  gene_axis_labels <- tibble(
    gene_id = links$gene_id,
    gene_name = gene_names[gene_id],
    highlight = gene_names[gene_id] %in% genes_to_highlight[[tf]],
    y = link_order
  ) %>%
    filter(highlight)
  
  if (nrow(gene_axis_labels) > 0) {
    repel_labels <- ggplot(gene_axis_labels, aes(0, -y)) +
      ggrepel::geom_text_repel(
        aes(label=gene_name),
        force        = 0.5,
        nudge_x      = 0.1,
        direction    = "y",
        hjust        = 0,
        segment.size = 0.2,
        min.segment.length = 0
      ) + 
      scale_y_continuous(breaks=NULL, limits=c(min(-link_order), max(-link_order)), expand = c(0,0))+
      scale_x_continuous(breaks=NULL, limits=c(0,1), expand=c(0,0)) +
      theme_void() +
      theme(plot.margin = unit(c(5.5, 5.5, 5.5, 0), "pt"))
  } else {
      repel_labels <- patchwork::plot_spacer()
  }
  
  rna_plot <- ggplot(rna_frame, aes(column, -row)) +
    geom_raster(aes(fill=value)) + 
    scale_x_discrete(breaks=NULL) +
    scale_y_continuous(breaks=-gene_axis_labels$y, 
                       labels=NULL, #gene_axis_labels$gene_name, 
                       position = "right", 
                       expand = c(0,0)) +
    scale_fill_gradientn(colors=ArchR::paletteContinuous("blueYellow")) +
    theme_void() +
    theme(plot.margin = unit(c(5.5, 0, 5.5, 5.5), "pt"),
          #axis.text.y.right = element_text(),
          axis.ticks.y = element_line(),
          axis.ticks.length.y = unit(3.5, "pt")) +
    labs(fill="RNA", x=NULL, y=NULL)
  
  atac_plot <- ggplot(peak_frame, aes(column, -row)) +
    geom_raster(aes(fill=value)) + 
    scale_x_discrete(breaks=NULL) +
    scale_y_continuous(breaks=NULL, expand = c(0,0)) +
    scale_fill_gradientn(colors=ArchR::paletteContinuous("solarExtra")) +
    theme_void() +
    theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
          axis.text.y.right = element_text(),
          axis.ticks.y = element_line(),
          axis.ticks.length.y = unit(3.5, "pt")) +
    labs(fill="ATAC", x=NULL, y=NULL)
  
  heatmaps[[tf]] <- adt_plot + adt_plot + patchwork::plot_spacer() + 
    atac_plot + rna_plot + repel_labels +
    patchwork::plot_layout(guides="collect", heights=c(0.1,1), widths=c(1,1, 0.6), nrow=2) +
    patchwork::plot_annotation(title = tf)
  
  tf_label <- ggplot(NULL, aes(0,0,label=!!tf)) + 
    geom_text(fontface="bold", size=4) + 
    theme_void()
  
  if (is.null(shared_plot)) {
    shared_plot <- adt_plot + adt_plot + tf_label + 
      atac_plot + rna_plot + repel_labels
  } else {
    shared_plot <- shared_plot + adt_plot + adt_plot + tf_label + 
      atac_plot + rna_plot + repel_labels
  }
}

tfs <- c("RORgT", "GATA3", "Tbet")
heatmap_heights <- top_links %>% as_tibble() %>%
  filter(adt %in% tfs) %>%
  group_by(adt) %>%
  tally() %>%
  pull(n, name=adt) %>%
  {.[tfs]} %>%
  {as.vector(rbind(20, .))}

pdf(output_path, useDingbats = FALSE, height=21)
shared_plot + patchwork::plot_layout(
  guides="collect", 
  widths=c(1,1,0.6), 
  heights=heatmap_heights,
  nrow=length(tfs)*2)
dev.off()
