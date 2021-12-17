# About 3.5 minutes to run, and 10GB of RAM
suppressPackageStartupMessages({
  library(ArchR)

  library(IRanges)
  library(tidyverse)
  library(data.table)
  library(dtplyr)
  library(ggrepel)
  
  library(GenomicRanges)
})

input_data_path <- "geo_downlaod"
output_path <- "revisions/outputs"

source("fig3_correlation_analysis/code/correlation_utils.R")
source("fig3_correlation_analysis/code/trackplot_utils.R")
source("code_utils/fragmentsS4.R")


tf_ensg_ids <- c(
  "RORgT" = "ENSG00000143365",
  "Tbet" = "ENSG00000073861",
  "FOXP3" = "ENSG00000049768",
  "GATA3" = "ENSG00000107485",
  "Helios" = "ENSG00000030419"
)

proj <- readRDS("fig2_CD4_Tcells/data/ArchR_HTOsinglets_CD4only_25XADT.rds")
proj_all <- readRDS("fig2_CD4_Tcells/data/ArchR_HTOsinglets_CD4only.rds")

rna_counts <- readRDS(file.path(input_data_path, "GSM5396333_CD4_RNA_counts.rds"))
peak_counts <- readRDS(file.path(input_data_path, "GSM5396336_CD4_Peak_matrix.rds"))

rna <- rna_counts %>%
  t() %>% {log10(10000 * . / rowMeans(.) + 1)} %>% as("dgCMatrix") %>%
  .[,colSums(.) > 0]
peaks <- peak_counts / proj_all$ReadsInTSS

rm(rna_counts, peak_counts)

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


gene_names <- read_tsv(
  file.path(input_data_path, "GSM5396333_lane1_features.tsv.gz"),
  col_names=c("gene_id", "gene_name", "feature_type", "chr", "start", "end"),
  col_types="ccccii"
) %>%
  filter(feature_type == "Gene Expression") %>%
  pull(gene_name, name=gene_id)

gene_ids <- names(gene_names) %>% set_names(gene_names)

bootstrap_ci <- function(x, q, iters=200) {
  sim <- map_dbl(seq_len(iters), ~ mean(sample(x, length(x), replace=TRUE)))
  quantile(sim, q)
}

make_plots <- function(gene_name, tf, peak_id) {
  gene_id <- gene_ids[gene_name]
  
  bulk_group_levels = c("0-20%", "20-40%", "40-60%", "60-80%", "80-100%")
  bulk_groups <- as.integer((rank(adt[proj$cellNames, tf])-1)*5/length(proj$cellNames)) + 1
  bulk_groups <- bulk_group_levels[bulk_groups]
  bulk_palette <- c(RColorBrewer::brewer.pal(6, "BuPu")[-1])
  names(bulk_palette) <- bulk_group_levels
  
  quintile_data <- tibble(
    group = bulk_groups,
    rna = rna[proj$cellNames, gene_id]/log10(2),
    peak = log2(1+median(proj$nFrags)*peaks[proj$cellNames, peak_id]),
    adt = adt[,tf]/log10(2)
  ) %>%
    group_by(group) %>%
    summarize(
      adt_mean = mean(adt),
      rna_mean = mean(rna),
      ci_low_rna = bootstrap_ci(rna, pnorm(1)), #rna_mean - sd(rna)/sqrt(n()),
      ci_high_rna = bootstrap_ci(rna, pnorm(-1)), #rna_mean + sd(rna)/sqrt(n()),
      peak_mean=mean(peak),
      ci_low_peak = bootstrap_ci(peak, pnorm(1)), #peak_mean - sd(peak)/sqrt(n()),
      ci_high_peak = bootstrap_ci(peak, pnorm(-1)),# peak_mean + sd(peak)/sqrt(n()),
      n=n()
    )
  sc_data <- tibble(
    rna = rna[proj$cellNames, gene_id]/log10(2),
    peak = log2(1+median(proj$nFrags)*peaks[proj$cellNames, peak_id]),
    adt = adt[,tf]/log10(2)
  )
  
  group_sizes <- pull(quintile_data, n, name=group)
  
  rna_plot <- ggplot(quintile_data, aes(adt_mean, rna_mean, ymin=ci_low_rna, ymax=ci_high_rna, color=group)) +
    geom_pointrange() +
    scale_color_manual(values=bulk_palette) +
    theme_bw() + theme(panel.grid = element_blank()) +
    labs(x=sprintf("%s log2-normalized ADT (a.u.)", tf),
         y=sprintf("%s log2-normalized RNA", gene_name),
         color=sprintf("%s quintile", tf),
         title=sprintf("%s ADT vs. %s RNA", tf, gene_name),
         subtitle=paste0("n=", paste0(group_sizes[bulk_group_levels], collapse=",")))
  
  atac_plot <- ggplot(quintile_data, aes(adt_mean, peak_mean, ymin=ci_low_peak, ymax=ci_high_peak, color=group)) +
    geom_pointrange() +
    scale_color_manual(values=bulk_palette) +
    theme_bw() + theme(panel.grid = element_blank()) +
    labs(x=sprintf("%s log2-normalized ADT (a.u.)", tf),
         y=sprintf("%s\nlog2-normalized accessibility (a.u)", peak_id),
         color=sprintf("%s quintile", tf),
         title=sprintf("%s ADT vs. Peak accessibility", tf),
         subtitle=peak_id) 
  
  sc_rna <- ggplot(sc_data, aes(adt, rna)) +
    stat_bin_hex(bins=100) + scale_fill_viridis_c(trans="log10") +
    stat_smooth(method = "loess", color="black") +
    labs(x=sprintf("%s log2-normalized ADT (a.u.)", tf),
         y=sprintf("%s log2-normalized RNA", gene_name),
         title=sprintf("%s ADT vs. %s RNA", tf, gene_name),
         subtitle=sprintf("Spearman rho=%.2f, p=%.2e", 
                          cor(sc_data$adt, sc_data$rna, method="spearman"),
                          cor.test(sc_data$adt, sc_data$rna, method="spearman")$p.value))
  
  sc_atac <- ggplot(sc_data, aes(adt, peak)) +
    stat_bin_hex(bins=100) + scale_fill_viridis_c(trans="log10") +
    stat_smooth(method = "loess", color="black") +
    labs(x=sprintf("%s log2-normalized ADT (a.u.)", tf),
         y=sprintf("%s\nlog2-normalized accessibility (a.u)", peak_id),
         color=sprintf("%s quintile", tf),
         title=sprintf("%s ADT vs. Peak accessibility", tf),
         subtitle=sprintf("%s\nSpearman rho=%.2f, p=%.2e", 
                          peak_id,
                          cor(sc_data$adt, sc_data$peak, method="spearman"),
                          cor.test(sc_data$adt, sc_data$peak, method="spearman")$p.value))
  
  detail_plot <- rna_plot + atac_plot + patchwork::plot_layout(guides="collect")
  return(
    (detail_plot / (sc_rna + sc_atac)) + 
      patchwork::plot_annotation(subtitle="Top: bootstrap 68% CI of mean [equivalent of +/- sem], n of ADT groups in RNA subtitle\nBottom: single cell plots with loess regression smoothing line")
  )
}


dir.create(output_path, showWarnings=FALSE)

# Specific plots for paper
pdf(file.path(output_path, "correlation_plots_ungrouped.pdf"), useDingbats = FALSE, width=12, height=12)

make_plots("CCR6", "RORgT", "chr6:167054110-167054610")
make_plots("CCR6", "RORgT", "chr6:167114436-167114936")

# Add a vertical line showing the SNP location in the TSEN54 trackplot
make_plots("TSEN54", "GATA3", "chr17:75518153-75518653")

make_trackplot("CCR4","GATA3", "chr3:32998996-32999496")


dev.off()


