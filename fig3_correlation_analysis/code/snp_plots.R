# Takes about 30 seconds to run, and 1GB of memory
suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
  library(patchwork)
})

source("fig3_correlation_analysis/code/correlation_utils.R")

input_data_path <- "geo_download"
output_path = "fig3_correlation_analysis/outputs/"

proj <- readRDS("fig2_CD4_Tcells/data/ArchR_HTOsinglets_CD4only_25XADT.rds")
proj_all <- readRDS("fig2_CD4_Tcells/data/ArchR_HTOsinglets_CD4only.rds")

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
  .[proj_all$cellNames,]

impute_mat <- get_impute_matrix(getReducedDims(proj_all))

snp_fragments <- bind_rows(
    lane1 = read_tsv("fig3_correlation_analysis/outputs/lane1_fragments.tsv",
                    col_names=c("chr", "start", "end", "barcode", "snp", "qual")),
    lane2 = read_tsv("fig3_correlation_analysis/outputs/lane2_fragments.tsv",
                    col_names=c("chr", "start", "end", "barcode", "snp", "qual")),
    .id="lane"
  ) %>%
  as_tibble() %>%
  mutate(cell_id = str_c(lane, "#", barcode)) %>%
  group_by(start, end, cell_id, snp) %>%
  summarize(qual=max(qual), duplicates=n(), .groups="drop")
  

nucleotide_colors <- c('A'='#109648','T'='#D62839','U'='#D62839', 'C'='#255C99', 'G'='#F7B32B', 'N'='#000000')

z_transform <- function(x) {(x-mean(x))/sd(x)}
quant_filter <- function(x, q, below=FALSE) {names(x)[below != (x >= quantile(x, q))]}

gata3 <- adt[,"GATA3"]
gata3[proj$cellNames] <- z_transform(gata3[proj$cellNames])
gata3[!(names(gata3) %in% proj$cellNames)] <- z_transform(gata3[!(names(gata3) %in% proj$cellNames)])
gata3_impute <- impute_data(gata3, impute_mat)
names(gata3_impute) <- proj_all$cellNames

bottom_50 <- quant_filter(gata3_impute, .5, below=TRUE)
top_10 <- quant_filter(gata3_impute, .9, below=FALSE)
gata3_quantiles <- c(rep_along(bottom_50, "<50%"), rep_along(top_10, ">90%")) %>%
  set_names(c(bottom_50, top_10))

  

plot_snp <- function(snp_fragments, cell_groups, comparison_base) {
  norm_factors <- tapply(
    proj_all$ReadsInTSS[match(names(cell_groups), proj_all$cellNames)], 
    cell_groups, sum)/1e6
  
  data <- snp_fragments %>%
    filter(cell_id %in% names(cell_groups)) %>%
    mutate(
      adt_group = cell_groups[cell_id]
    ) %>% 
    group_by(adt_group) %>%
    summarize(
      a = sum(snp=="A"),
      g = sum(snp=="G"),
      a.g=a/g, 
      A_norm = a / norm_factors[as.character(adt_group)][1],
      G_norm = g / norm_factors[as.character(adt_group)][1],
      .groups="drop"
    ) %>%
    as.data.frame()
  rownames(data) <- data$adt_group
  p1 <- data %>% 
    pivot_longer(ends_with("_norm"), names_to="base", values_to="insertions") %>%
    mutate(
      base=str_remove(base, "_norm"),
      base=c("A"="T", "G"="C", "T"="A", "C"="G")[base]
    ) %>%
    ggplot(aes(adt_group, insertions)) +
    geom_col(aes(fill=base)) +
    scale_fill_manual(values=nucleotide_colors) +
    cowplot::theme_minimal_hgrid() +
    scale_y_continuous(expand=c(0,0))+
    labs(y="Normalized insertions", x=NULL)
    
  base_rate <- data %>%
    filter(adt_group == comparison_base) %>%
    {.$g / (.$a + .$g)}
  p_for_ggpubr <- data %>%
    filter(adt_group != comparison_base) %>%
    transmute(
      group1 = adt_group,
      group2 = comparison_base,
      p_val = pbinom(g, a+g, base_rate),
      label = sprintf("p = %.2e", p_val),
      y.position=a.g + 1
    ) %>%
    filter(p_val < 0.05)
  
  p2 <- data %>%
    ggplot(aes(adt_group, a.g)) +
    geom_col() +
    scale_color_brewer(palette="Set1") +
    geom_label(aes(label=sprintf("%d T / %d C", a, g))) +
    labs(y="T:C ratio", x = "") +
    ggpubr::stat_pvalue_manual(p_for_ggpubr, label="label") +
    cowplot::theme_minimal_hgrid()
  
  p1 + p2 + patchwork::plot_annotation(title="TSEN54 Allele-specific Accessibility")
}



pdf("fig3_correlation_analysis/outputs/snp_plots.pdf", width=14, useDingbats=FALSE)

plot_snp(snp_fragments, gata3_quantiles, "<50%")

dev.off()