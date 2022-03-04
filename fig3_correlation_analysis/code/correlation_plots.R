# Running time 1.5 minutes, <1GB RAM
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(ggrastr)
})

input_data_path <- "geo_download"
output_dir <- "fig3_correlation_analysis/outputs"

options(ggrastr.default.dpi=600)


tf_ensg_ids <- c(
  "RORgT" = "ENSG00000143365",
  "Tbet" = "ENSG00000073861",
  "FOXP3" = "ENSG00000049768",
  "GATA3" = "ENSG00000107485",
  "Helios" = "ENSG00000030419"
)

peak_corr <- read_csv(file.path(input_data_path, "GSM5396330_corr_peak_adt.csv.gz"))
rna_corr <- read_csv(file.path(input_data_path, "GSM5396330_corr_gene_adt.csv.gz"))

gene_labels <- list(
  "RORgT" = c("RORA", "ABCB1", "CCR6", "MAP3K4", "RBMS1", "KIF5C", "PTPN13", "RUNX2", "RORC"),
  "GATA3" = c("GATA3", "NIBAN1", "SNED1", "RUNX3", "NR3C1", "IL4R"),
  "Helios" = c("IKZF2", "RTKN2", "AC093865.1", "TTN", "SGMS1", "STAM", "TOX", "IL12RB2", "CCDC141", "SESN3", "HPGD", "CTLA4", "FOXP3"),
  "FOXP3" = c("IKZF2", "RTKN2", "NBAN1", "HPGD", "AC093865.1", "IL12RB2", "SGMS1", "CTLA4", "ACTN4", "TTN", "LAYN", "VAV3", "FOXP3", "IL2RA", "STAM"),
  "Tbet" = c("IFNG-AS1", "CCL5", "LRIG1", "INPP5D", "PARP8", "LUZP1", "EDARADD", "TPRG1", "TBX21")
)

gene_adt_corr_plots <- list()
for (tf in names(tf_ensg_ids)) {
  gene_adt_corr_plots[[tf]] <- rna_corr %>%
    filter(adt == tf) %>%
    mutate(rank = row_number(desc(spearman))) %>% {
      ggplot(., aes(rank, spearman)) +
        ggrastr::geom_point_rast(data=filter(., p.ttest < 0.05 | row_number(spearman) %% 50 == 0)) +
        geom_point(data=filter(., gene_id == tf_ensg_ids[!!tf]), color="firebrick") +
        geom_hline(yintercept = c(-1,1) * min(abs(.$spearman[.$p.ttest.adj < 0.05])), linetype="dashed") +
        ggrepel::geom_label_repel(aes(label=ifelse(gene_name %in% gene_labels[[!!tf]] , gene_name, NA)), max.overlaps=200) +
        ggrepel::geom_label_repel(aes(label=ifelse(spearman < 0 & p.ttest.adj < 0.05, gene_name, NA)), max.overlaps=20) +
        cowplot::theme_minimal_hgrid() +
        labs(subtitle = sprintf("%s genes up, %s genes down (p.adj < 0.05)", 
                                sum(.$spearman > 0 & .$p.ttest.adj < 0.05),
                                sum(.$spearman < 0 & .$p.ttest.adj < 0.05)),
             x="Gene Rank", y="Spearman correlation", title=sprintf("%s ADT-Gene correlation", tf))
    }
}

pdf(file.path(output_dir, "gene_adt_corr_plots.pdf"), useDingbats = FALSE)
gene_adt_corr_plots
dev.off()


peak_adt_corr_plots <- list()
for (tf in names(tf_ensg_ids)) {
  peak_adt_corr_plots[[tf]] <- peak_corr %>%
    filter(adt==tf) %>%
    mutate(rank = row_number(desc(spearman))) %>% {
      ggplot(., aes(rank, spearman, color=has_motif)) +
        ggrastr::geom_point_rast(alpha=0.1) +
        geom_hline(yintercept = c(-1,1) * min(abs(.$spearman[.$p.ttest.adj < 0.05])), linetype="dashed") +
        scale_color_manual(values=c(`TRUE`="red", `FALSE`="grey")) +
        scale_x_continuous(labels=NULL) +
        guides(alpha=FALSE) +
        cowplot::theme_minimal_hgrid() +
        labs(x="Peak Rank", y="Spearman correlation", title=sprintf("%s ADT-Peak Correlation", tf),
             subtitle=sprintf("%s peaks up, %s peaks down (p-adj < 0.05)\nMotif enrichment in up=%.2f (p=%.2e)", 
                              sum(.$spearman > 0 & .$p.ttest.adj < 0.05),
                              sum(.$spearman < 0 & .$p.ttest.adj < 0.05), 
                              mean(.$has_motif[.$spearman>0 & .$p.ttest.adj < 0.05])/mean(.$has_motif),
                              phyper(sum(.$has_motif & .$spearman > 0 & .$p.ttest.adj < 0.05)-1, 
                                     sum(.$has_motif),
                                     sum(!.$has_motif),
                                     sum(.$spearman > 0 & .$p.ttest.adj < 0.05),
                                     lower.tail=FALSE))
             )
    }
}

pdf(file.path(output_dir, "peak_adt_corr_plots.pdf"), useDingbats = FALSE)
peak_adt_corr_plots
dev.off()

