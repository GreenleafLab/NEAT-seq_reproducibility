# About 3.5 minutes to run, and 10GB of RAM
suppressPackageStartupMessages({
   library(ArchR)
  
   library(chromVARmotifs)
   library(TFBSTools)
   library(BSgenome.Hsapiens.UCSC.hg38)
   library(motifmatchr)
  
   library(IRanges)
   library(tidyverse)
   library(data.table)
   library(dtplyr)
   library(ggrepel)
  
   library(GenomicRanges)
})

input_data_path <- "geo_download"
output_path <- "fig3_correlation_analysis/outputs"

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

refGenes <- read_gtf("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/genes/hg38.refGene.gtf.gz", 
        c("gene_id", "transcript_id", "gene_name")) %>%
    keepStandardChromosomes(species="Homo_sapiens", pruning.mode = "coarse") %>%
    .[.$feature %in% c("transcript", "exon")]



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
    file.path(input_data_path, "GSM5396333_lane1/features.tsv.gz"),
    col_names=c("gene_id", "gene_name", "feature_type", "chr", "start", "end"),
    col_types="ccccii"
  ) %>%
  filter(feature_type == "Gene Expression") %>%
  pull(gene_name, name=gene_id)

gene_ids <- names(gene_names) %>% set_names(gene_names)


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
  
all_peak_coords <- peak_motif_matches %>% {unique(.$peak_id)} %>% as("GRanges")
all_peak_coords$group <- "All Peaks"
adt_peak_coords <- peak_motif_matches %>%
  filter(has_motif) %>% {
    x <- as(.$peak_id, "GRanges")
    x$group <- .$adt
    x
  }
peak_coords <- c(all_peak_coords, adt_peak_coords)


links <- read_csv(file.path(input_data_path, "GSM5396330_putative_targets.csv.gz"))

link_coords <- tibble(
  chr = as.vector(seqnames(gene_coords[links$gene_id])),
  start = start(resize(gene_coords[links$gene_id], 1, fix="start")),
  end = start(resize(as(links$peak_id, "GRanges"), 1, fix="center")),
  correlation = links$spearman.peak_gene,
  pval = -log10(links$p.ttest.adj.peak_gene),
  gene_id = links$gene_id,
  adt = links$adt
) %>% 
  mutate(
    tmp = pmin(start,end),
    end = pmax(start,end),
    start = tmp
  ) %>%
  as("GRanges")

# 3.7GB ram in use before running this...
gc()

lane1_fragments <- fread(file.path(input_data_path, "GSM5396332_lane1_atac_fragments.tsv.gz"), 
                       col.names=c("chr", "start", "end", "barcode", "duplicates"),
                       stringsAsFactors = TRUE) %>% lazy_dt() %>%
        mutate(cell_id = str_c("lane1#", barcode)) %>%
        filter(cell_id %in% proj$cellNames) %>%
        collect() %>%
        makeGRangesFromDataFrame(keep.extra.columns=TRUE)

lane2_fragments <- fread(file.path(input_data_path, "GSM5396336_lane2_atac_fragments.tsv.gz"), 
                       col.names=c("chr", "start", "end", "barcode", "duplicates"),
                       stringsAsFactors = TRUE) %>% lazy_dt() %>%
        mutate(cell_id = str_c("lane1#", barcode)) %>%
        filter(cell_id %in% proj_all$cellNames) %>%
        collect() %>%
        makeGRangesFromDataFrame(keep.extra.columns=TRUE)

insertions <- Insertions(
    c(lane1_fragments, lane2_fragments), 
    cell_ids = c(lane1_fragments$cell_id, lane2_fragments$cell_id),
    cell_levels = proj$cellNames
)

rm(lane1_fragments, lane2_fragments)
gc()

make_trackplot <- function(gene_name, tf, peak_id, min_peak_gene=0.01, zoom_region=FALSE, bin_width=200) {
  gene_id <- gene_ids[gene_name]
  
  
  if (zoom_region) {
    gene_start = start(resize(gene_coords[gene_id], 1, fix="start"))
    coords <- GRanges(
      seqnames=seqnames(gene_coords[gene_id]),
      ranges = IRanges(
        start = min(gene_start, start(as(peak_id, "GRanges"))),
        end = max(gene_start, end(as(peak_id, "GRanges")))
      )
    ) %>%
      resize(2*width(.), fix="center") %>%
      `strand<-`(value="*")
  } else {
    coords <- GRanges(
      seqnames=seqnames(gene_coords[gene_id]),
      ranges = IRanges(
        start = min(start(gene_coords[gene_id]), start(as(peak_id, "GRanges"))),
        end = max(end(gene_coords[gene_id]), end(as(peak_id, "GRanges")))
      )
    ) %>%
      resize(width(.) + 2*25000, fix="center") %>%
      `strand<-`(value="*")
  }
  
  bulk_group_levels = c("0-20%", "20-40%", "40-60%", "60-80%", "80-100%")
  bulk_groups <- as.integer((rank(adt[proj$cellNames, tf])-1)*5/length(proj$cellNames)) + 1
  bulk_groups <- bulk_group_levels[bulk_groups]
  bulk_palette <- c(RColorBrewer::brewer.pal(6, "BuPu")[-1])
  names(bulk_palette) <- bulk_group_levels
  
  group_norms <- tapply(proj$ReadsInTSS, bulk_groups, sum) / 1e6
  bulk_track <- trackplot_bulk(coords, insertions[proj$cellNames], bulk_groups, group_norms, bin_width = bin_width,
                               palette=bulk_palette, legend_label = sprintf("%s ADT quintile", tf))
  
  gene_track <- trackplot_gene(coords, refGenes, refGenes$transcript_id, 
                               labels= ifelse(refGenes$feature=="transcript", refGenes$gene_name, NA),
                               widths= ifelse(refGenes$feature=="transcript", 1, 3)) +
    labs(color="Gene strand")
  
  has_motif_label <- str_c("Peak w/\n", tf, " motif")
  peak_palette <- c("All Peaks"="darkgrey", "highlighted"="firebrick")
  peak_palette[has_motif_label] <- bulk_palette[length(bulk_palette)]
  display_peaks <- c(
    peak_coords[peak_coords$group %in% c("All Peaks", tf)],
    as(peak_id, "GRanges")
  )
  display_peaks$group[length(display_peaks)] <- "highlighted"
  display_peak_groups <- factor(
    ifelse(display_peaks$group == tf, has_motif_label, display_peaks$group),
    levels=c("highlighted", "All Peaks", has_motif_label)
  )
  peak_track_motif <- trackplot_peak(coords, display_peaks, display_peak_groups, peak_palette) +
    labs(color="Peak type")
  
  keeper_links <- link_coords$pval > -log10(min_peak_gene) & link_coords$gene_id == gene_id &
    link_coords$adt == tf
  loop_track <- trackplot_loop(coords, 
                               link_coords[keeper_links], link_coords[keeper_links]$correlation, 
                               viridisLite::plasma(20)) +
    labs(color="Correlation")
  
  
  trackplot <- draw_trackplot_grid(
    bulk_track, peak_track_motif, gene_track, loop_track,
    labels=c("Bulk", "Peaks", "Genes", sprintf("Correlation\nto %s", gene_name)),
    heights=c(1, 1.5, 1.5, 1),
    label_width = 0.1,
    label_style = list(fontface="bold", size=4)
    ) +
    patchwork::plot_annotation(title=sprintf("%s regulation by %s", gene_name, tf), 
                               subtitle=sprintf("%s:%d-%d", seqnames(coords), start(coords), end(coords)))
  
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
      ci_low_rna = rna_mean - sd(rna)/sqrt(n()),
      ci_high_rna = rna_mean + sd(rna)/sqrt(n()),
      peak_mean=mean(peak),
      ci_low_peak = peak_mean - sd(peak)/sqrt(n()),
      ci_high_peak = peak_mean + sd(peak)/sqrt(n())
    )
  
  rna_plot <- ggplot(quintile_data, aes(adt_mean, rna_mean, ymin=ci_low_rna, ymax=ci_high_rna, color=group)) +
    geom_pointrange() +
    scale_color_manual(values=bulk_palette) +
    theme_bw() + theme(panel.grid = element_blank()) +
    labs(x=sprintf("%s log2-normalized ADT (a.u.)", tf),
         y=sprintf("%s log2-normalized RNA", gene_name),
         color=sprintf("%s quintile", tf),
         title=sprintf("%s ADT vs. %s RNA", tf, gene_name)) 
  
  atac_plot <- ggplot(quintile_data, aes(adt_mean, peak_mean, ymin=ci_low_peak, ymax=ci_high_peak, color=group)) +
    geom_pointrange() +
    scale_color_manual(values=bulk_palette) +
    theme_bw() + theme(panel.grid = element_blank()) +
    labs(x=sprintf("%s log2-normalized ADT (a.u.)", tf),
         y=sprintf("%s\nlog2-normalized accessibility (a.u)", peak_id),
         color=sprintf("%s quintile", tf),
         title=sprintf("%s ADT vs. Peak accessibility", tf),
         subtitle=peak_id) 
  
  detail_plot <- rna_plot + atac_plot + patchwork::plot_layout(guides="collect")
  return(
    list(trackplot, detail_plot)
  )
}


dir.create(output_path, showWarnings=FALSE)

# Specific plots for paper
pdf(file.path(output_path, "trackplots.pdf"), useDingbats = FALSE, width=12, height=6)

make_trackplot("CCR6", "RORgT", "chr6:167054110-167054610", zoom_region=TRUE)
make_trackplot("CCR6", "RORgT", "chr6:167114436-167114936", zoom_region=TRUE, bin_width=20)

# Add a vertical line showing the SNP location in the TSEN54 trackplot
x <- make_trackplot("TSEN54", "GATA3", "chr17:75518153-75518653", zoom_region=TRUE, bin_width=20)
for (i in 9:16) {
  x[[1]][[i]] <- x[[1]][[i]] + geom_vline(xintercept = 75518526)
}
x

make_trackplot("CCR4","GATA3", "chr3:32998996-32999496", zoom_region=TRUE)


dev.off()


