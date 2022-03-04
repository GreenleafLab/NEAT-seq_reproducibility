suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(dtplyr)
  library(GenomicRanges)
  library(patchwork)
  library(Matrix)
})
source("code_utils/fragmentsS4.R")
source("code_utils/plotting_utils.R")
source("fig1_species_mixing/code/species_mixing_utils.R")

input_data_path <- "geo_download"

cell_calls <- fread("fig1_species_mixing/outputs/cell_calls.tsv")

human_cells <- cell_calls %>% filter(consensus_species == "human") %>%
  pull(cell_id)

#### RNA COMPARISON ####
# Looks like the chimeric alignment is about as good as using the just human-aligned data, so 
# we'll do the chimeric alignment since it's already included in the data release
#rna <- read_10x_gene_h5("/oak/stanford/groups/wjg/bparks/TF_ADTs_Amy/01_raw_data/Dec2020_multiome_hg38/raw_feature_bc_matrix.h5") %>%
rna <- read_10x_gene_h5(file.path(input_data_path, "GSM5396329_barnyard_raw_feature_bc_matrix.h5")) %>%
  .[,human_cells] %>%
  rowSums() %>%
  {1e6 *. / sum(.)}


# From ENCSR000CPS K562 PolyA RNA nuclear fraction
rna_rep1 <- fread("https://www.encodeproject.org/files/ENCFF501IXI/@@download/ENCFF501IXI.tsv") %>%
  mutate(gene_id = str_remove(gene_id, "\\.[0-9]+$")) %>%
  as_tibble()

rna_rep2 <- fread("https://www.encodeproject.org/files/ENCFF142LZQ/@@download/ENCFF142LZQ.tsv") %>%
  mutate(gene_id = str_remove(gene_id, "\\.[0-9]+$")) %>%
  as_tibble() 

# From ENCSR530NHO K562 PolyA RNA nuclear fraction (older experiment w/o quantification spikeins)
rna_rep3 <- fread("https://www.encodeproject.org/files/ENCFF631TDY/@@download/ENCFF631TDY.tsv") %>%
  mutate(gene_id = str_remove(gene_id, "\\.[0-9]+$")) %>%
  as_tibble()

rna_rep4 <- fread("https://www.encodeproject.org/files/ENCFF826OCU/@@download/ENCFF826OCU.tsv") %>%
  mutate(gene_id = str_remove(gene_id, "\\.[0-9]+$")) %>%
  as_tibble()

shared_genes <- intersect(names(rna), rna_rep1$gene_id) %>%
  intersect(rna_rep3$gene_id)

rna_tpm <- tibble(
  "NEATseq" = rna[shared_genes],
  "K562_bulk_rep1" = rna_rep1$FPKM[match(shared_genes, rna_rep1$gene_id)],
  "K562_bulk_rep2" = rna_rep2$FPKM[match(shared_genes, rna_rep2$gene_id)],
  "K562_bulk_rep3" = rna_rep3$FPKM[match(shared_genes, rna_rep3$gene_id)],
  "K562_bulk_rep4" = rna_rep4$FPKM[match(shared_genes, rna_rep3$gene_id)],
  
) %>%
  filter(NEATseq != 0, K562_bulk_rep1 != 0, K562_bulk_rep2 != 0) %>%
  mutate(across(everything(), ~ log10(1 + 1e6 * .x / sum(.x))))


rna_GM12878 <- read_10x_gene_h5(file.path(input_data_path, "GSM5916384_filtered_feature_bc_matrix.h5")) %>%
  rowSums() %>%
  {1e6 *. / sum(.)}
gm_rna_rep1 <- fread("https://www.encodeproject.org/files/ENCFF387YXX/@@download/ENCFF387YXX.tsv") %>%
  mutate(gene_id = str_remove(gene_id, "\\.[0-9]+$")) %>%
  as_tibble()
gm_rna_rep2 <- fread("https://www.encodeproject.org/files/ENCFF089PEL/@@download/ENCFF089PEL.tsv") %>%
  mutate(gene_id = str_remove(gene_id, "\\.[0-9]+$")) %>%
  as_tibble()

shared_genes_gm <- intersect(names(rna_GM12878), gm_rna_rep1$gene_id)
rna_tpm_gm <- tibble(
  "10x_Multiome_GM" = rna_GM12878[shared_genes_gm],
  "GM_bulk_rep1" = gm_rna_rep1$FPKM[match(shared_genes_gm, gm_rna_rep1$gene_id)],
  "GM_bulk_rep2" = gm_rna_rep2$FPKM[match(shared_genes_gm, gm_rna_rep2$gene_id)]
) %>%
  filter(`10x_Multiome_GM` != 0, GM_bulk_rep1 != 0, GM_bulk_rep2 != 0) %>%
  mutate(across(everything(), ~ log10(1 + 1e6 * .x / sum(.x))))



p1 <- ggplot(rna_tpm, aes(K562_bulk_rep1, K562_bulk_rep3)) +
  stat_scatter_density() + scale_color_viridis_c() +
  ggpubr::stat_cor(aes(label=..r.label..)) 
  
p2 <- ggplot(rna_tpm, aes(NEATseq, K562_bulk_rep1)) +
  stat_scatter_density() + scale_color_viridis_c() +
  ggpubr::stat_cor(aes(label=..r.label..)) 

p3 <- ggplot(rna_tpm, aes(NEATseq, K562_bulk_rep3)) +
  stat_scatter_density() + scale_color_viridis_c() +
  ggpubr::stat_cor(aes(label=..r.label..)) 

p4 <- ggplot(rna_tpm_gm, aes(`10x_Multiome_GM`, GM_bulk_rep1)) +
  stat_scatter_density() + scale_color_viridis_c() +
  ggpubr::stat_cor(aes(label=..r.label..)) 




#### ATAC COMPARISON ####
insertion_from_bed <- function(path, valid_cells, blacklist_path="https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz") {
  fragments_table <- fread(path, col.names=c("chr", "start", "end", "barcode", "duplicates"),
                         stringsAsFactors = TRUE)
  fragments_granges <- makeGRangesFromDataFrame(fragments_table)
  blacklist_granges <- fread(blacklist_path,
                             sep="\t", col.names=c("chr", "start", "end", "reason")) %>%
    makeGRangesFromDataFrame()
  blacklist_overlaps <- countOverlaps(fragments_granges, blacklist_granges)
  
  cell_insertions <- fragments_table %>% lazy_dt() %>%
    filter(blacklist_overlaps == 0, barcode %in% valid_cells) %>%
    collect() %>% {
      Insertions(makeGRangesFromDataFrame(.), cell_ids = .$barcode, cell_levels = valid_cells)
    }
  return(cell_insertions)
}

cell_insertions <- insertion_from_bed(
  file.path(input_data_path, "GSM5396328_barnyard_atac_fragments.tsv.gz"),
  blacklist_path = "fig1_species_mixing/data/combined_blacklist.bed",
  valid_cells = human_cells
)

fix_chr_names_hg38 <- function(l) {
  l <- l[str_detect(names(l), "hg38.")]
  names(l) <- str_remove(names(l), "hg38.")
  return(l)
}

cell_insertions@iranges <- fix_chr_names_hg38(cell_insertions@iranges)
cell_insertions@cell_ids <- fix_chr_names_hg38(cell_insertions@cell_ids)

K562_rep1_insertions <- insertion_from_bed(
  "Supplementary_figures/outputs/K562_rep1/possorted.bed.gz",
  valid_cells= "K562_rep1"
)

# Take replicated peaks from the combination of rep1 + rep2
encode_K562_peaks <- fread("https://www.encodeproject.org/files/ENCFF558BLC/@@download/ENCFF558BLC.bed.gz",
                           col.names = c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")) %>%
  collect()

atac_counts <- tibble(
  NEATseq = peakMatrix(cell_insertions, GRanges(encode_K562_peaks)) %>% 
    rowSums(),
  bulk_rep1 =  peakMatrix(K562_rep1_insertions, GRanges(encode_K562_peaks)) %>%
    rowSums()
)

pdf("Supplementary_figures/outputs/correlation_plots.pdf", width=5, height=5, useDingbats = FALSE)

((p2 + labs(x="NEATseq K562")) + p4) * guides(color="none") * theme_classic() * coord_fixed() +
  patchwork::plot_annotation(title="Bulk vs. NEAT-seq RNA-seq", subtitle="K562=ENCSR000CPS, GM=ENCSR000CPO\nlog10(1 + TMP or FPKM)")

atac_counts %>%
  filter(NEATseq > 0, bulk_rep1 > 0) %>%
  ggplot(aes(log10(1+NEATseq), log10(1+bulk_rep1))) +
  stat_scatter_density() +
  scale_color_viridis_c() +
  ggpubr::stat_cor(aes(label=..r.label..)) +
  labs(x="log10(1 + NEAT-seq counts)", y="log10(1 + bulk K562 counts)", fill="Peaks in bin",
       title="ATAC-seq signal correlation",
       subtitle="Peaks from ENCFF558BLC, bulk data from ENCFF512VEZ") +
  coord_fixed() +
  theme_classic()

dev.off()