# Running time on Sherlock cluster: 1 CPU, 4.5 min, max 9.5GB RAM
# Most RAM usage comes from the ATAC-qc analysis, which loads all fragments
# into memory
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(dtplyr)
  library(glue)
  library(cowplot)
  library(GenomicRanges)
  library(Matrix)
  library(hdf5r)
  library(ggrastr)
})
source("fig1_species_mixing/code/species_mixing_utils.R")
source("code_utils/fragmentsS4.R")
source("code_utils/plotting_utils.R")

species_colors_full <- c("mouse" = "#e41a1c","human" = "#377eb8", "mix"="#984ea3", "low yield"="#cccccc", "hash mix"="#4daf4a", "hash empty"="#ff7f00")
species_colors <- species_colors_full[c("mouse", "human")]

input_data_path <- "geo_download"

output_path <- "fig1_species_mixing/outputs/cell_calls.tsv"
plots_output <- "fig1_species_mixing/outputs/plots.pdf"
adt_counts_output <- "fig1_species_mixing/outputs/adt_counts.csv.gz"
dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)

# Outputs a table of cell species + doublet calls based on RNA, ATAC, and antibody info
# Output columns: 
# - cell_id, 
# - ATAC_pass, RNA_pass, ADT_pass (whether droplets pass min reads in modality) 
# - ATAC_species, RNA_species (species call for the modality)
# - ADT_hash_count (count of positive hash oligos for ADT)
# - consensus_pass (whether droplet is a pass in min reads for all modalities)
# - consensus_species (human, mouse, or mix. Mix if ATAC, RNA, or ADT calls a mix)
get_purity_data <- function(hg38_counts, mm10_counts, cell_ids, hg38_purity, mm10_purity) {
  tibble(
    cell_id = cell_ids,
    hg38 = hg38_counts,
    mm10 = mm10_counts,
    reads = hg38 + mm10,
    purity = pmax(hg38, mm10)/(hg38 + mm10),
    nearest_species = ifelse(hg38 > mm10, "human", "mouse"),
    species = case_when(
      hg38/reads > hg38_purity ~ "human",
      mm10/reads > mm10_purity ~ "mouse",
      TRUE ~ "mix"
    )
  )
}

#### RNA Calls ####
rna_mix <- read_10x_gene_h5(file.path(input_data_path, "GSM5396329_barnyard_raw_feature_bc_matrix.h5"))
rna_mix_purity <- get_purity_data(
  hg38_counts = Matrix::colSums(rna_mix[str_detect(rownames(rna_mix), "^ENSG"),]),
  mm10_counts = Matrix::colSums(rna_mix[str_detect(rownames(rna_mix), "^ENSMUSG"),]),
  cell_ids = colnames(rna_mix),
  hg38_purity = .7,
  mm10_purity = .95
) %>%
  filter(reads > 10)

rna_mix_genes <- tibble(
  hg38_genes = Matrix::colSums(rna_mix[str_detect(rownames(rna_mix), "^ENSG"),] > 0),
  mm10_genes = Matrix::colSums(rna_mix[str_detect(rownames(rna_mix), "^ENSMUSG"),] > 0),
  cell_id = colnames(rna_mix)
) %>%
  inner_join(rna_mix_purity, by="cell_id") %>%
  filter(species != "mix", reads > 7500) %>%
  mutate(
    umis = if_else(species == "mouse", mm10, hg38),
    genes = if_else(species == "mouse", mm10_genes, hg38_genes)
  )

rna_calls <- tibble(
  cell_id = rna_mix_purity$cell_id,
  RNA_pass = rna_mix_purity$reads > 7500,
  RNA_species = rna_mix_purity$species,
  hg38_genes = Matrix::colSums(rna_mix[str_detect(rownames(rna_mix), "^ENSG"),cell_id] > 0),
  mm10_genes = Matrix::colSums(rna_mix[str_detect(rownames(rna_mix), "^ENSMUSG"),cell_id] > 0),
  hg38_umis = rna_mix_purity$hg38,
  mm10_umis = rna_mix_purity$mm10
)



#### ATAC Calls ####
fragments_mix_raw <- fread(file.path(input_data_path, "GSM5396328_barnyard_atac_fragments.tsv.gz"), 
                       col.names=c("chr", "start", "end", "barcode", "duplicates"),
                       stringsAsFactors = TRUE)
fragments_granges_raw <- makeGRangesFromDataFrame(fragments_mix_raw)
blacklist_granges <- fread("fig1_species_mixing/data/combined_blacklist.bed",
                   sep="\t", col.names=c("chr", "start", "end", "reason")) %>%
  makeGRangesFromDataFrame()

blacklist_overlaps <- countOverlaps(fragments_granges_raw, blacklist_granges)

atac_genome_counts <- fragments_mix_raw %>% lazy_dt() %>%
  filter(blacklist_overlaps == 0) %>%
  group_by(barcode) %>%
  summarize(
    hg38_counts = sum(str_detect(chr, "^hg38\\.chr[0-9X]")),
    mm10_counts = sum(str_detect(chr, "^mm10\\.chr[0-9X]"))
  ) %>%
  collect()


multiome_atac_purity <- get_purity_data(
  hg38_counts = atac_genome_counts$hg38_counts,
  mm10_counts = atac_genome_counts$mm10_counts,
  cell_ids = atac_genome_counts$barcode,
  hg38_purity = 0.95,
  mm10_purity = 0.95
) %>%
  filter(reads > 10)



atac_calls <- tibble(
  cell_id = multiome_atac_purity$cell_id,
  ATAC_pass = multiome_atac_purity$reads > 10000,
  ATAC_species = multiome_atac_purity$species,
  hg38_fragments = multiome_atac_purity$hg38,
  mm10_fragments = multiome_atac_purity$mm10
)


#### HTO Calls ####
adt_clr <- read_csv("fig1_species_mixing/data/ADT_table.csv.gz") %>%
  mutate(cell = str_c(cell, "-1")) %>%
  column_to_rownames("cell") %>%
  select(GATA1, SOX2, OCT4, starts_with("NPC")) %>%
  as.matrix() %>% t() %>% clr_normalize(base=2) %>% t() %>%
  as_tibble(rownames = "cell")
hto_matrix_raw <- read_csv("fig1_species_mixing/data/ADT_table.csv.gz") %>%
  filter(reads > 400) %>%
  mutate(cell = str_c(cell, "-1")) %>%
  column_to_rownames("cell") %>%
  select(starts_with("NPC")) %>%
  rename_with(~str_replace(.x, "NPC", "HTO"), starts_with("NPC")) %>%
  as.matrix() %>%
  t()


hto_classifications <- clr_normalize(hto_matrix_raw, base=2) %>%
  HTO_demultiplex()
hto_calls <- tibble(
  cell_id = names(hto_classifications$hto_count),
  ADT_hash_count = hto_classifications$hto_count,
  ADT_hash_assignment = hto_classifications$hto_assignment,
  ADT_pass = TRUE
)




### ATAC QC METRICS ###
cell_insertions <- fragments_mix_raw %>% lazy_dt() %>%
  filter(blacklist_overlaps == 0, barcode %in% multiome_atac_purity$cell_id) %>%
  collect() %>% {
    Insertions(makeGRangesFromDataFrame(.), cell_ids = .$barcode, cell_levels = multiome_atac_purity$cell_id)
  }

gene_coords_table <- fread("fig1_species_mixing/data/genes.gtf.gz",
                           sep="\t", col.names = c("seqname", "source", "feature", "start", "end",
                                                   "score", "strand", "frame", "attributes")) %>%
  as_tibble() %>%
  filter(feature == "gene") %>%
  mutate(
    gene_id = str_match(attributes, 'gene_id "([^"]*)"')[,2],
    gene_type = str_match(attributes, 'gene_type "([^"]*)"')[,2]
  )

tss_sites <- gene_coords_table %>% 
  makeGRangesFromDataFrame() %>%
  resize(1, fix="start") %>%
  .[str_detect(seqnames(.), "(mm10|hg38)\\.chr[0-9X]+")] %>%
  keepSeqlevels(., seqlevelsInUse(.))

tss_background_flank <- c(
  #Positive Flank
  GRanges(seqnames(tss_sites), IRanges(end(tss_sites) + 1901, end(tss_sites) + 2000)),
  #Negative Flank
  GRanges(seqnames(tss_sites), IRanges(start(tss_sites) - 2000, start(tss_sites) - 1901))
) 

tss_background_count <- peakMatrix(cell_insertions, tss_background_flank) %>% colSums()
tss_window_count <- peakMatrix(cell_insertions, resize(tss_sites, 101, fix="center")) %>% colSums()

tss_enrichment <- tss_window_count / 101 / pmax(tss_background_count/200, 0.1)

ATAC_qc <- tibble(
  cell_id = names(tss_enrichment),
  "TSS_enrichment" = tss_enrichment
)

cell_calls <- full_join(atac_calls, rna_calls, by="cell_id") %>%
  full_join(hto_calls, by="cell_id") %>%
  full_join(ATAC_qc, by="cell_id") %>%
  mutate(
    ATAC_pass = replace_na(ATAC_pass & (TSS_enrichment > 10), FALSE),
    RNA_pass = replace_na(RNA_pass, FALSE),
    ADT_pass = replace_na(ADT_pass, FALSE),
    consensus_pass = ATAC_pass & RNA_pass & ADT_pass,
    consensus_species = case_when(
      ATAC_species != RNA_species ~ "mix",
      ADT_hash_count != 1 ~ "mix",
      TRUE ~ ATAC_species
    )
  )
cell_calls %>%
  filter(ATAC_pass | RNA_pass) %>%
  select(cell_id, consensus_pass, consensus_species, ATAC_pass, ATAC_species, 
         RNA_pass, RNA_species, ADT_pass, ADT_hash_assignment, ADT_hash_count,
         TSS_enrichment, hg38_fragments, mm10_fragments, hg38_genes, mm10_genes, hg38_umis, mm10_umis) %>%
  write_tsv(output_path)

read_csv("fig1_species_mixing/data/ADT_table.csv.gz") %>%
  mutate(cell_id = str_c(cell, "-1")) %>%
  select(cell_id, GATA1, SOX2, OCT4, starts_with("NPC")) %>%
  filter(cell_id %in% cell_calls$cell_id[cell_calls$ATAC_pass | cell_calls$RNA_pass]) %>%
  write_csv(adt_counts_output)


cell_groups <- cell_calls %>% 
  mutate(species = case_when(consensus_pass ~ consensus_species, TRUE ~ "low yield")) %>%
  pull(species, name=cell_id) %>% 
  .[as.character(cell_insertions@cell_levels)]
tss_profile <- insertionProfile(cell_insertions, resize(tss_sites, 4001, fix="center"), 
                                cell_groups = cell_groups)


rm(cell_insertions)
gc()

length_distribution <- fragments_mix_raw %>% lazy_dt() %>%
  transmute(width = end-start + 1, species = cell_groups[as.character(barcode)]) %>%
  filter(species %in% c("human", "mouse")) %>%
  group_by(species, width) %>%
  summarize(count = n()) %>%
  group_by(species) %>%
  mutate(fraction = count/sum(count)) %>%
  collect()


### Plots ###
plot_hto_cutoffs <- function(raw_counts, cutoffs) {
  raw_counts %>% clr_normalize(base=2) %>% t() %>%
    as_tibble(rownames="cell_id") %>%
    pivot_longer(rownames(raw_counts)) %>%
    ggplot(aes(value, color=name)) +
    geom_density(key_glyph=draw_key_path) +
    geom_vline(data=cutoffs, aes(xintercept=cutoff, color=hto), 
               linetype="dashed",key_glyph=draw_key_path) +
    theme_cowplot()
}


plot_barnyard <- function(purity_data, min_reads) {
  stats <- purity_data %>%
    filter(reads > min_reads) %>%
    summarize(
      mouse_cells = sum(species == "mouse"),
      human_cells = sum(species == "human"),
      mix_cells = sum(species == "mix")
    ) %>%
    mutate(
      mouse_human_doublet_rate = mix_cells/(mouse_cells + human_cells + mix_cells),
      doublet_rate = mouse_human_doublet_rate * 
        (mouse_cells^2 + 2*mouse_cells*human_cells + human_cells^2) /
        (2*mouse_cells*human_cells)
    )
  
  purity_data %>% 
    filter(reads > min_reads) %>%
    ggplot(aes(hg38, mm10, color=species)) + geom_point() +
    labs(
      subtitle=sprintf("%d mouse, %d human, %d mix\nInferred doublet rate: %.1f%%", stats$mouse_cells, stats$human_cells, stats$mix_cells, stats$doublet_rate*100),
      y="mm10 reads",
      x="hg38 reads",
      color="species"
    ) +
    scale_color_manual(values=species_colors_full[unique(purity_data$species)]) +
    theme_cowplot() 
}

plot_purity_cutoff_histogram <- function(purity_data, min_reads, hg38_cutoff, mm10_cutoff) {
  purity_data %>%
    filter(reads > min_reads) %>%
    ggplot(aes(hg38/reads, color=nearest_species, fill=nearest_species)) + 
    geom_density(alpha=0.5) +
    labs(y="Density", x="Fraction hg38 reads") +
    scale_x_continuous(limits=c(0,1), labels=scales::percent_format()) +  
    scale_color_manual(values=species_colors) +
    scale_fill_manual(values=species_colors) +
    geom_vline(xintercept = hg38_cutoff, linetype="dashed") + 
    geom_vline(xintercept = 1-mm10_cutoff, linetype="dashed") + 
    theme_cowplot()
}

median_annotation <- function(value, species) {
  sprintf("Median: %s human, %s mouse",
          scales::comma(median(value[species == "human"])),
          scales::comma(median(value[species == "mouse"])))
}

plot_adt <- function(marker) {
  adt_clr %>% inner_join(cell_calls, by=c("cell"="cell_id")) %>%
    filter(consensus_pass, consensus_species %in% c("mouse", "human")) %>%
    ggplot(aes_string(marker, color="consensus_species", fill="consensus_species")) +
    geom_density(alpha=0.5, adjust=1.25) +
    scale_color_manual(values=species_colors) +
    scale_fill_manual(values=species_colors) +
    theme_cowplot() +
    labs(x = sprintf("%s CLR", marker), color="", fill="")
}

pdf(plots_output, useDingbats = FALSE)


plot_hto_cutoffs(hto_matrix_raw, hto_classifications$cutoffs) +
  scale_color_brewer(palette="Dark2") +
  labs(x="CLR", y="density", color="")

plot_adt("GATA1")
plot_adt("SOX2")
plot_adt("OCT4")

plot_purity_cutoff_histogram(rna_mix_purity, min_reads=7500 , hg38_cutoff=.7, mm10_cutoff=.95) +
  ggtitle("RNA species purity")

ggplot(NULL, aes(rank(-rna_mix_purity$reads), rna_mix_purity$reads)) +
  scale_x_log10() + scale_y_log10() +
  ggrastr::geom_point_rast() +
  geom_hline(yintercept = 7500) +
  ggtitle("RNA minimum reads cutoff")

rna_mix_genes %>%
  inner_join(cell_calls, by="cell_id") %>%
  filter(consensus_pass & consensus_species != "mix") %>% {
  ggplot(., aes(consensus_species, umis, fill=consensus_species, color=consensus_species)) + 
  geom_violin(alpha = 0.5) +
  geom_boxplot(width=0.5, outlier.shape = NA, alpha=0.5, key_glyph=draw_key_blank) +
  scale_fill_manual(values=species_colors) +
  scale_color_manual(values=species_colors) +
  scale_y_log10() +
  theme_cowplot() +
  labs(title="RNA Mix UMIs per cell", x="species", fill="species", color="species",
       subtitle=median_annotation(.$umis, .$consensus_species))
  }

rna_mix_genes %>%
  inner_join(cell_calls, by="cell_id") %>%
  filter(consensus_pass & consensus_species != "mix") %>% {
  ggplot(., aes(consensus_species, genes, fill=consensus_species, color=consensus_species)) + 
  geom_violin(alpha=0.5) +
  geom_boxplot(width=0.5, outlier.shape = NA, alpha=0.5, key_glyph=draw_key_blank) +
  scale_fill_manual(values=species_colors) +
  scale_color_manual(values=species_colors) +
  scale_y_log10() +
  theme_cowplot() +
  labs(title="RNA Mix Genes per cell", x="species", fill="species", color="species",
       subtitle=median_annotation(.$genes, .$consensus_species))
  }


plot_purity_cutoff_histogram(multiome_atac_purity, min_reads=10000 , hg38_cutoff=.95, mm10_cutoff=.95) +
  ggtitle("ATAC species purity")

ggplot(NULL, aes(rank(-multiome_atac_purity$reads), multiome_atac_purity$reads)) +
  scale_x_log10() + scale_y_log10() +
  ggrastr::geom_point_rast() +
  geom_hline(yintercept = 10000) +
  ggtitle("ATAC minimum reads cutoff")

atac_genome_counts %>%
  inner_join(cell_calls, by=c("barcode" = "cell_id")) %>%
  filter(consensus_pass & consensus_species != "mix") %>%
  mutate(fragments = ifelse(consensus_species == "human", hg38_counts, mm10_counts)) %>% {
  ggplot(., aes(consensus_species, fragments, color=consensus_species, fill=consensus_species)) + 
  geom_violin(alpha=0.5) +
  geom_boxplot(width=0.5, outlier.shape = NA, alpha=0.5, key_glyph=draw_key_blank) +
  scale_fill_manual(values=species_colors) +
  scale_color_manual(values=species_colors) +
  scale_y_log10() +
  theme_cowplot() +
  labs(title="ATAC Mix Fragments per cell", x="species", fill="species", color="species",
       subtitle=median_annotation(.$fragments, .$consensus_species))
  }


ggplot(length_distribution, aes(width, fraction, color=species)) +
  geom_line() +
  scale_y_continuous(labels=scales::label_percent()) +
  scale_color_manual(values=species_colors) +
  xlim(0,500) +
  theme_cowplot() +
  labs(x="Fragment length", y="Percent of reads", title="Fragment length distribution")

tss_tibble_species <- tibble(
  count = as.vector(tss_profile),
  position = rep(-2000:2000, ncol(tss_profile)),
  group = rep(colnames(tss_profile), each=nrow(tss_profile))
) %>%
  group_by(group) %>%
  mutate(enrichment = count*200 / sum(tss_background_count[cell_groups == cur_group()$group]))

tss_tibble_species %>%
  filter(group %in% c("human", "mouse")) %>%
  ggplot(aes(position, enrichment, color=group)) +
  geom_line() +
  scale_color_manual(values=species_colors) +
  theme_cowplot() +
  labs(x="Distance from TSS (bp)", y="Normalized Insertions", color="species", title="TSS profile")


ggplot(NULL, aes(log10(multiome_atac_purity$reads), tss_enrichment, color=cell_groups)) +
  scale_x_continuous(limits=c(log10(100),NA), labels=scales::label_math(), breaks=2:5) +
  stat_scatter_density() +
  scale_color_gradientn(colors=ArchR::paletteContinuous("sambaNight")) +
  guides(color=FALSE) +
  geom_vline(xintercept = log10(10000), linetype="dashed") + 
  geom_hline(yintercept = 10, linetype="dashed") +
  theme_cowplot() + 
  labs(x="Unique Fragments", y = "TSS Enrichment", title="TSS Enrichment cutoff",
       subtitle=sprintf("Median TSS enrichment: %.2f", 
                        median(tss_enrichment[multiome_atac_purity$reads > 10e4 &
                                                tss_enrichment > 10])))

multiome_atac_purity %>%
  inner_join(cell_calls) %>%
  inner_join(hto_calls) %>%
  mutate(tss_enrichment = tss_enrichment[as.character(cell_id)]) %>%
  filter(consensus_pass, tss_enrichment > 10) %>%
  mutate(
    species = ATAC_species
  ) %>%
  plot_barnyard(10000) + ggtitle("ATAC Mix barnyard (no HTO filtering)")

multiome_atac_purity %>%
  inner_join(cell_calls) %>%
  inner_join(hto_calls) %>%
  mutate(tss_enrichment = tss_enrichment[as.character(cell_id)]) %>%
  filter(consensus_pass, tss_enrichment > 10) %>%
  mutate(
    species = case_when(ADT_hash_count > 1 ~ "hash mix", TRUE ~ consensus_species)
  ) %>%
  plot_barnyard(10000) + ggtitle("ATAC Mix barnyard (HTO labeled and removed)")

multiome_atac_purity %>%
  inner_join(cell_calls) %>%
  inner_join(hto_calls) %>%
  mutate(tss_enrichment = tss_enrichment[as.character(cell_id)]) %>%
  filter(consensus_pass, tss_enrichment > 10) %>%
  mutate(
    species = case_when(ADT_hash_count > 1 ~ "hash mix", TRUE ~ consensus_species)
  ) %>%
  filter(species != "hash mix") %>%
  plot_barnyard(10000) + ggtitle("ATAC Mix barnyard (filtered)")


rna_mix_purity %>%
  inner_join(cell_calls) %>%
  inner_join(hto_calls) %>%
  mutate(tss_enrichment = tss_enrichment[as.character(cell_id)]) %>%
  mutate(
    species = RNA_species
  ) %>%
  filter(consensus_pass, tss_enrichment > 10, !is.na(consensus_species)) %>%
  plot_barnyard(7500) + ggtitle("RNA Mix barnyard (no filtering)")


rna_mix_purity %>%
  inner_join(cell_calls) %>%
  inner_join(hto_calls) %>%
  mutate(tss_enrichment = tss_enrichment[as.character(cell_id)]) %>%
  mutate(
    species = case_when(ADT_hash_count > 1 ~ "hash mix", is.na(consensus_species) ~ "empty", TRUE ~ consensus_species)
  ) %>%
  filter(consensus_pass, tss_enrichment > 10, consensus_species != "empty") %>%
  plot_barnyard(7500) + ggtitle("RNA Mix barnyard (HTO labeled and removed)")

rna_mix_purity %>%
  inner_join(cell_calls) %>%
  inner_join(hto_calls) %>%
  mutate(tss_enrichment = tss_enrichment[as.character(cell_id)]) %>%
  mutate(
    species = case_when(ADT_hash_count > 1 ~ "hash mix", is.na(consensus_species) ~ "empty", TRUE ~ consensus_species)
  ) %>%
  filter(consensus_pass, tss_enrichment > 10, consensus_species != "empty", species != "hash mix") %>%
  plot_barnyard(7500) + ggtitle("RNA Mix barnyard")


dev.off()