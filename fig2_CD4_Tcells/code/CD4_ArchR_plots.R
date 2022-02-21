# Need to run CD4_HTO_singlets_ADT_counts.R first

suppressPackageStartupMessages({
  library(tidyverse)
  library(ArchR)
  library(readr)
  library(dplyr)
  library(motifmatchr)
  library(chromVARmotifs)
})

set.seed(1)

addArchRGenome("hg38")

# create ArchR project

fragment_files <- c(
  "lane1"="geo_download/lane1_atac_fragments.tsv.gz",
  "lane2"="geo_download/lane2_atac_fragments.tsv.gz"
)

arrows_lane1 <- createArrowFiles(
  inputFiles = fragment_files["lane1"],
  sampleNames = "lane1",
  QCDir = "fig2_CD4_Tcells/outputs/QualityControl"
)
arrows_lane2 <- createArrowFiles(
  inputFiles = fragment_files["lane2"],
  sampleNames = "lane2",
  QCDir = "fig2_CD4_Tcells/outputs/QualityControl"
)
proj <- ArchRProject(
  ArrowFiles = c(arrows_lane1, arrows_lane2),
  outputDirectory = "fig2_CD4_Tcells/outputs/ArchR_HTOsinglets_CD4only/",
  copyArrows = FALSE
)

proj <- saveArchRProject(proj, "fig2_CD4_Tcells/outputs/ArchR_HTOsinglets_CD4only/")


# filter for singlets and high quality cells

singlets <- read.csv("fig2_CD4_Tcells/outputs/TF_ADT_counts_singlets_from_NPCandHTO.csv")

singlets$cell2 <- paste0(singlets$lane, "#", singlets$cell, "-1")

proj_cells <- as.data.frame(proj$cellNames)
names(proj_cells)[1] <- "cell_id"

singlets_inproj <- filter(singlets, cell2 %in% proj_cells$cell_id)

singlet_cellNames <- as.character(singlets_inproj$cell2)


proj <- subsetCells(proj, cellNames = singlet_cellNames)

keeper_cells <- proj$cellNames[proj$TSSEnrichment > 10 & proj$nFrags > 1000]

proj <- subsetCells(proj, cellNames = keeper_cells)





# filter for CD4 memory T cells (i.e remove contaminating CD8 T cells)

classification <- read_tsv("fig2_CD4_Tcells/data/MPAL_alignment.tsv")


Tcells <- filter(classification, cell_type == "22_CD4.M")

proj <- subsetCells(proj, cellNames = Tcells$cell_id)



proj <- saveArchRProject(proj)

# Cluster and plot UMAP

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 4)
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")
proj <- addImputeWeights(proj)
proj <- saveArchRProject(proj)

umap_plot <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

plotPDF(umap_plot,
        name = "Plot-UMAP.pdf",
        ArchRProj = proj,
        addDOC = FALSE, width = 5, height = 5)


# Add peak and motif deviations matrices

proj <- addGroupCoverages(proj, groupBy="Clusters")
proj <- addReproduciblePeakSet(proj, groupBy="Clusters")
proj <- saveArchRProject(proj)

proj <- addPeakMatrix(proj)


proj <- addMotifAnnotations(proj)
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(proj, peakAnnotation="Motif")
proj <- saveArchRProject(proj)


# Plot gene accessibility scores

markerGenes  <- c(
  "FOXP3",
  "GATA3",
  "TBX21",
  "IKZF2",
  "RORC",
  "TCF7",
  "LEF1",
  "IFNG",
  "IL4",
  "IL10",
  "IL17A"
  )

p_GS <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "GeneScoreMatrix",
  name = markerGenes,
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = getImputeWeights(proj)
)

plotPDF(plotList = p_GS,
        name = "Plot-UMAP-Marker-Genes-with-Imputation.pdf",
        ArchRProj = proj,
        addDOC = FALSE, width = 5, height = 5)



# plot motif deviations

motifs <- c("GATA3",
            "FOXP3",
            "TBX21",
            "RORC",
            "TCF7",
            "LEF1",
            "NFATC1")

markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs <- grep("z:", markerMotifs, value = TRUE)

p_motifs <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "MotifMatrix",
  name = sort(markerMotifs),
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj)
)

plotPDF(plotList = p_motifs,
        name = "Plot-UMAP-chromvar-Zscores-with-Imputation.pdf",
        ArchRProj = proj,
        addDOC = FALSE, width = 5, height = 5)

# add and plot Helios motif deviations

helios_pwms <- readRDS("fig2_CD4_Tcells/data/helios_pwms.Rds")

motif_matches <- motifmatchr::matchMotifs(
  pwms = helios_pwms,
  subject = getPeakSet(proj),
  genome = BSgenome.Hsapiens.UCSC.hg38
)

proj <- addDeviationsMatrix(proj, matches=motif_matches, matrixName="HeliosDeviations", threads=1)

mat <- getMatrixFromProject(proj, useMatrix = "HeliosDeviations")

for(n in rownames(mat)) {
  proj <- addCellColData(proj, data=assay(mat, "z")[n,], name=str_c("z_", n), cells = colnames(mat))
}

# PFM_1 is ChIP-seq motif

Helios <- plotEmbedding(proj, name="z_Helios_PFM_1", imputeWeights = getImputeWeights(proj))

plotPDF(Helios,
        name = "Plot-UMAP-ChromVAR-Helios_ChIP-seq.pdf",
        ArchRProj = proj,
        addDOC = FALSE, width = 5, height = 5)


# plot enriched motifs in clusters

# ID marker peaks
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)


heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 5, transpose = TRUE, cutOff = 5)

plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)


# Add and plot RNA

RNA_counts <- readRDS("geo_download/CD4_RNA_counts.rds")


gene_IDs <- c(
  "RORC" = "ENSG00000143365",
  "TBX21" = "ENSG00000073861",
  "FOXP3" = "ENSG00000049768",
  "GATA3" = "ENSG00000107485",
  "IKZF2" = "ENSG00000030419",
  "IL10" = "ENSG00000136634",
  "IFNG" = "ENSG00000111537",
  "IL4" = "ENSG00000113520",
  "IL17A" = "ENSG00000112115",
  "CCR7" = "ENSG00000126353",
  "NFATC1" = "ENSG00000131196"
)



RNA <- tibble(
  total_umis = colSums(RNA_counts),
  total_genes = colSums(RNA_counts != 0),
  cell_id = colnames(RNA_counts)
)

for (gene in names(gene_IDs)) {
  RNA[[gene]] <- RNA_counts[gene_IDs[gene], ]
  name <- paste0(gene, "_10knorm")
  RNA[[name]] <- log(((10000*(RNA[[gene]]/RNA$total_umis))+1), 2)
}



RNA$cell_id2 <- gsub("_", "#", RNA$cell_id)

RNA <- filter(RNA, cell_id2 %in% proj$cellNames)
RNA_cells <- as.character(RNA$cell_id2)

for (gene in names(gene_IDs)) {
  col_name <- paste0(gene, "_10knorm")
  out_name <- paste0(gene, "_RNA")
  proj <- addCellColData(proj, data = RNA[[col_name]], name = out_name, cells = RNA_cells, force = TRUE)
}

proj <- saveArchRProject(proj)

idxPass <- which(is.na(proj$RORC_RNA) == FALSE)
proj <- proj[idxPass,]

proj <- addImputeWeights(proj)

tf_RNA_names <- c()

for (gene in names(gene_IDs)) {
  out_name <- paste0(gene, "_RNA")
  tf_RNA_names <- c(tf_RNA_names, out_name)
}

plot_RNA <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = tf_RNA_names,
  embedding = "UMAP",
  plotAs = "hex",
  quantCut = c(0.02, 0.98),
  imputeWeights = getImputeWeights(proj)
)

plotPDF(plotList = plot_RNA,
        name = "Plot-UMAP-RNA-with-Imputation.pdf",
        ArchRProj = proj,
        addDOC = FALSE, width = 5, height = 5)


