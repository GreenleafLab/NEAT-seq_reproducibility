# Run CD4_HTO_singlet_ADT_counts.R and ArchR_CD4cells_add25xADT.R first

library(dplyr)
library(ggplot2)
library(ArchR)
library(cowplot)

dir.create(dirname("Supplementary_figures/outputs/ADT_norm_test.pdf"), showWarnings = FALSE, recursive = TRUE)



ADTs <- read.csv("fig2_CD4_Tcells/outputs/TF_ADT_counts_singlets_from_NPCandHTO.csv")

cell_info <- read.csv('fig2_CD4_Tcells/outputs/CD4cells_RNA_25XADT.csv')

ADTs$cell2 <- paste0(ADTs$lane, "#", ADTs$cell, "-1")

merged <- filter(ADTs, cell2 %in% rownames(cell_info))

# RORgT

merged$RORgT_NPCnorm <- log((250*(merged$RORgT/merged$NPC2) + 1), 2)
merged$RORgT_totalnorm <- log((250*(merged$RORgT/merged$total_counts) + 1), 2)


# FOXP3
merged$FOXP3_NPCnorm <- log((250*(merged$FOXP3/merged$NPC2) + 1), 2)
merged$FOXP3_totalnorm <- log((250*(merged$FOXP3/merged$total_counts) + 1), 2)

# Helios
merged$Helios_NPCnorm <- log((250*(merged$Helios/merged$NPC2) + 1), 2)
merged$Helios_totalnorm <- log((250*(merged$Helios/merged$total_counts) + 1), 2)

# GATA3
merged$GATA3_NPCnorm <- log((250*(merged$GATA3/merged$NPC2) + 1), 2)
merged$GATA3_totalnorm <- log((250*(merged$GATA3/merged$total_counts) + 1), 2)

# Tbet
merged$Tbet_NPCnorm <- log((250*(merged$Tbet/merged$NPC2) + 1), 2)
merged$Tbet_totalnorm <- log((250*(merged$Tbet/merged$total_counts) + 1), 2)

proj <- loadArchRProject("fig2_CD4_Tcells/outputs/ArchR_HTOsinglets_CD4only_ADT/")


TFs <- c("RORgT",
         "FOXP3",
         "Helios",
         "GATA3",
         "Tbet")



for (gene in TFs) {
  proj <-addCellColData(proj, data = merged[,paste0(gene, "_NPCnorm")], 
                        name = paste0(gene, "_NPCnorm"), 
                        cells = as.character(merged$cell2))
  
  proj <-addCellColData(proj, data = merged[,paste0(gene, "_totalnorm")], 
                        name = paste0(gene, "_totalnorm"), 
                        cells = as.character(merged$cell2))

  proj <-addCellColData(proj, data = log((merged[,gene]+1), 2), 
                        name = paste0(gene, "_log"), 
                        cells = as.character(merged$cell2) )
}



col <- c()

for (gene in TFs) {
  col <- c(col, paste0(gene, "_log"), paste0(gene, "_NPCnorm"), paste0(gene, "_totalnorm"))
}


p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "CellColData", 
  name = col, 
  embedding = "UMAP",
  quantCut = c(0.02, 0.98),
  imputeWeights = NULL,
  plotAs = "points",
  size = 2
)


pdf(file = "Supplementary_figures/outputs/ADT_norm_test.pdf", width = 5, height = 5)
print(p)
dev.off()



