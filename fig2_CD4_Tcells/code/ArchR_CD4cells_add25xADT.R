# Need to run ADT_normalization.R and CD4_ArchR_plots.R first

suppressPackageStartupMessages({
  library(tidyverse)
  library(ArchR)
  library(glue)
})


proj <- loadArchRProject("fig2_CD4_Tcells/outputs/ArchR_HTOsinglets_CD4only/")


ADT_25x <- read.csv("fig2_CD4_Tcells/outputs/ADT_25x_NPCnorm.csv")

ADT_25x$cell2 <- paste0(ADT_25x$lane, "#", ADT_25x$cell, "-1")

ADT_25x <- filter(ADT_25x, cell2 %in% proj$cellNames)

# Add TF ADT data
proj_25X <- addCellColData(proj, data = ADT_25x[['Helios_norm']], name = "Helios_ADT", cells = ADT_25x[["cell2"]])
proj_25X <- addCellColData(proj_25X, data = ADT_25x[['FOXP3_norm']], name = "FOXP3_ADT", cells = ADT_25x[["cell2"]])
proj_25X <- addCellColData(proj_25X, data = ADT_25x[['GATA3_norm']], name = "GATA3_ADT", cells = ADT_25x[["cell2"]])
proj_25X <- addCellColData(proj_25X, data = ADT_25x[['Tbet_norm']], name = "Tbet_ADT", cells = ADT_25x[["cell2"]])
proj_25X <- addCellColData(proj_25X, data = ADT_25x[['RORgT_norm']], name = "RORgT_ADT", cells = ADT_25x[["cell2"]])

# Need to remove NAs before allowed to plot
idxPass <- which(is.na(proj_25X$Helios_ADT) == FALSE)
proj_25X <- proj_25X[idxPass,]

proj_25X <- saveArchRProject(proj_25X, outputDirectory = "fig2_CD4_Tcells/outputs/ArchR_HTOsinglets_CD4only_ADT")


data_exp <- getCellColData(proj_25X)

write.table(data_exp, file = "fig2_CD4_Tcells/outputs/CD4cells_RNA_25XADT.csv", sep = ",", quote = F, row.names = T)


plot_ADT_25X <- plotEmbedding(
  ArchRProj = proj_25X,
  colorBy = "cellColData",
  name = c("Helios_ADT", "Tbet_ADT", "GATA3_ADT", "FOXP3_ADT", "RORgT_ADT"),
  embedding = "UMAP",
  plotAs = "points",
  quantCut = c(0.02, 0.98),
  imputeWeights = NULL
)


plotPDF(plotList = plot_ADT_25X,
        name = "Plot-UMAP-conc_25X_NPCnorm-no-Imputation.pdf",
        ArchRProj = proj_25X,
        addDOC = FALSE, width = 5, height = 5)


proj_25X <- addImputeWeights(proj_25X)


plot_ADT_25X_imp <- plotEmbedding(
  ArchRProj = proj_25X,
  colorBy = "cellColData",
  name = c("Helios_ADT", "Tbet_ADT", "GATA3_ADT", "FOXP3_ADT", "RORgT_ADT"),
  embedding = "UMAP",
  plotAs = "hex",
  quantCut = c(0.02, 0.98),
  imputeWeights = getImputeWeights(proj_25X)
)


plotPDF(plotList = plot_ADT_25X_imp,
        name = "Plot-UMAP-conc_25X_NPCnorm_with-Imputation.pdf",
        ArchRProj = proj_25X,
        addDOC = FALSE, width = 5, height = 5)



plot_hist_25x <- plotGroups(proj_25X,
                            groupBy = "Clusters",
                            name =  c("Helios_ADT", "Tbet_ADT", "GATA3_ADT", "FOXP3_ADT", "RORgT_ADT"),
                            colorBy = "cellColData",
                            plotAs = "ridges")


plotPDF(plot_hist_25x,
        name = "Plot-ADTs-bycluster-NPCnorm.pdf",
        ArchRProj = proj_25X,
        addDOC = FALSE, width = 5, height = 5)
