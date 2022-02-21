# Need to run ADT_normalization.R and CD4_ArchR_plots.R first

suppressPackageStartupMessages({
  library(tidyverse)
  library(ArchR)
  library(glue)
})


proj <- loadArchRProject("fig2_CD4_Tcells/outputs/ArchR_HTOsinglets_CD4only/")


ADT_5x <- read.csv("fig2_CD4_Tcells/outputs/ADT_5x_NPCnorm.csv")

ADT_5x$cell2 <- paste0(ADT_5x$lane, "#", ADT_5x$cell, "-1")

ADT_5x <- filter(ADT_5x, cell2 %in% proj$cellNames)

# Add TF ADT data
proj_5X <- addCellColData(proj, data = ADT_5x[['Helios_norm']], name = "Helios_ADT", cells = ADT_5x[["cell2"]])
proj_5X <- addCellColData(proj_5X, data = ADT_5x[['FOXP3_norm']], name = "FOXP3_ADT", cells = ADT_5x[["cell2"]])
proj_5X <- addCellColData(proj_5X, data = ADT_5x[['GATA3_norm']], name = "GATA3_ADT", cells = ADT_5x[["cell2"]])
proj_5X <- addCellColData(proj_5X, data = ADT_5x[['Tbet_norm']], name = "Tbet_ADT", cells = ADT_5x[["cell2"]])
proj_5X <- addCellColData(proj_5X, data = ADT_5x[['RORgT_norm']], name = "RORgT_ADT", cells = ADT_5x[["cell2"]])

# Need to remove NAs before allowed to plot
idxPass <- which(is.na(proj_5X$Helios_ADT) == FALSE)
proj_5X <- proj_5X[idxPass,]

data_exp <- getCellColData(proj_5X)

write.table(data_exp, file = "Supplementary_figures/outputs/CD4cells_RNA_5XADT.csv", sep = ",", quote = F, row.names = T)


plot_hist_5x <- plotGroups(proj_5X,
                            groupBy = "Clusters",
                            name =  c("Helios_ADT", "Tbet_ADT", "GATA3_ADT", "FOXP3_ADT", "RORgT_ADT"),
                            colorBy = "cellColData",
                            plotAs = "ridges")

plotPDF(plot_hist_5x,
        name = "Plot-5X_ADTs-bycluster-NPCnorm.pdf",
        ArchRProj = proj_5X,
        addDOC = FALSE, width = 5, height = 5)
