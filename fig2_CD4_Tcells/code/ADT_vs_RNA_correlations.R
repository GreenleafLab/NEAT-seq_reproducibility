# Need to run ArchR_CD4cells_add25xADT.R first

library(tidyverse)
library(ggpubr)


data <- read.csv("fig2_CD4_Tcells/outputs/CD4cells_RNA_25XADT.csv")


FOXP3_plot <- ggscatterhist(data, 
                             x = "FOXP3_ADT", y = "FOXP3_RNA", 
                             xlab= "ADT", ylab = "RNA", title = "\nFOXP3", 
                             margin.params = list(fill = "#d9d9d9"), margin.plot = "density",
                             group = "Clusters", 
                             palette = c("#d9d9d9","#272E6A","#d9d9d9","#d9d9d9","#d9d9d9", "#d9d9d9","#d9d9d9"),
                             shape = 19, size = 2, alpha = 0.4)

GATA3_plot <- ggscatterhist(data, 
                            x = "GATA3_ADT", y = "GATA3_RNA", 
                            xlab= "ADT", ylab = "RNA", title = "\nGATA3", 
                            margin.params = list(fill = "lightgray"), margin.plot = "density",
                            group = "Clusters", 
                            palette = c("#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9","#F47D2B", "#d9d9d9","#d9d9d9"),
                            shape = 19, size = 2, alpha = 0.4)

Helios_plot <- ggscatterhist(data, 
                             x = "Helios_ADT", y = "IKZF2_RNA", 
                             xlab= "ADT", ylab = "RNA", title = "\nHelios", 
                             margin.params = list(fill = "lightgray"), margin.plot = "density",
                             group = "Clusters", 
                             palette = c("#d9d9d9","#272E6A","#d9d9d9","#d9d9d9","#d9d9d9", "#d9d9d9","#d9d9d9"),
                             shape = 19, size = 2, alpha = 0.4)

Tbet_plot <- ggscatterhist(data, 
                           x = "Tbet_ADT", y = "TBX21_RNA", 
                           xlab= "ADT", ylab = "RNA", title = "\nTbet", 
                           margin.params = list(fill = "lightgray"), margin.plot = "density",
                           group = "Clusters", 
                           palette = c("#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9", "#FEE500","#d9d9d9"),
                           shape = 19, size = 2, alpha = 0.4)

RORgT_plot <- ggscatterhist(data, 
                            x = "RORgT_ADT", y = "RORC_RNA", 
                            xlab= "ADT", ylab = "RNA", title = "\nRORgT", 
                            margin.params = list(fill = "lightgray"), margin.plot = "density",
                            group = "Clusters", 
                            palette = c("#d9d9d9","#d9d9d9","#208A42","#d9d9d9","#d9d9d9", "#d9d9d9","#d9d9d9"),
                            shape = 19, size = 2, alpha = 0.4)


pdf(file = "fig2_CD4_Tcells/outputs/ADT_vs_RNA_counts.pdf", width = 7, height = 7, useDingbats = FALSE)

print(FOXP3_plot)
print(Helios_plot)
print(GATA3_plot)
print(Tbet_plot)
print(RORgT_plot)

dev.off()
