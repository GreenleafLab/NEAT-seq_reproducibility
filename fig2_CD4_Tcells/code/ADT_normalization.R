# Need to run CD4_HTO_singlet_ADT_counts.R first

library(tidyverse)

ADTs <- read.csv("fig2_CD4_Tcells/outputs/TF_ADT_counts_singlets_from_NPCandHTO.csv")


# split into different tables based on titration, then normalize to NPC
ADT_5x <- filter(ADTs, NPC1_pos == 1)
ADT_25x <- filter(ADTs, NPC2_pos == 1)

ADT_5x$FOXP3_norm <- log((250*(ADT_5x$FOXP3/ADT_5x$NPC1) + 1), 2)
ADT_5x$Helios_norm <- log((250*(ADT_5x$Helios/ADT_5x$NPC1) + 1), 2)
ADT_5x$Tbet_norm <- log((250*(ADT_5x$Tbet/ADT_5x$NPC1) + 1), 2)
ADT_5x$GATA3_norm <- log((250*(ADT_5x$GATA3/ADT_5x$NPC1) + 1), 2)
ADT_5x$RORgT_norm <- log((250*(ADT_5x$RORgT/ADT_5x$NPC1) + 1), 2)

ADT_5x_norm <- ADT_5x[, c(1,9,17:21)]

ADT_25x$FOXP3_norm <- log((250*(ADT_25x$FOXP3/ADT_25x$NPC2) + 1), 2)
ADT_25x$Helios_norm <- log((250*(ADT_25x$Helios/ADT_25x$NPC2) + 1), 2)
ADT_25x$Tbet_norm <- log((250*(ADT_25x$Tbet/ADT_25x$NPC2) + 1), 2)
ADT_25x$GATA3_norm <- log((250*(ADT_25x$GATA3/ADT_25x$NPC2) + 1), 2)
ADT_25x$RORgT_norm <- log((250*(ADT_25x$RORgT/ADT_25x$NPC2) + 1), 2)

ADT_25x_norm <- ADT_25x[,c(1,9, 17:21)]

write.table(ADT_5x_norm, file = "fig2_CD4_Tcells/outputs/ADT_5x_NPCnorm.csv", sep = ",", quote = F, row.names = F)

write.table(ADT_25x_norm, file = "fig2_CD4_Tcells/outputs/ADT_25x_NPCnorm.csv", sep = ",", quote = F, row.names = F)



