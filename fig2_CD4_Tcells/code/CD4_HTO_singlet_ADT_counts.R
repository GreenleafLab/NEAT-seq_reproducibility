library(tidyverse)

## Identify doublets based on first round of hashing oligos

dir.create(dirname('fig2_CD4_Tcells/outputs/TF_ADT_counts_singlets_from_NPCandHTO.csv'), showWarnings = FALSE, recursive = TRUE)

HTO_lane1 <- read.csv('geo_download/HTO_counts_lane1.csv')
HTO_lane2 <- read.csv('geo_download/HTO_counts_lane2.csv')

HTO_lane1$lane <- "lane1"
HTO_lane2$lane <- "lane2"


HTO_compiled <- rbind(HTO_lane1, HTO_lane2)

HTO_compiled$total_counts <- rowSums(HTO_compiled[,2:6])


# filter for total counts greater than 75

HTO_filtered <- filter(HTO_compiled, total_counts > 75)


# CLR normalize and then mark as positive/negative for each HTO

HTO_filtered$geom_mean <- (log(HTO_filtered$HTO2+1, 2) + 
  log(HTO_filtered$HTO3+1, 2) + 
  log(HTO_filtered$HTO4+1, 2) + 
  log(HTO_filtered$HTO5+1, 2) + 
  log(HTO_filtered$HTO6+1, 2))/5


HTO_filtered$HTO2_CLR <- log((HTO_filtered$HTO2 + 1)/(HTO_filtered$geom_mean), 2)
HTO_filtered$HTO3_CLR <- log((HTO_filtered$HTO3 + 1)/(HTO_filtered$geom_mean), 2)
HTO_filtered$HTO4_CLR <- log((HTO_filtered$HTO4 + 1)/(HTO_filtered$geom_mean), 2)
HTO_filtered$HTO5_CLR <- log((HTO_filtered$HTO5 + 1)/(HTO_filtered$geom_mean), 2)
HTO_filtered$HTO6_CLR <- log((HTO_filtered$HTO6 + 1)/(HTO_filtered$geom_mean), 2)


# cutoffs for positive:
# HTO2: 3
# HTO3: 3.5
# HTO4: 3
# HTO5: 3
# HTO6: 3.5


HTO_CLRnorm <- HTO_filtered[,c(1,8,10:14)]

hashed <- mutate(HTO_CLRnorm, HTO2_pos = ifelse(HTO2_CLR < 3, 0, 1), HTO3_pos = ifelse(HTO3_CLR < 3.5, 0, 1), 
                       HTO4_pos = ifelse(HTO4_CLR < 3, 0, 1), HTO5_pos = ifelse(HTO5_CLR < 3, 0, 1), 
                       HTO6_pos = ifelse(HTO6_CLR < 3.5, 0, 1))

hashed$Hash_sum <- rowSums(hashed[,8:12])

singlets <- filter(hashed, Hash_sum == 1)


## load ADT and titration hash oligos and find doublets based on titration hash oligos

ADT_lane1 <- read.csv("geo_download/ADT_counts_lane1.csv")
ADT_lane2 <- read.csv("geo_download/ADT_counts_lane2.csv")


ADT_lane1$lane <- "lane1"
ADT_lane2$lane <- "lane2"

ADT_compiled <- rbind(ADT_lane1, ADT_lane2)

ADT_compiled$total_counts <- rowSums(ADT_compiled[,2:8])

ADT_filtered <- filter(ADT_compiled, total_counts > 100)


# calculate CLR for NPCs

ADT_filtered$geom_mean_NPC <- (log(ADT_filtered$NPC1+1, 2) + 
                             log(ADT_filtered$NPC2+1, 2))/2

ADT_filtered$NPC1_CLR <- log((ADT_filtered$NPC1 + 1)/(ADT_filtered$geom_mean_NPC), 2)
ADT_filtered$NPC2_CLR <- log((ADT_filtered$NPC2 + 1)/(ADT_filtered$geom_mean_NPC), 2)


# cutoffs for positive:
# NPC1: 4
# NPC2: 3.5

hashed_NPC <- mutate(ADT_filtered, NPC1_pos = ifelse(NPC1_CLR < 4, 0, 1), NPC2_pos = ifelse(NPC2_CLR < 3.5, 0, 1))

hashed_NPC$Hash_sum <- rowSums(hashed_NPC[,14:15])

singlets_NPC <- filter(hashed_NPC, Hash_sum == 1)

singlets_NPC_and_HTO <- semi_join(singlets_NPC, singlets, by = c("cell" = "cell"))

write.table(singlets_NPC_and_HTO, file = "fig2_CD4_Tcells/outputs/TF_ADT_counts_singlets_from_NPCandHTO.csv", sep = ",", quote = F, row.names = F)
