# Need to run Seurat_markers.R and ArchR_CD4cells_add25xADT.R first

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


combined <- readRDS("fig2_CD4_Tcells/outputs/Tmem_allCD4cells.rds")


# load normalized ADT data

ADT <- read.csv("fig2_CD4_Tcells/outputs/CD4cells_RNA_25XADT.csv")

ADT <- tibble::rownames_to_column(ADT, "cell_id")

ADT$cell_id2 <- gsub("#", "_", ADT$cell_id)

ADT_only <- ADT[,c('Helios_ADT', "RORgT_ADT", "FOXP3_ADT", "Tbet_ADT", "GATA3_ADT", "cell_id2")]

ADT_only_filtered <- filter(ADT_only, cell_id2 %in% Cells(combined))

rownames(ADT_only_filtered) <- ADT_only_filtered[,c('cell_id2')]

ADT_only_filtered[,c('cell_id2')] <- NULL

# need matrix to have cells as column names, so transpose

ADT_transposed <- t(as.matrix(ADT_only_filtered))

ADT_assay <- CreateAssayObject(data = as.matrix(ADT_transposed))


# need to subset Seurat object to only have cells in ADT table before adding ADT matrix

combined <- subset(combined, cells = rownames(ADT_only_filtered))


# add ADT matrix to Seurat object and save RDS
combined[["ADT"]] <- ADT_assay

saveRDS(combined, file = "fig2_CD4_Tcells/outputs/Tmem_CD4cells_withADT.rds")



# fetch GATA3 data as data frame to set cutoffs for "high" and "low" GATA3

Gata3_dat <- FetchData(combined, vars = c('GATA3', 'GATA3-ADT', 'ident'))

# filter for high GATA3 RNA > 2.25

Gata3_highRNA <- filter(Gata3_dat, GATA3 > 2.25)


# filter for high GATA3 RNA and high GATA3 protein > 6.121

Gata3_highRNA_highADT <- filter(Gata3_highRNA, `adt_GATA3-ADT` > 6.121)


# get the lowest GATA3 protein cells within high GATA3 RNA population, matching num of cells in high protein population

Gata3_highRNA_ordered <- Gata3_highRNA[order(Gata3_highRNA$`adt_GATA3-ADT`),]
Gata3_highRNA_lowADT <- Gata3_highRNA_ordered[c(1:nrow(Gata3_highRNA_highADT)),]

# note that low GATA3 protein cutoff is 4.911566

# extract cell barcodes for each group

highRNA_highADT_cells <- row.names(Gata3_highRNA_highADT)

highRNA_lowADT_cells <- row.names(Gata3_highRNA_lowADT)

 

library(MASS)

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

Gata3_dat$dens <- get_density(Gata3_dat$GATA3, Gata3_dat$'adt_GATA3-ADT', n = 100)
p1 <- ggplot(Gata3_dat) + geom_point(aes(x = GATA3, y = `adt_GATA3-ADT`, color = dens)) + scale_color_viridis_c() +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(aspect.ratio=1) + xlab("RNA") + ylab("ADT") +
  geom_vline(xintercept = 2.25, linetype = "dashed", size =1) +
  geom_hline(yintercept = 6.121, linetype = "dashed", size =1) +
  geom_hline(yintercept = 4.911566, linetype = "dashed", size =1)

pdf(file = "fig2_CD4_Tcells/outputs/GATA3_cutoffs.pdf", width = 6, height = 6)

print(p1)

dev.off()


# differential expression

diff_allgenes <- FindMarkers(combined, slot = "data", ident.1 = highRNA_highADT_cells, ident.2 = highRNA_lowADT_cells, logfc.threshold = 0)
diff_allgenes$Log2FC <- log(exp(diff_allgenes$avg_logFC), 2)
diff_allgenes$adj_pval_BH <- p.adjust(diff_allgenes$p_val, method = "BH")

# add column for plotting colors
diff_allgenes$ann <- "blank"
diff_allgenes$ann[abs(diff_allgenes$Log2FC) > 0.5 & diff_allgenes$adj_pval_BH < 0.05] <- "sigdiff"
diff_allgenes$ann[abs(diff_allgenes$Log2FC) > 0.5 & diff_allgenes$adj_pval_BH > 0.05] <- "diff"
diff_allgenes$ann[abs(diff_allgenes$Log2FC) < 0.5 & diff_allgenes$adj_pval_BH < 0.05] <- "sig"


genes <- c("CAMK4", "GAB2", "EEF1G", "SNX9", "NELL2", 
           "PABPC4", "NIBAN1", "SNED1", "BTBD11", "RPL18A")

diff_allgenes <- tibble::rownames_to_column(diff_allgenes, "Gene")

# add labels
diff_allgenes$labels <- NA
diff_allgenes$labels[diff_allgenes$Gene %in% genes] <- diff_allgenes$Gene[diff_allgenes$Gene %in% genes]

vplot <- ggplot(diff_allgenes) + geom_point(aes(x = Log2FC, y = -log(adj_pval_BH, 10), color = ann)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(aspect.ratio=1) + xlab("Log2 Fold Change (High/Low)") + ylab("-log10 padj") +
  geom_vline(xintercept = 0.5, linetype = "dashed", size =1) +
  geom_vline(xintercept = -0.5, linetype = "dashed", size =1) +
  geom_hline(yintercept = -log(0.05, 10), linetype = "dashed", size =1) + 
  ggrepel::geom_text_repel(aes(x = Log2FC, y = -log(adj_pval_BH, 10), 
                               label = labels), size = 6, point.padding = 0.25, 
                           min.segment.length = 0.2) + 
  scale_color_manual(values=c("black", "#00B26E", "#1F7CEA", "#D60000"))


pdf(file = "fig2_CD4_Tcells/outputs/volcano_diff_wilcox_BH.pdf", width = 8, height = 8, useDingbats = FALSE)

print(vplot)

dev.off()


