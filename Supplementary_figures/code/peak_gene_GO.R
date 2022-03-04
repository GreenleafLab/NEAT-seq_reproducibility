library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

targets <- read.csv("geo_download/GSM5396330_putative_targets.csv.gz")

# get background gene set - all genes with at least one read across cells

all_genes <- readRDS("geo_download/GSM5396333_CD4_RNA_counts.rds")

all_genes_mat <- as.matrix(all_genes)


exp_genes <- all_genes_mat[rowSums(all_genes_mat[]) > 0,]

exp_gene_names <- as.data.frame(dimnames(exp_genes)[[1]])
colnames(exp_gene_names) <- "gene_id"


# RORgT

RORgT <- targets %>% filter(adt == "RORgT")


GO_RORgT <- enrichGO(gene         = RORgT$gene_id,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff  = 0.05,
                     universe = exp_gene_names$gene_id)

p_RORgT <- dotplot(GO_RORgT, showCategory = 15, title = "RORgT")


# GATA3

GATA3 <- targets %>% filter(adt == "GATA3")

GO_GATA3 <- enrichGO(gene         = GATA3$gene_id,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff  = 0.05,
                     universe = exp_gene_names$gene_id)

p_GATA3 <- dotplot(GO_GATA3, showCategory = 15, title = "GATA3")


# Tbet

Tbet <- targets %>% filter(adt == "Tbet")

GO_Tbet <- enrichGO(gene         = Tbet$gene_id,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    qvalueCutoff  = 0.05,
                    universe = exp_gene_names$gene_id)

p_Tbet <- dotplot(GO_Tbet, showCategory = 15, title = "Tbet")

pdf("Supplementary_figures/outputs/Peak_Gene_GOanalysis.pdf", useDingbats = FALSE, width = 10, height = 10)

print(p_RORgT)
print(p_GATA3)
print(p_Tbet)

dev.off()

