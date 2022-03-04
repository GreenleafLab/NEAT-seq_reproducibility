## Need to remove "lane1_" and "lane2_" prefix from file names in filered_feature_bc_matrix folder before running

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


set.seed(100)

lane1 <- Read10X('geo_download/GSM5396333_lane1/')


# Create Seurat objects from each lane, then merge into one object
# Multiome RNA data also includes Peaks matrix so need to specify Gene Expression matrix when creating object
lane1_obj <- CreateSeuratObject(counts = lane1$`Gene Expression`, project = "lane1_RNA")


lane2 <- Read10X('geo_download/GSM5396337_lane2/')

lane2_obj <- CreateSeuratObject(counts = lane2$`Gene Expression`, project = "lane2_RNA")

combined <- merge(lane1_obj, y = lane2_obj, add.cell.ids = c("lane1", "lane2"), project = "RNA_all")


# load table of all CD4 T cells from ArchR project

cell_table <- read.csv("geo_download/GSM5396332_CD4cells.csv.gz")
cell_table <- tibble::rownames_to_column(cell_table, "cell_id")

# Need to reformat ArchR cell IDs so that # is changed to _ to match Seurat cell IDs
cell_table$cell_id2 <- gsub("#", "_", cell_table$cell_id)

# subset Seurat object to filter for cells in ArchR object

combined <- subset(combined, cells = cell_table$cell_id2)

combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)


# Find top 2000 most variable features that can be used for downstream analysis
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)

# scale data for dim reduction analyses

all.genes <- rownames(combined)
combined <- ScaleData(combined, features = all.genes)

# Run PCA

combined <- RunPCA(combined, features = VariableFeatures(object = combined))

ElbowPlot(combined)

# Cluster cells

combined <- FindNeighbors(combined, dims = 1:15)
combined <- FindClusters(combined, resolution = 0.6)

# plot UMAP

combined <- RunUMAP(combined, dims = 1:15)

saveRDS(combined, file = "fig2_CD4_Tcells/outputs/Tmem_allCD4cells.rds")

# Add ATAC cluster info

cells_filtered <- filter(cell_table, cell_id2 %in% Cells(combined))

ATAC_clusterID <- cells_filtered[,c('cell_id2', 'Clusters')]

ATAC_clusterID <- tibble::column_to_rownames(ATAC_clusterID, "cell_id2")

combined <- AddMetaData(combined, ATAC_clusterID, col.name = "ATAC_clusters")

umap2 <- DimPlot(combined, reduction = "umap", group.by = "ATAC_clusters",
                cols = c("C6"="#FEE500","C4"="#89288F","C3"="#208A42","C5"="#F47D2B","C2"="#272E6A","C7"="#8A9FD1", "C1"="#D51F26")) + theme(aspect.ratio=1)

pdf(file = "fig2_CD4_Tcells/outputs/Seurat_UMAP_ATAC_cluster_coloring.pdf", width = 8, height = 8, useDingbats = FALSE)

print(umap2)

dev.off()

markers <- c("CXCR3", "TBX21", #Th1
             "CCR4", "GATA3", #Th2
             "CCR6", "RORC", #Th17
             "IL2RA", "CTLA4", "FOXP3", "IKZF2", #Treg
             "CD38", "CD69" #activated T cells
)

p <- DotPlot(combined, features = markers, 
             group.by = "ATAC_clusters",
             cols =  c("#3F00FF","#F8FA0D"),
             dot.min = 0.005,
             scale.by = "size") + scale_size(range = c(2, 8)) + RotatedAxis()

pdf(file = "fig2_CD4_Tcells/outputs/Markerexp_ATAC_clusters.pdf", width = 8, height = 8, useDingbats = FALSE)
print(p)
dev.off()

