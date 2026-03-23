# Libraries
.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")
library(Seurat)
library(SeuratDisk)
library(Signac)
library(GenomicRanges)
library(future)
library(ggplot2)
library(harmony)
library(ComplexHeatmap)
library(biovizBase)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)

# read in validated data and new data
multiomes_integrated <- readRDS("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/RNA_ATAC/multiomes_integrated_20240813_ASD_no_annotation.rds")
allen_brain <- readRDS("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/RNA_ATAC/allen_brain_map/HumanMultipleCorticalAreasSMARTseq/seurat_obj.rds")
DefaultAssay(multiomes_integrated) <- "RNA"

allen_brain$cell_type <- allen_brain$subclass_label
allen_brain$cell_type[allen_brain$class_label == "Glutamatergic"] <- "EN"
allen_brain$cell_type[allen_brain$class_label == "GABAergic"] <- "IN"
table(allen_brain$cell_type)

# transfer cell type labels
anchor <- FindTransferAnchors(reference = allen_brain, query = multiomes_integrated, dims = 1:30)
prediction <- TransferData(anchorset = anchor, refdata = allen_brain$cell_type, dims = 1:30)

multiomes_integrated <- AddMetaData(multiomes_integrated, metadata = prediction[,1, drop=F])

multiomes_integrated$predicted.id
sum(multiomes_integrated$predicted.id == "")
multiomes_integrated$predicted.id[multiomes_integrated$predicted.id == ""] <- "Unknown"
table(multiomes_integrated$predicted.id, multiomes_integrated$seurat_clusters)

# set each cluster as their most dominant cell type
multiomes_integrated[[c("predicted.id", "seurat_clusters")]] %>% 
  group_by(seurat_clusters, predicted.id) %>% 
  summarise("count" = n()) %>%
  group_by(seurat_clusters) %>%
  filter(count == max(count)) -> cluster_celltype

Idents(multiomes_integrated)

multiomes_integrated$celltype_transferred <- 
  factor(multiomes_integrated$seurat_clusters, levels = cluster_celltype$seurat_clusters, labels = cluster_celltype$predicted.id) 

table(multiomes_integrated$celltype_transferred, multiomes_integrated$seurat_clusters)

# Plot with transferred cell type annotations
DimPlot(multiomes_integrated, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("RNA")

png("figures/RNAseq/Umap_Integrated_original_groupby_celltype_transferred_20240813_ASD.png", res = 300, height = 2000, width = 2200)
DimPlot(multiomes_integrated, reduction = "umap.rna", group.by = "celltype_transferred", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("RNA")
dev.off()
p1 <- DimPlot(multiomes_integrated, reduction = "umap.rna", group.by = "celltype_transferred", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(multiomes_integrated, reduction = "umap.atac", group.by = "celltype_transferred", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(multiomes_integrated, reduction = "wnn.umap", group.by = "celltype_transferred", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("WNN")
png("figures/WSNN/Umap_WNNIntegrated_original_groupby_celltype_transferred_20240813_ASD.png", res = 300, height = 2200, width = 6000)
ggpubr::ggarrange(p1, p2, p3, nrow = 1, common.legend = T, legend = "bottom")
dev.off()

# Manually check
# cell type markers
markers <- c(
  "SYT1", "SNAP25", "NRGN", # All neurons
  "SATB2", "SLC17A6", "SLC17A7", # EN
  "GAD1", "GAD2", # IN
  "VIP", "SST", "PVALB", "RELN", # IN subtypes
  "SLC1A2", "SLC1A3", "AQP4", # Astrocyte
  "APBB1IP", # Microglia
  "SOX6", # OPC?
  "PDGFRA", # OPC
  "MOG", "PLP1" # Oligodendrocyte
)

DotPlot(multiomes_integrated, features = markers, col.min = 0, scale.by = 'size') + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 3))
ggsave("./figures/WSNN/cell_type/ASD_marker_expression_20240813.pdf", width = 10, height = 10)

multiomes_integrated$celltype <- as.character(multiomes_integrated$celltype_transferred)
multiomes_integrated$celltype[multiomes_integrated$celltype_transferred == "Unknown"] <- "OPC"
multiomes_integrated$celltype[multiomes_integrated$seurat_clusters == "6"] <- "Oligodendrocyte"

table(multiomes_integrated$celltype)

# Plot with new cell type annotations
p1 <- DimPlot(multiomes_integrated, reduction = "umap.rna", group.by = "celltype") + ggtitle("RNA")
p2 <- DimPlot(multiomes_integrated, reduction = "umap.atac", group.by = "celltype") + ggtitle("ATAC")
p3 <- DimPlot(multiomes_integrated, reduction = "wnn.umap", group.by = "celltype") + ggtitle("WNN")
ggpubr::ggarrange(p1, p2, p3, nrow = 1, common.legend = T, legend = "bottom")
ggsave("figures/WSNN/Umap_WNNIntegrated_original_groupbycelltype_manually_checked_20240813_ASD.pdf", width=12)


saveRDS(multiomes_integrated, "/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/RNA_ATAC/multiomes_integrated_20240813_ASD_annotated.rds")

write.table(multiomes_integrated[[c("dataset", "cell.id", "celltype")]], 
            "./results/celltype_annotation/annotation_table_20240813_ASD.txt", sep = "\t",
            quote = F, col.names = F, row.names = F)

multiomes_integrated[[c("dataset", "cell.id", "seurat_clusters", "celltype")]] %>% 
  group_by(seurat_clusters, celltype) %>% summarise() -> celltype_cluster_table

write.table(celltype_cluster_table, 
            "./results/celltype_annotation/annotation_cluster_celltype_20240813_ASD.txt", sep = "\t",
            quote = F, col.names = T, row.names = F)

