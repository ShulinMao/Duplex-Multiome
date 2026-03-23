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

cell_type <- read.table("./results/celltype_annotation/annotation_xuyu_20240119.txt")

# read in validated data and new data
multiomes_integrated <- readRDS("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/RNA_ATAC/multiomes_integrated_20240512_infant_adult_no_annotation.rds")
multiomes_integrated_old <- readRDS("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/RNA_ATAC/multiomes_integrated_20240119_infant_adult_xuyu_annotated.rds")

table(multiomes_integrated_old$celltype)

# transfer cell type labels
library(stringr)
multiomes_integrated_old[["cell_type_no_id"]] <- factor(multiomes_integrated_old$celltype, 
                                                        levels = levels(multiomes_integrated_old$celltype),
                                                        labels = str_remove(levels(multiomes_integrated_old$celltype), "[-]{0,1}[0-9]$"))

table(multiomes_integrated_old[["cell_type_no_id"]])

anchor <- FindTransferAnchors(reference = multiomes_integrated_old, query = multiomes_integrated, dims = 1:30)
prediction <- TransferData(anchorset = anchor, refdata = multiomes_integrated_old$cell_type_no_id, dims = 1:30)

multiomes_integrated <- AddMetaData(multiomes_integrated, metadata = prediction[,1, drop=F])
multiomes_integrated$predicted.id
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

# Plot with transferred cell type annotations
png("figures/RNAseq/Umap_Integrated_original_groupby_celltype_transferred_20240512_infant_adult.png", res = 300, height = 2000, width = 2200)
DimPlot(multiomes_integrated, reduction = "umap.rna", group.by = "celltype_transferred", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("RNA")
dev.off()
p1 <- DimPlot(multiomes_integrated, reduction = "umap.rna", group.by = "celltype_transferred", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(multiomes_integrated, reduction = "umap.atac", group.by = "celltype_transferred", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(multiomes_integrated, reduction = "wnn.umap", group.by = "celltype_transferred", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("WNN")
png("figures/WSNN/Umap_WNNIntegrated_original_groupbycelltype_transferred_20240512.png", res = 300, height = 2200, width = 6000)
ggpubr::ggarrange(p1, p2, p3, nrow = 1, common.legend = T, legend = "bottom")
dev.off()

saveRDS(multiomes_integrated, "/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/RNA_ATAC/multiomes_integrated_20240512_infant_adult_transfered_annotation.rds")

write.table(cluster_celltype[,1:2], "/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/results/celltype_annotation/annotation_transferring_20240512.txt", 
            append = F, quote = F, col.names = T, row.names = F)

# Manually check
multiomes_integrated <- readRDS("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/RNA_ATAC/multiomes_integrated_20240512_infant_adult_transfered_annotation.rds")
cluster_celltype_checked <- read.table("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/results/celltype_annotation/annotation_final_20240512.txt", header = T)
multiomes_integrated$celltype <- 
  factor(multiomes_integrated$seurat_clusters, levels = cluster_celltype_checked$seurat_clusters, labels = cluster_celltype_checked$predicted.id)
table(multiomes_integrated$celltype)

multiomes_integrated <- subset(multiomes_integrated, cell = colnames(multiomes_integrated)[multiomes_integrated$celltype != "remove"])
table(multiomes_integrated$celltype)

# Plot with new cell type annotations
p1 <- DimPlot(multiomes_integrated, reduction = "umap.rna", group.by = "celltype", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(multiomes_integrated, reduction = "umap.atac", group.by = "celltype", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(multiomes_integrated, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("WNN")
png("figures/WSNN/Umap_WNNIntegrated_original_groupbycelltype_manually_checked_20240512.png", res = 300, height = 2200, width = 6000)
ggpubr::ggarrange(p1, p2, p3, nrow = 1, common.legend = T, legend = "bottom")
dev.off()

p1 <- DimPlot(multiomes_integrated, reduction = "umap.rna", group.by = "dataset") + ggtitle("RNA")
p2 <- DimPlot(multiomes_integrated, reduction = "umap.atac", group.by = "dataset") + ggtitle("ATAC")
p3 <- DimPlot(multiomes_integrated, reduction = "wnn.umap", group.by = "dataset") + ggtitle("WNN")
ggpubr::ggarrange(p1, p2, p3, nrow = 1, ncol = 3, common.legend = T, legend = "bottom")
ggsave("figures/WSNN/Umap_WNNIntegrated_original_groupbysample_manually_checked_20240512.pdf", width = 9, height = 4)

saveRDS(multiomes_integrated, "/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/RNA_ATAC/multiomes_integrated_20240512_infant_adult_annotated.rds")

write.table(multiomes_integrated[[c("dataset", "cell.id", "celltype")]], 
            "./results/celltype_annotation/annotation_table_20240512.txt", sep = "\t",
            quote = F, col.names = F, row.names = F)

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
  "MOG", "PLP1", # Oligodendrocyte
  "COL1A1",
  "PECAM1"
)

multiomes_integrated[,!str_detect(multiomes_integrated$celltype, "doublet")] -> multiomes_integrated
multiomes_integrated[["cell_type"]] <- as.character(multiomes_integrated$celltype)
multiomes_integrated$cell_type[str_detect(multiomes_integrated$celltype, "EN")] <- "Excitatory\nNeuron"
multiomes_integrated$cell_type[str_detect(multiomes_integrated$celltype, "IN")] <- "Inhibitory\nNeuron"

DotPlot(multiomes_integrated, features = markers, col.min = 0, 
        scale.by = 'size', group.by = "cell_type") + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 3)) +
  labs(x = "Gene", y = "Cell type")
ggsave("./figures/WSNN/cell_type/control_brain_marker_expression_20240512.pdf", width = 10, height = 5)

