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
library(stringr)

library(httpgd)
#hgd()

data_processing_pipe <- function(seurat_obj) {
  # RNA analysis
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- NormalizeData(object = seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(object = seurat_obj, selection.method = 'vst', nfeatures = 2000)

  # ATAC
  DefaultAssay(seurat_obj) <- "ATAC"
  seurat_obj <- RunTFIDF(seurat_obj)
  seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = 20)
  seurat_obj <- RunSVD(seurat_obj)

  # Run Harmony for scATAC integration
  seurat_obj <- RunHarmony(
    object = seurat_obj,
    group.by.vars = 'dataset',
    reduction = 'lsi',
    assay.use = 'ATAC',
    reduction.save = "harmony_lsi",
    max.iter.harmony = 5e3,
    project.dim = FALSE
  )

  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- SCTransform(seurat_obj, vars.to.regress = c("percent.mt"), verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, npcs = 100)

  # Run Harmony for scRNA integration
  seurat_obj <- RunHarmony(
    object = seurat_obj,
    group.by.vars = 'dataset',
    reduction = 'pca',
    assay.use = 'SCT',
    reduction.save = "harmony_pca",
    max.iter.harmony = 1e3,
    project.dim = FALSE
  )

  seurat_obj <- FindMultiModalNeighbors(seurat_obj, reduction.list = list("harmony_pca", "harmony_lsi"), dims.list = list(1:50, 1:50))
  seurat_obj <- FindClusters(seurat_obj, graph.name = "wsnn", algorithm = 3, verbose = FALSE, res = 1.2)

  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, reduction = "harmony_pca", reduction.name = "umap.rna", reduction.key = "rnaUMAP_")
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, reduction = "harmony_lsi", reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  seurat_obj <- RunUMAP(seurat_obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
}

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

# ------------------------------------------------------------------------------
# EN subtype annotation
multiomes_integrated <- 
  readRDS("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/RNA_ATAC/multiomes_integrated_20240512_infant_adult_annotated.rds")
EN <- multiomes_integrated[, str_starts(multiomes_integrated$celltype, "EN")]
rm(multiomes_integrated)
gc()

DotPlot(EN, features = markers, col.min = 0, scale.by = 'size') + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 3))

#EN <- multiomes_integrated[, as.integer(as.character(multiomes_integrated$wsnn_res.0.8)) %in% c(5, 8, 15, 16, 20, 24)]

# RNA & ATAC analysis
EN <- data_processing_pipe(EN)

p1 <- DimPlot(EN, reduction = "umap.rna", group.by = "dataset") + ggtitle("RNA")
p2 <- DimPlot(EN, reduction = "umap.atac", group.by = "dataset") + ggtitle("ATAC")
p3 <- DimPlot(EN, reduction = "wnn.umap", group.by = "dataset") + ggtitle("WNN")
ggpubr::ggarrange(p1, p2, p3, common.legend = T, nrow = 1, legend = "bottom", widths = c(1, 1, 1))

p1 <- DimPlot(EN, reduction = "umap.rna", group.by = "seurat_clusters", label.size = 5.5, repel = TRUE, label = T) + ggtitle("RNA")
p2 <- DimPlot(EN, reduction = "umap.atac", group.by = "seurat_clusters", label.size = 5.5, repel = TRUE, label = T) + ggtitle("ATAC")
p3 <- DimPlot(EN, reduction = "wnn.umap", group.by = "seurat_clusters", label.size = 5.5, repel = TRUE, label = T) + ggtitle("WNN")
ggpubr::ggarrange(p1, p2, p3, common.legend = T, nrow = 1, legend = "bottom", widths = c(1, 1, 1))

DotPlot(EN, features = markers, col.min = 0, scale.by = 'size') + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 3))
ggsave("./figures/WSNN/cell_type/EN_marker_expression_20240715.pdf", width = 12, height = 12)

# Remove potential doublets
EN <- EN[, !(as.integer(as.character(EN$seurat_clusters)) %in% c(11, 12, 18, 21, 25, 26))]

# RNA & ATAC analysis
EN <- data_processing_pipe(EN)

p1 <- DimPlot(EN, reduction = "umap.rna", group.by = "seurat_clusters", label.size = 5.5, repel = TRUE, label = T) + ggtitle("RNA")
p2 <- DimPlot(EN, reduction = "umap.atac", group.by = "seurat_clusters", label.size = 5.5, repel = TRUE, label = T) + ggtitle("ATAC")
p3 <- DimPlot(EN, reduction = "wnn.umap", group.by = "seurat_clusters", label.size = 5.5, repel = TRUE, label = T) + ggtitle("WNN")
ggpubr::ggarrange(p1, p2, p3, common.legend = T, nrow = 1, legend = "bottom", widths = c(1, 1, 1))
ggsave("./figures/WSNN/cell_type/EN_UMAP_cluster_remove_doublets_20240715.pdf", width = 12, height = 5)


# Label transfer
allen_brain_EN <- readRDS("./data/RNA_ATAC/allen_brain_map/MTG_SMARTseq/exc_forannot_processed_June10_2020.rds")

DefaultAssay(EN) <- "SCT"
EN_anchors <- FindTransferAnchors(reference = allen_brain_EN, query = EN, dims = 1:30,
                                        reference.reduction = "pca")
predictions <- TransferData(anchorset = EN_anchors, refdata = allen_brain_EN$cluster, dims = 1:30)
EN <- AddMetaData(EN, metadata = predictions)

p1 <- DimPlot(EN, reduction = "umap.rna", group.by = "predicted.id", label.size = 5.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(EN, reduction = "umap.atac", group.by = "predicted.id", label.size = 5.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(EN, reduction = "wnn.umap", group.by = "predicted.id", label.size = 5.5, repel = TRUE) + ggtitle("WNN")
ggpubr::ggarrange(p1, p2, p3, common.legend = T, nrow = 1, legend = "bottom", widths = c(1, 1, 1))
ggsave("./figures/WSNN/cell_type/EN_UMAP_subcelltype_transfer_20240715.pdf", width = 12.5, height = 6)

# annotate upper and deep layer
# upper layer: L2 + L2-3 + L2-4
# deep layer: others
EN[[c("predicted.id", "seurat_clusters")]] %>% 
  group_by(seurat_clusters, predicted.id) %>% 
  summarise("count" = n()) -> cluster_celltype

EN[["cortex_layer"]] <- ifelse(EN$wsnn_res.1.2 %in% c(0, 1, 2, 3, 4, 5, 14, 9, 25, 27),
                              "upper_layer", "deep_layer")

p1 <- DimPlot(EN, reduction = "umap.rna", group.by = "cortex_layer", label.size = 5.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(EN, reduction = "umap.atac", group.by = "cortex_layer", label.size = 5.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(EN, reduction = "wnn.umap", group.by = "cortex_layer", label.size = 5.5, repel = TRUE) + ggtitle("WNN")
ggpubr::ggarrange(p1, p2, p3, common.legend = T, nrow = 1, legend = "bottom", widths = c(1, 1, 1))
ggsave("./figures/WSNN/cell_type/EN_UMAP_subcelltype_manual_20240715.pdf", width = 12, height = 5)


EN_markers <- markers <- c(
  "CUX1", "POU3F2", "RORB"
)

DotPlot(EN, features = EN_markers, col.min = 0, 
        scale.by = 'size', group.by = "cortex_layer") + 
  theme(panel.grid = element_blank(), panel.border = element_blank()) +
  labs(x = "Gene", y = "Subtype")
ggsave("./figures/WSNN/cell_type/control_brain_EN_marker_expression_20240715.pdf", width = 6, height = 3.5)

write.table(EN[[c("dataset", "cell.id", "predicted.id", "cortex_layer")]], 
            "./results/celltype_annotation/EN_annotation_transfer_20240715.txt", sep = "\t",
            quote = F, col.names = F, row.names = F)

saveRDS(EN, "./data/RNA_ATAC/cell_type/EN_20240715.rds")
rm(EN_anchors)
rm(predictions)
rm(allen_brain_EN)

# ------------------------------------------------------------------------------
# IN subtype annotation
multiomes_integrated <- 
  readRDS("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/RNA_ATAC/multiomes_integrated_20240512_infant_adult_annotated.rds")
IN <- multiomes_integrated[, str_starts(multiomes_integrated$celltype, "IN")]
rm(multiomes_integrated)
gc()

DotPlot(IN, features = markers, col.min = 0, scale.by = 'size') + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 3))

# RNA & ATAC analysis
IN <- data_processing_pipe(IN)

p1 <- DimPlot(IN, reduction = "umap.rna", group.by = "dataset", label.size = 5.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(IN, reduction = "umap.atac", group.by = "dataset", label.size = 5.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(IN, reduction = "wnn.umap", group.by = "dataset", label.size = 5.5, repel = TRUE) + ggtitle("WNN")
ggpubr::ggarrange(p1, p2, p3, common.legend = T, nrow = 1, legend = "bottom", widths = c(1, 1, 1))
#ggsave("./figures/WSNN/cell_type/IN_UMAP_sample_20240715.pdf", width = 12, height = 5)

p1 <- DimPlot(IN, reduction = "umap.rna", group.by = "seurat_clusters", label.size = 5.5, repel = TRUE, label = T) + ggtitle("RNA")
p2 <- DimPlot(IN, reduction = "umap.atac", group.by = "seurat_clusters", label.size = 5.5, repel = TRUE, label = T) + ggtitle("ATAC")
p3 <- DimPlot(IN, reduction = "wnn.umap", group.by = "seurat_clusters", label.size = 5.5, repel = TRUE, label = T) + ggtitle("WNN")
ggpubr::ggarrange(p1, p2, p3, common.legend = T, nrow = 1, legend = "bottom", widths = c(1, 1, 1))
ggsave("./figures/WSNN/cell_type/IN_UMAP_cluster_before_doublets_20240715.pdf", width = 12, height = 5)

DotPlot(IN, features = markers, col.min = 0, scale.by = 'size') + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 3))
ggsave("./figures/WSNN/cell_type/IN_marker_expression_20240715.pdf", width = 12, height = 14)

IN <- IN[, !(as.integer(as.character(IN$seurat_clusters)) %in% c(20, 21, 26, 27, 30, 34))]

# RNA & ATAC analysis
IN <- data_processing_pipe(IN)

p1 <- DimPlot(IN, reduction = "umap.rna", group.by = "seurat_clusters", label.size = 5.5, repel = TRUE, label = T) + ggtitle("RNA")
p2 <- DimPlot(IN, reduction = "umap.atac", group.by = "seurat_clusters", label.size = 5.5, repel = TRUE, label = T) + ggtitle("ATAC")
p3 <- DimPlot(IN, reduction = "wnn.umap", group.by = "seurat_clusters", label.size = 5.5, repel = TRUE, label = T) + ggtitle("WNN")
ggpubr::ggarrange(p1, p2, p3, common.legend = T, nrow = 1, legend = "bottom", widths = c(1, 1, 1))
ggsave("./figures/WSNN/cell_type/IN_UMAP_cluster_remove_doublets_20240715.pdf", width = 12, height = 5)

DotPlot(IN, features = markers, col.min = 0, scale.by = 'size') + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 3))

allen_brain_IN <- readRDS("./data/RNA_ATAC/allen_brain_map/MTG_SMARTseq/inh_forannot_processed_June10_2020.rds")

DefaultAssay(IN) <- "SCT"
IN_anchors <- FindTransferAnchors(reference = allen_brain_IN, query = IN, dims = 1:30,
                                  reference.reduction = "pca")
predictions <- TransferData(anchorset = IN_anchors, refdata = allen_brain_IN$cluster, dims = 1:30)
IN <- AddMetaData(IN, metadata = predictions)

# annotate subtypes
IN[[c("predicted.id", "seurat_clusters")]] %>% 
  group_by(seurat_clusters, predicted.id) %>% 
  summarise("count" = n()) -> cluster_celltype

DotPlot(IN, features = c("GAD1", "ADARB2", "LAMP5", "PAX6", "VIP", "LHX6", "SST", "PVALB"), col.min = 0, scale.by = 'size') + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 3))

IN <- IN[, (as.integer(as.character(IN$seurat_clusters)) != 0)]
IN[["subtype"]] <- ifelse(IN$wsnn_res.1.2 %in% c(4, 5, 6, 11, 15, 18, 19, 21, 22, 23, 25, 27, 28, 29, 30, 12, 13, 14, 16),
                               "CGE", "MGE")

IN[[c("subtype", "predicted.id")]] %>% 
  group_by(predicted.id, subtype) %>% 
  summarise("count" = n()) -> cluster_subcelltype

p1 <- DimPlot(IN, reduction = "umap.rna", group.by = "predicted.id", label.size = 5.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(IN, reduction = "umap.atac", group.by = "predicted.id", label.size = 5.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(IN, reduction = "wnn.umap", group.by = "predicted.id", label.size = 5.5, repel = TRUE) + ggtitle("WNN")
ggpubr::ggarrange(p1, p2, p3, common.legend = T, nrow = 1, legend = "bottom", widths = c(1, 1, 1))
ggsave("./figures/WSNN/cell_type/IN_UMAP_subcelltype_transfer_20240715.pdf", width = 12, height = 6)


p1 <- DimPlot(IN, reduction = "umap.rna", group.by = "subtype", label.size = 5.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(IN, reduction = "umap.atac", group.by = "subtype", label.size = 5.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(IN, reduction = "wnn.umap", group.by = "subtype", label.size = 5.5, repel = TRUE) + ggtitle("WNN")
ggpubr::ggarrange(p1, p2, p3, common.legend = T, nrow = 1, legend = "bottom", widths = c(1, 1, 1))
ggsave("./figures/WSNN/cell_type/IN_UMAP_subcelltype_manual_20240715.pdf", width = 12, height = 5)

IN_markers <- markers <- c(
  "ADARB2", "LHX6",
  "VIP", "SST", "PVALB", "RELN" # IN subtypes
)

DotPlot(IN, features = IN_markers, col.min = 0, 
        scale.by = 'size', group.by = "subtype") + 
  theme(panel.grid = element_blank(), panel.border = element_blank()) +
  labs(x = "Gene", y = "Subtype")
ggsave("./figures/WSNN/cell_type/control_brain_IN_marker_expression_20240715.pdf", width = 8, height = 3.5)


write.table(IN[[c("dataset", "cell.id", "predicted.id", "subtype")]], 
            "./results/celltype_annotation/IN_annotation_transfer_20240715.txt", sep = "\t",
            quote = F, col.names = F, row.names = F)

saveRDS(IN, "./data/RNA_ATAC/cell_type/IN_20240715.rds")

