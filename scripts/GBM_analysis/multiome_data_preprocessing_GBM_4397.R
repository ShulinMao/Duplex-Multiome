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

source("./scripts/multiome_data_analysis/multiome_data_processing_functions/multiome_data_processing_functions.R")

# Read the folder names
data_path <- "/lab-share/Gene-Walsh-e2/Public/indData/AK/Multiome-strand-02212022/"
folder_names <- c(
  "Novogene_05152025/Mult_4397T/outs",
  "Novogene_05152025/Mult_4397N/outs"
)
sample_names <- c("UMB4397_tumor", "UMB4397_normal")


########## 1. Creating a common peak set ##########
# If the peaks were identified independently in each experiment then they will likely not overlap perfectly. We can merge peaks from all the datasets to create a common peak set, and quantify this peak set in each experiment prior to merging the objects.
# First we’ll load the peak coordinates for each experiment and convert them to genomic ranges, the use the GenomicRanges::reduce function to create a common set of peaks to quantify in each dataset.
# read in peak sets
peaks_list <- list()
for (i in 1:length(folder_names)) {
  folder.i <- folder_names[i]
  peak.i <- read.table(file = paste(data_path, folder.i, "/atac_peaks.bed", sep = ""), col.names = c("chr", "start", "end"))
  peaks_list[[sample_names[i]]] <- peak.i
}

# convert to genomic ranges
gr_list <- list()
for (i in 1:length(peaks_list)) {
  sample.i <- names(peaks_list)[i]
  gr_list[[sample.i]] <- makeGRangesFromDataFrame(peaks_list[[sample.i]])
}

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(
  gr_list$UMB4397_tumor,
  gr_list$UMB4397_normal))

# # Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths < 10000 & peakwidths > 20]
combined.peaks

# ########## 2. Create Fragment objects ##########
# # To quantify our combined set of peaks we’ll need to create a Fragment object for each experiment. The Fragment class is a specialized class defined in Signac to hold all the information related to a single fragment file.
# # First we’ll load the cell metadata for each experiment so that we know what cell barcodes are contained in each file, then we can create Fragment objects using the CreateFragmentObject function. The CreateFragmentObject function performs some checks to ensure that the file is present on disk and that it is compressed and indexed, computes the MD5 sum for the file and the tabix index so that we can tell if the file is modified at any point, and checks that the expected cells are present in the file.
# # load metadata
md_list <- list()
for (i in 1:length(folder_names)) {
  folder.i <- folder_names[i]
  md.i <- read.table(
    file = paste(data_path, folder.i, "/per_barcode_metrics.csv", sep = ""),
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )
  md_list[[sample_names[i]]] <- md.i
}
lapply(md_list, dim)

# # create fragment objects
fragobj_list <- list()
for (i in 1:length(folder_names)) {
  folder.i <- folder_names[i]
  counts.i <- Read10X_h5(paste(data_path, folder.i, "/filtered_feature_bc_matrix.h5", sep = ""))
  cat("Number of cells for", sample_names[i], ":", ncol(counts.i$Peaks), identical(colnames(counts.i$`Gene Expression`), colnames(counts.i$Peaks)), "\n")
  frags.i <- CreateFragmentObject(
    path = paste(data_path, folder.i, "/atac_fragments.tsv.gz", sep = ""),
    cells = colnames(counts.i$Peaks)
  )
  fragobj_list[[sample_names[i]]] <- frags.i
  md_list[[sample_names[i]]] <- md_list[[sample_names[i]]][colnames(counts.i$Peaks), ] # Update md metadata with only cells of interest
}
lapply(md_list, dim)

# ########## 3. Quantify peaks in each dataset ##########
# # We can now create a matrix of peaks x cell for each sample using the FeatureMatrix function. This function is parallelized using the future package. See the parallelization vignette for more information about using future.
# # Skip this step if FeatureMatrix is already created
counts_list <- list()
for (i in 1:length(folder_names)) {
  folder.i <- folder_names[i]
  counts.i <- FeatureMatrix(
    fragments = fragobj_list[[sample_names[i]]],
    features = combined.peaks,
    cells = rownames(md_list[[sample_names[i]]]),
    verbose = FALSE
  )
  counts_list[[sample_names[i]]] <- counts.i
}
lapply(counts_list, dim)

########## 4. Create the objects ##########
# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
# We will now use the quantified matrices to create a Seurat object for each dataset, storing the Fragment object for each dataset in the assay.
# counts_list <- readRDS("IntermediateResults/Peaks_FeatureMatrix_4samples.rds")
multiomes_list <- list()
for (i in 1:length(folder_names)) {
  folder_name <- folder_names[i]
  sample_name <- sample_names[i]
  print(sample_name)
  
  # load the ATAC data
  counts <- Read10X_h5(paste(data_path, folder_name, "/filtered_feature_bc_matrix.h5", sep = ""))
  rna_counts <- counts$`Gene Expression`
  
  # Create Seurat object
  multiomes <- CreateSeuratObject(counts = rna_counts, assay = "RNA")
  multiomes[["percent.mt"]] <- PercentageFeatureSet(multiomes, pattern = "^MT-")
  
  # Load scATAC data
  atac_assay <- CreateChromatinAssay(
    counts_list[[sample_name]], 
    sep = c(":", "-"), 
    min.cells = 10, 
    fragments = fragobj_list[[sample_name]],
    annotation = annotation
  )
  multiomes[["ATAC"]] <- atac_assay
  
  # QC violin plot before filtering
  v1 <- VlnPlot(
    object = multiomes,
    features = c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC", "percent.mt"),
    ncol = 5,
    pt.size = 0.005
  )
  png(paste("figures/ATAC/QC/Vlnplot_", sample_name, "_prefiltering_20250604_GBMQC.png", sep = ""), res = 300, height = 1400, width = 3000)
  print(v1)
  dev.off()
  
  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  plot1 <- FeatureScatter(multiomes, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(multiomes, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2
  
  # QC filtering scRNA
  rna_assay <- CreateSeuratObject(counts = rna_counts, assay = "RNA", min.cells = 3, min.features = 200)
  identifiedDoublets <- isDoublet(rna_assay) # Identifying doublets
  rna_assay <- rna_assay[, !identifiedDoublets]
  #rna_assay <- scQC(rna_assay, mtThreshold = 0.1, minLSize = 500)
  multiomes <- multiomes[, colnames(rna_assay)]
  
  # QC filtering scATAC
  multiomes <- subset(
    x = multiomes,
    subset = nCount_ATAC > 1000 &
      percent.mt < 10 &
      nFeature_RNA > 200
  )
  
  # QC violin plot after filtering
  v2 <- VlnPlot(
    object = multiomes,
    features = c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC", "percent.mt"),
    ncol = 5,
    pt.size = 0.005
  )
  png(paste("figures/ATAC/QC/Vlnplot_", sample_name, "_afterfiltering_20250604_GBMQC.png", sep = ""), res = 300, height = 1400, width = 3000)
  print(v2)
  dev.off()
  
  # RNA analysis
  DefaultAssay(multiomes) <- "RNA"
  multiomes <- NormalizeData(object = multiomes, normalization.method = "LogNormalize", scale.factor = 10000)
  multiomes <- FindVariableFeatures(object = multiomes, selection.method = 'vst', nfeatures = 2000)
  
  # ATAC analysis
  DefaultAssay(multiomes) <- "ATAC"
  multiomes <- RunTFIDF(multiomes)
  multiomes <- FindTopFeatures(multiomes, min.cutoff = 20)
  multiomes <- RunSVD(multiomes)
  
  # Update multiomes_list
  multiomes$cell.id <- colnames(multiomes)
  multiomes_list[[sample_name]] <- multiomes
}

########## 5. Merge objects ##########
# Now that the objects each contain an assay with the same set of features, we can use the standard merge function to merge the objects. This will also merge all the fragment objects so that we retain the fragment information for each cell in the final merged object.
# add information to identify dataset of origin
sample_list = list()
for (sample_id in names(multiomes_list)){
  temp <- multiomes_list[[sample_id]]
  temp$dataset <- sample_id
  sample_list[[sample_id]] <- temp
}

# Merge all datasets
combined <- merge(
  x = sample_list[[1]],
  y = sample_list[2:length(sample_list)],
  add.cell.ids = names(sample_list)
)

########## 6. Integrate scATAC (harmony) ##########
# scATAC analysis for merged Seurat object
DefaultAssay(combined) <- "ATAC"
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)

# Run Harmony for scATAC integration
multiomes_integrated <- RunHarmony(
  object = combined,
  group.by.vars = 'dataset',
  reduction = 'lsi',
  assay.use = 'ATAC',
  reduction.save = "harmony_lsi",
  max.iter.harmony = 5e3,
  project.dim = FALSE
)

# re-compute the UMAP using corrected LSI embeddings for scATAC
multiomes_integrated <- RunUMAP(multiomes_integrated, dims = 1:30, reduction = "harmony_lsi", reduction.name = "umap.atac", reduction.key = "atacUMAP_")
p <- DimPlot(multiomes_integrated, group.by = "dataset", reduction = "umap.atac", pt.size = 0.2)
png("figures/ATAC/UMAP/umap_scATAC_harmonyintegrated_20250604_GBMQC.png", res = 300, height = 1600, width = 1600)
p + ggtitle("Harmony Integrated")
dev.off()

########## 7. Integrate scRNA (harmony) ##########
# scRNASeq analysis
DefaultAssay(multiomes_integrated) <- "RNA"
multiomes_integrated <- SCTransform(multiomes_integrated, vars.to.regress = c("percent.mt"), verbose = FALSE)
multiomes_integrated <- RunPCA(multiomes_integrated, npcs = 100)

# Run Harmony for scRNA integration
multiomes_integrated <- RunHarmony(
  object = multiomes_integrated,
  group.by.vars = 'dataset',
  reduction = 'pca',
  assay.use = 'SCT',
  reduction.save = "harmony_pca",
  max.iter.harmony = 1e3,
  project.dim = FALSE
)

# re-compute the UMAP using corrected PCA embeddings for scRNA
multiomes_integrated <- RunUMAP(multiomes_integrated, dims = 1:30, reduction = "harmony_pca", reduction.name = "umap.rna", reduction.key = "rnaUMAP_")
p <- DimPlot(multiomes_integrated, group.by = "dataset", reduction = "umap.rna", pt.size = 0.2)
png("figures/RNAseq/UMAP/umap_scRNASeq_harmonyintegrated_20250604_GBMQC.png", res = 300, height = 1600, width = 1600)
p + ggtitle("Harmony Integrated")
dev.off()

########## 8. WNN for integrating scATAC and scRNA embedding ##########
# WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. We use this graph for UMAP visualization and clustering
multiomes_integrated <- FindMultiModalNeighbors(multiomes_integrated, reduction.list = list("harmony_pca", "harmony_lsi"), dims.list = list(1:50, 1:50))
multiomes_integrated <- RunUMAP(multiomes_integrated, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
multiomes_integrated <- FindClusters(multiomes_integrated, graph.name = "wsnn", algorithm = 3, verbose = FALSE, res = 0.8)

# Visualize clustering based on gene expression, ATAC-seq, or WNN analysis
p1 <- DimPlot(multiomes_integrated, reduction = "umap.rna", group.by = "dataset", label.size = 5.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(multiomes_integrated, reduction = "umap.atac", group.by = "dataset", label.size = 5.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(multiomes_integrated, reduction = "wnn.umap", group.by = "dataset", label.size = 5.5, repel = TRUE) + ggtitle("WNN")
png(paste("figures/RNAseq/UMAP/Umap_WNNIntegrated_original_groupbysamples_20250604_GBMQC.png", sep = ""), res = 300, height = 1600, width = 4800)
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

########## 9. Find cluster markers ##########
markers <- FindAllMarkers(multiomes_integrated, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(markers, "./results/celltype_annotation/markers_20250604_GBM.csv")

png("figures/RNAseq/Umap_Integrated_original_groupbyclusters_20250604_GBM.png", res = 300, height = 2000, width = 2200)
DimPlot(multiomes_integrated, reduction = "umap.rna", group.by = "wsnn_res.0.8", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("RNA")
dev.off()
p1 <- DimPlot(multiomes_integrated, reduction = "umap.rna", group.by = "dataset") + ggtitle("RNA")
p2 <- DimPlot(multiomes_integrated, reduction = "umap.atac", group.by = "dataset") + ggtitle("ATAC")
p3 <- DimPlot(multiomes_integrated, reduction = "wnn.umap", group.by = "dataset") + ggtitle("WNN")
png("figures/WSNN/Umap_WNNIntegrated_original_groupbydataset_20250604_GBM.png", res = 300, height = 2200, width = 6000)
ggpubr::ggarrange(p1, p2, p3, nrow = 1, common.legend = T, legend = "bottom")
dev.off()
p1 <- DimPlot(multiomes_integrated, reduction = "umap.rna", group.by = "wsnn_res.0.8", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(multiomes_integrated, reduction = "umap.atac", group.by = "wsnn_res.0.8", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(multiomes_integrated, reduction = "wnn.umap", group.by = "wsnn_res.0.8", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("WNN")
png("figures/WSNN/Umap_WNNIntegrated_original_groupbycluster_20250604_GBM.png", res = 300, height = 2500, width = 6000)
ggpubr::ggarrange(p1, p2, p3, nrow = 1, common.legend = T, legend = "bottom")
dev.off()

saveRDS(multiomes_integrated, "./data/RNA_ATAC/multiomes_integrated_20250604_GBM_no_annotation.rds")

########## 9. Cell type annotation ##########
multiomes_integrated <- readRDS("./data/RNA_ATAC/multiomes_integrated_20250604_GBM_no_annotation.rds")

# cell type markers
markers <- c(
  "EGFR", "PTPRZ1", # GBM
  "SYT1", "SNAP25", "NRGN", # All neurons
  "SATB2", "SLC17A6", "SLC17A7", # EN
  "GAD1", "GAD2", # IN
  "VIP", "SST", "PVALB", "RELN", # IN subtypes
  "SLC1A2", "SLC1A3", "AQP4", # Astrocyte
  "APBB1IP", # Microglia
  "SOX6", "PDGFRA", # OPC
  "MOG", "PLP1", # Oligodendrocyte
  "COL1A1",
  "PECAM1"
)

DotPlot(multiomes_integrated, features = markers, col.min = 0, 
        scale.by = 'size', group.by = "seurat_clusters") + 
  theme(panel.grid = element_blank(), panel.border = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = 3)) +
  labs(x = "Gene", y = "Cluster")

cell_type_annotation <- c(
  "0" = "Oligodendorcyte",
  "1" = "Oligodendorcyte",
  "2" = "Tumor",
  "3" = "Oligodendorcyte",
  "4" = "Tumor",
  "5" = "Oligodendorcyte",
  "6" = "Tumor",
  "7" = "Microglia",
  "8" = "Oligodendorcyte",
  "9" = "Astrocyte",
  "10" = "Microglia",
  "11" = "Tumor",
  "12" = "Oligodendorcyte",
  "13" = "Tumor",
  "14" = "Astrocyte",
  "15" = "Unknown",
  "16" = "Oligodendorcyte",
  "17" = "Endothelial",
  "18" = "Neuron"
)

multiomes_integrated <- RenameIdents(multiomes_integrated, cell_type_annotation)
multiomes_integrated$celltype <- Idents(multiomes_integrated)
table(Idents(multiomes_integrated), multiomes_integrated$seurat_clusters)

write.table(multiomes_integrated[[c("dataset", "cell.id", "celltype")]], 
            "./results/celltype_annotation/annotation_table_GBM_20250604.txt", sep = "\t",
            quote = F, col.names = F, row.names = F)
saveRDS(multiomes_integrated, "./data/RNA_ATAC/multiomes_integrated_20250604_GBM_annotated.rds")

multiomes_integrated <- readRDS("./data/RNA_ATAC/multiomes_integrated_20250604_GBM_annotated.rds")

tumor <- multiomes_integrated[,multiomes_integrated$celltype == "Tumor"]

tumor <- RunPCA(object = tumor)
tumor <- FindNeighbors(object = tumor, dims = 1:30)
tumor <- RunUMAP(object = tumor, dims = 1:30, reduction.name = "umap.tumor.rna")

# meta-modula score
DefaultAssay(tumor) <- "RNA"

all.genes <- row.names(tumor)
ScaleData(tumor, assay = "RNA", do.scale = F, features = all.genes) -> tumor

Idents(tumor) <- "all"
gene_average_expression <- AverageExpression(tumor)$RNA
gene_average_expression <- data.frame(gene_average_expression)

gene_average_expression <- gene_average_expression[gene_average_expression$all > quantile(gene_average_expression$all, 0.1), ,drop=F]
gene_expression_bin <- gene_average_expression %>%
  mutate(expression_bin = ntile(all, 30))

detected_gene <- rownames(gene_expression_bin)

subtype_marker <- read.csv("./data/others/GBM_subtype_marker.csv")[,1:6]

NPC1_marker <- subtype_marker[["NPC1"]]
NPC1_marker <- NPC1_marker[NPC1_marker != ""]
NPC2_marker <- subtype_marker[["NPC2"]]
NPC2_marker <- NPC2_marker[NPC2_marker != ""]
NPC_marker <- c(NPC1_marker, NPC2_marker)

MES1_marker <- subtype_marker[["MES1"]]
MES1_marker <- MES1_marker[MES1_marker != ""]
MES2_marker <- subtype_marker[["MES2"]]
MES2_marker <- MES2_marker[MES2_marker != ""]
MES_marker <- c(MES1_marker, MES2_marker)

AC_marker <- subtype_marker[["AC"]]
AC_marker <- AC_marker[AC_marker != ""]

OPC_marker <- subtype_marker[["OPC"]]
OPC_marker <- OPC_marker[OPC_marker != ""]

marker_list <- list("NPC" = c(NPC1_marker, NPC2_marker), "OPC" = OPC_marker, 
                    "MES" = c(MES1_marker, MES2_marker), "AC" = AC_marker)

subtype_score <- list()
for(subtype in names(marker_list)){
  marker <- marker_list[[subtype]]
  marker <- marker[marker != ""]
  marker <- marker[marker %in% detected_gene]
  
  control <- c()
  for (gene in marker){
    bin <- gene_expression_bin[gene, "expression_bin"]
    random_selected_gene <- sample(detected_gene[gene_expression_bin$expression_bin == bin], 100, replace = F)
    control <- c(control, random_selected_gene)
  }
  
  subtype_score[[subtype]] <- apply(tumor@assays$RNA@scale.data[marker,], 2, mean) - apply(tumor@assays$RNA@scale.data[control,], 2, mean)
}
subtype_score <- data.frame(subtype_score)
subtype_score[colnames(tumor),] -> subtype_score

AddMetaData(tumor, subtype_score) -> tumor
DotPlot(tumor, features = c("MES1", "MES2", "AC", "OPC", "NPC1", "NPC2"), group.by = "seurat_clusters")
DotPlot(tumor, features = c("MES", "AC", "OPC", "NPC"), group.by = "seurat_clusters")

p1 <- FeaturePlot(tumor, reduction = "umap.tumor.rna", features = "MES")
p3 <- FeaturePlot(tumor, reduction = "umap.tumor.rna", features = "AC")
p4 <- FeaturePlot(tumor, reduction = "umap.tumor.rna", features = "OPC")
p5 <- FeaturePlot(tumor, reduction = "umap.tumor.rna", features = "NPC")
ggpubr::ggarrange(p1, p3, p4, p5, ncol = 4, nrow = 1)

tumor$tumor_seurat_clusters <- factor(tumor$tumor_seurat_clusters, levels = c(2, 4, 6, 11, 13),
                                      levels =  c("NPC-/OPC-like #1", "MES-like", "NPC-/OPC-like #2", "OPC-like", "AC-like"))

saveRDS(tumor, "./data/RNA_ATAC/multiomes_integrated_20250604_GBM_tumor_only_annotated.rds")
