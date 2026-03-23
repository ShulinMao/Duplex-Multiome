.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")
library(Seurat)

# build a seurat object
expr_matrix <- read.csv("./data/RNA_ATAC/allen_brain_map/HumanMultipleCorticalAreasSMARTseq/matrix.csv", 
                        header = T, row.names = 1)
expr_matrix <- t(expr_matrix)
metadata <- read.csv("./data/RNA_ATAC/allen_brain_map/HumanMultipleCorticalAreasSMARTseq/metadata.csv", 
                     header = T, row.names = 1)

allen_brain <- CreateSeuratObject(counts = expr_matrix, project = "Allen_brain_Map")
allen_brain@meta.data
allen_brain <- AddMetaData(object = allen_brain, metadata = metadata)
table(allen_brain$subclass_label)

# Perform data normalization and filtering
allen_brain <- NormalizeData(allen_brain)
allen_brain <- FindVariableFeatures(allen_brain, selection.method = "vst", nfeatures = 2000)
allen_brain <- ScaleData(allen_brain)

saveRDS(allen_brain, "./data/RNA_ATAC/allen_brain_map/HumanMultipleCorticalAreasSMARTseq/seurat_obj.rds")

cell_type_obj_list <- SplitObject(allen_brain, split.by = "class_label")

for (cell_type in names(cell_type_obj_list)){
  if (cell_type == ""){next}
  print(cell_type)
  
  cell_type_obj <- cell_type_obj_list[[cell_type]]
  cell_type_obj <- NormalizeData(cell_type_obj)
  cell_type_obj <- FindVariableFeatures(cell_type_obj)
  cell_type_obj <- ScaleData(cell_type_obj)
  cell_type_obj <- RunPCA(cell_type_obj)
  cell_type_obj <- FindNeighbors(cell_type_obj, dims = 1:30)
  cell_type_obj <- FindClusters(cell_type_obj)
  cell_type_obj <- RunUMAP(cell_type_obj, dims = 1:30)
  
  saveRDS(cell_type_obj, paste0("./data/RNA_ATAC/allen_brain_map/HumanMultipleCorticalAreasSMARTseq/", cell_type, ".rds"))
}

