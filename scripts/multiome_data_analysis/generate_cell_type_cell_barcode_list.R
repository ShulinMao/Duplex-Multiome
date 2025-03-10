# Libraries
.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(stringr)

# ------------------------------------------------------------------------------
# Healthy human samples
cell_type <- read.table("./results/celltype_annotation/annotation_table_20240512.txt",
                        sep = "\t")
colnames(cell_type) <- c("sample_id", "cell_id", "cell_type")
cell_type <- cell_type[!str_detect(cell_type$cell_type, "doublet"),]
table(cell_type$cell_type)

sample_list <- unique(cell_type$sample_id)

for (sample_id in sample_list){
  print(sample_id)
  print(table(cell_type$cell_type[cell_type$sample_id == sample_id]))
  
  if(!file.exists(paste0("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/", sample_id, "/cell_type"))){
    system(paste0("mkdir /lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/", sample_id, "/cell_type")) 
  }
  
  # only keep cell type with >200 cells for burden estimation
  EN <- cell_type[cell_type$sample_id == sample_id & str_detect(cell_type$cell_type, "EN"),]
  if(dim(EN)[1] > 200){
    write.table(EN, paste0("./data/", sample_id, "/cell_type/EN.txt"), sep = "\t", quote = F, col.names = F, row.names = F)
  }
  
  IN <- cell_type[cell_type$sample_id == sample_id & str_detect(cell_type$cell_type, "IN"),]
  if(dim(IN)[1] > 200){
    write.table(IN, paste0("./data/", sample_id, "/cell_type/IN.txt"), sep = "\t", quote = F, col.names = F, row.names = F)
  }
  
  Astrocyte <- cell_type[cell_type$sample_id == sample_id & str_detect(cell_type$cell_type, "Astrocyte"),]
  if(dim(Astrocyte)[1] > 200){
    write.table(Astrocyte, paste0("./data/", sample_id, "/cell_type/Astrocyte.txt"), sep = "\t", quote = F, col.names = F, row.names = F)
  }
  
  Microglia <- cell_type[cell_type$sample_id == sample_id & str_detect(cell_type$cell_type, "Microglia"),]
  if(dim(Microglia)[1] > 200){
    write.table(Microglia, paste0("./data/", sample_id, "/cell_type/Microglia.txt"), sep = "\t", quote = F, col.names = F, row.names = F)
  }
  
  Oligodendrocyte <- cell_type[cell_type$sample_id == sample_id & str_detect(cell_type$cell_type, "Oligodendrocyte"),]
  if(dim(Oligodendrocyte)[1] > 200){
    write.table(Oligodendrocyte, paste0("./data/", sample_id, "/cell_type/Oligodendrocyte.txt"), sep = "\t", quote = F, col.names = F, row.names = F)
  }
  
  OPC <- cell_type[cell_type$sample_id == sample_id & str_detect(cell_type$cell_type, "OPC"),]
  if(dim(OPC)[1] > 100){
    write.table(OPC, paste0("./data/", sample_id, "/cell_type/OPC.txt"), sep = "\t", quote = F, col.names = F, row.names = F)
  }
}

# ------------------------------------------------------------------------------
# subtype
# EN
subtype <- read.table("./results/celltype_annotation/EN_annotation_transfer_20240715.txt",
                        sep = "\t")
colnames(subtype) <- c("sample_id", "cell_id", "transferred_labels", "subtype")
table(subtype$sample_id, subtype$subtype)

sample_list <- unique(subtype$sample_id)

for (sample_id in sample_list){
  print(sample_id)
  print(table(subtype$subtype[subtype$sample_id == sample_id]))
  
  if(!file.exists(paste0("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/", sample_id, "/cell_type/subtype"))){
    system(paste0("mkdir /lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/", sample_id, "/cell_type/subtype")) 
  }
  
  # only keep cell type with >200 cells for burden estimation
  upper_layer <- subtype[subtype$sample_id == sample_id & str_detect(subtype$subtype, "upper_layer"),]
  if(dim(upper_layer)[1] > 200){
    write.table(upper_layer, paste0("./data/", sample_id, "/cell_type/subtype/EN_upper_layer.txt"), sep = "\t", quote = F, col.names = F, row.names = F)
  }
  
  deep_layer <- subtype[subtype$sample_id == sample_id & str_detect(subtype$subtype, "deep_layer"),]
  if(dim(deep_layer)[1] > 200){
    write.table(deep_layer, paste0("./data/", sample_id, "/cell_type/subtype/EN_deep_layer.txt"), sep = "\t", quote = F, col.names = F, row.names = F)
  }
}

# IN
subtype <- read.table("./results/celltype_annotation/IN_annotation_transfer_20240715.txt",
                      sep = "\t")
colnames(subtype) <- c("sample_id", "cell_id", "transferred_labels", "subtype")
table(subtype$sample_id, subtype$subtype)

sample_list <- unique(subtype$sample_id)

for (sample_id in sample_list){
  print(sample_id)
  print(table(subtype$subtype[subtype$sample_id == sample_id]))
  
  if(!file.exists(paste0("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/", sample_id, "/cell_type/subtype"))){
    system(paste0("mkdir /lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/", sample_id, "/cell_type/subtype")) 
  }
  
  # only keep cell type with >200 cells for burden estimation
  CGE <- subtype[subtype$sample_id == sample_id & str_detect(subtype$subtype, "CGE"),]
  if(dim(CGE)[1] > 200){
    write.table(CGE, paste0("./data/", sample_id, "/cell_type/subtype/IN_CGE.txt"), sep = "\t", quote = F, col.names = F, row.names = F)
  }
  
  MGE <- subtype[subtype$sample_id == sample_id & str_detect(subtype$subtype, "MGE"),]
  if(dim(MGE)[1] > 200){
    write.table(MGE, paste0("./data/", sample_id, "/cell_type/subtype/IN_MGE.txt"), sep = "\t", quote = F, col.names = F, row.names = F)
  }
}

# ------------------------------------------------------------------------------
# COLO829BLT50
cell_type <- read.table("./results/celltype_annotation/annotation_COLO829BLT50_20240523.txt",
                        sep = "\t")
colnames(cell_type) <- c("sample_id", "cell_id", "cell_type")

sample_list <- unique(cell_type$sample_id)

for (sample_id in sample_list){
  print(sample_id)
  print(table(cell_type$cell_type[cell_type$sample_id == sample_id]))
  
  if(!file.exists(paste0("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/", sample_id, "/cell_type"))){
    system(paste0("mkdir /lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/", sample_id, "/cell_type")) 
  }
  
  B_lymphoblast <- cell_type[cell_type$sample_id == sample_id & str_detect(cell_type$cell_type, "B_lymphoblast"),]
  write.table(B_lymphoblast, paste0("./data/", sample_id, "/cell_type/B_lymphoblast.txt"), sep = "\t", quote = F, col.names = F, row.names = F)

  melanoma_fibroblast <- cell_type[cell_type$sample_id == sample_id & str_detect(cell_type$cell_type, "melanoma_fibroblast"),]
  write.table(melanoma_fibroblast, paste0("./data/", sample_id, "/cell_type/melanoma_fibroblast.txt"), sep = "\t", quote = F, col.names = F, row.names = F)
}

# ------------------------------------------------------------------------------
# Healthy human samples ATAC recovered
cell_type <- read.table("./results/celltype_annotation/annotation_table_20240512_ATAC_recovered.txt",
                        sep = "\t")
colnames(cell_type) <- c("sample_id", "cell_id", "celltype_ATAC", "celltype_wsnn")
cell_type <- na.omit(cell_type)
table(cell_type$celltype_ATAC)
table(cell_type$celltype_wsnn)

sample_list <- unique(cell_type$sample_id)

for (sample_id in sample_list){
  print(sample_id)
  
  if(!file.exists(paste0("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/", sample_id, "/cell_type/celltype_ATAC"))){
    system(paste0("mkdir /lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/", sample_id, "/cell_type/celltype_ATAC")) 
  }
  
  if(!file.exists(paste0("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/", sample_id, "/cell_type/celltype_wsnn"))){
    system(paste0("mkdir /lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/", sample_id, "/cell_type/celltype_wsnn")) 
  }
  
  celltype_list <- unique(cell_type[cell_type$sample_id == sample_id, "celltype_ATAC"])
  print(celltype_list)
  for (celltype in celltype_list){
    write.table(cell_type[cell_type$sample_id == sample_id & cell_type$celltype_ATAC == celltype, ], 
                paste0("./data/", sample_id, "/cell_type/celltype_ATAC/", celltype, ".txt"), sep = "\t", quote = F, col.names = F, row.names = F)
  }
  
  celltype_list <- unique(cell_type[cell_type$sample_id == sample_id, "celltype_wsnn"])
  print(celltype_list)
  for (celltype in celltype_list){
    write.table(cell_type[cell_type$sample_id == sample_id & cell_type$celltype_wsnn == celltype, ], 
                paste0("./data/", sample_id, "/cell_type/celltype_wsnn/", celltype, ".txt"), sep = "\t", quote = F, col.names = F, row.names = F)
  }
  
}


