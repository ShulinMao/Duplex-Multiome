# Libraries
.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(dplyr)
library(tidyr)
library(stringr)

# ------------------------------------------------------------------------------
# Healthy cells
case_id_list <- c("Brain1901106_deeper", "UMB1278_deeper", "UMB1864_deeper", "UMB5823_deeper", "UMB4638", "UMB1465", "UMB5451", "UMB5657")

for (a in seq(2, 6, 2)){
  s <- a%/%2
  stringency <- paste0("a", a, "s", s)
  print(stringency)
  
  context_3mer_summary <- list()
  
  for(case_id in case_id_list){
    # cell type
    file_list <- dir(paste0("./data/", case_id, "/cell_type/"))
    cell_type_list <- str_remove(file_list, "[.]txt")
    
    for (cell_type in cell_type_list){
      read.table(paste0("./data/", case_id, "/cell_type/", cell_type, ".txt"), sep = "\t") -> selected_cells
      
      observed_count <- rep(0, 32)
      # read in normalized spectrum
      for (cell_id in selected_cells$V2){
        # read in 3-mer context
        context_3mer_cell <-
          read.table(paste0("./data/", case_id, "/single_cell/burden_estimation/Dan/", stringency, "/", stringency, "_filterGermline_burden.", cell_id, ".txt"), 
                     skip = 7, nrow = 32, comment.char = "", header = T)[,-1]
        observed_count <- observed_count + context_3mer_cell$ObservedCount
      }
      
      context_3mer_summary[[paste(case_id, cell_type, sep="_")]] <- observed_count
    }
  }
  context_3mer_summary <- data.frame(context_3mer_summary)
  saveRDS(context_3mer_summary, paste0("./results/context_3mer_ratio/", "a", a, "s", s, "_sample_celltype.rds"))
  context_3mer_summary <- readRDS(paste0("./results/context_3mer_ratio/", "a", a, "s", s, "_sample_celltype.rds"))
  
  context_3mer_celltype <- list()
  for (cell_type in cell_type_list){
    context_3mer <- context_3mer_summary[,str_detect(colnames(context_3mer_summary), paste0("_", cell_type, "$"))]
    context_3mer_celltype[[cell_type]] <- rowSums(context_3mer)
  }
  context_3mer_celltype <- data.frame(context_3mer_celltype)
  context_3mer_celltype <- apply(context_3mer_celltype, 2, function(x){x/sum(x)})
  
  # use this from the last cell since it's a constant for hg38
  context_3mer_genome <- context_3mer_cell$GenomeCount/sum(context_3mer_cell$GenomeCount)
  context_3mer_ratio <- context_3mer_celltype / context_3mer_genome
  rownames(context_3mer_ratio) <- context_3mer_cell$Context
  
  saveRDS(context_3mer_ratio, paste0("./results/context_3mer_ratio/", "a", a, "s", s, "_celltype.rds"))
}


# ------------------------------------------------------------------------------
# COLO829BLT50
case_id_list <- c("COLO829BLT50_rep1", "COLO829BLT50_rep2")

a <- 2
s <- a%/%2
stringency <- paste0("a", a, "s", s)

context_3mer_summary <- list()

for(case_id in case_id_list){
  # cell type
  file_list <- dir(paste0("./data/", case_id, "/cell_type/"))
  cell_type_list <- str_remove(file_list, "[.]txt")
  
  for (cell_type in cell_type_list){
    read.table(paste0("./data/", case_id, "/cell_type/", cell_type, ".txt"), sep = "\t") -> selected_cells
    
    observed_count <- rep(0, 32)
    # read in normalized spectrum
    for (cell_id in selected_cells$V2){
      # read in 3-mer context
      context_3mer_cell <-
        read.table(paste0("./data/", case_id, "/single_cell/burden_estimation/Dan/", stringency, "/", stringency, "_filterGermline_burden.", cell_id, ".txt"), 
                   skip = 7, nrow = 32, comment.char = "", header = T)[,-1]
      observed_count <- observed_count + context_3mer_cell$ObservedCount
    }
    
    context_3mer_summary[[paste(case_id, cell_type, sep="_")]] <- observed_count
  }
}
context_3mer_summary <- data.frame(context_3mer_summary)
saveRDS(context_3mer_summary, paste0("./results/context_3mer_ratio/COLO829BLT50_", "a", a, "s", s, "_sample_celltype.rds"))
context_3mer_summary <- readRDS(paste0("./results/context_3mer_ratio/COLO829BLT50_", "a", a, "s", s, "_sample_celltype.rds"))

context_3mer_celltype <- list()
for (cell_type in cell_type_list){
  context_3mer <- context_3mer_summary[,str_detect(colnames(context_3mer_summary), paste0("_", cell_type, "$"))]
  context_3mer_celltype[[cell_type]] <- rowSums(context_3mer)
}
context_3mer_celltype <- data.frame(context_3mer_celltype)
context_3mer_celltype <- apply(context_3mer_celltype, 2, function(x){x/sum(x)})

# use this from the last cell since it's a constant for hg38
context_3mer_genome <- context_3mer_cell$GenomeCount/sum(context_3mer_cell$GenomeCount)
context_3mer_ratio <- context_3mer_celltype / context_3mer_genome
rownames(context_3mer_ratio) <- context_3mer_cell$Context

saveRDS(context_3mer_ratio, paste0("./results/context_3mer_ratio/COLO829BLT50_", "a", a, "s", s, "_celltype.rds"))
