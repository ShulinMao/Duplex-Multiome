# Libraries
.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(stringr)
library(tidyr)
library(dplyr)

germline_VAF <- 0.1

# Healthy postmortem samples
case_id_list <- c("Brain1901106_deeper", "UMB1278_deeper", "UMB1864_deeper", 
                  "UMB5823_deeper", "UMB4638", "UMB1465", "UMB5451", "UMB5657")

clone <- c()
for (case_id in case_id_list){
  print(case_id)
    
  # read in all clonal mutations
  all_clonal_mutation <- read.table(paste0("./data/", case_id, "/single_cell/call_variants/filter_germline_mutation/a2s0_AF", germline_VAF, ".vcf"), header = F, sep = "\t")
  colnames(all_clonal_mutation) <- c("chr", "pos", "cell", "ref", "alt", "qual", "filter", "info", "fromat", "value")
  all_clonal_mutation[,"key"] <- unite(all_clonal_mutation, col = "key", c("chr", "pos", "ref", "alt"), sep = ":")[,"key"]
  all_clonal_mutation[,"case_id"] <- case_id
  all_clonal_mutation <- all_clonal_mutation[!(duplicated(str_c(all_clonal_mutation$cell, all_clonal_mutation$key, all_clonal_mutation$case_id))),]
  
  all_clonal_mutation_summary <- unite(all_clonal_mutation, col = "key", c("chr", "pos", "ref", "alt"), sep = ":") %>% group_by(key) %>% summarise("occurance" = n())
  clonal_mutation_list <- all_clonal_mutation_summary[all_clonal_mutation_summary$occurance > 2,]
  
  for (large_clonal_snv in clonal_mutation_list$key){
    clone <- rbind(clone, all_clonal_mutation[all_clonal_mutation$key == large_clonal_snv, ])
  }
}

cell_type <- read.table("./results/celltype_annotation/annotation_table_20240512.txt",
                        sep = "\t")
colnames(cell_type) <- c("case_id", "cell", "cell_type")
clone <- left_join(clone, cell_type, by = c("case_id", "cell"))

write.csv(clone, "./results/clonal_mutation/clonal_mutation_list_20240812.csv", row.names = F, quote = F)


# COLO829BLT50
germline_VAF <- 0.4
case_id_list <- c("COLO829BLT50_rep1", "COLO829BLT50_rep2")

clone <- c()
for (case_id in case_id_list){
  # read in all clonal mutations
  all_clonal_mutation <- read.table(paste0("./data/", case_id, "/single_cell/call_variants/filter_germline_mutation/a2s0_AF", germline_VAF, ".vcf"), header = F, sep = "\t")
  colnames(all_clonal_mutation) <- c("chr", "pos", "cell", "ref", "alt", "qual", "filter", "info", "fromat", "value")
  all_clonal_mutation[,"key"] <- unite(all_clonal_mutation, col = "key", c("chr", "pos", "ref", "alt"), sep = ":")[,"key"]
  all_clonal_mutation[,"case_id"] <- case_id
  all_clonal_mutation <- all_clonal_mutation[!(duplicated(str_c(all_clonal_mutation$cell, all_clonal_mutation$key, all_clonal_mutation$case_id))),]
  
  all_clonal_mutation_summary <- unite(all_clonal_mutation, col = "key", c("chr", "pos", "ref", "alt"), sep = ":") %>% group_by(key) %>% summarise("occurance" = n())
  clonal_mutation_list <- all_clonal_mutation_summary[all_clonal_mutation_summary$occurance > 2,]
  
  for (large_clonal_snv in clonal_mutation_list$key){
    clone <- rbind(clone, all_clonal_mutation[all_clonal_mutation$key == large_clonal_snv, ])
  }
}

cell_type <- read.table("./results/celltype_annotation/annotation_COLO829BLT50_20240523.txt",
                        sep = "\t")
colnames(cell_type) <- c("case_id", "cell", "cell_type")
clone <- left_join(clone, cell_type, by = c("case_id", "cell"))

write.csv(clone, "./results/clonal_mutation/clonal_mutation_list_COLO829BLT50_20240523.csv", row.names = F, quote = F)



# COLO829BLT50 BL bulk
case_id_list <- c("COLO829BLT50_rep1", "COLO829BLT50_rep2")

clone <- c()
for (case_id in case_id_list){
  # read in all clonal mutations
  all_clonal_mutation <- read.table(paste0("./data/", case_id, "/single_cell/call_variants/filter_germline_mutation/a2s0_BL_bulk.vcf"), header = F, sep = "\t")
  colnames(all_clonal_mutation) <- c("chr", "pos", "cell", "ref", "alt", "qual", "filter", "info", "fromat", "value")
  all_clonal_mutation[,"key"] <- unite(all_clonal_mutation, col = "key", c("chr", "pos", "ref", "alt"), sep = ":")[,"key"]
  all_clonal_mutation[,"case_id"] <- case_id
  all_clonal_mutation <- all_clonal_mutation[!(duplicated(str_c(all_clonal_mutation$cell, all_clonal_mutation$key, all_clonal_mutation$case_id))),]
  
  all_clonal_mutation_summary <- unite(all_clonal_mutation, col = "key", c("chr", "pos", "ref", "alt"), sep = ":") %>% group_by(key) %>% summarise("occurance" = n())
  clonal_mutation_list <- all_clonal_mutation_summary[all_clonal_mutation_summary$occurance > 2,]
  
  for (large_clonal_snv in clonal_mutation_list$key){
    clone <- rbind(clone, all_clonal_mutation[all_clonal_mutation$key == large_clonal_snv, ])
  }
}

cell_type <- read.table("./results/celltype_annotation/annotation_COLO829BLT50_20240523.txt",
                        sep = "\t")
colnames(cell_type) <- c("case_id", "cell", "cell_type")
clone <- left_join(clone, cell_type, by = c("case_id", "cell"))

write.csv(clone, "./results/clonal_mutation/clonal_mutation_list_COLO829BLT50_BL_bulk_20240604.csv", row.names = F, quote = F)

# COLO829BLT50 tumor bulk
case_id_list <- c("COLO829BLT50_rep1", "COLO829BLT50_rep2")

clone <- c()
for (case_id in case_id_list){
  # read in all clonal mutations
  all_clonal_mutation <- read.table(paste0("./data/", case_id, "/single_cell/call_variants/filter_germline_mutation/a2s0_tumor_bulk.vcf"), header = F, sep = "\t")
  colnames(all_clonal_mutation) <- c("chr", "pos", "cell", "ref", "alt", "qual", "filter", "info", "fromat", "value")
  all_clonal_mutation[,"key"] <- unite(all_clonal_mutation, col = "key", c("chr", "pos", "ref", "alt"), sep = ":")[,"key"]
  all_clonal_mutation[,"case_id"] <- case_id
  all_clonal_mutation <- all_clonal_mutation[!(duplicated(str_c(all_clonal_mutation$cell, all_clonal_mutation$key, all_clonal_mutation$case_id))),]
  
  all_clonal_mutation_summary <- unite(all_clonal_mutation, col = "key", c("chr", "pos", "ref", "alt"), sep = ":") %>% group_by(key) %>% summarise("occurance" = n())
  clonal_mutation_list <- all_clonal_mutation_summary[all_clonal_mutation_summary$occurance > 2,]
  
  for (large_clonal_snv in clonal_mutation_list$key){
    clone <- rbind(clone, all_clonal_mutation[all_clonal_mutation$key == large_clonal_snv, ])
  }
}

cell_type <- read.table("./results/celltype_annotation/annotation_COLO829BLT50_20240523.txt",
                        sep = "\t")
colnames(cell_type) <- c("case_id", "cell", "cell_type")
clone <- left_join(clone, cell_type, by = c("case_id", "cell"))

write.csv(clone, "./results/clonal_mutation/clonal_mutation_list_COLO829BLT50_tumor_bulk_20240604.csv", row.names = F, quote = F)

# ASD
case_id_list <- c("ABNCR6D", "AN06365")

clone <- c()
for (case_id in case_id_list){
  print(case_id)
  
  # read in all clonal mutations
  all_clonal_mutation <- read.table(paste0("./data/", case_id, "/single_cell/call_variants/filter_germline_mutation/a2s0_AF", germline_VAF, ".vcf"), header = F, sep = "\t")
  colnames(all_clonal_mutation) <- c("chr", "pos", "cell", "ref", "alt", "qual", "filter", "info", "fromat", "value")
  all_clonal_mutation[,"key"] <- unite(all_clonal_mutation, col = "key", c("chr", "pos", "ref", "alt"), sep = ":")[,"key"]
  all_clonal_mutation[,"case_id"] <- case_id
  all_clonal_mutation <- all_clonal_mutation[!(duplicated(str_c(all_clonal_mutation$cell, all_clonal_mutation$key, all_clonal_mutation$case_id))),]
  
  all_clonal_mutation_summary <- unite(all_clonal_mutation, col = "key", c("chr", "pos", "ref", "alt"), sep = ":") %>% group_by(key) %>% summarise("occurance" = n())
  clonal_mutation_list <- all_clonal_mutation_summary[all_clonal_mutation_summary$occurance > 2,]
  
  for (large_clonal_snv in clonal_mutation_list$key){
    clone <- rbind(clone, all_clonal_mutation[all_clonal_mutation$key == large_clonal_snv, ])
  }
}

cell_type <- read.table("./results/celltype_annotation/annotation_table_20240813_ASD.txt",
                        sep = "\t")
colnames(cell_type) <- c("case_id", "cell", "cell_type")
clone <- left_join(clone, cell_type, by = c("case_id", "cell"))

write.csv(clone, "./results/clonal_mutation/clonal_mutation_list_ASD_20241028.csv", row.names = F, quote = F)

