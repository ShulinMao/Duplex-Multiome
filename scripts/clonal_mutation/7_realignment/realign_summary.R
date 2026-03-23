.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(dplyr)
library(stringr)

original_clonal_mutation <- read.csv("./results/clonal_mutation/clonal_mutation_list_ASD_filtered_20241028.csv")
output_dir <- "/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/results/clonal_mutation/realigned/ASD/"
bulk_WGS_mutations <- readRDS("./data/others/bulk_WGS_germline_mutations/bulk_WGS_mutations.rds")

# ------------------------------------------------------------------------------
# Realign with bwa mem
# INDEL panelty 4 (open), 4 (extend)
realign_clonal_mutation <- c()
for (i in seq(dim(original_clonal_mutation)[1])) {
  case_id <- original_clonal_mutation$case_id[i]
  variant <- original_clonal_mutation$key[i]
  cell <- original_clonal_mutation$cell[i]
  wd <- paste0("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/", case_id, "/single_cell/")
  RF <- str_split(original_clonal_mutation$value[i], ":", simplify = T)[1, 5]
  
  file_name <- paste0(output_dir, case_id, ".", cell, ".", variant, ".realigned.vcf")
  try_read_vcf <- tryCatch({read.table(file_name, sep = "\t", colClasses = c("character"))},
                           error = function(e){return(NULL)})
  
  if(!is.null(try_read_vcf)){
    try_read_vcf <- try_read_vcf[, c(1:2, 4:5)]
    colnames(try_read_vcf) <- c("chr", "pos", "ref", "alt")
    try_read_vcf <- cbind(try_read_vcf, case_id, cell)
  }

  realign_clonal_mutation <- rbind(realign_clonal_mutation, try_read_vcf)
}

realign_clonal_mutation$pos <- as.integer(realign_clonal_mutation$pos)
realign_clonal_mutation$realign <- T
left_join(original_clonal_mutation, realign_clonal_mutation) -> original_realign_clonal_mutation

realign_clonal_mutation_summary <- original_realign_clonal_mutation %>%
  group_by(case_id, key) %>% summarise("recalled" = sum(realign, na.rm = T), "count" = n())

realign_clonal_mutation_summary_pass <- 
  realign_clonal_mutation_summary[realign_clonal_mutation_summary$recalled >= realign_clonal_mutation_summary$count*0.5,]

realign_clonal_mutation_summary_filtered <- 
  realign_clonal_mutation_summary[realign_clonal_mutation_summary$recalled < realign_clonal_mutation_summary$count*0.5,]

# ------------------------------------------------------------------------------
# GATK haplotypecaller
haplotypecaller_clonal_mutation <- c()

for (i in seq(dim(original_clonal_mutation)[1])){
  case_id <- original_clonal_mutation$case_id[i]
  variant <- original_clonal_mutation$key[i]
  cell <- original_clonal_mutation$cell[i]
  wd <- paste0("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/", case_id, "/single_cell/")
  RF <- str_split(original_clonal_mutation$value[i], ":", simplify = T)[1, 5]
  
  file_name <- paste0(output_dir, case_id, "_", cell, "_", variant, ".vcf")
  try_read_vcf <- tryCatch({read.table(file_name, sep = "\t", colClasses = c("character"))},
                           error = function(e){return(NULL)})
  
  if(!is.null(try_read_vcf)){
    # remove variants if there are multiple variants (except germline variants) surrounding it
    try_read_vcf_variants <- paste(try_read_vcf$V1, try_read_vcf$V2, try_read_vcf$V4, try_read_vcf$V5, sep = ":")
    clustered_variants <- sum(!(setdiff(try_read_vcf_variants, str_replace_all(variant, "-", ":")) %in% bulk_WGS_mutations[[case_id]]))
    if(clustered_variants > 0){
      try_read_vcf <- NULL
    } else{
      try_read_vcf <- try_read_vcf[, c(1:2, 4:5)]
      colnames(try_read_vcf) <- c("chr", "pos", "ref", "alt")
      try_read_vcf <- cbind(try_read_vcf, case_id, cell)  
    }
  }
  
  haplotypecaller_clonal_mutation <- rbind(haplotypecaller_clonal_mutation, try_read_vcf)
}
haplotypecaller_clonal_mutation$pos <- as.integer(haplotypecaller_clonal_mutation$pos)
haplotypecaller_clonal_mutation$realign <- T
left_join(original_clonal_mutation, haplotypecaller_clonal_mutation) -> original_haplotypecaller_clonal_mutation

haplotypecaller_clonal_mutation_summary <- original_haplotypecaller_clonal_mutation %>%
  group_by(case_id, key) %>% summarise("recalled" = sum(realign, na.rm = T), "count" = n())

haplotypecaller_clonal_mutation_summary_pass <- 
  haplotypecaller_clonal_mutation_summary[haplotypecaller_clonal_mutation_summary$recalled >= haplotypecaller_clonal_mutation_summary$count*0.5,]

haplotypecaller_clonal_mutation_summary_filtered <- 
  haplotypecaller_clonal_mutation_summary[haplotypecaller_clonal_mutation_summary$recalled < haplotypecaller_clonal_mutation_summary$count*0.5,]

clonal_mutation_summary_pass <- haplotypecaller_clonal_mutation_summary_pass[haplotypecaller_clonal_mutation_summary_pass$key %in% realign_clonal_mutation_summary_pass$key, 1:2]
clonal_mutation_summary_pass_cell <- original_clonal_mutation[original_clonal_mutation$key %in% clonal_mutation_summary_pass$key,]

write.csv(clonal_mutation_summary_pass_cell, "./results/clonal_mutation/clonal_mutation_list_ASD_realigned_20241028.csv")

