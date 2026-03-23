.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")
library(stringr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(MASS)
library(truncnorm)
set.seed(42)

case_id_list <- c("ABNCR6D")
  #c("Brain1901106_deeper", "UMB1278_deeper", "UMB1864_deeper", "UMB5823_deeper", "UMB4638", "UMB1465", "UMB5451", "UMB5657")
  #c("COLO829BLT50_rep1", "COLO829BLT50_rep2")


for(case_id in case_id_list){
  print(case_id)
  
  # cell type
  file_list <- dir(paste0("./data/", case_id, "/cell_type/"))
  cell_type_list <- str_remove(file_list, "[.]txt")
  
  error_rate_incorrect_rf_cell_type <- list()
  for(cell_type in cell_type_list){
    print(cell_type)
    read.table(paste0("./data/", case_id, "/cell_type/", cell_type, ".txt"), sep = "\t") -> selected_cells
    
    # ------------------------------------------------------------------------------
    # Estimate the # of errors from incorrect read families
    ssnv_rate_strand12 <- c()
    ssnv_rate_strand2 <- c()
    for(a in seq(2, 10, 2)){
      print(paste0("a", a))
      
      # strand1 + strand2(strand1bc) + strand2
      ssnv_summary <- read.table(paste0("./data/", case_id, "/single_cell/call_variants/autosomal_variants/a", a, "s0_autosome_variants.txt"))
      ssnv_summary[, "read_family"] <- str_split(ssnv_summary$V10, pattern = ":", simplify = T)[,5]
      ssnv_summary[, "cell_id"] <- str_extract(str_split(ssnv_summary$V11, pattern = "[.]", simplify = T)[,1], "[ATCG]*-1$")
      
      ssnv_count <- c()
      callable_region_read_family_bc2_size <- c()
      for (cell_id in selected_cells$V2){
        # read in callable region 
        callable_region_file <- paste0("./data/", case_id, "/single_cell/call_variants/autosomal_callable_region/a", a, "s0/", cell_id, ".families.analysis.calledSites.bed")
        
        # if no callable region, skip the cell
        if (file.info(callable_region_file)$size == 0){
          callable_region_read_family_bc2_size[cell_id] <- 0
          ssnv_count[cell_id] <- 0
          next
        }
        
        callable_region <- read.table(callable_region_file)
        colnames(callable_region)[4] <- "read_family"
        
        # read in read family information
        read_family_info_bc2 <- read.table(paste0("./data/", case_id, "/single_cell/read_family_strand_barcode/all/", cell_id, "_read_family_strand_barcode.txt"))
        read_family_info_bc2[, "read_family"] <- str_split(read_family_info_bc2$V1, pattern = ":", simplify = T)[,3]
        read_family_info_bc2$read_family <- as.integer(read_family_info_bc2$read_family)
        read_family_info_bc2 <- read_family_info_bc2[read_family_info_bc2$V5 > 0 & read_family_info_bc2$V7 > 0,]
        
        callable_region_read_family_bc2 <- inner_join(callable_region, read_family_info_bc2, by = "read_family")
        
        if(dim(callable_region_read_family_bc2)[1]){
          callable_region_read_family_bc2_size[cell_id] <- 
            sum(ifelse(callable_region_read_family_bc2$V3.x - callable_region_read_family_bc2$V2.x > 22,
                       callable_region_read_family_bc2$V3.x - callable_region_read_family_bc2$V2.x, 0))
        }else{
          callable_region_read_family_bc2_size[cell_id] <- 0
        }
        
        # count the number of ssnv from the read families with only strand2 bc
        if(sum(ssnv_summary$cell_id == cell_id)){
          ssnv_count[cell_id] <- sum(ssnv_summary[ssnv_summary$cell_id == cell_id, "read_family"] %in% read_family_info_bc2$read_family)
        }else{
          ssnv_count[cell_id] <- 0
        }
        
        if(length(ssnv_count) %% 100 == 0){print(length(ssnv_count))}
      }
      
      ssnv_rate_strand12[paste0("a", a)] <- sum(ssnv_count)/sum(callable_region_read_family_bc2_size)
      
      # strand2(strand1bc) + strand2
      ssnv_summary <- read.table(paste0("./data/", case_id, "/single_cell/call_variants/autosomal_variants/a", a, "s0_autosome_variants.txt"))
      ssnv_summary[, "read_family"] <- str_split(ssnv_summary$V10, pattern = ":", simplify = T)[,5]
      ssnv_summary[, "cell_id"] <- str_extract(str_split(ssnv_summary$V11, pattern = "[.]", simplify = T)[,1], "[ATCG]*-1$")
      
      ssnv_count <- c()
      callable_region_read_family_bc2_size <- c()
      for (cell_id in selected_cells$V2){
        # read in callable region 
        callable_region_file <- paste0("./data/", case_id, "/single_cell/call_variants/autosomal_callable_region/a", a, "s0/", cell_id, ".families.analysis.calledSites.bed")
        
        # if no callable region, skip the cell
        if (file.info(callable_region_file)$size == 0){
          callable_region_read_family_bc2_size[cell_id] <- 0
          ssnv_count[cell_id] <- 0
          next
        }
        
        callable_region <- read.table(callable_region_file)
        colnames(callable_region)[4] <- "read_family"
        
        # read in read family information
        read_family_info_bc2 <- read.table(paste0("./data/", case_id, "/single_cell/read_family_strand_barcode/strand2_only/", cell_id, "_read_family_strand_barcode.txt"))
        read_family_info_bc2[, "read_family"] <- str_split(read_family_info_bc2$V1, pattern = ":", simplify = T)[,3]
        read_family_info_bc2$read_family <- as.integer(read_family_info_bc2$read_family)
        
        callable_region_read_family_bc2 <- inner_join(callable_region, read_family_info_bc2, by = "read_family")
        
        if(dim(callable_region_read_family_bc2)[1]){
          callable_region_read_family_bc2_size[cell_id] <- 
            sum(ifelse(callable_region_read_family_bc2$V3.x - callable_region_read_family_bc2$V2.x > 22,
                       callable_region_read_family_bc2$V3.x - callable_region_read_family_bc2$V2.x, 0))
        }else{
          callable_region_read_family_bc2_size[cell_id] <- 0
        }
        
        # count the number of ssnv from the read families with only strand2 bc
        if(sum(ssnv_summary$cell_id == cell_id)){
          ssnv_count[cell_id] <- sum(ssnv_summary[ssnv_summary$cell_id == cell_id, "read_family"] %in% read_family_info_bc2$read_family)
        }else{
          ssnv_count[cell_id] <- 0
        }
        
        if(length(ssnv_count) %% 100 == 0){print(length(ssnv_count))}
      }
      
      ssnv_rate_strand2[paste0("a", a)] <- sum(ssnv_count)/sum(callable_region_read_family_bc2_size)
    }
    
    error_rate_incorrect_rf <- cbind(ssnv_rate_strand2, ssnv_rate_strand12)
    print(error_rate_incorrect_rf)
    error_rate_incorrect_rf_cell_type[[cell_type]] <- error_rate_incorrect_rf
    
  }
  
  saveRDS(error_rate_incorrect_rf_cell_type, paste0("./results/error_rate_incorrect_rf/", case_id, ".rds"))
}

