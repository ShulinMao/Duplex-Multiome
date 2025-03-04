.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(stringr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(MASS)
library(truncnorm)
library(VariantAnnotation)
set.seed(42)

# read in sample & cell type from command lines
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

case_id <- args[1]
rm(args)

print(case_id)

# read in germline mutation file
germline_mutation <- readVcfAsVRanges(paste0("./data/", case_id, "/bulk/germline_mutation/", case_id, "_germline_het_variants.vcf.gz"), "hg38")

# cell type
file_list <- dir(paste0("./data/", case_id, "/cell_type/"))
cell_type_list <- str_remove(file_list, "[.]txt")

ssDNA_mut_rate_cell_type <- list()
ssnv_bc2_list <- list()
for(cell_type in cell_type_list){
  print(cell_type)
  read.table(paste0("./data/", case_id, "/cell_type/", cell_type, ".txt"), sep = "\t") -> selected_cells
  
  ssDNA_mut_rate <- c()
  # ------------------------------------------------------------------------------
  # Estimate single strand burden
  for(a in seq(2, 10, 2)){
    
    # strand2(strand1bc) + strand2
    ssnv_summary <- read.table(paste0("./data/", case_id, "/single_cell/call_variants/autosomal_variants/a2s0_autosome_variants.txt"))
    ssnv_summary[, "read_family"] <- str_split(ssnv_summary$V10, pattern = ":", simplify = T)[,5]
    ssnv_summary[, "cell_id"] <- str_extract(str_split(ssnv_summary$V11, pattern = "[.]", simplify = T)[,1], "[ATCG]*-1$")
    ssnv_summary[, "PS"] <- as.integer(str_split(ssnv_summary$V10, pattern = ":", simplify = T)[,3])
    ssnv_summary[, "MS"] <- as.integer(str_split(ssnv_summary$V10, pattern = ":", simplify = T)[,4])    
    ssnv_summary <- ssnv_summary[(ssnv_summary$PS + ssnv_summary$MS) > a, ] 
    
    ssnv_count <- c()
    total_germline_count <- c()
    detected_germline_mutation_count <- c()
    callable_region_read_family_bc2_size <- c()
    ssnv_bc2 <- c()
    for (cell_id in selected_cells$V2){
      # detected germline mutation
      germline_het_mutations_file <- paste0("./data/", case_id, "/single_cell/call_variants/germline_mutations/all_germline_het_mutations_a2s0/", cell_id, "_germline_het_mutations.vcf.gz")
      # if no detected germline mutationn, skip the cell
      germline_het_mutations_file_last_line <- system(paste("zcat", germline_het_mutations_file, "| tail -n 1"), intern = T)
      if (str_starts(germline_het_mutations_file_last_line, "#CHROM")){
        callable_region_read_family_bc2_size[cell_id] <- 0
        ssnv_count[cell_id] <- 0
        detected_germline_mutation_count[cell_id] <- 0
        total_germline_count[cell_id] <- 0
        next
      }
      
      detected_germline_mutation <- read.table(germline_het_mutations_file)
      detected_germline_mutation[, "read_family"] <- str_split(detected_germline_mutation$V10, pattern = ":", simplify = T)[,5]
      detected_germline_mutation[, "PS"] <- as.integer(str_split(detected_germline_mutation$V10, pattern = ":", simplify = T)[,3])
      detected_germline_mutation[, "MS"] <- as.integer(str_split(detected_germline_mutation$V10, pattern = ":", simplify = T)[,4])
      detected_germline_mutation_stringency <- detected_germline_mutation[(detected_germline_mutation$PS > a | detected_germline_mutation$MS > a),]
      
      # read in callable region 
      callable_region_file <- paste0("./data/", case_id, "/single_cell/call_variants/autosomal_callable_region/a2s0/", cell_id, ".families.analysis.calledSites.bed")
      
      # if no callable region, skip the cell
      if (file.info(callable_region_file)$size == 0){
        callable_region_read_family_bc2_size[cell_id] <- 0
        ssnv_count[cell_id] <- 0
        detected_germline_mutation_count[cell_id] <- 0
        total_germline_count[cell_id] <- 0
        next
      }
      
      callable_region <- read.table(callable_region_file)
      colnames(callable_region)[4] <- "read_family"
      callable_region <- callable_region[callable_region$V7 >= a*2 | callable_region$V8 >= a*2, ]
      
      
      # read in read family information
      read_family_info_bc2 <- read.table(paste0("./data/", case_id, "/single_cell/read_family_strand_barcode/strand2_only/", cell_id, "_read_family_strand_barcode.txt"))
      read_family_info_bc2 <- read_family_info_bc2[read_family_info_bc2$V7 > 1,]
      
      # if no read families, skip the cell
      if (dim(read_family_info_bc2)[1] == 0){
        callable_region_read_family_bc2_size[cell_id] <- 0
        ssnv_count[cell_id] <- 0
        next
      }
      
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
        ssnv_bc2 <- rbind(ssnv_bc2,
                          ssnv_summary[ssnv_summary$cell_id == cell_id & ssnv_summary$read_family %in% callable_region_read_family_bc2$read_family,])
        ssnv_count[cell_id] <- sum(ssnv_summary[ssnv_summary$cell_id == cell_id, "read_family"] %in% callable_region_read_family_bc2$read_family)
      }else{
        ssnv_count[cell_id] <- 0
      }
      
      # count total germline mutation within callable regions
      callable_region_read_family_bc2 <- 
        callable_region_read_family_bc2[callable_region_read_family_bc2$V3.x - callable_region_read_family_bc2$V2.x > 22,]
      callable_region_read_family_bc2_grange <- GRanges(seqnames = callable_region_read_family_bc2$V1.x,
                                                        ranges = IRanges(start = callable_region_read_family_bc2$V2.x+11,
                                                                         end = callable_region_read_family_bc2$V3.x-11))
      intersection <- findOverlaps(germline_mutation, callable_region_read_family_bc2_grange)
      intersected_ranges <- germline_mutation[queryHits(intersection)]
      total_germline_count[cell_id] <- length(intersected_ranges)
      
      # count detected germline mutation within callable regions
      if(dim(detected_germline_mutation_stringency)[1] > 0){
        #detected_germline_mutation_stringency[detected_germline_mutation_stringency$cell_id == cell_id & detected_germline_mutation_stringency$read_family %in% callable_region_read_family_bc2$read_family,]
        detected_germline_mutation_count[cell_id] <- 
          sum(detected_germline_mutation_stringency[, "read_family"] %in% callable_region_read_family_bc2$read_family)
      }else{
        detected_germline_mutation_count[cell_id] <- 0
      }
      
      if(length(ssnv_count) %% 100 == 0){print(length(ssnv_count))}
    }
    
    sensitivity <- sum(detected_germline_mutation_count)/(sum(total_germline_count)/2)
    ssDNA_mut_rate <- rbind(ssDNA_mut_rate, c(sum(ssnv_count)/sum(callable_region_read_family_bc2_size)/sensitivity, sum(ssnv_count), sum(callable_region_read_family_bc2_size), sensitivity, a))
    print(ssDNA_mut_rate)
    if(sum(ssnv_count) < 3){break}
    if(a == 2){ssnv_bc2_list[[cell_type]] <- ssnv_bc2}
  }
  ssDNA_mut_rate_cell_type[[cell_type]] <- ssDNA_mut_rate
}
saveRDS(ssDNA_mut_rate_cell_type, paste0("./results/ssDNA_mut/", case_id, ".rds"))
saveRDS(ssnv_bc2_list, paste0("./results/ssDNA_mut/", case_id, "_spectrum.rds"))

