.libPaths(c("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/"))
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(stringr)
library(reshape2)
library(dplyr)
library(VariantAnnotation)
library(GenomicRanges)
library(ggplot2)
set.seed(123)

case_id_list <- c(
  "ABNCR6D"
  #"AN06365")
  #c("Brain1901106_deeper", "UMB1278_deeper", "UMB1864_deeper", "UMB5823_deeper", "UMB4638", "UMB1465", "UMB5451", "UMB5657")
  #c("COLO829BLT50_rep1", "COLO829BLT50_rep2"
  )

# ------------------------------------------------------------------------------
# collect detected/covered
for(case_id in case_id_list){
  print(case_id)
  
  # read in germline mutation file
  germline_mutation <- readVcfAsVRanges(paste0("./data/", case_id, "/bulk/germline_mutation/", case_id, "_germline_het_variants.vcf.gz"), "hg38")
  germline_mutation <- germline_mutation[sapply(germline_mutation@alt, nchar) == 1 & sapply(germline_mutation@ref, nchar) == 1 ,] # remove indels
  germline_mutation_id <- str_c(germline_mutation@seqnames, germline_mutation@ranges@start, germline_mutation@ref, germline_mutation@alt, sep = "-")
  
  # store a map between cell ids and covered regions
  cell_cover_germline_mutation_site <- data.frame(row.names = germline_mutation_id)
  cell_cover_germline_mutation_site[, "detected"] <- 0
  cell_cover_germline_mutation_site[, "covered"] <- 0
  
  # read in covered region info 
  dir_path <- paste0("./data/", case_id, "/single_cell/read_families/a2s0_0minReadFamilyLength/")
  cell_callable_region_file_list <- dir(dir_path, pattern = ".*[.]families[.]analysis[.]calledSites[.]bed")
  
  count <- 1
  for(callable_region_file in cell_callable_region_file_list){
    cell_id <- str_extract(callable_region_file, "^[ATCG]+-1")
    
    # callable regions
    # if no callable region, skip the cell
    if (file.info(paste0(dir_path, callable_region_file))$size == 0){next}
    
    callable_region <- read.table(paste0(dir_path, callable_region_file))
    #callable_region <- callable_region[callable_region$V3 - callable_region$V2 > 22,]
    callable_region_grange <- GRanges(seqnames = callable_region$V1,
                                      ranges = IRanges(start = callable_region$V2,
                                                       end = callable_region$V3))
    
    # find cell with reads covering germline mutations
    intersection <- findOverlaps(germline_mutation, callable_region_grange)
    intersected_ranges <- germline_mutation[queryHits(intersection)]
    covered_germline_mutation_site <- 
      str_c(intersected_ranges@seqnames, intersected_ranges@ranges@start, intersected_ranges@ref, intersected_ranges@alt, sep = "-")
    cell_cover_germline_mutation_site[covered_germline_mutation_site, "covered"] <- 
      cell_cover_germline_mutation_site[covered_germline_mutation_site, "covered"] + 1


    # called germline mutations
    # if no callable region, skip the cell
    if (file.info(paste0("./data/", case_id, "/single_cell/call_variants/germline_mutations/all_germline_het_mutations_a2s0/", cell_id, "_germline_het_mutations.vcf.gz"))$size < 3000){next}

    detected_germline_mutation <- read.table(paste0("./data/", case_id, "/single_cell/call_variants/germline_mutations/all_germline_het_mutations_a2s0/", cell_id, "_germline_het_mutations.vcf.gz"))
    detected_germline_mutation_id <-
      str_c(detected_germline_mutation$V1, detected_germline_mutation$V2, detected_germline_mutation$V4, detected_germline_mutation$V5, sep = "-")
    cell_cover_germline_mutation_site[detected_germline_mutation_id, "detected"] <-
      cell_cover_germline_mutation_site[detected_germline_mutation_id, "detected"] + 1

    count <- count + 1
    if(count%%500 == 0){print(paste(case_id, count))}
  }
  sum(is.na(cell_cover_germline_mutation_site))
  #cell_cover_germline_mutation_site <- na.omit(cell_cover_germline_mutation_site)
  saveRDS(cell_cover_germline_mutation_site, paste0("./results/germline_mutation_VAF/", case_id, ".rds"))
}

