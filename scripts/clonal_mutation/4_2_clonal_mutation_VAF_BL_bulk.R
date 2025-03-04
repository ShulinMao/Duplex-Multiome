.libPaths(c("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/"))
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(stringr)
library(reshape2)
library(dplyr)
library(VariantAnnotation)
library(GenomicRanges)
library(ggplot2)

case_id_list <- c("COLO829BLT50_rep1", "COLO829BLT50_rep2")
cell_type <- "melanoma_fibroblast"

cell_type_list <- read.table("./results/celltype_annotation/annotation_COLO829BLT50_20240523.txt",
                        sep = "\t")
colnames(cell_type_list) <- c("case_id", "cell", "cell_type")
cell_type_list <- cell_type_list[cell_type_list$cell_type == cell_type,]

# ------------------------------------------------------------------------------
# collect detected/covered
for(case_id in case_id_list){
  print(case_id)
  
  # read in clonal mutation file
  clonal_mutation <- read.csv("./results/clonal_mutation/clonal_mutation_list_COLO829BLT50_BL_bulk_20240604.csv")
  clonal_mutation <- na.omit(clonal_mutation)
  clonal_mutation <- clonal_mutation[clonal_mutation$cell_type == cell_type,]
  clonal_mutation <- clonal_mutation[clonal_mutation$case_id == case_id,]
  clonal_mutation$key <- str_replace_all(clonal_mutation$key, pattern = "[:]", replacement = "-")
  
  # store a map between cell ids and covered regions
  cell_cover_clonal_mutation_site <- clonal_mutation %>% group_by(key) %>% summarise("detected" = n())
  cell_cover_clonal_mutation_site[, "covered"] <- 0
  cell_cover_clonal_mutation_site <- data.frame(cell_cover_clonal_mutation_site)
  row.names(cell_cover_clonal_mutation_site) <- cell_cover_clonal_mutation_site$key
  cell_cover_clonal_mutation_site <- cell_cover_clonal_mutation_site[,-1]
  
  clonal_mutation <- clonal_mutation[!duplicated(clonal_mutation$key),]
  clonal_mutation_grange <- GRanges(seqnames = clonal_mutation$chr,
                                    ranges = IRanges(start = clonal_mutation$pos, 
                                                     end = clonal_mutation$pos),
                                    ref = clonal_mutation$ref, 
                                    alt = clonal_mutation$alt)
  
  # read in covered region info 
  dir_path <- paste0("./data/", case_id, "/single_cell/read_families/a2s0_0minReadFamilyLength/")
  cell_callable_region_file_list <- str_c(cell_type_list[cell_type_list$case_id == case_id, "cell"],
                                          ".families.analysis.calledSites.bed")
  
  count <- 1
  for(callable_region_file in cell_callable_region_file_list){
    cell_id <- str_extract(callable_region_file, "^[ATCG]+-1")
    
    # callable regions
    # if no callable region, skip the cell
    if (file.info(paste0(dir_path, callable_region_file))$size == 0){next}
    
    callable_region <- read.table(paste0(dir_path, callable_region_file))
    callable_region_grange <- GRanges(seqnames = callable_region$V1,
                                      ranges = IRanges(start = callable_region$V2,
                                                       end = callable_region$V3))
    
    # find cell with reads covering germline mutations
    intersection <- findOverlaps(clonal_mutation_grange, callable_region_grange)
    intersected_ranges <- clonal_mutation_grange[queryHits(intersection)]
    if(length(intersected_ranges) > 0){
      covered_clonal_mutation_site <- 
        str_c(intersected_ranges@seqnames, intersected_ranges@ranges@start, intersected_ranges$ref, intersected_ranges$alt, sep = "-")
      cell_cover_clonal_mutation_site[covered_clonal_mutation_site, "covered"] <- 
        cell_cover_clonal_mutation_site[covered_clonal_mutation_site, "covered"] + 1
    }
    
    count <- count + 1
    if(count%%50 == 0){print(paste(case_id, count))}
  }
  sum(is.na(cell_cover_clonal_mutation_site))
  #cell_cover_clonal_mutation_site <- na.omit(cell_cover_clonal_mutation_site)
  saveRDS(cell_cover_clonal_mutation_site, paste0("./results/clonal_mutation/clonal_mutation_VAF_BL_bulk_melanoma_fibroblast_", case_id, ".rds"))
}

