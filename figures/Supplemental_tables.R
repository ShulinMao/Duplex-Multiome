# Libraries
.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(stringr)
library(WriteXLS)

# Neurotypical samples
duplex_metadata <- read.csv("./data/others/metadata_07012024.csv")

duplex_snv_candidate_a6s3 <- c()
for (sample_id in duplex_metadata$Sample){
  
  file_path <- paste0("./data/", sample_id, "/single_cell/call_variants/autosomal_variants/a6s3_autosome_variants.txt")
  if(file.info(file_path)$size==0){next}
  duplex_snv_candidate_a6s3_sample <- read.table(paste0("./data/", sample_id, "/single_cell/call_variants/autosomal_variants/a6s3_autosome_variants.txt"))
  colnames(duplex_snv_candidate_a6s3_sample) <- 
    c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
  if (sample_id == "Brain1901106_deeper"){sample_id <- "Brain190106"}
  duplex_snv_candidate_a6s3_sample[,"SAMPLE_ID"] <- str_remove(sample_id, "_deeper$")
  
  duplex_snv_candidate_a6s3 <- rbind(duplex_snv_candidate_a6s3, duplex_snv_candidate_a6s3_sample)
}

# COLO829BLT50
duplex_snv_candidate_a6s3_COLO829 <- c()
for (sample_id in c("COLO829BLT50_rep1", "COLO829BLT50_rep2")){
  
  file_path <- paste0("./data/", sample_id, "/single_cell/call_variants/autosomal_variants/a6s3_autosome_variants.txt")
  if(file.info(file_path)$size==0){next}
  duplex_snv_candidate_a6s3_sample <- read.table(paste0("./data/", sample_id, "/single_cell/call_variants/autosomal_variants/a6s3_autosome_variants.txt"))
  colnames(duplex_snv_candidate_a6s3_sample) <- 
    c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
  duplex_snv_candidate_a6s3_sample[,"SAMPLE_ID"] <- str_remove(sample_id, "_deeper$")
  
  duplex_snv_candidate_a6s3_COLO829 <- rbind(duplex_snv_candidate_a6s3_COLO829, duplex_snv_candidate_a6s3_sample)
}

# ASD
duplex_snv_candidate_a6s3_ASD <- c()
for (sample_id in c("AN06365", "ABNCR6D")){
  
  file_path <- paste0("./data/", sample_id, "/single_cell/call_variants/autosomal_variants/a6s3_autosome_variants.txt")
  if(file.info(file_path)$size==0){next}
  duplex_snv_candidate_a6s3_sample <- read.table(paste0("./data/", sample_id, "/single_cell/call_variants/autosomal_variants/a6s3_autosome_variants.txt"))
  colnames(duplex_snv_candidate_a6s3_sample) <- 
    c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
  duplex_snv_candidate_a6s3_sample[,"SAMPLE_ID"] <- str_remove(sample_id, "_deeper$")
  
  duplex_snv_candidate_a6s3_ASD <- rbind(duplex_snv_candidate_a6s3_ASD, duplex_snv_candidate_a6s3_sample)
}

write.csv(duplex_snv_candidate_a6s3, 
          "./results/supplemental_tables/neurotypical_sample_private_variants_a6s3.csv", quote = F, row.names = F)
write.csv(duplex_snv_candidate_a6s3_ASD, 
          "./results/supplemental_tables/ASD_sample_private_variants_a6s3.csv", quote = F, row.names = F)
write.csv(duplex_snv_candidate_a6s3_COLO829, 
          "./results/supplemental_tables/COLO829BLT50_sample_private_variants_a6s3.csv", quote = F, row.names = F)


