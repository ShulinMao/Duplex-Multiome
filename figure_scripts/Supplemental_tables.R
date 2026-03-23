# Libraries
.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(stringr)
library(WriteXLS)
library(tidyr)

# Neurotypical samples
duplex_metadata <- read.csv("./data/others/metadata_07012024.csv")

duplex_snv_candidate_a6s3 <- c()
for (sample_id in duplex_metadata$Sample){
  
  file_path <- paste0("./data/", sample_id, "/single_cell/call_variants/autosomal_variants/a6s3_autosome_variants.txt")
  if(file.info(file_path)$size==0){next}
  duplex_snv_candidate_a6s3_sample <- read.table(paste0("./data/", sample_id, "/single_cell/call_variants/autosomal_variants/a6s3_autosome_variants.txt"))
  colnames(duplex_snv_candidate_a6s3_sample) <- 
    c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "VALUE", "CELL_ID")
  if (sample_id == "Brain1901106_deeper"){sample_id <- "Brain190106"}
  duplex_snv_candidate_a6s3_sample[,"SAMPLE_ID"] <- str_remove(sample_id, "_deeper$")
  
  duplex_snv_candidate_a6s3 <- rbind(duplex_snv_candidate_a6s3, duplex_snv_candidate_a6s3_sample)
  
}
duplex_snv_candidate_a6s3 <- data.frame(duplex_snv_candidate_a6s3)
duplex_snv_candidate_a6s3 %>% separate(col = VALUE, 
                                               into = c("GENOTYPE", "DEPTH", 
                                                        "STRAND1_DEPTH", "STRAND2_DEPTH", 
                                                        "READ_FAMILY"), sep = ":") -> duplex_snv_candidate_a6s3
duplex_snv_candidate_a6s3[,c("CHROM", "POS", "ID", "REF",
                                     "ALT", "QUAL", "FILTER",
                                     "DEPTH", "STRAND1_DEPTH", "STRAND2_DEPTH", 
                                     "READ_FAMILY", "CELL_ID", "SAMPLE_ID")] -> duplex_snv_candidate_a6s3

# COLO829BLT50
duplex_snv_candidate_a6s3_COLO829 <- c()
for (sample_id in c("COLO829BLT50_rep1", "COLO829BLT50_rep2")){
  
  file_path <- paste0("./data/", sample_id, "/single_cell/call_variants/autosomal_variants/a6s3_autosome_variants.txt")
  if(file.info(file_path)$size==0){next}
  duplex_snv_candidate_a6s3_sample <- read.table(paste0("./data/", sample_id, "/single_cell/call_variants/autosomal_variants/a6s3_autosome_variants.txt"))
  colnames(duplex_snv_candidate_a6s3_sample) <- 
    c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "VALUE", "CELL_ID")
  duplex_snv_candidate_a6s3_sample[,"SAMPLE_ID"] <- str_remove(sample_id, "_deeper$")
  
  duplex_snv_candidate_a6s3_COLO829 <- rbind(duplex_snv_candidate_a6s3_COLO829, duplex_snv_candidate_a6s3_sample)
}
duplex_snv_candidate_a6s3_COLO829 <- data.frame(duplex_snv_candidate_a6s3_COLO829)
duplex_snv_candidate_a6s3_COLO829 %>% separate(col = VALUE, 
                                           into = c("GENOTYPE", "DEPTH", 
                                                    "STRAND1_DEPTH", "STRAND2_DEPTH", 
                                                    "READ_FAMILY"), sep = ":") -> duplex_snv_candidate_a6s3_COLO829
duplex_snv_candidate_a6s3_COLO829[,c("CHROM", "POS", "ID", "REF",
                                 "ALT", "QUAL", "FILTER",
                                 "DEPTH", "STRAND1_DEPTH", "STRAND2_DEPTH", 
                                 "READ_FAMILY", "CELL_ID", "SAMPLE_ID")] -> duplex_snv_candidate_a6s3_COLO829

# ASD
duplex_snv_candidate_a6s3_ASD <- c()
for (sample_id in c("AN06365", "ABNCR6D")){
  
  file_path <- paste0("./data/", sample_id, "/single_cell/call_variants/autosomal_variants/a6s3_autosome_variants.txt")
  if(file.info(file_path)$size==0){next}
  duplex_snv_candidate_a6s3_sample <- read.table(paste0("./data/", sample_id, "/single_cell/call_variants/autosomal_variants/a6s3_autosome_variants.txt"))
  colnames(duplex_snv_candidate_a6s3_sample) <- 
    c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "VALUE", "CELL_ID")
  duplex_snv_candidate_a6s3_sample[,"SAMPLE_ID"] <- str_remove(sample_id, "_deeper$")
  
  duplex_snv_candidate_a6s3_ASD <- rbind(duplex_snv_candidate_a6s3_ASD, duplex_snv_candidate_a6s3_sample)
}
duplex_snv_candidate_a6s3_ASD <- data.frame(duplex_snv_candidate_a6s3_ASD)
duplex_snv_candidate_a6s3_ASD %>% separate(col = VALUE, 
                                           into = c("GENOTYPE", "DEPTH", 
                                                    "STRAND1_DEPTH", "STRAND2_DEPTH", 
                                                    "READ_FAMILY"), sep = ":") -> duplex_snv_candidate_a6s3_ASD
duplex_snv_candidate_a6s3_ASD[,c("CHROM", "POS", "ID", "REF",
                                 "ALT", "QUAL", "FILTER",
                                 "DEPTH", "STRAND1_DEPTH", "STRAND2_DEPTH", 
                                 "READ_FAMILY", "CELL_ID", "SAMPLE_ID")] -> duplex_snv_candidate_a6s3_ASD

# GBM
duplex_snv_candidate_a6s3_GBM <- c()
for (sample_id in c("UMB4397_tumor", "UMB4397_normal")){
  
  file_path <- paste0("./data/", sample_id, "/single_cell/call_variants/autosomal_variants/a6s3_autosome_variants.txt")
  if(file.info(file_path)$size==0){next}
  duplex_snv_candidate_a6s3_sample <- read.table(paste0("./data/", sample_id, "/single_cell/call_variants/autosomal_variants/a6s3_autosome_variants.txt"))
  colnames(duplex_snv_candidate_a6s3_sample) <- 
    c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "VALUE", "CELL_ID")
  duplex_snv_candidate_a6s3_sample[,"SAMPLE_ID"] <- str_remove(sample_id, "_deeper$")
  
  duplex_snv_candidate_a6s3_GBM <- rbind(duplex_snv_candidate_a6s3_GBM, duplex_snv_candidate_a6s3_sample)
}
duplex_snv_candidate_a6s3_GBM <- data.frame(duplex_snv_candidate_a6s3_GBM)
duplex_snv_candidate_a6s3_GBM %>% separate(col = VALUE, 
                                           into = c("GENOTYPE", "DEPTH", 
                                                    "STRAND1_DEPTH", "STRAND2_DEPTH", 
                                                    "READ_FAMILY"), sep = ":") -> duplex_snv_candidate_a6s3_GBM
duplex_snv_candidate_a6s3_GBM[,c("CHROM", "POS", "ID", "REF",
                                 "ALT", "QUAL", "FILTER",
                                 "DEPTH", "STRAND1_DEPTH", "STRAND2_DEPTH", 
                                 "READ_FAMILY", "CELL_ID", "SAMPLE_ID")] -> duplex_snv_candidate_a6s3_GBM

write.csv(duplex_snv_candidate_a6s3, 
          "./results/supplemental_tables/neurotypical_sample_private_variants_a6s3.csv", quote = F, row.names = F)
write.csv(duplex_snv_candidate_a6s3_ASD, 
          "./results/supplemental_tables/ASD_sample_private_variants_a6s3.csv", quote = F, row.names = F)
write.csv(duplex_snv_candidate_a6s3_COLO829, 
          "./results/supplemental_tables/COLO829BLT50_sample_private_variants_a6s3.csv", quote = F, row.names = F)
write.csv(duplex_snv_candidate_a6s3_GBM, 
          "./results/supplemental_tables/GBM_sample_private_variants_a6s3.csv", quote = F, row.names = F)


