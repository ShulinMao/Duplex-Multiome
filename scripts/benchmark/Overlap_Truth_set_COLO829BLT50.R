.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(stringr)
library(tidyr)
library(dplyr)
library(ggpubr)

# Truth set
smaht_truthset <- read.table("./data/others/COLO829BLT50_truth_set/SMaHT_COLO829_SNV_truth_set_v1.0.vcf")
smaht_truthset <- str_c(smaht_truthset$V1, smaht_truthset$V2, smaht_truthset$V4, smaht_truthset$V5, sep = ":")

# Duplex Multiome results
# a6s3 calls 
case_id_list <- c("COLO829BLT50_rep1", "COLO829BLT50_rep2")
duplex_multiome_a6s3_mutation_BLT50_bulk <- c()
for (case_id in case_id_list){
  mutation <- read.table(paste0("./data/", case_id, "/single_cell/call_variants/autosomal_variants/a6s3_autosome_variants.txt"), header = F, sep = "\t")
  colnames(mutation) <- c("chr", "pos", "id", "ref", "alt", "qual", "filter", "info", "fromat", "value")
  mutation <- separate(mutation, col = "value", into = c("value", "cell"), sep = " ")
  mutation[,"key"] <- unite(mutation, col = "key", c("chr", "pos", "ref", "alt"), sep = ":")[,"key"]
  mutation[,"case_id"] <- case_id
  
  duplex_multiome_a6s3_mutation_BLT50_bulk <- rbind(duplex_multiome_a6s3_mutation_BLT50_bulk, mutation)
}
cell_type <- read.table("./results/celltype_annotation/annotation_COLO829BLT50_20240523.txt",
                        sep = "\t")
colnames(cell_type) <- c("case_id", "cell", "cell_type")
duplex_multiome_a6s3_mutation_BLT50_bulk <- left_join(duplex_multiome_a6s3_mutation_BLT50_bulk, cell_type, by = c("case_id", "cell"))
duplex_multiome_a6s3_mutation_BLT50_bulk <- duplex_multiome_a6s3_mutation_BLT50_bulk[!is.na(duplex_multiome_a6s3_mutation_BLT50_bulk$cell_type),]
duplex_multiome_a6s3_mutation_BLT50_bulk <- duplex_multiome_a6s3_mutation_BLT50_bulk[duplex_multiome_a6s3_mutation_BLT50_bulk$cell_type == "melanoma_fibroblast",]
duplex_multiome_a6s3_mutation_BLT50_bulk <- unique(duplex_multiome_a6s3_mutation_BLT50_bulk$key)

duplex_multiome_a6s3_mutation_BL_bulk <- c()
for (case_id in case_id_list){
  # read in all clonal mutations
  mutation <- read.table(paste0("./data/", case_id, "/single_cell/call_variants/filter_germline_mutation/a2s0_BL_bulk.vcf"), header = F, sep = "\t")
  colnames(mutation) <- c("chr", "pos", "cell", "ref", "alt", "qual", "filter", "info", "fromat", "value")
  mutation[,"key"] <- unite(mutation, col = "key", c("chr", "pos", "ref", "alt"), sep = ":")[,"key"]
  mutation[,"case_id"] <- case_id
  mutation <- mutation[!(duplicated(str_c(mutation$cell, mutation$key, mutation$case_id))),]
  mutation <- separate(mutation, col = "value", sep = ":", into = c("GT", "DP", "PS", "MS", "RF"), remove = F)
  
  mutation <- mutation[mutation$PS >= 3 & mutation$MS >= 3,]
  duplex_multiome_a6s3_mutation_BL_bulk <- rbind(duplex_multiome_a6s3_mutation_BL_bulk, mutation)
}
cell_type <- read.table("./results/celltype_annotation/annotation_COLO829BLT50_20240523.txt",
                        sep = "\t")
colnames(cell_type) <- c("case_id", "cell", "cell_type")
duplex_multiome_a6s3_mutation_BL_bulk <- left_join(duplex_multiome_a6s3_mutation_BL_bulk, cell_type, by = c("case_id", "cell"))
duplex_multiome_a6s3_mutation_BL_bulk <- duplex_multiome_a6s3_mutation_BL_bulk[!is.na(duplex_multiome_a6s3_mutation_BL_bulk$cell_type),]
duplex_multiome_a6s3_mutation_BL_bulk <- duplex_multiome_a6s3_mutation_BL_bulk[duplex_multiome_a6s3_mutation_BL_bulk$cell_type == "melanoma_fibroblast",]
duplex_multiome_a6s3_mutation_BL_bulk <- unique(duplex_multiome_a6s3_mutation_BL_bulk$key)

# clonal calls
duplex_multiome_clonal_mutation_vaf <- c()
for (case_id in case_id_list){
  mutation <- readRDS(paste0("./results/clonal_mutation/clonal_mutation_VAF_BL_bulk_melanoma_fibroblast_", case_id, ".rds"))
  mutation[,"key"] <- str_replace_all(row.names(mutation), pattern = "-", replacement = ":")
  
  duplex_multiome_clonal_mutation_vaf <- rbind(duplex_multiome_clonal_mutation_vaf, mutation)
}
duplex_multiome_clonal_mutation_vaf <- 
  duplex_multiome_clonal_mutation_vaf %>% group_by(key) %>% summarise("detected" = sum(detected),
                                                                    "covered" = sum(covered))
duplex_multiome_clonal_mutation_vaf[,"VAF"] <- duplex_multiome_clonal_mutation_vaf$detected/duplex_multiome_clonal_mutation_vaf$covered

duplex_multiome_clonal_mutation <- read.csv("results/clonal_mutation/clonal_mutation_list_COLO829BLT50_BL_bulk_realigned_20241007.csv")
duplex_multiome_clonal_mutation <- duplex_multiome_clonal_mutation[!is.na(duplex_multiome_clonal_mutation$cell_type),]
duplex_multiome_clonal_mutation <- duplex_multiome_clonal_mutation[duplex_multiome_clonal_mutation$cell_type == "melanoma_fibroblast",]
duplex_multiome_clonal_mutation <- unique(duplex_multiome_clonal_mutation$key)

duplex_multiome_clonal_mutation_vaf <- 
  duplex_multiome_clonal_mutation_vaf[duplex_multiome_clonal_mutation_vaf$key %in% duplex_multiome_clonal_mutation,]
hist(duplex_multiome_clonal_mutation_vaf$VAF[duplex_multiome_clonal_mutation_vaf$key %in% smaht_truthset])
hist(duplex_multiome_clonal_mutation_vaf$VAF[!(duplex_multiome_clonal_mutation_vaf$key %in% smaht_truthset)])

duplex_multiome_clonal_mutation_0.25 <- 
  duplex_multiome_clonal_mutation_vaf$key[duplex_multiome_clonal_mutation_vaf$VAF > 0.25*(0.32/0.5)]

sum(duplex_multiome_a6s3_mutation_BLT50_bulk %in% smaht_truthset)/length(duplex_multiome_a6s3_mutation_BLT50_bulk)
sum(duplex_multiome_a6s3_mutation_BL_bulk %in% smaht_truthset)/length(duplex_multiome_a6s3_mutation_BL_bulk)
sum(duplex_multiome_clonal_mutation %in% smaht_truthset)/length(duplex_multiome_clonal_mutation)
sum(duplex_multiome_clonal_mutation_0.25 %in% smaht_truthset)/length(duplex_multiome_clonal_mutation_0.25)

overlap_df <- data.frame("duplex_multiome_a6s3_mutation_BLT50_bulk" = c(sum(duplex_multiome_a6s3_mutation_BLT50_bulk %in% smaht_truthset),
                                                                        sum(!(duplex_multiome_a6s3_mutation_BLT50_bulk %in% smaht_truthset))),
                         "duplex_multiome_a6s3_mutation_BL_bulk" = c(sum(duplex_multiome_a6s3_mutation_BL_bulk %in% smaht_truthset),
                                                                     sum(!(duplex_multiome_a6s3_mutation_BL_bulk %in% smaht_truthset))),
                         "duplex_multiome_clonal_mutation" =  c(sum(duplex_multiome_clonal_mutation %in% smaht_truthset),
                                                                sum(!(duplex_multiome_clonal_mutation %in% smaht_truthset))),
                         "duplex_multiome_clonal_mutation_0.25" = c(sum(duplex_multiome_clonal_mutation_0.25 %in% smaht_truthset),
                                                                    sum(!(duplex_multiome_clonal_mutation_0.25 %in% smaht_truthset)))
)
overlap_df[,"group"] <- c("In the Truth set", "Not in the Truth set")
overlap_df <- reshape2::melt(overlap_df)
overlap_df <- overlap_df %>% 
  group_by(variable) %>%
  mutate(ypos = cumsum(value)- 0.5*value)
overlap_df$group <- factor(overlap_df$group, levels = c("Not in the Truth set", "In the Truth set"))
saveRDS(overlap_df, "./figures/manuscript_figures/Figure2/benchmark_truthset_COLO829BLT50.rds")

