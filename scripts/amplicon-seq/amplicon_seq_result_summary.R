.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(readxl)
library(stringr)
library(tidyr)
library(dplyr)

find_second_largest <- function(x){
  as.integer(sort(x, decreasing = T)[2])/sum(x)
}

# ------------------------------------------------------------------------------
# amplicon-seq result summary
base <- c("A", "G", "C", "T")

# read in the amplicon variant list
amplicon_result_dir <- "./data/amplicon-seq/"

amplicon_fastq_dir <- "/lab-share/Gene-Walsh-e2/Public/indData/AK/Multiome-strand-02212022/Clonal-Amplicon-seq-12042024/30-1115699648/00_fastq/"
amplicon_fastq <- dir(amplicon_fastq_dir,
                      pattern = "_R1_001.fastq.gz$", include.dirs = F)
amplicon_fastq <- str_remove(amplicon_fastq, "_R1_001.fastq.gz$")

amplicon_variant_list <- read_excel("data/amplicon-seq/Validation_12102024.xlsx", 
                                    col_types = c("text", "text", "skip", 
                                                  "skip", "skip", "skip", "skip", "skip", 
                                                  "skip"))
colnames(amplicon_variant_list) <- c("id", "mutation")

str_remove(amplicon_variant_list$id, "_*AK") -> amplicon_variant_list$id
separate(amplicon_variant_list, col = "mutation", 
         into = c("Chrom", "Pos", "Ref", "Alt"), sep = "-") -> amplicon_variant_list
amplicon_variant_list$Pos <- as.integer(amplicon_variant_list$Pos)

amplicon_result <- c()
for (i in seq(length(amplicon_fastq))){
  mosaic_hunter_result <- read.table(paste0(amplicon_result_dir, amplicon_fastq[i], "/NaiveCaller/", amplicon_fastq[i], ".NC.tsv"), header = T)
  
  amplicon_variant_list[amplicon_variant_list$id %in% str_split(amplicon_fastq[i], "-", simplify = T)[1,],] -> current_amplicon_variant_list
  amplicon_result <- rbind(amplicon_result,
                           inner_join(current_amplicon_variant_list, mosaic_hunter_result, c("Chrom", "Pos", "Ref")))
}

# calculate the count threshold for candidate positions
second_largest_allelic_count_threshold <- data.frame(second_largest_allelic_count_threshold)
colnames(second_largest_allelic_count_threshold) <- c("id", "Chrom", "Pos", "Threshold")
second_largest_allelic_count_threshold$Pos <- as.integer(second_largest_allelic_count_threshold$Pos)
second_largest_allelic_count_threshold$Threshold <- as.numeric(second_largest_allelic_count_threshold$Threshold)

amplicon_result <- inner_join(amplicon_result, second_largest_allelic_count_threshold)

amplicon_result$Filter_sum_following_two <- FALSE # Filter out the variants with less than the sum of the following two bases
amplicon_result$Filter_3x <- FALSE # Filter out the variants with less than 3x coverage of the non-variant and non-ref base
for (i in seq(dim(amplicon_result)[1])){
  non_variant_base <- setdiff(base, c(amplicon_result$Ref[i], amplicon_result$Alt[i]))
  non_variant_base_count <- max(amplicon_result[i, paste0("hQCnt_", non_variant_base)])
  amplicon_result$Filter_3x[i] <- (amplicon_result[i, paste0("hQCnt_", amplicon_result$Alt[i])] > non_variant_base_count*3)[1,1]

  non_variant_base_count <- sum(amplicon_result[i, paste0("hQCnt_", non_variant_base)])
  amplicon_result$Filter_sum_following_two[i] <- (amplicon_result[i, paste0("hQCnt_", amplicon_result$Alt[i])] > non_variant_base_count)[1,1]  
}

# write the result to a csv file
write.csv(amplicon_result, file = paste0(amplicon_result_dir, "healthy_brain_1.csv"),
          quote = F, row.names = F)


