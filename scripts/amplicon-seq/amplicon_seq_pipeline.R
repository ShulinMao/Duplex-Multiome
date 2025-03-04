.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(readxl)
library(stringr)
library(tidyr)

# set up the output directory
output_dir <- "/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/amplicon-seq/"

# read in the amplicon variant list
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
         into = c("chr", "pos", "ref", "alt"), sep = "-") -> amplicon_variant_list

# write the script to run the amplicon-seq pipeline for each amplicon
cat("echo \"start\"", file = "./scripts/amplicon-seq/run_amplicon_seq_calling.sh")
for (amplicon_list in amplicon_fastq){
  amplicons <- str_split(amplicon_list, "-", simplify = T)[1,]
  selected <- amplicon_variant_list$id %in% amplicons
  
  dir.create(paste0(output_dir, amplicon_list), showWarnings = F)
  
  amplicon_variant_list[selected, c("chr", "pos")] -> selected_amplicon_variant_list
  
  # write the candidate positions to a bed file
  write.table(selected_amplicon_variant_list[, c("chr", "pos", "pos")],
              paste0(output_dir, amplicon_list, "/candidate_pos.bed"), 
              append = F, quote = F, sep = "\t", row.names = F, col.names = F)
  
  as.integer(selected_amplicon_variant_list$pos) -> selected_amplicon_variant_list$pos
  amplicon_variant_surrounding_pos <- c()
  for (amplicon in seq(dim(selected_amplicon_variant_list)[1])){
    amplicon_variant_surrounding_pos <- rbind(amplicon_variant_surrounding_pos,
      cbind(rep(selected_amplicon_variant_list$chr[amplicon], 20),
            c(seq(selected_amplicon_variant_list$pos[amplicon]-10, selected_amplicon_variant_list$pos[amplicon]-1),
              seq(selected_amplicon_variant_list$pos[amplicon]+1, selected_amplicon_variant_list$pos[amplicon]+10)))
    )
  }
  
  # write the surrounding positions to a bed file
  write.table(amplicon_variant_surrounding_pos[, c(1, 2, 2)],
              paste0(output_dir, amplicon_list, "/surrounding_pos.bed"), 
              append = F, quote = F, sep = "\t", row.names = F, col.names = F)
  
  cmd <- paste("sbatch /lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/scripts/amplicon-seq/amplicon_seq_calling.sh",
               amplicon_fastq_dir, amplicon_list, paste0(output_dir, amplicon_list), "\n")
  cat(cmd, file = "./scripts/amplicon-seq/run_amplicon_seq_calling.sh", append = T)
}



