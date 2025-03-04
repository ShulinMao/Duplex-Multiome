# Libraries
.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(stringr)
library(tidyr)

ref_genome="BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = T)
chr_orders <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages_sig/")
library(MutationalPatterns)

# Tumor Truth set spectrum
smaht_truthset <- read_vcfs_as_granges("./data/others/COLO829BLT50_truth_set/SMaHT_COLO829_SNV_truth_set_v1.0.vcf", 
                     "Truth Set", ref_genome)
mut_mat <- mut_matrix(vcf_list = smaht_truthset, ref_genome = ref_genome, condensed = T)
plot_96_profile(mut_mat)  
ggsave("./figures/benchmark_COLO829BLT50/SMaHT_COLO829_SNV_truth_set_v1.0_spectrum.pdf",
       height = 3, width = 14)


# BL spectrum
COLO829BL_truthset <- read_vcfs_as_granges("./data/others/COLO829BLT50_truth_set/COLO829BL_somatic_07252024_snv_af_gt_25_dp_gt_100.vcf.gz",
                                           "Truth Set", ref_genome)

mut_mat <- mut_matrix(vcf_list = COLO829BL_truthset, ref_genome = ref_genome)
plot_96_profile(mut_mat, ymax = 0.08, condensed = T) + 
  scale_y_continuous(breaks = seq(0, 0.08, 0.02)) 
ggsave("./figures/benchmark_COLO829BLT50/SMaHT_COLO829_SNV_B_spectrum-VAF_25.pdf", 
       height = 3, width = 14)

COLO829BL_truthset <- read_vcfs_as_granges("./data/others/COLO829BLT50_truth_set/COLO829BL_somatic_07252024_snv_af1_dp1000.vcf.gz",
                                           "Truth Set", ref_genome)

mut_mat <- mut_matrix(vcf_list = COLO829BL_truthset, ref_genome = ref_genome)
plot_96_profile(mut_mat, ymax = 0.08, condensed = T) + 
  scale_y_continuous(breaks = seq(0, 0.08, 0.02)) 

