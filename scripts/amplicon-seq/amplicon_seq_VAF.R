.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(readxl)
library(ggplot2)

# read in the amplicon-seq VAF
ampliconseq_VAF <- read_excel("./data/amplicon-seq/healthy_brain_ampliconseq_VAF_01142025.xlsx")
ampliconseq_VAF_pass <- ampliconseq_VAF[ampliconseq_VAF$Filter_sum_following_two,] # Filter out the variants with less than the sum of the following two bases

# plot the VAF distribution
p1 <- ggplot(ampliconseq_VAF_pass, aes(x = `Amplicon-seq VAF`)) +
  geom_histogram(bins = 50) +
  labs(
    y = "Frequency"
  ) +
  theme_bw()

p2 <- ggplot(ampliconseq_VAF_pass, aes(x = `Duplex-Multiome VAF`)) +
  geom_histogram(bins = 50) +
  labs(
    y = "Frequency"
  ) +
  theme_bw()

ggpubr::ggarrange(p1, p2)
ggsave("./figures/amplicon_seq/amplicon_seq_VAF.pdf")

# read in the clonal mutation list
clonal_mutation <- read.csv("./results/clonal_mutation/clonal_mutation_list_healthy_brain_realigned_ATAC_celltype_20241028.csv")
clonal_mutation$vaf <- clonal_mutation$detected/clonal_mutation$covered

clonal_mutation <- clonal_mutation[!duplicated(clonal_mutation$key),]
# plot the VAF distribution
ggplot(clonal_mutation, aes(x = case_id, y = vaf, color = case_id)) +
  geom_violin() +
  geom_jitter() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x = "Sample", y = "Duplex-Multiome VAF") +
  guides(color = FALSE)
ggsave("./figures/amplicon_seq/duplex_multiome_VAF.pdf")

