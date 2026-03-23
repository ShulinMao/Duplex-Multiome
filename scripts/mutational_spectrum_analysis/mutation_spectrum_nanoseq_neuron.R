# Libraries
.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(pheatmap)
library(ggpubr)
library(scales)
library(reshape2)
library(stringi)

ref_genome="BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = T)
chr_orders <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages_sig/")
library(MutationalPatterns)

# font
font = "Helvetica"
font_size = 7
axis_font_size = 6

# Nanoseq
nanoseq_neuron <- read.csv("./data/others/nanoseq_trinucleotide_context.csv")
nanoseq_neuron <- nanoseq_neuron[!str_detect(nanoseq_neuron$SampleID, "AD"),] # remove AD cells
nanoseq_neuron <- nanoseq_neuron[,-1]
nanoseq_neuron <- colSums(nanoseq_neuron)

# Duples Multiome
spectrum_duplex_OL_EN <- readRDS("./results/spectrum/corrected_mutation_spectrum_a6s3_20240627.rds")
mut_mat_EN <- cbind(spectrum_duplex_OL_EN[["EN"]], nanoseq_neuron)
colnames(mut_mat_EN) <- c("Duplex-Multiome", "NanoSeq")

cos_sim_nanoseq <- round(cos_sim(mut_mat_EN[,1], mut_mat_EN[,2]), 2)

plot_96_profile(mut_mat_EN[, "NanoSeq", drop=F], ymax = 0.075, condensed = T) + 
  scale_y_continuous(breaks = seq(0, 0.075, 0.025)) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
  ) +
  labs(title = "Neuron (NanoSeq)", 
       caption = paste0("Cosine similarity to Duplex-Multiome excitatory neuron spectrum: ", cos_sim_nanoseq))
ggsave("./figures/mutation_spectrum/nanoseq_neuron.pdf", 
       width = 7, height = 2, units = "in")

