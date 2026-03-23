# Libraries
.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(stringr)
library(tidyr)
library(rngtools)
library(stringr)
library(ggpubr)

ref_genome="BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = T)
chr_orders <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages_sig/")
library(MutationalPatterns)

case_id <- "UMB1932"
a <- 2

read.table(paste0("./data/UMB1932/single_cell/per_barcode_metrics.cells"), sep = "\t") -> selected_cells

mut_mat_summary <- c()

# a2s1
s <- a%/%2
stringency <- paste0("a", a, "s", s)

normalized_spectrum <- rep(0, 96)
original_spectrum <- rep(0, 96)
# read in normalized spectrum
for (cell_id in selected_cells$V1){
  specturm_cell <-
    read.table(paste0("./data/", case_id, "/single_cell/burden_estimation/Dan/", stringency, "/", stringency, "_filterGermline_burden.", cell_id, ".txt"), 
               skip = 43, comment.char = "", header = T)
  specturm_cell$OriginalCount[is.na(specturm_cell$OriginalCount)] <- 0
  specturm_cell$Adjusted[is.na(specturm_cell$AdjustedCount)] <- 0
  original_spectrum <- original_spectrum + specturm_cell$OriginalCount
  normalized_spectrum <- normalized_spectrum + specturm_cell$AdjustedCount
  
  mut_mat <- matrix(c(normalized_spectrum, original_spectrum), ncol = 2, byrow = F)
  colnames(mut_mat) <- c(paste(case_id, "Normalized", stringency, sep="_"), paste(case_id, "Original", stringency, sep="_"))
  row.names(mut_mat) <- 
    paste0(substring(specturm_cell$Context, 1, 1), 
           "[", specturm_cell$Variant, "]", 
           substring(specturm_cell$Context, 3, 3))
}
mut_mat_summary <- cbind(mut_mat_summary, mut_mat)

# a2s0
s <- 0
stringency <- paste0("a", a, "s", s)

normalized_spectrum <- rep(0, 96)
original_spectrum <- rep(0, 96)
# read in normalized spectrum
for (cell_id in selected_cells$V1){
  specturm_cell <-
    read.table(paste0("./data/", case_id, "/single_cell/burden_estimation/Dan/", stringency, "/", stringency, "_filterGermline_burden.", cell_id, ".txt"), 
               skip = 43, comment.char = "", header = T)
  specturm_cell$OriginalCount[is.na(specturm_cell$OriginalCount)] <- 0
  specturm_cell$Adjusted[is.na(specturm_cell$AdjustedCount)] <- 0
  original_spectrum <- original_spectrum + specturm_cell$OriginalCount
  normalized_spectrum <- normalized_spectrum + specturm_cell$AdjustedCount
  
  mut_mat <- matrix(c(normalized_spectrum, original_spectrum), ncol = 2, byrow = F)
  colnames(mut_mat) <- c(paste(case_id, "Normalized", stringency, sep="_"), paste(case_id, "Original", stringency, sep="_"))
  row.names(mut_mat) <- 
    paste0(substring(specturm_cell$Context, 1, 1), 
           "[", specturm_cell$Variant, "]", 
           substring(specturm_cell$Context, 3, 3))
}
mut_mat_summary <- cbind(mut_mat_summary, mut_mat)
mut_mat_summary <- mut_mat_summary[, c("UMB1932_Original_a2s1", "UMB1932_Original_a2s0")]

saveRDS(mut_mat_summary, "./figures/manuscript_figures/Figure1/UMB1932_mutation_spectrum.rds")


get_known_signatures(incl_poss_artifacts = T)[,"SBS45", drop=F] -> SBS45
rownames(SBS45) <- row.names(mut_mat_summary)
colnames(SBS45) <- "COSMIC\n"
cos_sim_SBS45 <- cos_sim((mut_mat_summary[, "UMB1932_Original_a2s0"]/sum(mut_mat_summary[, "UMB1932_Original_a2s0"]) - 
                            mut_mat_summary[, "UMB1932_Original_a2s1"]/sum(mut_mat_summary[, "UMB1932_Original_a2s1"])), 
                         SBS45[,1])

# font
font = "Helvetica"
font_size = 7
axis_font_size = 6
plot_96_profile(SBS45, ymax = 0.25, condensed = T) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_blank(),
    axis.ticks.x=element_blank(), 
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    panel.spacing.y = unit(1, "lines")
  ) +
  labs(title = "SBS45")
ggsave("./figures/mutation_spectrum/SBS45.pdf", width = 3, height = 1.3, units = "in")
plot_compare_profiles(mut_mat_summary[, "UMB1932_Original_a2s0"],
                      mut_mat_summary[, "UMB1932_Original_a2s1"],
                      profile_names = c("w/ duplex\nconcensus", "w/o duplex\nconcensus"),
                      condensed = T, profile_ymax = 0.06) + 
  labs(title = case_id, 
       caption = paste0("Cosine similarity to COSMIC SBS45: ", round(cos_sim_SBS45, 2))) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_blank(),
    axis.ticks.x=element_blank(), 
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    panel.spacing.y = unit(0.2, "lines")
  ) 
ggsave("./figures/mutation_spectrum/UMB1932_SBS45.pdf", height = 2.7, width = 3, units = "in")

