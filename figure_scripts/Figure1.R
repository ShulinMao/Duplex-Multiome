# Libraries
.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(dplyr)
library(ggplot2)
library(ggpubr)

ref_genome="BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = T)
chr_orders <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages_sig/")
library(MutationalPatterns)

# font
font = "Helvetica"
font_size = 7
axis_font_size = 6

# UMB1932 specturm vs SBS45
readRDS("./figures/manuscript_figures/Figure1/UMB1932_mutation_spectrum.rds") -> mut_mat_summary

get_known_signatures(incl_poss_artifacts = T)[,"SBS45", drop=F] -> SBS45
rownames(SBS45) <- row.names(mut_mat_summary)
colnames(SBS45) <- "COSMIC\nSBS45"
cos_sim_SBS45 <- cos_sim((mut_mat_summary[, "UMB1932_Original_a2s0"]/sum(mut_mat_summary[, "UMB1932_Original_a2s0"]) - 
                            mut_mat_summary[, "UMB1932_Original_a2s1"]/sum(mut_mat_summary[, "UMB1932_Original_a2s1"])), 
                         SBS45[,1])

p1 <- plot_compare_profiles(mut_mat_summary[, "UMB1932_Original_a2s0"],
                      mut_mat_summary[, "UMB1932_Original_a2s1"],
                      profile_names = c("w/o duplex\nconcensus", "w/ duplex\nconcensus"),
                      condensed = T, profile_ymax = 0.06, diff_ylim = c(-0.015, 0.025)) + 
  labs(caption = paste0("Cosine similarity to COSMIC SBS45: ", round(cos_sim_SBS45, 2)),
       title = NULL, x = NULL, y = NULL) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_blank(),
    axis.ticks.x=element_blank(), 
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    panel.spacing.y = unit(0.2, "lines"),
    panel.border=element_blank(),
    plot.background=element_blank(),
    plot.margin=unit(c(0.1, 0.1, 0.1, 0), "cm")
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) 

p2 <- plot_96_profile(SBS45, ymax = 0.25, condensed = T) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), 
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    panel.spacing.y = unit(0.2, "lines"),
    plot.margin=unit(c(0.1, 0.1, 0.1, 0.1), "cm")
  ) +
  scale_y_continuous(breaks = seq(0, 0.2, 0.1), labels = scales::number_format(accuracy = 0.01)) +
  labs(y = NULL)
ggarrange(p1, p2, nrow = 2, heights = c(3, 1.3))

ggsave("./figures/manuscript_figures/Figure1/UMB1932_SBS45.pdf", height = 2.7, width = 3, units = "in")

# burden duplex vs non-duplex
burden_list_melt <- readRDS("./figures/manuscript_figures/Figure1/UMB1932_burden_duplex_vs_non_duplex.rds")
ggplot(burden_list_melt, aes(x=Var1, y=value, fill = duplex)) +
  geom_col() +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  scale_fill_simpsons() +
  guides(fill = "none") + 
  labs(x = "Stringency", y = "sSNV candidates per bp")
ggsave("./figures/manuscript_figures/Figure1/UMB1932_burden_duplex_vs_non_duplex.pdf", height = 2.25, width = 2, units = "in")

