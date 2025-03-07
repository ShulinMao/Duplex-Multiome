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

ref_genome="BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = T)
chr_orders <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages_sig/")
library(MutationalPatterns)


# font
font = "Helvetica"
font_size = 7
axis_font_size = 6

# ------------------------------------------------------------------------------
# Spectrum
# PTA
spectrum_PTA_OL_EN <- readRDS("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/figures/manuscript_figures/Figure4/spectrum_PTA_OL_EN.rds")
# Duples Multiome
spectrum_duplex <- readRDS("./results/spectrum/corrected_mutation_spectrum_a6s3_20240627.rds")

mut_mat_EN <- cbind(spectrum_PTA_OL_EN[, "EN_PTA"], spectrum_duplex[["EN"]])
mut_mat_OL <- cbind(spectrum_PTA_OL_EN[, "OL_PTA"], spectrum_duplex[["Oligodendrocyte"]])

colnames(mut_mat_EN) <- c("PTA", "Duplex-Multiome")
colnames(mut_mat_OL) <- c("PTA", "Duplex-Multiome")

p1 <- plot_96_profile(mut_mat_EN, ymax = 0.075, condensed = T) + 
  scale_y_continuous(breaks = seq(0, 0.075, 0.025))  +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_blank(),
    axis.ticks.x=element_blank(), 
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    panel.spacing.y = unit(1, "lines")
  )+
  labs(title = "Excitatory Neuron (EN)")

p2 <- plot_96_profile(mut_mat_OL, ymax = 0.075, condensed = T) + 
  scale_y_continuous(breaks = seq(0, 0.075, 0.025))  +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_blank(),
    axis.ticks.x=element_blank(), 
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    panel.spacing.y = unit(1, "lines")
  ) +
  labs(title = "Oligodendrocyte (OL)")

ggarrange(p1, p2, ncol = 2)
ggsave("./figures/manuscript_figures/Figure4/figure4_EN_OL_spectrum.pdf", 
       width = 7, height = 2.7, units = "in")

# similarity
cos_sim_EN <- round(cos_sim(mut_mat_EN[, 1], mut_mat_EN[, 2]), 2)
cos_sim_OL <- round(cos_sim(mut_mat_OL[, 1], mut_mat_OL[, 2]), 2)
df_cos_sim <- data.frame("cos_sim" = c(cos_sim_EN, cos_sim_OL),
           "cell_type" = c("EN", "OL"))

p3 <- ggplot(df_cos_sim, aes(x = cell_type, y = cos_sim, fill = cell_type)) + 
  geom_col(width = 0.5) +
  geom_text(aes(label = cos_sim, y = cos_sim + 0.05), size = 2.5) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("#7CAE00", "#00A9FF")) +
  labs(x = "", y = "Cosine Similarity")

# ------------------------------------------------------------------------------
# Signature
# Contribution
sbs96_selected <- get_known_signatures()[,c("SBS5", "SBS16", "SBS1", "SBS32", "SBS19")]

df_spectrum_duplex_OL_EN <- cbind(spectrum_duplex[["EN"]], spectrum_duplex[["Oligodendrocyte"]])
colnames(df_spectrum_duplex_OL_EN) <- c("EN", "OL")
fit <- fit_to_signatures(df_spectrum_duplex_OL_EN, sbs96_selected)
duplex_OL_EN_signature_contribution <- fit$contribution
colnames(duplex_OL_EN_signature_contribution) <- c("EN", "Oligodendrocyte")
p4 <- plot_contribution(fit$contribution, coord_flip = T, mode = "relative") +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 1, 0.5))  +
  guides(fill=guide_legend(nrow=3, byrow=F)) +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = 8, family = font),
    legend.position = "bottom",
    legend.box.spacing = unit(0, "in"),
    legend.key.size = unit(0.1, "in")
  ) +
  labs(fill = "")

# Aging trends
# EN + OL adjusted burdens
simulation_results_EN_Oligo <- readRDS("./results/burden_simulation/burden_summary_EN_OL_20241018.rds")
(simulation_results_EN_Oligo$sd/simulation_results_EN_Oligo$coverage)*5.75e9 -> simulation_results_EN_Oligo$burden_sd
simulation_results_EN_Oligo$burden_sd * 1.96 -> simulation_results_EN_Oligo$burden_95CI

burden_regression <- readRDS("./figures/manuscript_figures/Figure3/burden_regression_results.rds")

duplex_OL_EN_signature_contribution <- cbind(duplex_OL_EN_signature_contribution[,1]/colSums(duplex_OL_EN_signature_contribution)[1],
                                             duplex_OL_EN_signature_contribution[,2]/colSums(duplex_OL_EN_signature_contribution)[2])
colnames(duplex_OL_EN_signature_contribution) <- c("EN", "Oligodendrocyte")
duplex_OL_EN_signature_contribution <- melt(duplex_OL_EN_signature_contribution)
colnames(duplex_OL_EN_signature_contribution) <- c("signature", "cell_type", "contribution")

simulation_results_EN_Oligo_signature <- left_join(simulation_results_EN_Oligo, 
                                         duplex_OL_EN_signature_contribution, 
                                         "cell_type")
simulation_results_EN_Oligo_signature[, "signature_burden_corrected"] <-
  simulation_results_EN_Oligo_signature$contribution * simulation_results_EN_Oligo_signature$burden_corrected
simulation_results_EN_Oligo_signature[, "signature_burden_95CI"] <-
  simulation_results_EN_Oligo_signature$contribution * simulation_results_EN_Oligo_signature$burden_95CI


simulation_results_EN_Oligo_signature[,"burden_shade_mean"] <- NA
simulation_results_EN_Oligo_signature[simulation_results_EN_Oligo_signature$cell_type == "EN", "burden_shade_mean"] <-
  burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$EN[["mean_slope"]] * simulation_results_EN_Oligo_signature$Age[simulation_results_EN_Oligo_signature$cell_type == "EN"] + burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$EN[["mean_intercept"]]
simulation_results_EN_Oligo_signature[simulation_results_EN_Oligo_signature$cell_type == "Oligodendrocyte", "burden_shade_mean"] <-
  burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$Oligodendrocyte[["mean_slope"]] * simulation_results_EN_Oligo_signature$Age[simulation_results_EN_Oligo_signature$cell_type == "Oligodendrocyte"] + burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$Oligodendrocyte[["mean_intercept"]]
simulation_results_EN_Oligo_signature[,"signature_burden_shade_mean"] <-
  simulation_results_EN_Oligo_signature$burden_shade_mean * simulation_results_EN_Oligo_signature$contribution 

simulation_results_EN_Oligo_signature[,"burden_shade_sd"] <- NA
simulation_results_EN_Oligo_signature[simulation_results_EN_Oligo_signature$cell_type == "EN", "burden_shade_sd"] <-
  sqrt((burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$EN[["sd_slope"]] * simulation_results_EN_Oligo_signature$Age[simulation_results_EN_Oligo_signature$cell_type == "EN"])^2 + (burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$EN[["sd_intercept"]])^2) * 1.96
simulation_results_EN_Oligo_signature[simulation_results_EN_Oligo_signature$cell_type == "Oligodendrocyte", "burden_shade_sd"] <-
  sqrt((burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$Oligodendrocyte[["sd_slope"]] * simulation_results_EN_Oligo_signature$Age[simulation_results_EN_Oligo_signature$cell_type == "Oligodendrocyte"])^2 + (burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$Oligodendrocyte[["sd_intercept"]])^2) * 1.96
simulation_results_EN_Oligo_signature[,"signature_burden_shade_sd"] <-
  simulation_results_EN_Oligo_signature$burden_shade_sd * simulation_results_EN_Oligo_signature$contribution 

# signature enrichment
signature_enrichment_analysis_OL <- readRDS("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/figures/manuscript_figures/Figure4/signature_enrichment_analysis_Oligodendrocyte_a6s3.rds")
signature_enrichment_analysis_EN <- readRDS("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/figures/manuscript_figures/Figure4/signature_enrichment_analysis_EN_a6s3.rds")

context_3mer_ratio <- readRDS(paste0("./results/context_3mer_ratio/a6s3_celltype.rds"))
context_3mer_ratio <- context_3mer_ratio[c(1:16, 1:16, 1:16, 17:32, 17:32, 17:32),]

signature_enrichment_analysis <- cbind(signature_enrichment_analysis_OL$covered_region/context_3mer_ratio[,"Oligodendrocyte"],
                                       signature_enrichment_analysis_OL$all,
                                       signature_enrichment_analysis_EN$covered_region/context_3mer_ratio[,"EN"],
                                       signature_enrichment_analysis_EN$all)
colnames(signature_enrichment_analysis) <- c("OL_covered_region", "OL_all",
                                             "EN_covered_region", "EN_all")
signature_contribution <- 
  fit_to_signatures(signature_enrichment_analysis[, c(2, 4)], sbs96_selected)$contribution
signature_contribution[,1] <- signature_contribution[,1]/colSums(signature_enrichment_analysis)[2]
signature_contribution[,2] <- signature_contribution[,2]/colSums(signature_enrichment_analysis)[4]

# SBS5
signature_figure <- "SBS5"
signature_enrichment_analysis_contribution <- 
  fit_to_signatures(signature_enrichment_analysis, sbs96_selected[, signature_figure, drop = F])$contribution
signature_enrichment_analysis_contribution <- signature_enrichment_analysis_contribution/colSums(signature_enrichment_analysis)

signature_enrichment_ratio <- cbind(signature_enrichment_analysis_contribution[,"EN_covered_region"]/signature_enrichment_analysis_contribution[,"EN_all"],
                                    signature_enrichment_analysis_contribution[,"OL_covered_region"]/signature_enrichment_analysis_contribution[,"OL_all"])
colnames(signature_enrichment_ratio) <- c("EN", "OL")
row.names(signature_enrichment_ratio) <- signature_figure
signature_enrichment_ratio <- melt(signature_enrichment_ratio)
colnames(signature_enrichment_ratio) <- c("signature", "cell_type", "enrichment_ratio")

p5 <- ggplot(simulation_results_EN_Oligo_signature[simulation_results_EN_Oligo_signature$signature == signature_figure, ]) +
  geom_point(size = 0.2, aes(x = Age, y = signature_burden_corrected, color = cell_type)) +
  # geom_errorbar(aes(ymax = signature_burden_corrected + signature_burden_95CI,
  #                   ymin = signature_burden_corrected - signature_burden_95CI), 
  #               size = 0.2, width = 4) +
  # geom_smooth(method = "lm", se = F, show.legend = F, size = 0.5) +
  geom_abline(slope = burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$EN["mean_slope"]*
                duplex_OL_EN_signature_contribution$contribution[duplex_OL_EN_signature_contribution$signature == signature_figure & duplex_OL_EN_signature_contribution$cell_type == "EN"], 
              intercept = burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$EN["mean_intercept"]*
                duplex_OL_EN_signature_contribution$contribution[duplex_OL_EN_signature_contribution$signature == signature_figure & duplex_OL_EN_signature_contribution$cell_type == "EN"], 
              size = 0.5, color = "#7CAE00") +
  geom_abline(slope = burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$Oligodendrocyte["mean_slope"]*
                duplex_OL_EN_signature_contribution$contribution[duplex_OL_EN_signature_contribution$signature == signature_figure & duplex_OL_EN_signature_contribution$cell_type == "Oligodendrocyte"], 
              intercept = burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$Oligodendrocyte["mean_intercept"]*
                duplex_OL_EN_signature_contribution$contribution[duplex_OL_EN_signature_contribution$signature == signature_figure & duplex_OL_EN_signature_contribution$cell_type == "Oligodendrocyte"], 
              size = 0.5, color = "#00A9FF") +
  geom_ribbon(aes(x = Age,
                  ymin = signature_burden_shade_mean-signature_burden_shade_sd, 
                  ymax = signature_burden_shade_mean+signature_burden_shade_sd,
                  fill = cell_type),
              alpha = 0.1, inherit.aes = F) +
  geom_abline(slope = 14.5, intercept = 107*signature_contribution[signature_figure, "EN_all"],
              linetype = 2, size = 0.5, color = "#7CAE00") +
  geom_abline(slope = 22.7, intercept = 165*signature_contribution[signature_figure, "OL_all"],
              linetype = 2, size = 0.5, color = "#00A9FF") +
  scale_color_manual(values = c("#7CAE00", "#00A9FF")) +
  scale_fill_manual(values = c("#7CAE00", "#00A9FF")) +
  guides(color=guide_legend(nrow=2), override.aes = list(size=1)) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = 6, family = font),
    legend.position = "none"
  ) +
  labs(x = "", y = "Signature exposure\n(per cell)", subtitle = signature_figure)

p6 <- ggplot(signature_enrichment_ratio[signature_enrichment_ratio$signature == signature_figure,], 
       aes(x = cell_type, y = enrichment_ratio, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_classic() + 
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = 6, family = font),
    legend.position = "none"
  ) +
  scale_y_continuous(breaks = c(0, 0.4, 0.8, 1, 1.2), limits = c(0, 1.2)) +
  scale_fill_manual(values = c("#7CAE00", "#00A9FF")) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "Signature enrichment ratio\n(covered region / whole genome)") 

ggarrange(p5, p6, nrow = 2)

# SBS1
signature_figure <- "SBS1"
signature_enrichment_analysis_contribution <- 
  fit_to_signatures(signature_enrichment_analysis, sbs96_selected[, signature_figure, drop = F])$contribution
signature_enrichment_analysis_contribution <- signature_enrichment_analysis_contribution/colSums(signature_enrichment_analysis)

signature_enrichment_ratio <- cbind(signature_enrichment_analysis_contribution[,"EN_covered_region"]/signature_enrichment_analysis_contribution[,"EN_all"],
                                    signature_enrichment_analysis_contribution[,"OL_covered_region"]/signature_enrichment_analysis_contribution[,"OL_all"])
colnames(signature_enrichment_ratio) <- c("EN", "OL")
row.names(signature_enrichment_ratio) <- signature_figure
signature_enrichment_ratio <- melt(signature_enrichment_ratio)
colnames(signature_enrichment_ratio) <- c("signature", "cell_type", "enrichment_ratio")

p7 <- ggplot(simulation_results_EN_Oligo_signature[simulation_results_EN_Oligo_signature$signature == signature_figure, ], 
       aes(x = Age, y = signature_burden_corrected, color = cell_type)) +
  geom_point(size = 0.2) +
  # geom_errorbar(aes(ymax = signature_burden_corrected + signature_burden_95CI,
  #                   ymin = signature_burden_corrected - signature_burden_95CI), 
  #               size = 0.2, width = 4) +
  # geom_smooth(method = "lm", se = F, show.legend = F, size = 0.5) +
  geom_abline(slope = burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$EN["mean_slope"]*
                duplex_OL_EN_signature_contribution$contribution[duplex_OL_EN_signature_contribution$signature == signature_figure & duplex_OL_EN_signature_contribution$cell_type == "EN"], 
              intercept = burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$EN["mean_intercept"]*
                duplex_OL_EN_signature_contribution$contribution[duplex_OL_EN_signature_contribution$signature == signature_figure & duplex_OL_EN_signature_contribution$cell_type == "EN"], 
              size = 0.5, color = "#7CAE00") +
  geom_abline(slope = burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$Oligodendrocyte["mean_slope"]*
                duplex_OL_EN_signature_contribution$contribution[duplex_OL_EN_signature_contribution$signature == signature_figure & duplex_OL_EN_signature_contribution$cell_type == "Oligodendrocyte"], 
              intercept = burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$Oligodendrocyte["mean_intercept"]*
                duplex_OL_EN_signature_contribution$contribution[duplex_OL_EN_signature_contribution$signature == signature_figure & duplex_OL_EN_signature_contribution$cell_type == "Oligodendrocyte"], 
              size = 0.5, color = "#00A9FF") +
  geom_ribbon(aes(x = Age,
                  ymin = signature_burden_shade_mean-signature_burden_shade_sd, 
                  ymax = signature_burden_shade_mean+signature_burden_shade_sd,
                  fill = cell_type),
              alpha = 0.1, inherit.aes = F) +
  geom_abline(slope = 0.29, intercept = 107*signature_contribution[signature_figure, "EN_all"],
              linetype = 2, size = 0.5, color = "#7CAE00") +
  geom_abline(slope = 2.77, intercept = 165*signature_contribution[signature_figure, "OL_all"],
              linetype = 2, size = 0.5, color = "#00A9FF") +
  scale_color_manual(values = c("#7CAE00", "#00A9FF")) +
  scale_fill_manual(values = c("#7CAE00", "#00A9FF")) +
  guides(color=guide_legend(nrow=2), override.aes = list(size=1)) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = 6, family = font),
    legend.position = "none"
  ) +
  labs(x = "Age (yrs)", y = "", subtitle = signature_figure)

p8 <- ggplot(signature_enrichment_ratio[signature_enrichment_ratio$signature == signature_figure,], 
             aes(x = cell_type, y = enrichment_ratio, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_classic() + 
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = 6, family = font),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("#7CAE00", "#00A9FF")) +
  scale_y_continuous(breaks = c(0, 0.4, 0.8, 1, 1.2), limits = c(0, 1.2)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "") 

ggarrange(p7, p8, nrow = 2)

# SBS16
signature_figure <- "SBS16"
signature_enrichment_analysis_contribution <- 
  fit_to_signatures(signature_enrichment_analysis, sbs96_selected[, signature_figure, drop = F])$contribution
signature_enrichment_analysis_contribution <- signature_enrichment_analysis_contribution/colSums(signature_enrichment_analysis)

signature_enrichment_ratio <- cbind(signature_enrichment_analysis_contribution[,"EN_covered_region"]/signature_enrichment_analysis_contribution[,"EN_all"],
                                    #NA)
                                    signature_enrichment_analysis_contribution[,"OL_covered_region"]/signature_enrichment_analysis_contribution[,"OL_all"])
colnames(signature_enrichment_ratio) <- c("EN", "OL")
row.names(signature_enrichment_ratio) <- signature_figure
signature_enrichment_ratio <- melt(signature_enrichment_ratio)
colnames(signature_enrichment_ratio) <- c("signature", "cell_type", "enrichment_ratio")

p9 <- ggplot(simulation_results_EN_Oligo_signature[simulation_results_EN_Oligo_signature$signature == signature_figure, ], 
       aes(x = Age, y = signature_burden_corrected, color = cell_type)) +
  geom_point(size = 0.2) +
  # geom_errorbar(aes(ymax = signature_burden_corrected + signature_burden_95CI,
  #                   ymin = signature_burden_corrected - signature_burden_95CI),
  #               size = 0.2, width = 4) +
  # geom_smooth(method = "lm", se = F, show.legend = F, size = 0.5) +
  geom_abline(slope = burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$EN["mean_slope"]*
                duplex_OL_EN_signature_contribution$contribution[duplex_OL_EN_signature_contribution$signature == signature_figure & duplex_OL_EN_signature_contribution$cell_type == "EN"],
              intercept = burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$EN["mean_intercept"]*
                duplex_OL_EN_signature_contribution$contribution[duplex_OL_EN_signature_contribution$signature == signature_figure & duplex_OL_EN_signature_contribution$cell_type == "EN"],
              size = 0.5, color = "#7CAE00") +
  geom_abline(slope = burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$Oligodendrocyte["mean_slope"]*
                duplex_OL_EN_signature_contribution$contribution[duplex_OL_EN_signature_contribution$signature == signature_figure & duplex_OL_EN_signature_contribution$cell_type == "Oligodendrocyte"],
              intercept = burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$Oligodendrocyte["mean_intercept"]*
                duplex_OL_EN_signature_contribution$contribution[duplex_OL_EN_signature_contribution$signature == signature_figure & duplex_OL_EN_signature_contribution$cell_type == "Oligodendrocyte"],
              size = 0.5, color = "#00A9FF") +
  geom_ribbon(aes(x = Age,
                  ymin = signature_burden_shade_mean-signature_burden_shade_sd,
                  ymax = signature_burden_shade_mean+signature_burden_shade_sd,
                  fill = cell_type), alpha = 0.1, inherit.aes = F) +
  geom_abline(slope = 2.0, intercept = 107*signature_contribution[signature_figure, "EN_all"],
              linetype = 2, size = 0.5, color = "#7CAE00") +
  geom_abline(slope = 0.18, intercept = 165*signature_contribution[signature_figure, "OL_all"],
              linetype = 2, size = 0.5, color = "#00A9FF") +
  scale_color_manual(values = c("#7CAE00", "#00A9FF")) +
  scale_fill_manual(values = c("#7CAE00", "#00A9FF")) +
  guides(color=guide_legend(nrow=2), override.aes = list(size=1)) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = 6, family = font),
    legend.position = "none"
  ) +
  labs(x = "", y = "", subtitle = signature_figure)

p10 <- ggplot(signature_enrichment_ratio[signature_enrichment_ratio$signature == signature_figure,], 
             aes(x = cell_type, y = enrichment_ratio, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_classic() + 
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = 6, family = font),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("#7CAE00", "#00A9FF")) +
  scale_y_continuous(breaks = c(0, 0.4, 0.8, 1, 1.2), limits = c(0, 1.2))  +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "") 

ggarrange(p9, p10, nrow = 2)

# SBS32
signature_figure <- "SBS32"
signature_enrichment_analysis_contribution <- 
  fit_to_signatures(signature_enrichment_analysis, sbs96_selected[, signature_figure, drop = F])$contribution
signature_enrichment_analysis_contribution <- signature_enrichment_analysis_contribution/colSums(signature_enrichment_analysis)

signature_enrichment_ratio <- cbind(#signature_enrichment_analysis_contribution[,"EN_covered_region"]/signature_enrichment_analysis_contribution[,"EN_all"],
                                    0,
                                    signature_enrichment_analysis_contribution[,"OL_covered_region"]/signature_enrichment_analysis_contribution[,"OL_all"])
colnames(signature_enrichment_ratio) <- c("EN", "OL")
row.names(signature_enrichment_ratio) <- signature_figure
signature_enrichment_ratio <- melt(signature_enrichment_ratio)
colnames(signature_enrichment_ratio) <- c("signature", "cell_type", "enrichment_ratio")

p11 <- ggplot(simulation_results_EN_Oligo_signature[simulation_results_EN_Oligo_signature$signature == signature_figure, ], 
             aes(x = Age, y = signature_burden_corrected, color = cell_type)) +
  geom_point(size = 0.2) +
  # geom_errorbar(aes(ymax = signature_burden_corrected + signature_burden_95CI,
  #                   ymin = signature_burden_corrected - signature_burden_95CI), 
  #               size = 0.2, width = 4) +
  # geom_smooth(method = "lm", se = F, show.legend = F, size = 0.5) +
  geom_abline(slope = burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$EN["mean_slope"]*
                duplex_OL_EN_signature_contribution$contribution[duplex_OL_EN_signature_contribution$signature == signature_figure & duplex_OL_EN_signature_contribution$cell_type == "EN"],
              intercept = burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$EN["mean_intercept"]*
                duplex_OL_EN_signature_contribution$contribution[duplex_OL_EN_signature_contribution$signature == signature_figure & duplex_OL_EN_signature_contribution$cell_type == "EN"],
              size = 0.5, color = "#7CAE00") +
  geom_abline(slope = burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$Oligodendrocyte["mean_slope"]*
                duplex_OL_EN_signature_contribution$contribution[duplex_OL_EN_signature_contribution$signature == signature_figure & duplex_OL_EN_signature_contribution$cell_type == "Oligodendrocyte"],
              intercept = burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$Oligodendrocyte["mean_intercept"]*
                duplex_OL_EN_signature_contribution$contribution[duplex_OL_EN_signature_contribution$signature == signature_figure & duplex_OL_EN_signature_contribution$cell_type == "Oligodendrocyte"],
              size = 0.5, color = "#00A9FF") +
  geom_ribbon(aes(x = Age,
                  ymin = signature_burden_shade_mean-signature_burden_shade_sd,
                  ymax = signature_burden_shade_mean+signature_burden_shade_sd,
                  fill = cell_type), alpha = 0.1, inherit.aes = F) +
  geom_abline(slope = 0, intercept = 107*signature_contribution[signature_figure, "EN_all"],
              linetype = 2, size = 0.5, color = "#7CAE00") +
  geom_abline(slope = 2.5, intercept = 165*signature_contribution[signature_figure, "OL_all"],
              linetype = 2, size = 0.5, color = "#00A9FF") +
  scale_color_manual(values = c("#7CAE00", "#00A9FF")) +
  guides(color=guide_legend(nrow=2), override.aes = list(size=1)) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = 6, family = font),
    legend.position = "none"
  ) +
  labs(x = "Age (yrs)", y = "", subtitle = signature_figure)

p12 <- ggplot(signature_enrichment_ratio[signature_enrichment_ratio$signature == signature_figure,], 
             aes(x = cell_type, y = enrichment_ratio, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_classic() + 
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = 6, family = font),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("#7CAE00", "#00A9FF")) +
  scale_y_continuous(breaks = seq(0, 1.2, 0.2), limits = c(0, 1.2)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "") 

ggarrange(p11, p12, nrow = 2)

ggarrange(plotlist = list(p3, p5, p7, p9, p4, p6, p8, p10), 
                              nrow = 2, ncol = 4, widths = c(1.05, 1.11, 1, 1))
ggsave("./figures/manuscript_figures/Figure4/figure4_signature.pdf", 
       width = 7, height = 3, units = "in")

# ------------------------------------------------------------------------------
# Spectrum (Astrocyte & Microglia)
spectrum_duplex <- readRDS("./results/spectrum/corrected_mutation_spectrum_a6s3_20240627.rds")

mut_mat_Astro <- spectrum_duplex[["Astrocyte"]]
mut_mat_Microglia <- spectrum_duplex[["Microglia"]]
mut_mat_IN <- spectrum_duplex[["IN"]]
mut_mat_OPC <- spectrum_duplex[["OPC"]]

colnames(mut_mat_Astro) <- c("Duplex-Multiome")
colnames(mut_mat_Microglia) <- c("Duplex-Multiome")
colnames(mut_mat_IN) <- c("Duplex-Multiome")
colnames(mut_mat_OPC) <- c("Duplex-Multiome")

p13 <- plot_96_profile(mut_mat_Astro, ymax = 0.1, condensed = T) + 
  scale_y_continuous(breaks = seq(0, 0.1, 0.05))  +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_blank(),
    axis.ticks.x=element_blank(), 
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    panel.spacing.y = unit(1, "lines")
  )+
  labs(title = "Astrocyte")

p14 <- plot_96_profile(mut_mat_Microglia, ymax = 0.1, condensed = T) + 
  scale_y_continuous(breaks = seq(0, 0.1, 0.05))  +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
  ) +
  labs(title = "Microglia")

p15 <- plot_96_profile(mut_mat_IN, ymax = 0.075, condensed = T) + 
  scale_y_continuous(breaks = seq(0, 0.075, 0.025))  +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_blank(),
    axis.ticks.x=element_blank(), 
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    panel.spacing.y = unit(1, "lines")
  ) +
  labs(title = "Inhibitory Neuron")

p16 <- plot_96_profile(mut_mat_OPC, ymax = 0.12, condensed = T) + 
  scale_y_continuous(breaks = seq(0, 0.12, 0.04))  +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
  ) +
  labs(title = "Oligodendrocyte precursor cell (OPC)")

ggarrange(p15, p13, ncol = 2)
ggsave("./figures/manuscript_figures/Figure4/figure4_spectrum_Astrocyte_IN.pdf", 
       width = 7, height = 1.6, units = "in")

ggarrange(p14, p16, nrow = 2)
ggsave("./figures/manuscript_figures/Figure4/figure4_spectrum_Microglia_OPC.pdf", 
       width = 7, height = 3.7, units = "in")

