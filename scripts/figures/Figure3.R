# Libraries
.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(stringr)
library(tidyr)
library(Seurat)
library(SeuratDisk)
library(Signac)
library(pheatmap)
library(igraph)
library(ggpubr)
library(scales)

# font
font = "Helvetica"
font_size = 7
axis_font_size = 6

# Custom function to format scientific notation with superscripts and handle zero
scientific_10 <- function(x) {
  ifelse(x == 0, "0", parse(text = gsub("e", " %*% 10^", scientific_format()(x))))
}

colors_cell_type <- c("#2E2585", "#337538", "#7CAE00", "#E69F00", "#DCCD7D", "#00A9FF", "#CC79A7", "#7E2954")

# ------------------------------------------------------------------------------
# UMAP of brain samples
brain <- readRDS("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/RNA_ATAC/multiomes_integrated_20240512_infant_adult_annotated.rds")

brain[,!str_detect(brain$celltype, "doublet")] -> brain
brain[["cell_type"]] <- as.character(brain$celltype)
brain$cell_type[str_detect(brain$celltype, "EN")] <- "Excitatory\nNeuron"
brain$cell_type[str_detect(brain$celltype, "IN")] <- "Inhibitory\nNeuron"

DimPlot(brain, reduction = "wnn.umap", 
        group.by = "cell_type", pt.size = 1e-10, 
        cols = colors_cell_type) + 
  labs(title = NULL, color = "Cell type") +
  guides(color=guide_legend(nrow=4, byrow=T, override.aes = list(size=2))) +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "bottom",
    #legend.key.size = unit(0, 'in'),
    legend.title = element_text(angle = 90, hjust = 0.5, vjust = 2.5)
  )
ggsave("./figures/manuscript_figures/Figure3/Brain_UMAP.pdf", height = 3.5, width = 2.3, units = "in")

p1 <- DimPlot(brain, reduction = "umap.rna", group.by = "cell_type", pt.size = 1e-10, 
              cols = colors_cell_type) + 
  labs(title = NULL, color = "Cell type") +
  guides(color=guide_legend(override.aes = list(size=2))) +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
  )

p2 <- DimPlot(brain, reduction = "umap.atac", group.by = "cell_type", pt.size = 1e-10, 
              cols = colors_cell_type) + 
  labs(title = NULL, color = "Cell type") +
  guides(color=guide_legend(override.aes = list(size=2))) +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
  )
ggarrange(p1,p2, ncol = 2, common.legend = T, legend = "bottom")
ggsave("./figures/manuscript_figures/Figure3/Brain_RNA_ATAC.pdf", height = 3, width = 5, units = "in")

# ------------------------------------------------------------------------------
# EN + OL adjusted burdens
simulation_results_EN_Oligo <- readRDS("./results/burden_simulation/burden_summary_EN_OL_20241018.rds")
(simulation_results_EN_Oligo$sd/simulation_results_EN_Oligo$coverage)*5.75e9 -> simulation_results_EN_Oligo$burden_sd
simulation_results_EN_Oligo$burden_sd * 1.96 -> simulation_results_EN_Oligo$burden_95CI
simulation_results_EN_Oligo$cell_type[simulation_results_EN_Oligo$cell_type == "EN"] <- "Excitatory Neuron"

burden_regression <- readRDS("./figures/manuscript_figures/Figure3/burden_regression_results.rds")

ggplot(simulation_results_EN_Oligo, aes(x = Age, y = burden_corrected, color = cell_type)) +
  geom_point(size = 0.2) +
  geom_errorbar(aes(ymax = burden_corrected + burden_95CI, ymin = burden_corrected - burden_95CI), size = 0.2, width = 4) +
  geom_abline(slope = burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$EN["mean_slope"], 
              intercept = burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$EN["mean_intercept"], 
              size = 0.5, color = "#7CAE00") +
  geom_abline(slope = burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$Oligodendrocyte["mean_slope"], 
              intercept = burden_regression$EN_OLigodendrocyte_corrected_burden_regression_results$Oligodendrocyte["mean_intercept"], 
              size = 0.5, color = "#00A9FF") +
  geom_abline(slope = 16, intercept = 107, linetype = 2, size = 0.5, color = "#7CAE00") +
  geom_abline(slope = 29, intercept = 165, linetype = 2, size = 0.5, color = "#00A9FF") +
  #annotate(geom="text", x=20, y=2750, parse = T, label = "'s=32.5' * 'cell'^{-1} * 'year'^{-1}", color = "#00A9FF", family = font, size = 2) +
  #annotate(geom="text", x=20, y=3000, parse = T, label = "'s=15.3' * 'cell'^{-1} * 'year'^{-1}", color = "#7CAE00", family = font, size = 2) +
  annotate(geom="text", x=20, y=3000, parse = T, 
           label = paste0("italic(P) == ", scientific_10(burden_regression$EN_OLigodendrocyte_corrected_burden_p_value)), 
           color = "black", family = font, size = 2) +
  scale_color_manual(values = c("#7CAE00", "#00A9FF")) +
  guides(color=guide_legend(nrow=2), override.aes = list(size=1)) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = 6, family = font),
    legend.position = c(0.7, 0.12),
    legend.key.size = unit(0.1, 'in'),
    legend.title = element_blank(),
    legend.background = element_rect(fill='transparent'),
  ) +
  scale_y_continuous(labels = scientific_10) +
  labs(x = "Age (yrs)", y = "Mutational burden (per cell)")
ggsave("./figures/manuscript_figures/Figure3/OL_EN_burden_adjusted_enrichment.pdf", width = 2, height = 2, units = "in")

# ------------------------------------------------------------------------------
# EN + OL enrichment
EN_enrichment <- readRDS("./results/open_region_enrichment/EN.rds")

EN_enrichment %>% group_by(stringency) %>% 
  summarise("mean_enrichment_ratio" = mean(enrichment_ratio),
            "sd_enrichment_ratio"= sd(enrichment_ratio)) -> EN_enrichment_summary
EN_enrichment_summary["celltype"] <- "Excitatory Neuron"

OL_enrichment <- readRDS("./results/open_region_enrichment/Oligodendrocyte.rds")

OL_enrichment %>% group_by(stringency) %>% 
  summarise("mean_enrichment_ratio" = mean(enrichment_ratio),
            "sd_enrichment_ratio"= sd(enrichment_ratio)) -> OL_enrichment_summary
OL_enrichment_summary["celltype"] <- "Oligodendrocyte"

enrichment_summary <- rbind(EN_enrichment_summary, OL_enrichment_summary)
enrichment_summary$sd_enrichment_ratio * 1.96 -> enrichment_summary$CI95_enrichment_ratio

ggplot(enrichment_summary, aes(x = stringency, y = mean_enrichment_ratio, 
                               fill = celltype, group = celltype)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = mean_enrichment_ratio - CI95_enrichment_ratio,
                    ymax = mean_enrichment_ratio + CI95_enrichment_ratio),
                width = 0.25,
                position = position_dodge(width = 0.9)) +
  geom_hline(yintercept=1, linetype="dashed", color="red") +
  #facet_wrap(~celltype,  ncol=1, strip.position = "right") +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = 6, family = font),
    legend.position = "bottom",
    #axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1)
  ) +
  scale_fill_manual(values = c("#7CAE00", "#00A9FF")) +
  labs(x = "Stringency", y = "Mutation enrichment ratio \n (covered region / whole genome)")
ggsave("./figures/manuscript_figures/Figure3/OL_EN_enrichment.pdf", width = 2.3, height = 2.3, units = "in")

# ------------------------------------------------------------------------------
# All cell type
simulation_results <- 
  readRDS("./results/burden_simulation/burden_summary_all_cell_type_20241018.rds")

neuron_simulation_results <- simulation_results[simulation_results$cell_type %in% c("EN", "IN"),]
neuron_simulation_results$cell_type_figure <- 
  factor(neuron_simulation_results$cell_type, levels = c("EN", "IN"),
         labels = c("Excitatory Neuron", "Inhibitory Neuron"))
neuron_simulation_results$rate_sd * 1.96 -> neuron_simulation_results$rate_95CI

# EN + IN
ggplot(neuron_simulation_results, aes(x = Age, y = rate, color = cell_type_figure)) +
  geom_point(size = 0.2) +
  #geom_smooth(method = "lm", se = F, size = 0.5) +
  #stat_cor() +
  geom_abline(slope = burden_regression$raw_rate_regression_results$EN["mean_slope"], 
              intercept = burden_regression$raw_rate_regression_results$EN["mean_intercept"], 
              size = 0.5, color = "#7CAE00") +
  geom_abline(slope = burden_regression$raw_rate_regression_results$IN["mean_slope"], 
              intercept = burden_regression$raw_rate_regression_results$IN["mean_intercept"], 
              size = 0.5, color = "#E69F00") +
  annotate(geom="text", x=20, y=4e-7, parse = T, 
           label = paste0("italic(P) == ", round(burden_regression$raw_rate_adj_p_value_list["EN_vs_IN"], 3)), 
           color = "black", family = font, size = 2) +
  theme_classic() +  
  geom_errorbar(aes(ymax = rate + rate_95CI, ymin = rate - rate_95CI), size = 0.2, width = 4) +
  guides(color=guide_legend(nrow=2), override.aes = list(size=1)) +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = 6, family = font),
    legend.position = c(0.7, 0.12),
    legend.key.size = unit(0.1, 'in'),
    legend.title = element_blank(),
    legend.background = element_rect(fill='transparent'),
  ) +
  scale_y_continuous(labels = scientific_10) +
  scale_color_manual(values = c("#7CAE00", "#E69F00")) +
  labs(x = "Age (yrs)", y = "Mutational Burden (sSNVs/bp)")
ggsave("./figures/manuscript_figures/Figure3/IN_EN_burden_raw.pdf", width = 2, height = 2, units = "in")

# Glial cell types
OL_OPC_simulation_results <- simulation_results[simulation_results$cell_type %in% c("Oligodendrocyte", "OPC"),]
OL_OPC_simulation_results$rate_sd * 1.96 -> OL_OPC_simulation_results$rate_95CI

ggplot(OL_OPC_simulation_results, aes(x = Age, y = rate, color = cell_type)) +
  geom_point(size = 0.2) +
  geom_abline(slope = burden_regression$raw_rate_regression_results$Oligodendrocyte["mean_slope"], 
              intercept = burden_regression$raw_rate_regression_results$Oligodendrocyte["mean_intercept"], 
              size = 0.5, color = "#00A9FF") +
  geom_abline(slope = burden_regression$raw_rate_regression_results$OPC["mean_slope"], 
              intercept = burden_regression$raw_rate_regression_results$OPC["mean_intercept"], 
              size = 0.5, color = "#CC79A7") +
  annotate(geom="text", x=20, y=5.5e-7, parse = T, 
           label = paste0("italic(P) == ", round(burden_regression$raw_rate_adj_p_value_list["OPC_vs_Oligodendrocyte"], 3)), 
           color = "black", family = font, size = 2) +  
  theme_classic() +  
  geom_errorbar(aes(ymax = rate + rate_95CI, ymin = rate - rate_95CI), size = 0.2, width = 4) +
  guides(color=guide_legend(nrow=2), override.aes = list(size=1)) +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = 6, family = font),
    legend.position = c(0.7, 0.12),
    legend.key.size = unit(0.1, 'in'),
    legend.title = element_blank(),
    legend.background = element_rect(fill='transparent'),
  ) +
  scale_y_continuous(labels = scientific_10) +
  scale_color_manual(values = c("#00A9FF", "#CC79A7")) +
  labs(x = "Age (yrs)", y = "Mutational Burden (sSNVs/bp)")
ggsave("./figures/manuscript_figures/Figure3/OL_OPC_burden_raw.pdf", width = 2, height = 2, units = "in")

Astro_Micro_simulation_results <- simulation_results[simulation_results$cell_type %in% c("Astrocyte", "Microglia"),]
Astro_Micro_simulation_results$rate_sd * 1.96 -> Astro_Micro_simulation_results$rate_95CI

ggplot(Astro_Micro_simulation_results, aes(x = Age, y = rate, color = cell_type)) +
  geom_point(size = 0.2) +
  #geom_smooth(method = "lm", se = F, size = 0.5) +
  #stat_cor() +
  geom_abline(slope = burden_regression$raw_rate_regression_results$Astrocyte["mean_slope"], 
              intercept = burden_regression$raw_rate_regression_results$Astrocyte["mean_intercept"], 
              size = 0.5, color = "#2E2585") +
  geom_abline(slope = burden_regression$raw_rate_regression_results$Microglia["mean_slope"], 
              intercept = burden_regression$raw_rate_regression_results$Microglia["mean_intercept"], 
              size = 0.5, color = "#DCCD7D") +
  annotate(geom="text", x=20, y=5.5e-7, parse = T, 
           label = paste0("italic(P) == ", round(burden_regression$raw_rate_adj_p_value_list["Astrocyte_vs_Microglia"], 3)), 
           color = "black", family = font, size = 2) +  
  theme_classic() +  
  geom_errorbar(aes(ymax = rate + rate_95CI, ymin = rate - rate_95CI), size = 0.2, width = 4) +
  guides(color=guide_legend(nrow=4), override.aes = list(size=1)) +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = 6, family = font),
    legend.position = c(0.7, 0.12),
    legend.key.size = unit(0.1, 'in'),
    legend.title = element_blank(),
    legend.background = element_rect(fill='transparent'),
  ) +
  scale_y_continuous(labels = scientific_10) +
  scale_color_manual(values = c("#2E2585", "#DCCD7D")) +
  labs(x = "Age (yrs)", y = "Mutational Burden (sSNVs/bp)")
ggsave("./figures/manuscript_figures/Figure3/Astro_Micro_burden_raw.pdf", width = 2, height = 2, units = "in")


# ------------------------------------------------------------------------------
# EN sub type
subtype_simulation_results <- 
  readRDS("./results/burden_simulation/burden_summary_subtype_20241018.rds")
subtype_simulation_results <- subtype_simulation_results[str_detect(subtype_simulation_results$cell_type, "EN"),]
subtype_simulation_results$cell_type_figure <- 
  factor(subtype_simulation_results$cell_type, levels = c("EN_upper_layer", "EN_deep_layer"),
         labels = c("Upper-layer EN", "Deep-layer EN"))
subtype_simulation_results$rate_sd * 1.96 -> subtype_simulation_results$rate_95CI

burden_regression_subtype <- readRDS("./figures/manuscript_figures/Figure3/burden_regression_results_subtype.rds")

ggplot(subtype_simulation_results, aes(x = Age, y = rate, color = cell_type_figure)) +
  geom_point(size = 0.2) +
  geom_abline(slope = burden_regression_subtype$raw_rate_regression_results$EN_deep_layer["mean_slope"], 
              intercept = burden_regression_subtype$raw_rate_regression_results$EN_deep_layer["mean_intercept"], 
              size = 0.5, color = "#0F65A1") +
  geom_abline(slope = burden_regression_subtype$raw_rate_regression_results$EN_upper_layer["mean_slope"], 
              intercept = burden_regression_subtype$raw_rate_regression_results$EN_upper_layer["mean_intercept"],
              size = 0.5, color = "#6A4A3C") +
  annotate(geom="text", x=20, y=6e-7, parse = T, 
           label = paste0("italic(P) == ", round(burden_regression_subtype$raw_rate_p_value, 3)), 
           color = "black", family = font, size = 2) +
  theme_classic() +  
  geom_errorbar(aes(ymax = rate + rate_95CI, ymin = rate - rate_95CI), size = 0.2, width = 4) +
  guides(color=guide_legend(nrow=2), override.aes = list(size=1)) +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = 6, family = font),
    legend.position = c(0.7, 0.09),
    legend.key.size = unit(0.1, 'in'),
    legend.title = element_blank(),
    legend.background = element_rect(fill='transparent'),
  ) +
  scale_y_continuous(labels = scientific_10) +
  scale_color_manual(values = c("#6A4A3C", "#0F65A1")) +
  labs(x = "Age (yrs)", y = "Mutation Burden (sSNVs/bp)")
ggsave("./figures/manuscript_figures/Figure3/EN_subtype_burden_raw.pdf", width = 2, height = 2, units = "in")

