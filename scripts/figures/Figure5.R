# Libraries
.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(pheatmap)
library(ggpubr)
library(scales)
library(reshape2)
library(patchwork)

ref_genome="BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = T)
chr_orders <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages_sig/")
library(MutationalPatterns)


# font
font = "Helvetica"
font_size = 7
axis_font_size = 6

# Custom function to format scientific notation with superscripts and handle zero
scientific_10 <- function(x) {
  ifelse(x == 0, "0", parse(text = gsub("e[\\+]*", " %*% 10^", label_scientific()(x))))
}


# scRNAseq
brain <- readRDS("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/RNA_ATAC/multiomes_integrated_20240512_infant_adult_no_QC.rds")

brain[["cell_type"]] <- as.character(brain$celltype_transferred_ATAC)
brain$cell_type[str_detect(brain$celltype_transferred_ATAC, "EN")] <- "Excitatory\nNeuron"
brain$cell_type[str_detect(brain$celltype_transferred_ATAC, "IN")] <- "Inhibitory\nNeuron"

brain <- subset(
  x = brain,
  subset = nCount_ATAC > 1000
)
gc()
# brain <- readRDS("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/RNA_ATAC/multiomes_integrated_20240512_infant_adult_annotated.rds")
# 
# brain[,!str_detect(brain$celltype, "doublet")] -> brain
# brain[["cell_type"]] <- as.character(brain$celltype)
# brain$cell_type[str_detect(brain$celltype, "EN")] <- "Excitatory\nNeuron"
# brain$cell_type[str_detect(brain$celltype, "IN")] <- "Inhibitory\nNeuron"

brain_mutation_map <- cbind(brain@reductions$umap.atac@cell.embeddings, brain[["cell_type"]])
brain_mutation_map <- brain_mutation_map[str_starts(row.names(brain_mutation_map), case_id),]

ggplot(brain_mutation_map, aes(x = atacUMAP_1, y = atacUMAP_2))+
  geom_point(aes(color = cell_type), size = 0.01) +
  scale_color_manual(values = c("#2E2585", "#337538", "#7CAE00", "#E69F00", "#DCCD7D", "#00A9FF", "#CC79A7", "#7E2954")) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "bottom",
    legend.title = element_blank()
  )
ggsave("./figures/manuscript_figures/Figure5/clonal_mutation_UMAP.pdf", width = 2.1, height = 2.8, unit = "in")


# clonal mutation
clonal_mutations <- read.csv("./results/clonal_mutation/clonal_mutation_list_healthy_brain_realigned_20241028.csv")

# Brain1901106 Chr7 155297830 T>C
case_id <- "Brain1901106_deeper"
mutation <- "chr7-155297830-T-C"
detected <- clonal_mutations$cell[clonal_mutations$case_id == case_id & clonal_mutations$key == mutation]
covered <- readRDS(paste0("./figures/manuscript_figures/Figure5/", mutation, "_covered_cell.rds"))
vaf <- round(length(detected)/length(covered), 3)

detected <- paste0(case_id, "_", detected)
detected <- detected[detected %in% colnames(brain)]
covered <- paste0(case_id, "_", covered)
covered <- covered[covered %in% colnames(brain)]

brain[["Coverage"]] <- "No coverage"
brain$Coverage[colnames(brain) %in% covered] <- "Covered"
brain$Coverage[colnames(brain) %in% detected] <- "Alt allele detected"

brain_mutation_map <- cbind(brain@reductions$umap.atac@cell.embeddings, brain[["Coverage"]])
brain_mutation_map <- brain_mutation_map[str_starts(row.names(brain_mutation_map), case_id),]

p1 <- ggplot(brain_mutation_map %>% arrange(desc(Coverage)), aes(x = atacUMAP_1, y = atacUMAP_2))+
  geom_point(aes(shape = Coverage, size = Coverage, color = Coverage)) +
  annotate(geom="text", x = -12, y=18, label = paste0("Duplex-Multiome VAF = ", vaf), 
           color = "black", family = font, size = 2.5, hjust = 0) +
  annotate(geom="text", x = -12, y=16, label = paste0("Bulk WGS VAF = not detected"), 
           color = "black", family = font, size = 2.5, hjust = 0) +
  scale_color_manual(values = c("red3", "gold1", "grey")) +
  scale_shape_manual(values = c(8, 16, 16)) +
  scale_size_manual(values = c(2, 1, 1)) +
  labs(subtitle = paste(mutation, "(Brain1901106)")) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# UMB5823 chr4:7208281 C>T
case_id <- "UMB5823_deeper"
mutation <- "chr4-7208281-C-T"
detected <- clonal_mutations$cell[clonal_mutations$case_id == case_id & clonal_mutations$key == mutation]
covered <- readRDS(paste0("./figures/manuscript_figures/Figure5/", mutation, "_covered_cell.rds"))
vaf <- round(length(detected)/length(covered), 3)

detected <- paste0(case_id, "_", detected)
detected <- detected[detected %in% colnames(brain)]
covered <- paste0(case_id, "_", covered)
covered <- covered[covered %in% colnames(brain)]

brain[["Coverage"]] <- "No coverage"
brain$Coverage[colnames(brain) %in% covered] <- "Covered"
brain$Coverage[colnames(brain) %in% detected] <- "Alt allele detected"

brain_mutation_map <- cbind(brain@reductions$umap.atac@cell.embeddings, brain[["Coverage"]])
brain_mutation_map <- brain_mutation_map[str_starts(row.names(brain_mutation_map), case_id),]

p2 <- ggplot(brain_mutation_map %>% arrange(desc(Coverage)), aes(x = atacUMAP_1, y = atacUMAP_2))+
  geom_point(aes(shape = Coverage, size = Coverage, color = Coverage)) +
  annotate(geom="text", x = -12, y=18, label = paste0("Duplex-Multiome VAF = ", vaf), 
           color = "black", family = font, size = 2.5, hjust = 0) +
  annotate(geom="text", x = -12, y=16, label = paste0("Bulk WGS VAF = 0.077"), 
           color = "black", family = font, size = 2.5, hjust = 0) +
  scale_color_manual(values = c("red3", "gold1", "grey")) +
  scale_shape_manual(values = c(8, 16, 16)) +
  scale_size_manual(values = c(2, 1, 1)) +
  labs(subtitle = paste0(mutation, " (UMB5823)")) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# UMB5823 chr9:89311105 G>A

case_id <- "UMB5823_deeper"
mutation <- "chr9-89311105-G-A"
detected <- clonal_mutations$cell[clonal_mutations$case_id == case_id & clonal_mutations$key == mutation]
covered <- readRDS(paste0("./figures/manuscript_figures/Figure5/", mutation, "_covered_cell.rds"))
vaf <- round(length(detected)/length(covered), 3)

detected <- paste0(case_id, "_", detected)
detected <- detected[detected %in% colnames(brain)]
covered <- paste0(case_id, "_", covered)
covered <- covered[covered %in% colnames(brain)]


brain[["Coverage"]] <- "No coverage"
brain$Coverage[colnames(brain) %in% covered] <- "Covered"
brain$Coverage[colnames(brain) %in% detected] <- "Alt allele detected"

brain_mutation_map <- cbind(brain@reductions$umap.atac@cell.embeddings, brain[["Coverage"]])
brain_mutation_map <- brain_mutation_map[str_starts(row.names(brain_mutation_map), case_id),]

p3 <- ggplot(brain_mutation_map %>% arrange(desc(Coverage)), aes(x = atacUMAP_1, y = atacUMAP_2))+
  geom_point(aes(shape = Coverage, size = Coverage, color = Coverage)) +
  annotate(geom="text", x = -12, y=18, label = paste0("Duplex-Multiome VAF = ", vaf), 
           color = "black", family = font, size = 2.5, hjust = 0) +
  annotate(geom="text", x = -12, y=16, label = paste0("Bulk WGS VAF = not detected"), 
           color = "black", family = font, size = 2.5, hjust = 0) +
  scale_color_manual(values = c("red3", "gold1", "grey")) +
  scale_shape_manual(values = c(8, 16, 16)) +
  scale_size_manual(values = c(2, 1, 1)) +
  labs(subtitle = paste(mutation, "(UMB5823)")) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "bottom",
    legend.title = element_blank()
  )

ggarrange(p1, p2, p3, ncol = 3, common.legend = T, legend = "bottom")
ggsave("./figures/manuscript_figures/Figure5/clonal_mutation_examples.pdf", width = 6.3, height = 2.6, unit = "in")

p4 <- DimPlot(brain, reduction = "umap.rna", group.by = "celltype_transferred_ATAC", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("RNA")
p5 <- DimPlot(brain, reduction = "umap.atac", group.by = "celltype_transferred_ATAC", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("ATAC")
p6 <- DimPlot(brain, reduction = "wnn.umap", group.by = "celltype_transferred_ATAC", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("WNN")
ggpubr::ggarrange(p4, p5, p6, nrow = 1, common.legend = T, legend = "bottom")

# Coverage VS clonal mutation number
clonal_mutation_count_coverage <- readRDS("./figures/manuscript_figures/Figure5/clonal_mutation_count_coverage.rds")
summary(lm(count ~ genomecov_clone, clonal_mutation_count_coverage)) -> clonal_mutation_count_coverage_lm

ggplot(clonal_mutation_count_coverage, aes(x = genomecov_clone, y = count)) +
  geom_point() +
  geom_smooth(method = 'lm', formula = y~x) +
  annotate(geom="text", x=2e10, y=100, parse = T, 
           label = paste0("italic(P) == ", round(clonal_mutation_count_coverage_lm$coefficients["genomecov_clone", "Pr(>|t|)"], 3)), 
           color = "black", family = font, size = 2) +
  labs(x = "Callable region size", y = "Number of clonal sSNVs detected") +
  scale_x_continuous(labels = scientific_10) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font)
  )
ggsave("./figures/manuscript_figures/Figure5/clonal_mutation_vs_callable_region_size.pdf", 
       width = 2, height = 2, units = "in")

# clonal mutation spectrum
type_occurrences <- readRDS("./figures/manuscript_figures/Figure5/clonal_mutation_type_occurrences.rds")
plot_spectrum(type_occurrences, CT = TRUE) + 
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font)
  )
ggsave("./figures/manuscript_figures/Figure5/clonal_mutation_type_occurrence.pdf", 
       width = 3, height = 2.2, units = "in")


# Glial vs Neuron
readRDS("./figures/manuscript_figures/Figure5/clonal_mutation_count_coverage_cell_type.rds") -> clonal_mutation_count_coverage_celltype

wilcox.test(x = clonal_mutation_count_coverage_celltype$rate[clonal_mutation_count_coverage_celltype$cell_type == "neuronal"],
            y = clonal_mutation_count_coverage_celltype$rate[clonal_mutation_count_coverage_celltype$cell_type == "glial"], 
            paired = T, alternative = "less") -> wilcox_result

# Main
p1 <- ggplot(clonal_mutation_count_coverage_celltype, aes(x = cell_type, y = rate, group = Sample)) +
  geom_point(size = 1) +
  geom_line() +
  labs(x = "Cell type", y = "Detected clonal sSNVs per bp") +
  theme_classic() +
  theme(    
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font)) +
  scale_x_discrete(expand = expansion(mult = c(1.1, 1.1)), labels = str_to_title) + 
  coord_cartesian(ylim = c(0, 5e-10)) +
  scale_y_continuous(breaks = seq(0, 5e-10, 2.5e-10), labels = scientific_10)

# Outlier
p2 <- ggplot(clonal_mutation_count_coverage_celltype, aes(x = cell_type, y = rate, group = Sample)) +
  geom_point(size = 1) +
  geom_line() +
  labs(x = "Callable region size", y = NULL) +
  theme_classic() +
  coord_cartesian(ylim = c(1.4e-9, 1.5e-9)) +  # "Upper" segment
  scale_y_continuous(breaks = seq(1.4e-9, 1.5e-9, 1e-10), labels = scientific_10) +
  scale_x_discrete(expand = expansion(mult = c(1.1, 1.1))) + 
  annotate(geom="text", x=1, y=1.45e-9, parse = T, 
           label = paste0("italic(P) == ", round(wilcox_result$p.value, 3)), 
           color = "black", size = 2, family = font) + 
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
  )

# Stack them vertically
combined_plot <- p2 / p1 + plot_layout(heights = c(1, 5))
combined_plot

ggsave("./figures/manuscript_figures/Figure5/glial_vs_neuron.pdf", 
       width = 2, height = 2, units = "in")
