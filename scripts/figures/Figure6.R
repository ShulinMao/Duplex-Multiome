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
library(GenomicRanges)
library(Seurat)
library(SeuratDisk)
library(Signac)

# font
font = "Helvetica"
font_size = 7
axis_font_size = 6


# Custom function to format scientific notation with superscripts and handle zero
scientific_10 <- function(x) {
  ifelse(x == 0, "0", parse(text = gsub("e[\\+]*", " %*% 10^", label_scientific()(x))))
}

colors_cell_type <- c("#2E2585", "#337538", "#7CAE00", "#E69F00", "#DCCD7D", "#00A9FF", "#CC79A7", "#7E2954")

# ------------------------------------------------------------------------------
# UMAP of AN06365
brain <- readRDS("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/RNA_ATAC/multiomes_integrated_20240813_ASD_annotated.rds")
Idents(brain) <- brain$celltype
AN06365 <- brain[, brain$dataset == "AN06365"]

# cell type
AN06365[["cell_type"]] <- as.character(AN06365$celltype)
AN06365$cell_type[str_detect(AN06365$celltype, "EN")] <- "Excitatory\nNeuron"
AN06365$cell_type[str_detect(AN06365$celltype, "IN")] <- "Inhibitory\nNeuron"

DimPlot(AN06365, reduction = "wnn.umap", group.by = "cell_type", 
        pt.size = 1e-10, cols = colors_cell_type) + 
  labs(title = NULL, color = "Cell type") +
  guides(color=guide_legend(nrow=3, byrow=T, override.aes = list(size=2))) +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "bottom",
    legend.key.size = unit(0, 'in'),
    legend.title = element_text(angle = 90, hjust = 0.3, vjust = 2.5)
  )
ggsave("./figures/manuscript_figures/Figure6/AN06365_UMAP_cell_type.pdf", height = 3, width = 2.5, units = "in")

# mutation
# Chr1 205284910 T>C
detected <- read.table("./results/celltype_annotation/AN06365_chr1_205284910_T_C.txt")[,1, drop=T]
detected <- paste0("AN06365_", detected)
detected <- detected[detected %in% colnames(AN06365)]
covered <- read.table("./data/AN06365/single_cell/chr1_205284910_cov_CB_RF.txt")[,1, drop=T]
covered <- str_replace(covered, "^CB:Z:", "AN06365_")
covered <- unique(covered)
covered <- covered[covered %in% colnames(AN06365)]

vaf <- round(length(detected)/length(covered), 3)

AN06365[["Coverage"]] <- "No coverage"
AN06365$Coverage[colnames(AN06365) %in% covered] <- "Covered"
AN06365$Coverage[colnames(AN06365) %in% detected] <- "Alt allele detected"

AN06365_mutation_map <- cbind(AN06365@reductions$wnn.umap@cell.embeddings, AN06365[[c("Coverage", "cell_type")]])
AN06365_mutation_map$Coverage <- factor(AN06365_mutation_map$Coverage)

vaf_cell_type <- AN06365_mutation_map %>% count(cell_type, Coverage, name = "count", .drop = F)
vaf_cell_type <- vaf_cell_type[vaf_cell_type$Coverage != "No coverage",]
vaf_cell_type %>% group_by(cell_type) %>% mutate("vaf" = count/sum(count)) -> vaf_cell_type
vaf_cell_type <- vaf_cell_type[vaf_cell_type$Coverage == "Alt allele detected",]
vaf_cell_type

ggplot(AN06365_mutation_map %>% arrange(desc(Coverage)), aes(x = wnnUMAP_1, y = wnnUMAP_2))+
  geom_point(aes(shape = Coverage, size = Coverage, color = Coverage)) +
  annotate(geom="text", x = -9, y=15, label = paste0("Duplex-Multiome VAF = ", vaf), 
           color = "black", family = font, size = 2, hjust = 0) +
  annotate(geom="text", x = -9, y=13, label = paste0("Bulk WGS VAF = 0.053"), 
           color = "black", family = font, size = 2, hjust = 0) +
  scale_color_manual(values = c("red3", "gold1", "grey")) +
  scale_shape_manual(values = c(8, 16, 16)) +
  scale_size_manual(values = c(2, 1, 1)) +
  labs(title = NULL) +
  guides(size = guide_legend(nrow=3, byrow=T, override.aes = list(size=2))) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.direction = "vertical"
  )
ggsave("./figures/manuscript_figures/Figure6/AN06365_UMAP_chr1_205284910_T_C.pdf", 
       height = 3.2, width = 2.3, units = "in")

# ------------------------------------------------------------------------------
# Mutation effects
DefaultAssay(brain) <- "ATAC" 
DefaultAssay(AN06365) <- "ATAC"

Idents(AN06365) <- AN06365$cell_type
ranges.show <- StringToGRanges("chr1-205284730-205285079")
CoveragePlot(
  object = AN06365, idents = c("Astrocyte", "Microglia", "Excitatory\nNeuron", "Inhibitory\nNeuron", "Oligodendrocyte"),
  region = c("chr1-205284910-205284910"),
  ranges = ranges.show,
  ranges.title = "ENCODE cCREs",
  annotation = F,
  extend.upstream = 500, extend.downstream = 500
) +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
  )
ggsave("./figures/manuscript_figures/Figure6/AN06365_chr1_205284910_500_window.pdf",
       width = 5, height = 3, units = "in")


Idents(brain) <- brain$dataset
ranges.show <- StringToGRanges("chr1-205284730-205285079")
CoveragePlot(
  object = brain, idents = "AN06365",
  region = c("chr1-205284910-205284910"),
  ranges = ranges.show,
  peaks = F,
  extend.upstream = 500000, extend.downstream = 500000
)+
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
  )
ggsave("./figures/manuscript_figures/Figure6/AN06365_chr1_205284910_ATAC_1mb_window.pdf",
       width = 10, height = 10, units = "in")

# ------------------------------------------------------------------------------
# Gene expression
NUSCK1_expression_summary <- readRDS("./figures/manuscript_figures/Figure6/AN06365_NUSCK1_expr.rds")
scientific_10(as.numeric(NUSCK1_expression_summary[2, "p"])) -> NUSCK1_expression_summary[2, "p"]
rownames(NUSCK1_expression_summary)[2] <- "Controls"
NUSCK1_expression_summary$condition <- rownames(NUSCK1_expression_summary)
NUSCK1_expression_summary$group2 <- "Controls"
ggplot(NUSCK1_expression_summary, aes(x = condition, y = mean)) +
  geom_point(size = 0.5) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15) +
  stat_pvalue_manual(NUSCK1_expression_summary[2,], 
                     label = "P = {p}", tip.length = 0.01, 
                     y.position = c(0.275)) +
  theme_classic() +
  labs(x = NULL, y = "Expression level", subtitle = "NUCKS1") +
  ylim(c(0.23, 0.28)) +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )
ggsave("./figures/manuscript_figures/Figure6/AN06365_NUSCK1_expr.pdf", 
       width = 1.6, height = 1.6, unit = "in")


# DEG
DEG_data <- readRDS("./figures/manuscript_figures/Figure6/AN06365_DEG.rds")
log2fc <- DEG_data$log2fc
pvalue <- DEG_data$pvalue

genes <- c("NFASC", "CNTN2", "RBBP5", "TMEM81", "DSTYK", "TMCC2", "NUAK2", "KLHDC8A", "LEMD1", 
           "BLACAT1", "RP11-576D8.4", "CDK18", "MFSD4A", "ELK4", "SLC45A3", "NUCKS1","RAB29")
genes <- intersect(genes, rownames(log2fc))

log2fc <- log2fc[genes,]
pvalue <- pvalue[genes,]

colnames(log2fc) <- c("Astrocyte", "Microglia", "Excitatory Neuron", "Inhibitory Neuron", "Oligodendrocyte")
colnames(pvalue) <- c("Astrocyte", "Microglia", "Excitatory Neuron", "Inhibitory Neuron", "Oligodendrocyte")

cutoff.distance <- 0
makeColorRampPalette <- function(colors, cutoff.fraction, num.colors.in.palette)
{
  stopifnot(length(colors) == 4)
  ramp1 <- colorRampPalette(colors[1:2])(num.colors.in.palette * cutoff.fraction)
  ramp2 <- colorRampPalette(colors[3:4])(num.colors.in.palette * (1 - cutoff.fraction))
  return(c(ramp1, ramp2))
}
cols <- makeColorRampPalette(c("royalblue3", "white",    # distances 0 to 3 colored from white to red
                               "white", "orange"), # distances 3 to max(distmat) colored from green to black
                             abs(cutoff.distance - min(log2fc))/ (max(log2fc) - min(log2fc)),
                             100)
annotation <- data.frame(row.names(log2fc))

pheatmap(t(log2fc), scale = "none", color = cols, 
         cluster_rows = F, cluster_cols = F, 
         display_numbers = t(pvalue), number_color = "red",
         cellheight = 10, cellwidth = 10, 
         fontsize = axis_font_size)

pheatmap(t(log2fc), scale = "none", color = cols, 
         cluster_rows = F, cluster_cols = F, 
         display_numbers = t(pvalue), 
         number_color = "red", fontsize_number = font_size,
         cellheight = 10, cellwidth = 15, 
         fontsize = axis_font_size, angle_col = 45,
         filename = "./figures/manuscript_figures/Figure6/AN06365_DEG_nearby.pdf",
         width = 10, height = 5)

