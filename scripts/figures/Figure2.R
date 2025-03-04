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
library(Polychrome)
library(pheatmap)
library(igraph)

# font
font = "Helvetica"
font_size = 7
axis_font_size = 6

# ------------------------------------------------------------------------------
# Schematic of COLO829-BLT50 cell line mixture experiment

multiomes_integrated <- readRDS("./data/RNA_ATAC/multiomes_integrated_20240523_COLO829BLT50.rds")
multiomes_integrated$celltype_figure <- factor(multiomes_integrated$celltype, 
                                               levels = c("melanoma_fibroblast", "B_lymphoblast"),
                                               labels = c("Melanoma", "B-cell"))


p1 <- DimPlot(multiomes_integrated, reduction = "umap.rna", group.by = "celltype_figure", 
              #label = TRUE, label.size = 3.5, repel = TRUE, 
              pt.size = 0.05) + 
  #ggtitle("RNA") + 
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "bottom"
  ) +
  scale_color_manual(values = c("#A34C49", "lightskyblue1"))+
  labs(title = NULL)
# p2 <- DimPlot(multiomes_integrated, reduction = "umap.atac", group.by = "celltype_figure", 
#               #label = TRUE, label.size = 3.5, repel = TRUE, 
#               pt.size = 0.05) + 
#   ggtitle("ATAC") + 
#   theme(
#     text = element_text(color = "black", size = font_size, family = font),
#     axis.text = element_text(color = "black", size = axis_font_size, family = font),
#     legend.text = element_text(color = "black", size = axis_font_size, family = font),
#     #legend.position = "none"
#   ) +
#   scale_color_manual(values = c("#A34C49", "lightskyblue1"))
# ggpubr::ggarrange(p1, p2, nrow = 1, common.legend = T, legend = "bottom")
ggsave("./figures/manuscript_figures/Figure2/UMAP_COLO829BLT50.pdf", p1, 
       width = 1.6, height = 2, units = "in")

table(multiomes_integrated$celltype_figure)/dim(multiomes_integrated)[2]
ggplot(multiomes_integrated[[c("cell.id", "celltype_figure")]], 
       aes(x = celltype_figure, fill = celltype_figure)) +
  geom_bar() +
  annotate(geom="text", x=1, y=600, label = "128\n(2.5%)", 
           color = "black", family = font, size = 2) +
  annotate(geom="text", x=2, y=5500, label = "4989\n(97.5%)", 
           color = "black", family = font, size = 2) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  ylim(0, 5700) +
  scale_fill_manual(values = c("#A34C49", "lightskyblue1")) +
  labs(x = "", y = "Number of cells")
ggsave("./figures/manuscript_figures/Figure2/Cell_composition_COLO829BLT50.pdf", 
       width = 1.2, height = 2, units = "in")

# ------------------------------------------------------------------------------
# Lineage
load("./figures/manuscript_figures/Figure2/lineage_COLO829BLT50.RData")

# triangle vertex shape
mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
# clips as a circle
add_shape("triangle", clip=shapes("circle")$clip,
          plot=mytriangle)


label_prop_color <- palette36.colors(15)
V(g)$color <- as.character(factor(label_prop_cluster, labels = label_prop_color))
pdf("./figures/manuscript_figures/Figure2/lineage_COLO829BLT50.pdf", width = 10, height = 10)
plot(g, layout=layout, edge.width=0.5, vertex.size=4, vertex.frame.edge=0)
dev.off()

df_label_prop_cluster_count$celltype_figure <- factor(df_label_prop_cluster_count$cell_type, 
                           levels = c("melanoma_fibroblast", "filtered", "B_lymphoblast"),
                           labels = c("Melanoma", "Unable to determine\n(Low RNA quality)", "B-cell"))
ggplot(df_label_prop_cluster_count, aes(x=label_prop_cluster, fill=celltype_figure, 
                                        group=celltype_figure)) +
  
  geom_bar(stat="count") +
  scale_x_continuous(breaks=seq(1, length(unique(df_label_prop_cluster$label_prop_cluster))),
                     labels=seq(1, length(unique(df_label_prop_cluster$label_prop_cluster)))) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    #legend.position = "bottom",
    legend.key.size = unit(0.15, 'in'),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  scale_fill_manual(values=c("#A34C49", "gray83", "lightskyblue1")) +
  #coord_flip() +
  labs(x = "Clusters", y = "Number of cells", fill = "Cell type")
ggsave("./figures/manuscript_figures/Figure2/label_prop_cluster_lineage_COLO829BLT50.pdf", 
       width = 3.8, height = 2.7, units = "in")

# ------------------------------------------------------------------------------
# TP rates
overlap_df <- readRDS("./figures/manuscript_figures/Figure2/benchmark_truthset_COLO829BLT50.rds")
ggplot(overlap_df[overlap_df$variable == "duplex_multiome_clonal_mutation_0.25",],
       aes(x = "", y = value, fill = group)) +
  geom_bar(stat="identity", width=1) +
  annotate(geom="text", x=-1, y=0, label = "Clonal mutation\n(VAF > 0.25,\nsurviving B-cell\nbulk filter)", 
           color = "black", family = font, size = 2.5) +
  coord_polar("y", start=0) +
  theme_void() +
  theme(
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "bottom",
    legend.direction = "horizontal", 
    legend.margin=margin(),
  ) +
  geom_text(aes(y = ypos, label = value), color = "white", size=3) +
  labs(fill = "")
ggsave("./figures/manuscript_figures/Figure2/overlap_duplex_multiome_clonal_mutation_VAF25_BL_bulk_COLO829BLT50.pdf", 
       width = 2.7, height = 2.7, units = "in")

# ------------------------------------------------------------------------------
# Spectrum

# Truthset spectrum
ref_genome="BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = T)
chr_orders <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages_sig/")
library(MutationalPatterns)

# Tumor Truth set spectrum
smaht_truthset <- read_vcfs_as_granges("./data/others/COLO829BLT50_truth_set/SMaHT_COLO829_SNV_truth_set_v1.0.vcf", 
                                       "Melanoma", ref_genome)
mut_mat_tumor_truthset <- mut_matrix(vcf_list = smaht_truthset, ref_genome = ref_genome)

# BL truth set spectrum
COLO829BL_truthset <- read_vcfs_as_granges("./data/others/COLO829BLT50_truth_set/COLO829BL_somatic_07252024_snv_af_gt_25_dp_gt_100.vcf.gz",
                                           "B-cell", ref_genome)
mut_mat_COLO829BL_truthset <- mut_matrix(vcf_list = COLO829BL_truthset, ref_genome = ref_genome)

# COLO829BLT clonal spectrum
clonal_mutation_list_spectrum_mat <- readRDS("./figures/manuscript_figures/Figure2/spectrum_clonal_mutation_COLO829BLT50.rds")

# COLO829BLT private spectrum
private_mutation_list_spectrum_mat <- readRDS("./figures/manuscript_figures/Figure2/spectrum_private_a6s3_COLO829BLT50.rds")

mut_mat_tumor <- cbind(mut_mat_tumor_truthset,
                       clonal_mutation_list_spectrum_mat$melanoma_fibroblast_clonal_mutation_list_BL_bulk,
                       private_mutation_list_spectrum_mat$melanoma_fibroblast,
                       get_known_signatures()[, "SBS7a", drop=T])
colnames(mut_mat_tumor) <- c("Truth set", "B-cell bulk", "BLT 50 bulk", "SBS7a")

mut_mat_BL <- cbind(mut_mat_COLO829BL_truthset,
                       clonal_mutation_list_spectrum_mat$B_lymphoblast_clonal_mutation_list_tumor_bulk,
                       private_mutation_list_spectrum_mat$B_lymphoblast,
                    get_known_signatures()[, "SBS18", drop=T])
colnames(mut_mat_BL) <- c("Truth set", "Melanoma bulk", "BLT 50 bulk", "SBS18")

plot_96_profile(mut_mat_tumor, condensed = T, ymax = 0.35) + 
  #scale_y_continuous(breaks = seq(0, 0.05, 0.01))  +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_blank(),
    axis.ticks.x=element_blank(), 
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    panel.spacing.y = unit(1, "lines"),
    strip.text.y = element_blank()
  ) 
ggsave("./figures/manuscript_figures/Figure2/SMaHT_COLO829_SNV_Tumor_spectrum.pdf", 
       height = 4, width = 4, unit = "in")

plot_96_profile(mut_mat_BL, ymax = 0.075, condensed = T) + 
  scale_y_continuous(breaks = seq(0, 0.075, 0.025))  +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_blank(),
    axis.ticks.x=element_blank(), 
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    panel.spacing.y = unit(1, "lines"),
    strip.text.y = element_blank()
  )
ggsave("./figures/manuscript_figures/Figure2/SMaHT_COLO829_SNV_B_spectrum.pdf", 
       height = 4, width = 4, unit = "in")

plot_96_profile(mut_mat_BL, ymax = 0.12, condensed = T) + 
  scale_y_continuous(breaks = seq(0, 0.1, 0.05))  +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_blank(),
    axis.ticks.x=element_blank(), 
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    panel.spacing.y = unit(1, "lines"),
    strip.text.y = element_blank()
  )
ggsave("./figures/manuscript_figures/Figure2/SMaHT_COLO829_SNV_B_spectrum_SBS18.pdf", 
       height = 4, width = 4, unit = "in")
