# Libraries
.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(Seurat)
library(ggplot2)
library(scales)
library(ggsci)
library(dplyr)
library(stringr)

# font
font = "Helvetica"
font_size = 7
axis_font_size = 6

# Custom function to format scientific notation with superscripts and handle zero
scientific_10 <- function(x) {
  ifelse(x == 0, "0", parse(text = gsub("e\\+?", " %*% 10^", scientific_format()(x))))
}

# ------------------------------------------------------------------------------
# UMAP 4397
multiomes_integrated <- readRDS("./data/RNA_ATAC/multiomes_integrated_20250604_GBM_annotated.rds")
multiomes_integrated$cell_type <- as.character(multiomes_integrated$cell_type)
multiomes_integrated$cell_type[multiomes_integrated$cell_type == "Oligodendorcyte"] <- "Oligodendrocyte"

DimPlot(multiomes_integrated, reduction = "umap.rna", 
        group.by = "cell_type", pt.size = 1e-10
) + 
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
ggsave("./figures/Figure6/GBM_all_UMAP.pdf", height = 3.3, width = 2.3, units = "in")

tumor <- readRDS("./data/RNA_ATAC/multiomes_integrated_20250604_GBM_tumor_only_annotated.rds")
DimPlot(tumor, reduction = "umap.tumor.rna", 
        group.by = "tumor_seurat_clusters", pt.size = 1e-10
) + 
  labs(title = NULL, color = "Tumor cell\nsubcluster") +
  guides(color=guide_legend(nrow=4, byrow=T, override.aes = list(size=2))) +
  scale_color_d3() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "bottom",
    #legend.key.size = unit(0, 'in'),
    #legend.title = element_text(angle = 90, hjust = 0.5, vjust = 2.5)
  )
ggsave("./figures/Figure6/GBM_tumor_UMAP.pdf", height = 3.3, width = 2.3, units = "in")

# Cell type composition
brain_cell_type <- read.table("./results/celltype_annotation/annotation_table_GBM_20250604.txt", sep = "\t")
colnames(brain_cell_type) <- c("case_id", "cell_id", "cell_type")
brain_cell_type$case_id <- factor(brain_cell_type$case_id, 
                                  levels = c("UMB4397_tumor", "UMB4397_normal"), 
                                  labels = c("More infiltrated\nby tumor cells", "Less infiltrated\nby tumor cells"))

ggplot(brain_cell_type, aes(x = case_id)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 2) +
  geom_bar(aes(fill = cell_type)) +
  theme_classic() +
  labs(x = "Brain region", y = "Cell number", fill = "Cell type") +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.key.size = unit(0.1, 'in')
  )
ggsave("./figures/Figure6/GBM_tumor_cell_type_composition.pdf", height = 2, width = 2.2, units = "in")


# ------------------------------------------------------------------------------
# burden
load("./results/burden_simulation/GBM/UMB4397_burden_UMAP.Rdata")
burden_umap <- cbind(tumor@reductions$umap.tumor.rna@cell.embeddings, tumor[[c("burden_smooth", "burden_smooth_top")]])
ggplot(burden_umap, aes(x = UMAP_1, y = UMAP_2))+
  geom_point(aes(color = burden_smooth), size = 0.3) +
  labs(color = "Mutational burden\n(per base)") +
  scale_color_gradientn(colours = c('#fdf4af', '#f9b64b', '#a51a49'), labels = scientific_10, breaks = breaks_extended(3)) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "bottom",
    legend.title = element_text(color = "black", size = axis_font_size, family = font)
  )
ggsave("./figures/Figure6/GBM_tumor_burden_UMAP.pdf", height = 3, width = 2.3, units = "in")

median(burden_umap$burden_smooth)
q95 <- quantile(burden_umap$burden_smooth, 0.95)
dens <- density(burden_umap$burden_smooth)
dens_df <- data.frame(x = dens$x, y = dens$y)
ggplot(burden_umap, aes(x = burden_smooth)) +
  geom_density() +
  #scale_y_continuous(position = "right") + 
  geom_area(data = subset(dens_df, x > q95),
            aes(x, y),
            fill = "darkred",
            alpha = 0.6) +
  geom_vline(xintercept = q95, linetype = "dashed", color = "darkred", size = 0.5) +
  annotate("text", x = q95*1.06, y = max(dens_df$y), label = "95%",
           hjust = 0.05, color = "darkred", size = 2) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 90, vjust = 0.5),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "bottom",
    legend.title = element_text(color = "black", size = axis_font_size, family = font)
  ) + 
  labs(x = "Mutational burden (per base)", y = "Density") +
  coord_flip() +
  scale_y_reverse(labels = scientific_10) +
  scale_x_continuous(position = "top", labels = scientific_10)
ggsave("./figures/Figure6/GBM_tumor_burden_density.pdf", height = 2.3, width = 1, units = "in")

ggplot(burden_umap, aes(x = UMAP_1, y = UMAP_2))+
  geom_point(aes(color = burden_smooth_top), size = 0.3) +
  labs(color = "Mutational burden Top 5%") +
  scale_color_manual(values = c("snow3", "darkred")) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "bottom",
    legend.title = element_text(color = "black", size = axis_font_size, family = font)
  )
ggsave("./figures/Figure6/GBM_tumor_burden_top5_UMAP.pdf", height = 2.9, width = 2.3, units = "in")

# ------------------------------------------------------------------------------
# high-burden cluster DEGs
DEG_burden_smooth_top[DEG_burden_smooth_top$p_val_adj < 0.05 & DEG_burden_smooth_top$avg_log2FC > 0.25,] -> DEG_up
DEG_burden_smooth_top[DEG_burden_smooth_top$p_val_adj < 0.05 & abs(DEG_burden_smooth_top$avg_log2FC) > 0.25,]
write.csv(DEG_burden_smooth_top, "./DEG_high_burden_OPC_NPC_to_Luis.csv", quote = F)

go_bp_results <- enrichGO(
  gene          = rownames(DEG_up),
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP", # Biological Process
  pAdjustMethod = "BH", # Benjamini-Hochberg for FDR
  pvalueCutoff  = 0.05
)

go_bp_results <- simplify(go_bp_results, cutoff=0.5, by="p.adjust", select_fun=min)
go_bp_results <- go_bp_results@result[go_bp_results@result$p.adjust < 0.05,]
go_bp_results$Description <- factor(go_bp_results$Description, 
                                    levels = go_bp_results$Description[order(go_bp_results$p.adjust, decreasing = T)])
ggplot(go_bp_results, aes(x = -log10(p.adjust), y = Description)) +
  # remove epithelium cell related terms since they share the same processes with ECs
  geom_col(fill = "steelblue1") +
  theme_classic() +
  geom_vline(xintercept = -log10(0.05), linetype="dashed", 
             color = "black", size=0.5) +
  labs(title = "GO Biological Processes (up-regulated genes)", x = "-log10(p.adj)", y = NULL) +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 90, vjust = 0.5),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(color = "black", size = axis_font_size, family = font, hjust = 1.5)
  )
ggsave("./figures/Figure6/GBM_tumor_burden_top5_DEG_GO.pdf", height = 2.8, width = 3, units = "in")


# ------------------------------------------------------------------------------
# Functional analysis and mutation spatial distribution analysis
# Please find the source code under "scripts/GBM_analysis"
