# Libraries
#.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
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
library(foreach)
library(doParallel)

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


# ------------------------------------------------------------------------------
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
  coord_cartesian(ylim = c(6e-10, 1.01e-9)) +  # "Upper" segment
  scale_y_continuous(breaks = seq(7.5e-10, 1e-9, 2.5e-10), labels = scientific_10) +
  scale_x_discrete(expand = expansion(mult = c(1.1, 1.1))) + 
  annotate(geom="text", x=1, y=1e-9, parse = T, 
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


# ------------------------------------------------------------------------------
# UMB4638 chr2 19990244 G A
# scRNAseq
brain <- readRDS("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/RNA_ATAC/multiomes_integrated_20240512_infant_adult_annotated.rds")

brain[,!str_detect(brain$celltype, "doublet")] -> brain
brain[["cell_type"]] <- as.character(brain$celltype)
brain$cell_type[str_detect(brain$celltype, "EN")] <- "Excitatory\nNeuron"
brain$cell_type[str_detect(brain$celltype, "IN")] <- "Inhibitory\nNeuron"

# UMAP of UMB4638
UMB4638 <- brain[, brain$dataset == "UMB4638"]

DimPlot(UMB4638, reduction = "wnn.umap", group.by = "cell_type", pt.size = 1e-10) + 
  labs(title = NULL, color = "Cell type") +
  guides(color=guide_legend(nrow=3, byrow=T, override.aes = list(size=2))) +
  labs(title = "UMB4638") +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "bottom",
    legend.key.size = unit(0, 'in'),
    legend.title = element_text(angle = 90, hjust = 0.3, vjust = 2.5)
  )
ggsave("./figures/Figure5/UMB4638_UMAP_cell_type.pdf", height = 3.7, width = 3, units = "in")

# clonal mutation
clonal_mutations <- read.csv("./results/clonal_mutation/clonal_mutation_list_healthy_brain_realigned_20241028.csv")

# UMAP of the mutation
# UMB4638 chr2:19990244 G>A
case_id <- "UMB4638"
mutation <- "chr2-19990244-G-A"
detected <- clonal_mutations$cell[clonal_mutations$case_id == case_id & clonal_mutations$key == mutation]
detected <- paste0(case_id, "_", detected)
detected <- detected[detected %in% colnames(brain)]
covered <- readRDS(paste0("./figures/manuscript_figures/Figure5/", mutation, "_covered_cell.rds"))
covered <- paste0(case_id, "_", covered)
covered <- covered[covered %in% colnames(brain)]

vaf <- round(length(detected)/length(covered), 3)

brain[["Coverage"]] <- "No coverage"
brain$Coverage[colnames(brain) %in% covered] <- "Covered"
brain$Coverage[colnames(brain) %in% detected] <- "Alt allele detected"

brain_mutation_map <- cbind(brain@reductions$wnn.umap@cell.embeddings, brain[["Coverage"]])
brain_mutation_map <- brain_mutation_map[str_starts(row.names(brain_mutation_map), case_id),]

ggplot(brain_mutation_map %>% arrange(desc(Coverage)), aes(x = wnnUMAP_1, y = wnnUMAP_2))+
  geom_point(aes(shape = Coverage, size = Coverage, color = Coverage)) +
  annotate(geom="text", x = -9, y=18, label = paste0("Duplex-Multiome VAF = ", vaf), 
           color = "black", family = font, size = 2.5, hjust = 0) +
  annotate(geom="text", x = -9, y=16, label = paste0("Bulk WGS VAF = 0.098"), 
           color = "black", family = font, size = 2.5, hjust = 0) +
  scale_color_manual(values = c("red3", "gold1", "grey")) +
  scale_shape_manual(values = c(8, 16, 16)) +
  scale_size_manual(values = c(2, 1, 1)) +
  labs(subtitle = paste(mutation, "(UMB4638)")) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "bottom",
    legend.title = element_blank()
  )
ggsave("./figures/Figure5/UMB4638_UMAP_chr2_19990244_G_A.pdf", 
       height = 3.5, width = 3, units = "in")

# ------------------------------------------------------------------------------
# UMB4638 chr2:19990244 G>A and LAPTM4A expression
# read in gtf gene annotation
gene_annotation <- read.table("../ref/gene_annotation/genes.gtf", sep = "\t")
colnames(gene_annotation) <- c("chr", "pipeline", "type", "start", "end", "V6", "V7", "V8", "info")
gene_name <- str_extract(gene_annotation$info, "gene_name [^;]*")
gene_name <- str_split(gene_name, " ", simplify = T)[,2]
gene_annotation[, "gene_name"] <- gene_name

gene_annotation <- gene_annotation %>% group_by(gene_name) %>% summarise("chr" = last(chr),
                                                                         "start" = mean(start),
                                                                         "end" = max(end))

# Read in clonal mutations after the germline filter
clonal_mutation_list <- 
  read.csv("./results/clonal_mutation/clonal_mutation_list_filtered_ATAC_celltype_20240812.csv")
clonal_mutation_list <- clonal_mutation_list[!is.na(clonal_mutation_list$cell_type),]
clonal_mutation_list[(duplicated(str_c(clonal_mutation_list$cell, clonal_mutation_list$key, clonal_mutation_list$case_id))),]
selected_clonal_mutation <- "chr2-19990244-G-A"

# read in RNA-seq data
multiomes_integrated <- readRDS("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/RNA_ATAC/multiomes_integrated_20240512_infant_adult_annotated.rds")
DefaultAssay(multiomes_integrated) <- "RNA"
multiomes_integrated <- ScaleData(multiomes_integrated, features = rownames(multiomes_integrated))
gc()
multiomes_integrated$celltype <- as.character(multiomes_integrated$celltype)
multiomes_integrated$celltype[str_detect(multiomes_integrated$celltype, "EN")] <- "EN"
multiomes_integrated$celltype[str_detect(multiomes_integrated$celltype, "IN")] <- "IN"

case_id <- "UMB4638"
print(selected_clonal_mutation)
# clone
clone <- str_c(case_id,
               clonal_mutation_list[clonal_mutation_list$key == selected_clonal_mutation, "cell"],
               sep = "_")
print(multiomes_integrated$celltype[colnames(multiomes_integrated) %in% clone])


cell_not_clone <- colnames(multiomes_integrated)[!(colnames(multiomes_integrated) %in% clone) & 
                                                   str_detect(multiomes_integrated$dataset, case_id)]
cell_not_clone_celltype <- multiomes_integrated$celltype[cell_not_clone]
clone_size <- length(clone)

# nearby genes
mutation_chr <- str_split(selected_clonal_mutation, "-", simplify = T)[1,1]
mutation_pos <- as.integer(str_split(selected_clonal_mutation, "-", simplify = T)[1,2])

nearby_gene_range <- 5e5 # 100 kb
nearby_gene <- gene_annotation[(gene_annotation$chr == mutation_chr & 
                                  ((gene_annotation$end >= mutation_pos - nearby_gene_range & gene_annotation$end <= mutation_pos) |
                                     (gene_annotation$start <= mutation_pos + nearby_gene_range & gene_annotation$start >= mutation_pos) |
                                     (gene_annotation$end >= mutation_pos & gene_annotation$start <= mutation_pos))),
                               "gene_name", drop = T]

# express in > 10% of cells
(rowSums(multiomes_integrated@assays$RNA@counts[nearby_gene, c(clone, cell_not_clone)] > 0)/length(c(clone, cell_not_clone)))
expressed_nearby_gene <- (rowSums(multiomes_integrated@assays$RNA@counts[nearby_gene, c(clone, cell_not_clone)] > 0)/length(c(clone, cell_not_clone))) > 0.1
nearby_gene[expressed_nearby_gene] -> nearby_gene
nearby_gene_expr <- multiomes_integrated@assays$RNA@data[nearby_gene, c(clone, cell_not_clone), drop = F]

# sampling + build empirical distribution
clone_expr <- apply(nearby_gene_expr[, clone, drop = F], 1, sum)
pseudo_clone <- c()
sampling_time <- 10000

# Set up a parallel backend
cl <- makeCluster(10)  # Create a cluster with 10 cores
registerDoParallel(cl)  # Register the cluster as the parallel backend

sample_expr <- function(cell_not_clone, clone_size){
  cell_not_clone_sample <- c(sample(cell_not_clone[cell_not_clone_celltype == "EN"], 5), 
                             sample(cell_not_clone[cell_not_clone_celltype == "IN"], 3),
                             sample(cell_not_clone[cell_not_clone_celltype == "Oligodendrocyte"], 1),
                             sample(cell_not_clone[cell_not_clone_celltype == "OPC"], 1),
                             sample(cell_not_clone[cell_not_clone_celltype == "Microglia"], 5))
  
  return(cell_not_clone_sample)
}

pseudo_clone <- foreach(i = 1:sampling_time, .combine = cbind) %dopar% {
  sample_expr(cell_not_clone, clone_size)
}

# Clean up the parallel backend
stopCluster(cl)

pseudo_clone_expr <- apply(pseudo_clone, 2, function(X){apply(nearby_gene_expr[, X, drop = F], 1, sum)})

percentile_gene <- c()
if(length(nearby_gene)==1){
  percentile_gene[nearby_gene] <- (rank(c(pseudo_clone_expr, clone_expr[nearby_gene]))[sampling_time+1]-1)/sampling_time
} else {
  for (gene in nearby_gene){
    percentile_gene[gene] <- (rank(c(pseudo_clone_expr[gene,], clone_expr[gene]))[sampling_time+1]-1)/sampling_time
  }
}

if(sum(percentile_gene >= 0.95 | percentile_gene <= 0.05)){
  print(percentile_gene[percentile_gene >= 0.95 | percentile_gene <= 0.05])
}
print(selected_clonal_mutation)
print(percentile_gene)
p.adjust(ifelse(percentile_gene > 0.5, 1-percentile_gene, percentile_gene), method = "BH")

pseudo_clone_expr <- data.frame(t(pseudo_clone_expr/15))
max_density_y <- max(density(pseudo_clone_expr$LAPTM4A)$y)

# arrow (mutant cell expr)
arrow_start_x <- clone_expr["LAPTM4A"]/15 # there are 15 cells
arrow_start_y <- max_density_y * 0.2
arrow_end_x <- clone_expr["LAPTM4A"]/15
arrow_end_y <- 0 

percentile_95 <- quantile(pseudo_clone_expr[,"LAPTM4A"], 0.95)
ggplot(pseudo_clone_expr, aes(x = LAPTM4A)) +
  geom_density(fill = "lightblue", color = "blue", alpha = 0.7) +
  geom_vline(xintercept = percentile_95, linetype = "dashed", color = "red", size = 0.5) +
  annotate("text", x = percentile_95, y = max_density_y * 0.9,
           label = paste0("95th Percentile: ", round(percentile_95, 2)),
           color = "red", hjust = -0.05, vjust = -0.5, size = 2) +
  geom_segment(aes(x = arrow_start_x, y = arrow_start_y,
                   xend = arrow_end_x, yend = arrow_end_y),
               arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
               color = "darkgreen", size = 0.5) +
  annotate("text", x = arrow_start_x, y = max_density_y * 0.25,
           label = paste0("Mutant cells:", round(arrow_start_x, 2)),
           color = "darkgreen", hjust = 0.5, vjust = -0.5, size = 2) +
  labs(
    x = "LAPTM4A expression\n(pseudoclone)",
    y = "Density"
  ) +
  theme_bw() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
  )
ggsave("./figures/Figure5/UMB4638_clonal_mutation_LAPTM4A.pdf", width = 6, height = 4, units = "cm")
