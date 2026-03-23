setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/")
library(Seurat)
library(stringr)
library(ggplot2)
library(dplyr)

# font
font = "Helvetica"
font_size = 7
axis_font_size = 6

all_covered_cell <- readRDS("./results/clonal_mutation/clonal_mutation_site_covered_cell_GBM_08152025.rds")
all_covered_cell <- data.frame(all_covered_cell)
colnames(all_covered_cell) <- c("dataset", "mutation", "cell.id")
all_covered_cell$cell.id <- str_c(all_covered_cell$dataset, all_covered_cell$cell.id, sep = "_")

clonal_mutation_list <- 
  read.csv("./results/clonal_mutation/clonal_mutation_list_realigned_GBM_20250604.csv")
clonal_mutation_list[(duplicated(str_c(clonal_mutation_list$cell, clonal_mutation_list$key, clonal_mutation_list$case_id))),]
clonal_mutation_list$cell <- str_c(clonal_mutation_list$case_id, clonal_mutation_list$cell, sep="_")
clonal_mutation_list[is.na(clonal_mutation_list$cell_type), "cell_type"] <- "filtered"
unique(clonal_mutation_list$key) -> clonal_mutations

multiomes_integrated <- readRDS("./data/RNA_ATAC/multiomes_integrated_20250604_GBM_annotated.rds")
tumor <- readRDS("./data/RNA_ATAC/multiomes_integrated_20250604_GBM_tumor_only_annotated.rds")

for (mutation in clonal_mutations){
  covered_cell <- all_covered_cell[all_covered_cell$mutation == mutation, "cell.id"]
  detected_cell <- clonal_mutation_list[clonal_mutation_list$key == mutation, "cell"]
  
  # all cells
  multiomes_integrated[["Coverage"]] <- "No coverage"
  multiomes_integrated$Coverage[colnames(multiomes_integrated) %in% covered_cell] <- "Covered"
  multiomes_integrated$Coverage[colnames(multiomes_integrated) %in% detected_cell] <- "Alt allele detected"
  
  mutation_map <- cbind(multiomes_integrated@reductions$umap.rna@cell.embeddings, multiomes_integrated[["Coverage"]])
  
  ggplot(mutation_map %>% arrange(desc(Coverage)), aes(x = rnaUMAP_1, y = rnaUMAP_2))+
    geom_point(aes(shape = Coverage, size = Coverage, color = Coverage)) +
    scale_color_manual(values = c("red3", "gold1", "grey")) +
    scale_shape_manual(values = c(8, 16, 16)) +
    scale_size_manual(values = c(0.5, 0.65, 0.1)) +
    guides(color=guide_legend(nrow=4, byrow=T, override.aes = list(size=2))) +
    labs(subtitle = paste(mutation)) +
    theme_classic() +
    theme(
      text = element_text(color = "black", size = font_size, family = font),
      axis.text = element_text(color = "black", size = axis_font_size, family = font),
      legend.text = element_text(color = "black", size = axis_font_size, family = font),
      legend.position = "bottom"
    )
  ggsave(paste0("./figures/GBM_clonal_variants/all_cell_types_4397/", mutation, ".pdf"), height = 3.6, width = 2.3, units = "in")
  
  # tumor
  tumor[["Coverage"]] <- "No coverage"
  tumor$Coverage[colnames(tumor) %in% covered_cell] <- "Covered"
  tumor$Coverage[colnames(tumor) %in% detected_cell] <- "Alt allele detected"

  if(sum(tumor$Coverage == "Alt allele detected") == 0){next}

  mutation_map <- cbind(tumor@reductions$umap.tumor.rna@cell.embeddings, tumor[["Coverage"]])

  ggplot(mutation_map %>% arrange(desc(Coverage)), aes(x = UMAP_1, y = UMAP_2))+
    geom_point(aes(shape = Coverage, size = Coverage, color = Coverage)) +
    scale_color_manual(values = c("red3", "gold1", "grey")) +
    scale_shape_manual(values = c(8, 16, 16)) +
    scale_size_manual(values = c(0.5, 0.65, 0.1)) +
    labs(subtitle = paste(mutation)) +
    guides(color=guide_legend(nrow=4, byrow=T, override.aes = list(size=2))) +
    theme_classic() +
    theme(
      text = element_text(color = "black", size = font_size, family = font),
      axis.text = element_text(color = "black", size = axis_font_size, family = font),
      legend.text = element_text(color = "black", size = axis_font_size, family = font),
      legend.position = "bottom"
    )
  ggsave(paste0("./figures/GBM_clonal_variants/tumor_only_4397/", mutation, ".pdf"), height = 3.6, width = 2.3, units = "in")
}
