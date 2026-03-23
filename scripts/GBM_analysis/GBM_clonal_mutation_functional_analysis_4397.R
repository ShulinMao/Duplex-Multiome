# Libraries
.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(dplyr)
library(stringr)
library(tidyr)
library(Seurat)
library(SeuratDisk)
library(Signac)
library(GenomicRanges)
library(future)
library(ggplot2)
library(harmony)
library(ComplexHeatmap)
library(biovizBase)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(foreach)
library(doParallel)

# font
font = "Helvetica"
font_size = 7
axis_font_size = 6

# Custom function to format scientific notation with superscripts and handle zero
scientific_10 <- function(x) {
  ifelse(x == 0, "0", parse(text = gsub("e", " %*% 10^", scientific_format()(x))))
}


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
  read.csv("./results/clonal_mutation/clonal_mutation_list_realigned_GBM_20250604.csv")
clonal_mutation_list[is.na(clonal_mutation_list$cell_type), "cell_type"] <- "filtered"

# cell type mutations
cell_type <- "Tumor"
clonal_mutation_list <- clonal_mutation_list[clonal_mutation_list$cell_type == "Tumor",]

clonal_mutation_count <- clonal_mutation_list %>% group_by(key) %>% summarise("count" = n())
selected_clonal_mutation <- clonal_mutation_count[clonal_mutation_count$count > 5, "key", drop = T]
clonal_mutation_list <- clonal_mutation_list[clonal_mutation_list$key %in% selected_clonal_mutation,]

# read in RNA-seq data
multiomes_integrated <- readRDS("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/RNA_ATAC/multiomes_integrated_20250604_GBM_annotated.rds")
DefaultAssay(multiomes_integrated) <- "RNA"
multiomes_integrated <- multiomes_integrated[,multiomes_integrated$celltype == "Tumor"]
multiomes_integrated <- ScaleData(multiomes_integrated, features = rownames(multiomes_integrated))
gc()

percentile_gene_summary <- list()
for (mutation in selected_clonal_mutation){
  print(mutation)
  # clone
  clone <- str_c(clonal_mutation_list[clonal_mutation_list$key == mutation, "case_id"],
                 clonal_mutation_list[clonal_mutation_list$key == mutation, "cell"],
                 sep = "_")
  cell_not_clone <- colnames(multiomes_integrated)[!(colnames(multiomes_integrated) %in% clone) & str_detect(multiomes_integrated$celltype, cell_type)]
  clone_size <- length(clone)
  
  # nearby genes
  mutation_chr <- str_split(mutation, "-", simplify = T)[1,1]
  mutation_pos <- as.integer(str_split(mutation, "-", simplify = T)[1,2])
  
  nearby_gene_range <- 5e5 # 1mb
  nearby_gene <- gene_annotation[(gene_annotation$chr == mutation_chr & 
                                    ((gene_annotation$end >= mutation_pos - nearby_gene_range & gene_annotation$end <= mutation_pos) |
                                       (gene_annotation$start <= mutation_pos + nearby_gene_range & gene_annotation$start >= mutation_pos) |
                                       (gene_annotation$end >= mutation_pos & gene_annotation$start <= mutation_pos))),
                                 "gene_name", drop = T]
  # express in > 10% of cells
  if(length(nearby_gene) == 0){
    next
  }else if(length(nearby_gene) == 1){
    expressed_nearby_gene <- (sum(multiomes_integrated@assays$RNA@counts[nearby_gene, c(clone, cell_not_clone)] > 0)/length(c(clone, cell_not_clone))) > 0.1
  }else{
    expressed_nearby_gene <- (rowSums(multiomes_integrated@assays$RNA@counts[nearby_gene, c(clone, cell_not_clone)] > 0)/length(c(clone, cell_not_clone))) > 0.1
  }
  nearby_gene[expressed_nearby_gene] -> nearby_gene
 
  if(length(nearby_gene) == 0){
    next
  }
  
  nearby_gene_expr <- multiomes_integrated@assays$RNA@data[nearby_gene,, drop = F]
  
  
  # sampling + build empirical distribution
  clone_expr <- apply(nearby_gene_expr[, clone, drop = F], 1, sum)
  pseudo_clone <- c()
  sampling_time <- 50000
  
  # Set up a parallel backend
  cl <- makeCluster(10)  # Create a cluster with 10 cores
  registerDoParallel(cl)  # Register the cluster as the parallel backend
  
  sample_expr <- function(cell_not_clone, clone_size){
    cell_not_clone_sample <- sample(cell_not_clone, clone_size)
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
  
  percentile_gene_summary[[mutation]] <- percentile_gene
  
  if (length(nearby_gene) == 1){
    pseudo_clone_expr <- data.frame(pseudo_clone_expr/clone_size)
    colnames(pseudo_clone_expr) <- nearby_gene
  }else{
    pseudo_clone_expr <- data.frame(t(pseudo_clone_expr/clone_size))    
  }
  colnames(pseudo_clone_expr) <- nearby_gene
  
  # Plot
  for (gene in nearby_gene){
    density_y <- density(pseudo_clone_expr[, gene])
    max_density_y <- max(density_y$y)
    
    # End point (head of arrow, pointing to the x-axis)
    arrow_end_x <- clone_expr[gene]/clone_size
    arrow_end_y <- density_y$y[which.min(abs(density_y$x - arrow_end_x))]  + 0.02 * max_density_y
    
    # Start point (tail of arrow)
    arrow_start_x <- clone_expr[gene]/clone_size
    arrow_start_y <- arrow_end_y + 0.15 * max_density_y
    
    percentile_5 <- quantile(pseudo_clone_expr[,gene], 0.05)
    percentile_95 <- quantile(pseudo_clone_expr[,gene], 0.95)
    
    ggplot(pseudo_clone_expr, aes(x = !!sym(gene))) +
      # Draw the density curve
      geom_density(fill = "lightblue", color = "blue", alpha = 0.7) +
      # Mark the 5th and 95th percentile with a vertical line
      geom_vline(xintercept = percentile_5, linetype = "dashed", color = "red", size = 0.5) +
      geom_vline(xintercept = percentile_95, linetype = "dashed", color = "red", size = 0.5) +
      # Add label for the 5th and 95th percentile line
      annotate("text", x = percentile_5, y = max_density_y * 0.9,
               label = paste0("5th: ", round(percentile_5, 2)),
               color = "red", hjust = -0.05, vjust = -0.5, size = 2) +
      annotate("text", x = percentile_95, y = max_density_y * 0.9,
               label = paste0("95th: ", round(percentile_95, 2)),
               color = "red", hjust = -0.05, vjust = -0.5, size = 2) +
      geom_segment(aes(x = arrow_start_x, y = arrow_start_y,
                       xend = arrow_end_x, yend = arrow_end_y),
                   arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
                   color = "darkgreen", size = 0.5) +
      annotate("text", x = arrow_start_x, y = arrow_start_y + 0.01 * max_density_y, # Place label slightly above arrow tail
               label = paste0("Mutant cells:", round(arrow_start_x, 2)),
               color = "darkgreen", hjust = 0.5, vjust = -0.5, size = 2) +
      labs(
        x = paste0(gene, " expression\n(pseudoclone)"),
        y = "Density"
      ) +
      theme_bw() +
      theme(
        text = element_text(color = "black", size = font_size, family = font),
        axis.text.x = element_text(color = "black", size = axis_font_size, family = font),
        axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
        legend.position = "none",
      )
    ggsave(paste0("./figures/functional_analysis/GBM/", mutation, "_", gene, ".pdf"), width = 6, height = 4, units = "cm")
  }
}

saveRDS(percentile_gene_summary, "./results/functional_analysis/GBM_4397_tumor_01012026.rds")
percentile_gene_summary <- readRDS("./results/functional_analysis/GBM_4397_tumor_01012026.rds")

percentile_gene_summary_df <- c()
for (mutation in names(percentile_gene_summary)){
  percentile_gene <- percentile_gene_summary[[mutation]]
  percentile_gene_summary_df <- rbind(percentile_gene_summary_df,
                                      cbind(data.frame(percentile_gene), mutation))
  if(sum(percentile_gene > 0.975 | percentile_gene < 0.025)){
    gene <- names(percentile_gene[percentile_gene > 0.975 | percentile_gene < 0.025])
    if(sum(p.adjust(ifelse(percentile_gene < 0.5, percentile_gene, 1-percentile_gene), method = "BH") < 0.05)){
      print(mutation)
      
      print(percentile_gene[p.adjust(ifelse(percentile_gene < 0.5, percentile_gene, 1-percentile_gene), method = "BH") < 0.05])
    }
  }
}


