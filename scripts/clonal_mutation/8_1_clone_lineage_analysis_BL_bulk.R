.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(dplyr)
library(stringr)
library(Polychrome)
library(ggplot2)
library(pheatmap)
library(igraph)

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


#-------------------------------------------------------------------------------
# Read in clonal mutations after the germline filter
clonal_mutation_list <- read.csv("results/clonal_mutation/clonal_mutation_list_COLO829BLT50_BL_bulk_realigned_20241007.csv")
clonal_mutation_list$key <- str_replace_all(clonal_mutation_list$key, pattern = ":", replacement = "-")
clonal_mutation_list$cell <- str_c(clonal_mutation_list$case_id, clonal_mutation_list$cell, sep = "_")

clonal_mutation_count <- clonal_mutation_list %>% group_by(key) %>% summarise("count" = n())
selected_clonal_mutation <- clonal_mutation_count[clonal_mutation_count$count >= 5, "key", drop = T]
#selected_clonal_mutation <- clonal_mutation_count[clonal_mutation_count$count > 3, "key", drop = T]
clonal_mutation_list <- clonal_mutation_list[clonal_mutation_list$key %in% selected_clonal_mutation,]
clonal_mutation_list[is.na(clonal_mutation_list$cell_type), "cell_type"] <- "filtered"

#-------------------------------------------------------------------------------
adj_matrix <- matrix(0, nrow = length(unique(clonal_mutation_list$cell)), ncol = length(unique(clonal_mutation_list$cell)))
rownames(adj_matrix) <- unique(clonal_mutation_list$cell)
colnames(adj_matrix) <- unique(clonal_mutation_list$cell)
for (cell_i in rownames(adj_matrix)){
  cell_i_mutation <- clonal_mutation_list[clonal_mutation_list$cell == cell_i, "key"]
  for (cell_j in colnames(adj_matrix)){
    cell_j_mutation <- clonal_mutation_list[clonal_mutation_list$cell == cell_j, "key"]
    adj_matrix[cell_i, cell_j] <- sum(cell_i_mutation %in% cell_j_mutation)
  }
}
diag(adj_matrix) <- 0

# Function to perform DFS and check if the depth exceeds 4
dfs <- function(adj_matrix, start_node){
  n <- ncol(adj_matrix)
  visited <- rep(FALSE, n) # Vector to track visited nodes
  
  dfs_util <- function(node, current_depth, visited) {
    visited[node] <- TRUE
    if (current_depth >= 5) {
      return(TRUE)
    }
    
    for (neighbor in 1:n) {
      if (adj_matrix[node, neighbor] > 0 & !visited[neighbor]) {
        if (dfs_util(neighbor, current_depth + 1, visited)) {
          return(TRUE)
        }
      }
    }
    return(FALSE)
  }
  return(dfs_util(start_node, 0, visited))
}

connected_cell_count <- c()
for (cell in 1:nrow(adj_matrix)){
  connected_cell_count[cell] <- dfs(adj_matrix, cell)
}
adj_matrix[connected_cell_count, connected_cell_count] -> adj_matrix

cell_type <- c()
case <- c()
for (cell in rownames(adj_matrix)){
  cell_type[cell] <- clonal_mutation_list[(clonal_mutation_list$cell == cell), "cell_type"][1]
  case[cell] <- clonal_mutation_list[(clonal_mutation_list$cell == cell), "case_id"][1]
}

g <- graph_from_adjacency_matrix(adj_matrix, mode="undirected", weighted = T)
V(g)$label <- ""
V(g)$size <- 4
layout <- layout_with_fr(g, dim = 2, niter = 2000, start.temp = 100)

V(g)$shape <- as.character(factor(cell_type, levels = c("B_lymphoblast", "melanoma_fibroblast", "filtered"), labels = c("square", "circle", "triangle")))
V(g)$color <- as.character(factor(cell_type, levels = c("B_lymphoblast", "melanoma_fibroblast", "filtered"), labels = c("cadetblue1", "coral1", "whitesmoke")))
pdf("./figures/lineage_tracing/graph/cell_type_BL_bulk_20241007.pdf", width = 10, height = 10)
plot(g, layout=layout)
dev.off()

V(g)$color <- as.character(factor(case, levels = c("COLO829BLT50_rep1", "COLO829BLT50_rep2"), labels = c("blue", "white")))
pdf("./figures/lineage_tracing/graph/rep_BL_bulk_20241007.pdf", width = 10, height = 10)
plot(g, layout=layout)
dev.off()


# # Use the Louvain method to detect communities
# louvain_communities <- cluster_louvain(g)
# table(cell_type, membership(louvain_communities)[names(cell_type)])
# louvain_cluster <- membership(louvain_communities)[names(cell_type)]
# V(g)$color <- as.character(factor(louvain_cluster, labels = createPalette(length(unique(louvain_cluster)), c("#ff0000", "#00ff00", "#0000ff"))))
# plot(g, layout=layout, main="Louvain")
# 
# # Use the edge betweenness method to detect communities
# edge_betweenness_communities <- cluster_edge_betweenness(g)
# table(cell_type, membership(edge_betweenness_communities)[names(cell_type)])
# edge_betweenness_cluster <- membership(edge_betweenness_communities)[names(cell_type)]
# V(g)$color <- as.character(factor(edge_betweenness_cluster, labels = createPalette(length(unique(edge_betweenness_cluster)), c("#ff0000", "#00ff00", "#0000ff"))))
# plot(g, layout=layout, main="Edge Betweenness")

# Use the label propagation method to detect communities
label_prop_communities <- cluster_label_prop(g)
table(cell_type, membership(label_prop_communities)[names(cell_type)])
label_prop_cluster <- membership(label_prop_communities)[names(cell_type)]
label_prop_color <- palette36.colors(15)
V(g)$color <- as.character(factor(label_prop_cluster, labels = label_prop_color))
pdf("./figures/lineage_tracing/graph/label_prop_clustering_BL_bulk_20241007.pdf", width = 10, height = 10)
plot(g, layout=layout, main="Label Propagation")
dev.off()

# Use the Walktrap method to detect communities
# walktrap_communities <- cluster_walktrap(g)
# table(cell_type, membership(walktrap_communities)[names(cell_type)])
# walktrap_cluster <- membership(walktrap_communities)[names(cell_type)]
# V(g)$color <- as.character(factor(walktrap_cluster, labels = createPalette(length(unique(walktrap_cluster)), c("#ff0000", "#00ff00", "#0000ff"))))
# plot(g, layout=layout, main="Walktrap")

#-------------------------------------------------------------------------------
df_label_prop_cluster <- data.frame(label_prop_cluster)
df_label_prop_cluster[,"cell"] <- row.names(df_label_prop_cluster)
df_label_prop_cluster <- inner_join(clonal_mutation_list, df_label_prop_cluster, by = "cell")
df_label_prop_cluster$cell_type <- factor(df_label_prop_cluster$cell_type, levels = c("B_lymphoblast", "melanoma_fibroblast", "filtered"))

df_label_prop_cluster_count <- df_label_prop_cluster[!(duplicated(str_c(df_label_prop_cluster$cell, df_label_prop_cluster$label_prop_cluster, df_label_prop_cluster$case_id))),]
ggplot(df_label_prop_cluster_count, aes(x=label_prop_cluster, fill=cell_type, group=cell_type)) +
  
  geom_bar(stat="count") +
  scale_x_continuous(breaks=seq(1, length(unique(df_label_prop_cluster$label_prop_cluster))),
                   labels=seq(1, length(unique(df_label_prop_cluster$label_prop_cluster)))) +
  theme_classic() +
  scale_fill_manual(values=c("cadetblue1", "coral1", "whitesmoke")) +
  labs(x = "Clusters", y = "# of cells", fill = "Cell Type")
ggsave("./figures/lineage_tracing/label_prop_cluster_BL_bulk_20241007.pdf", width = 8, height = 4)

clonal_mutation_matrix <- matrix(0, nrow = length(unique(df_label_prop_cluster$cell)), ncol = length(unique(df_label_prop_cluster$key)))
rownames(clonal_mutation_matrix) <- unique(df_label_prop_cluster$cell)
colnames(clonal_mutation_matrix) <- unique(df_label_prop_cluster$key)
for (cell in rownames(clonal_mutation_matrix)){
  cell_mutation <- df_label_prop_cluster[df_label_prop_cluster$cell == cell, "key"]
  clonal_mutation_matrix[cell,] <- colnames(clonal_mutation_matrix) %in% cell_mutation
}
cell_order <- unique(df_label_prop_cluster$cell[order(df_label_prop_cluster$label_prop_cluster, decreasing = F)])
clonal_mutation_matrix <- clonal_mutation_matrix[cell_order,]

annotation_col <- df_label_prop_cluster[, c("cell", "label_prop_cluster", "cell_type")]
annotation_col <- annotation_col[!duplicated(annotation_col$cell),]
row.names(annotation_col) <- annotation_col$cell
annotation_col <- annotation_col[cell_order, c("label_prop_cluster", "cell_type")]
annotation_col$label_prop_cluster <- factor(annotation_col$label_prop_cluster,
                                            levels = seq(1, length(unique(df_label_prop_cluster$label_prop_cluster))))
annotation_colors <- list(cell_type = c("B_lymphoblast" = "cadetblue1", 
                                        "melanoma_fibroblast" = "coral1", 
                                        "filtered" = "whitesmoke"),
                          label_prop_cluster = label_prop_color)
names(annotation_colors[["label_prop_cluster"]]) <- seq(1, length(unique(df_label_prop_cluster$label_prop_cluster)))

pheatmap(t(clonal_mutation_matrix), scale = "none", show_colnames = F,
         cluster_rows = T, cluster_cols = F, legend = F,
         annotation_col = annotation_col, 
         annotation_colors = annotation_colors,
         #border_color = NA,
         cellwidth = 1, 
         cellheight = 10)
         #filename = "figures/lineage_tracing/heatmap_label_prop_cluster_BL_bulk_20240604.pdf", height = 10)

#-------------------------------------------------------------------------------
# Read in clonal mutations after the germline filter (tumor bulk)
# clonal_mutation_list_tumor_bulk <- read.csv("results/clonal_mutation/clonal_mutation_list_COLO829BLT50_tumor_bulk_20240604.csv")
# clonal_mutation_list_tumor_bulk$key <- str_replace_all(clonal_mutation_list_tumor_bulk$key, pattern = ":", replacement = "-")
# clonal_mutation_list_tumor_bulk[clonal_mutation_list_tumor_bulk$cell_type == "B_lymphoblast",] -> clonal_mutation_list_tumor_bulk
# clonal_mutation_list_tumor_bulk$cell <- str_c(clonal_mutation_list_tumor_bulk$case_id, clonal_mutation_list_tumor_bulk$cell, sep = "_")
# 
# clonal_mutation_list_tumor_bulk_cell_selected <- clonal_mutation_list_tumor_bulk[clonal_mutation_list_tumor_bulk$cell %in% rownames(clonal_mutation_matrix),]
# 
# clonal_mutation_matrix <- matrix(0, nrow = length(unique(df_label_prop_cluster$cell)), ncol = length(unique(clonal_mutation_list_tumor_bulk_cell_selected$key)))
# rownames(clonal_mutation_matrix) <- unique(df_label_prop_cluster$cell)
# colnames(clonal_mutation_matrix) <- unique(clonal_mutation_list_tumor_bulk_cell_selected$key)
# for (cell in rownames(clonal_mutation_matrix)){
#   cell_mutation <- clonal_mutation_list_tumor_bulk_cell_selected[clonal_mutation_list_tumor_bulk_cell_selected$cell == cell, "key"]
#   clonal_mutation_matrix[cell,] <- colnames(clonal_mutation_matrix) %in% cell_mutation
# }
# cell_order <- unique(df_label_prop_cluster$cell[order(df_label_prop_cluster$label_prop_cluster, decreasing = F)])
# clonal_mutation_matrix <- clonal_mutation_matrix[cell_order,]
# 
# annotation_col <- df_label_prop_cluster[, c("cell", "label_prop_cluster", "cell_type")]
# annotation_col <- annotation_col[!duplicated(annotation_col$cell),]
# row.names(annotation_col) <- annotation_col$cell
# annotation_col <- annotation_col[cell_order, c("label_prop_cluster", "cell_type")]
# annotation_col$label_prop_cluster <- factor(annotation_col$label_prop_cluster,
#                                             levels = seq(1, length(unique(df_label_prop_cluster$label_prop_cluster))))
# annotation_colors <- list(cell_type = c("B_lymphoblast" = "cadetblue1", 
#                                         "melanoma_fibroblast" = "coral1", 
#                                         "filtered" = "whitesmoke"),
#                           label_prop_cluster = label_prop_color)
# names(annotation_colors[["label_prop_cluster"]]) <- seq(1, length(unique(df_label_prop_cluster$label_prop_cluster)))
# 
# pheatmap(t(clonal_mutation_matrix), scale = "none", show_colnames = F,
#          cluster_rows = T, cluster_cols = F, legend = F,
#          annotation_col = annotation_col, 
#          annotation_colors = annotation_colors)
# dev.off()

# ------------------------------------------------------------------------------
# collect detected/covered
case_id_list <- unique(df_label_prop_cluster$case_id)
cell_list <- unique(df_label_prop_cluster$cell)

clonal_mutation_unique <- df_label_prop_cluster[!duplicated(df_label_prop_cluster$key),]
clonal_mutation_grange <- GRanges(seqnames = clonal_mutation_unique$chr,
                                  ranges = IRanges(start = clonal_mutation_unique$pos,
                                                   end = clonal_mutation_unique$pos),
                                  ref = clonal_mutation_unique$ref,
                                  alt = clonal_mutation_unique$alt)

cell_cover_clonal_mutation_site <- matrix(0, nrow = length(unique(df_label_prop_cluster$cell)), ncol = length(unique(df_label_prop_cluster$key)))
rownames(cell_cover_clonal_mutation_site) <- unique(df_label_prop_cluster$cell)
colnames(cell_cover_clonal_mutation_site) <- unique(df_label_prop_cluster$key)

count <- 1
for(cell_id in cell_list){
  case_id <- str_extract(cell_id, "COLO829BLT50_rep[12]")
  # read in covered region info 
  dir_path <- paste0("./data/", case_id, "/single_cell/read_families/a2s0_0minReadFamilyLength/")
  
  
  callable_region_file <- paste0(str_remove(cell_id, "COLO829BLT50_rep[12]_"), ".families.analysis.calledSites.bed")
  
  # callable regions
  # if no callable region, skip the cell
  if (is.na(file.info(paste0(dir_path, callable_region_file))$size)){next}
  if (file.info(paste0(dir_path, callable_region_file))$size == 0){next}
  
  callable_region <- read.table(paste0(dir_path, callable_region_file))
  callable_region_grange <- GRanges(seqnames = callable_region$V1,
                                    ranges = IRanges(start = callable_region$V2,
                                                     end = callable_region$V3))
  
  # find cell with reads covering germline mutations
  intersection <- findOverlaps(clonal_mutation_grange, callable_region_grange)
  intersected_ranges <- clonal_mutation_grange[queryHits(intersection)]
  if(length(intersected_ranges) > 0){
    covered_clonal_mutation_site <- 
      str_c(intersected_ranges@seqnames, intersected_ranges@ranges@start, intersected_ranges$ref, intersected_ranges$alt, sep = "-")
    cell_cover_clonal_mutation_site[cell_id, covered_clonal_mutation_site] <- 1
  }
  
  count <- count + 1
  if(count%%50 == 0){print(paste(case_id, count))}
}

cell_cover_clonal_mutation_site[row.names(clonal_mutation_matrix), colnames(clonal_mutation_matrix)] -> cell_cover_clonal_mutation_site
colSums(clonal_mutation_matrix == 1)/colSums(cell_cover_clonal_mutation_site == 1)
pheatmap(t(clonal_mutation_matrix+cell_cover_clonal_mutation_site), scale = "none", show_colnames = F,
         cluster_rows = T, cluster_cols = F, legend = F,
         annotation_col = annotation_col, 
         annotation_colors = annotation_colors,
         border_color = NA)

cbind(colSums(clonal_mutation_matrix), colSums(cell_cover_clonal_mutation_site)) -> clonal_mutation_detected_covered
pvalue <- pbinom(clonal_mutation_detected_covered[,1], clonal_mutation_detected_covered[,2], prob = 0.5)
padj <- p.adjust(pvalue, method = "BH")
padj_sig <- padj < 0.05
sum(padj < 0.05)
annotation_row <- data.frame(padj_sig/1)
colnames(annotation_row) <- "p.adj<0.05"
annotation_row <- annotation_row[colnames(clonal_mutation_matrix), ,drop=F]

annotation_col <- df_label_prop_cluster[, c("cell", "label_prop_cluster", "cell_type")]
annotation_col <- annotation_col[!duplicated(annotation_col$cell),]
row.names(annotation_col) <- annotation_col$cell
annotation_col <- annotation_col[cell_order, c("label_prop_cluster", "cell_type")]
annotation_col$label_prop_cluster <- factor(annotation_col$label_prop_cluster,
                                            levels = seq(1, length(unique(df_label_prop_cluster$label_prop_cluster))))
annotation_col$cell_type <- factor(annotation_col$cell_type, 
                                                      levels = c("melanoma_fibroblast", "B_lymphoblast", "filtered"),
                                                      labels = c("Melanoma", "B-cell", "Unable to determine (Low RNA quality)"))
colnames(annotation_col) <- c("Cluster", "Cell type")
annotation_colors <- list(`Cell type` = c("B-cell" = "cadetblue1", 
                                        "Melanoma" = "coral1", 
                                        "Unable to determine (Low RNA quality)" = "whitesmoke"),
                          Cluster = label_prop_color)
names(annotation_colors[["Cluster"]]) <- seq(1, length(unique(df_label_prop_cluster$label_prop_cluster)))

pheatmap(t(clonal_mutation_matrix+cell_cover_clonal_mutation_site), scale = "none", show_colnames = F,
         cluster_rows = T, cluster_cols = F, legend = F,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         #annotation_row = annotation_row,
         border_color = NA,
         cellwidth = 1, 
         cellheight = 10)
dev.off()
pheatmap(t(clonal_mutation_matrix+cell_cover_clonal_mutation_site), scale = "none", show_colnames = F,
         cluster_rows = T, cluster_cols = F, legend = F,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         #annotation_row = annotation_row,
         border_color = NA,
         cellwidth = 1, 
         cellheight = 10, 
         filename = "figures/lineage_tracing/heatmap_label_prop_cluster_BL_bulk_coverage_20240604.pdf", height = 10)

save.image("./figures/manuscript_figures/Figure2/lineage_COLO829BLT50.RData")
