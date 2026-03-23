setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/")
library(randomForest)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# font
font = "Helvetica"
font_size = 7
axis_font_size = 6


tumor <- readRDS("./data/RNA_ATAC/multiomes_integrated_20250604_GBM_annotated.rds")
seurat_cluster <- tumor[[c("cell_type")]]
seurat_cluster[,"cell.id"] <- row.names(seurat_cluster)

potential_tumor_doublet <- readRDS("./results/celltype_annotation/GBM_potential_tumor_doublet.rds")

clonal_mutation_list <- 
  read.csv("./results/clonal_mutation/clonal_mutation_list_realigned_GBM_20250604.csv")
clonal_mutation_list[(duplicated(str_c(clonal_mutation_list$cell, clonal_mutation_list$key, clonal_mutation_list$case_id))),]
clonal_mutation_list$cell <- str_c(clonal_mutation_list$case_id, clonal_mutation_list$cell, sep="_")
clonal_mutation_list[!is.na(clonal_mutation_list$cell_type), ]-> clonal_mutation_list
clonal_mutation_list <- clonal_mutation_list[!clonal_mutation_list$cell %in% potential_tumor_doublet,]

#clonal_mutation_list <- clonal_mutation_list[clonal_mutation_list$cell_type == "Tumor",]
#clonal_mutation_list <- left_join(clonal_mutation_list, seurat_cluster, by = c("cell" = "cell.id"))
#clonal_mutation_list$cell_type <- clonal_mutation_list$tumor_seurat_clusters

clonal_mutation_count <- clonal_mutation_list %>% group_by(key) %>% summarise("count" = n())
selected_clonal_mutation <- clonal_mutation_count[clonal_mutation_count$count >= 2, "key", drop = T]
clonal_mutation_list <- clonal_mutation_list[clonal_mutation_list$key %in% selected_clonal_mutation,]

write.csv(clonal_mutation_list, 
          "./results/clonal_mutation/clonal_mutation_list_realigned_GBM_ATAC_celltype_remove_doublet_20250604_new.csv", quote = F, row.names = F)

clonal_mutation_list$cell_type <- as.factor(clonal_mutation_list$cell_type)
clonal_mutation_list %>% group_by(key, cell_type, .drop=F) %>% count() -> clonal_mutation_detected_cell_count
colnames(clonal_mutation_detected_cell_count)[3] <- "detected"

all_covered_cell <- readRDS("./results/clonal_mutation/clonal_mutation_site_covered_cell_GBM_08152025.rds")
all_covered_cell <- data.frame(all_covered_cell)
colnames(all_covered_cell) <- c("dataset", "mutation", "cell.id")
all_covered_cell$cell.id <- str_c(all_covered_cell$dataset, all_covered_cell$cell.id, sep = "_")
all_covered_cell <- inner_join(all_covered_cell, seurat_cluster)
all_covered_cell <- all_covered_cell[!all_covered_cell$cell.id %in% potential_tumor_doublet,]

all_covered_cell$cell_type <- as.factor(all_covered_cell$cell_type)
all_covered_cell %>% group_by(mutation, cell_type, .drop=F) %>% count() -> covered_cell_count
colnames(covered_cell_count) <- c("key", "cell_type", "covered")

clonal_mutation_detected_covered_cell_count <- left_join(clonal_mutation_detected_cell_count, covered_cell_count, by = c("key", "cell_type"))

clonal_mutation_detected_covered_cell_count$detected/clonal_mutation_detected_covered_cell_count$covered -> clonal_mutation_detected_covered_cell_count$vaf
clonal_mutation_detected_covered_cell_count$vaf[clonal_mutation_detected_covered_cell_count$covered==0] <- 0
clonal_mutation_detected_covered_cell_count$vaf[clonal_mutation_detected_covered_cell_count$vaf > 1] <- 1
clonal_mutation_detected_covered_cell_count[, c(1,2,5)] %>%
  pivot_wider(
    names_from = cell_type,
    values_from = vaf
  ) -> clonal_mutation_vaf_subtype

# Order variants by VAF in Tumor cells
variant_order <- clonal_mutation_detected_covered_cell_count %>%
  select(key, cell_type, detected) %>%
  pivot_wider(names_from = cell_type, values_from = detected) %>%
  arrange(desc(Tumor),
          desc(Oligodendorcyte),
          desc(Microglia),
          desc(Astrocyte)
          ) %>%
  pull(key)

clonal_mutation_detected_covered_cell_count <- clonal_mutation_detected_covered_cell_count %>%
  mutate(key = factor(key, levels = variant_order))

key_plot <- str_replace(variant_order, "[-]", ":")
key_plot <- str_replace(key_plot, "[-]", " ")
key_plot <- str_replace(key_plot, "[-]", ">")

clonal_mutation_detected_covered_cell_count <- clonal_mutation_detected_covered_cell_count %>%
  mutate(key_plot = factor(key, levels = variant_order, labels = key_plot))

clonal_mutation_detected_covered_cell_count_plot <- clonal_mutation_detected_covered_cell_count %>% filter(detected > 0)
clonal_mutation_detected_covered_cell_count_plot <- clonal_mutation_detected_covered_cell_count_plot[clonal_mutation_detected_covered_cell_count_plot$cell_type %in% c("Astrocyte", "Microglia", "Oligodendorcyte", "Tumor"),]
clonal_mutation_detected_covered_cell_count_plot$cell_type <- as.character(clonal_mutation_detected_covered_cell_count_plot$cell_type)
clonal_mutation_detected_covered_cell_count_plot$cell_type[clonal_mutation_detected_covered_cell_count_plot$cell_type == "Oligodendorcyte"] <- "Oligodendrocyte"

ggplot(clonal_mutation_detected_covered_cell_count_plot, aes(x = cell_type, y = key_plot)) +
  geom_point(aes(size = detected, color = vaf)) +
  scale_size_continuous(range = c(0, 2)) +   # adjust dot size range
  scale_color_viridis_c(option = "plasma", begin = 0, end = 1) +   # nicer color scale
  theme_minimal() +
  labs(x = "Cell type", y = "Clonal sSNV", 
       size = "Number of cells\n(Alt allele detected)", 
       color = "alt allele detected/locus covered") +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, hjust = 1, vjust = 1.1),
    axis.text.y = element_text(color = "black", size = 4, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
  )
ggsave("./figures/GBM_clonal_variants/clonal_SNV_VAF_clean.pdf", width = 8.5, height = 18, units = "cm")

# 1. Identify all unique mutations in the "Tumor" group
tumor_mutations <- clonal_mutation_detected_covered_cell_count_plot %>%
  filter(cell_type == "Tumor") %>%
  pull(key) %>%
  unique()

# 2. Calculate the number of shared mutations for each non-Tumor cell type
shared_counts_df <- clonal_mutation_detected_covered_cell_count_plot %>%
  # Exclude the "Tumor" group itself from the summary table
  filter(cell_type != "Tumor") %>%
  # Group by the remaining cell types
  group_by(cell_type) %>%
  # Summarize the count of shared mutations
  summarize(
    shared_mutation_count = sum(key %in% tumor_mutations),
    .groups = 'drop'
  ) %>%
  # Order the bars by count (optional)
  arrange(desc(shared_mutation_count))

# Convert cell_type to a factor with the desired display order (descending count)
shared_counts_df$cell_type <- factor(
  shared_counts_df$cell_type,
  levels = shared_counts_df$cell_type
)

# Create the bar plot
ggplot(shared_counts_df, aes(x = cell_type, y = shared_mutation_count)) +
  # Add the bars
  geom_col(fill = "#1f78b4") +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "bottom"
  ) +
  labs(
    x = "Non-tumor cell type",
    y = paste0("Number of mutations\nshared with tumor cells \n(Total Tumor Mutations:", length(tumor_mutations), ")")
  )
ggsave("./figures/GBM_clonal_variants/clonal_SNV_shared_with_tumor.pdf", width = 3, height = 4, units = "cm")


# ------------------------------------------------------------------------------
only_oligo <- clonal_mutation_detected_covered_cell_count %>%
  group_by(key) %>%
  # keep variants where:
  # - detected in Oligodendrocyte (n_cells > 0)
  # - AND not detected in any other cell type
  filter(any(cell_type == "Oligodendorcyte" & detected > 0),
         all((cell_type == "Oligodendorcyte") | (detected == 0))
         ) %>%
  ungroup() %>%
  filter(cell_type == "Oligodendorcyte")

all_covered_cell$dataset <- as.factor(all_covered_cell$dataset)
all_covered_cell %>% group_by(mutation, cell_type, dataset, .drop=F) %>% count() -> covered_cell_count
colnames(covered_cell_count) <- c("key", "cell_type", "case_id", "covered")
clonal_mutation_list$case_id <- as.factor(clonal_mutation_list$case_id)
clonal_mutation_list %>% group_by(key, cell_type, case_id, .drop=F) %>% count() -> clonal_mutation_detected_cell_count
colnames(clonal_mutation_detected_cell_count) <- c("key", "cell_type", "case_id", "detected")

clonal_mutation_detected_cell_count <- left_join(clonal_mutation_detected_cell_count, covered_cell_count)
clonal_mutation_detected_cell_count_oligo_only <- clonal_mutation_detected_cell_count[clonal_mutation_detected_cell_count$key %in% only_oligo$key,]
clonal_mutation_detected_cell_count_oligo_only <- clonal_mutation_detected_cell_count_oligo_only[clonal_mutation_detected_cell_count_oligo_only$cell_type == "Oligodendorcyte",]

clonal_mutation_detected_cell_count_oligo_only_plot <- data.frame(clonal_mutation_detected_cell_count_oligo_only) %>%
  mutate(key = as.character(key),
         case_id = as.character(case_id)) %>%
  complete(key, case_id, fill = list(detected = 0, covered = 0)) %>%
  # Value to color tiles; only color when detected > 0
  mutate(
    rate  = ifelse(covered > 0, detected / covered, NA_real_),
    fill_value = ifelse(detected > 0, rate, NA_real_),      # NA -> blank cell
    label = sprintf("%d/%d", detected, covered)
  )

clonal_mutation_detected_cell_count_oligo_only_plot$case_id <- factor(clonal_mutation_detected_cell_count_oligo_only_plot$case_id, 
                                                                      levels = c("UMB4397_tumor", "UMB4397_normal"), 
                                                                      labels = c("More infiltrated\nby tumor cellsr", "Less infiltrated\nby tumor cells"))
clonal_mutation_detected_cell_count_oligo_only_plot$key_plot <- as.character(clonal_mutation_detected_cell_count_oligo_only_plot$key)
clonal_mutation_detected_cell_count_oligo_only_plot$key_plot <- str_replace(clonal_mutation_detected_cell_count_oligo_only_plot$key_plot, "[-]", ":")
clonal_mutation_detected_cell_count_oligo_only_plot$key_plot <- str_replace(clonal_mutation_detected_cell_count_oligo_only_plot$key_plot, "[-]", " ")
clonal_mutation_detected_cell_count_oligo_only_plot$key_plot <- str_replace(clonal_mutation_detected_cell_count_oligo_only_plot$key_plot, "[-]", ">")

ggplot(clonal_mutation_detected_cell_count_oligo_only_plot, aes(x = key_plot, y = case_id)) +
  # Heatmap tiles; NA fill shows as white
  geom_tile(aes(fill = fill_value), color = "grey90") +
  # Text context inside each cell
  geom_text(aes(label = label, color = detected > 0), size = 2) +
  # Color scales 
  scale_fill_gradient2(
    low = "#e6f0ff",
    high = "#0033cc",
    na.value = "white",
    name = "alt allele detected/locus covered",
  ) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  coord_fixed() +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "bottom"
  ) +
  labs(x = "Clonal sSNVs (Oligodendrocyte-only)", y = "Brain region")
ggsave("./figures/GBM_clonal_variants/OL_only_clonal_SNV_VAF_clean.pdf", width = 10, height = 6, units = "cm")

clonal_mutation_detected_cell_count_tumor <- clonal_mutation_detected_cell_count[clonal_mutation_detected_cell_count$cell_type == "Tumor",]

# Step 1: collapse per key which locs are present
key_status <- clonal_mutation_detected_cell_count_tumor %>%
  group_by(key, case_id) %>%
  summarise(det = sum(detected), .groups = "drop") %>%
  mutate(present = det > 0) %>%
  group_by(key) %>%
  summarise(
    has_A = any(case_id == "UMB4397_tumor" & present),
    has_B = any(case_id == "UMB4397_normal" & present),
    .groups = "drop"
  ) %>%
  mutate(category = case_when(
    has_A & has_B ~ "Both regions",
    has_A & !has_B ~ "Only region more infiltrated by tumor cells",
    !has_A & has_B ~ "Only region less infiltrated by tumor cells",
    TRUE ~ "None"
  ))

# Step 2: count how many keys fall into each category
counts <- key_status %>%
  count(category)
counts <- counts[counts$category != "None",]

# Step 3: pie chart
ggplot(counts, aes(x = "", y = n, fill = category)) +
  geom_col(color = "white") +
  coord_polar(theta = "y") +
  theme_void() +
  ggsci::scale_fill_npg() +
  labs(fill = "Tumor clonal sSNV detections", y = NULL, x = NULL) +
  geom_text(aes(label = n),
            position = position_stack(vjust = 0.5),
            color = "black", size = 2) +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.key.size = unit(0.3, "cm")
  )
ggsave("./figures/GBM_clonal_variants/Tumor_clonal_SNV_region_clean.pdf", width = 7, height = 3, units = "cm")

