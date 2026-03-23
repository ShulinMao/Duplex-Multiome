# Libraries
#.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(scales)
library(rtracklayer)

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

case_id_list <- c("UMB4428", "UMB1278_deeper", "UMB1864_deeper", "UMB4638", "UMB1465", 
                  "UMB5451", "Brain1901106_deeper", "UMB5657", "UMB5823_deeper")
COLO829BLT_case_id_list <- c("COLO829BLT50_rep1", "COLO829BLT50_rep2")
brain_COLO829BLT_case_id_list <- c(case_id_list, COLO829BLT_case_id_list)


# ------------------------------------------------------------------------------
# depth of mutation calls
duplex_snv_candidate_depth <- c()
for (sample_id in brain_COLO829BLT_case_id_list){
  for (a in c(2, 4, 6)){
    # with duplex
    s <- a %/% 2
    stringency <- str_c("a", a, "s", s)
    file_path <- paste0("./data/", sample_id, "/single_cell/call_variants/autosomal_variants/", stringency, "_autosome_variants.txt")
    if(file.info(file_path)$size==0){
      duplex_snv_candidate_count <- rbind(duplex_snv_candidate_count, c(sample_id, stringency, 0))
      next
    }
    duplex_snv_candidate_sample <- read.table(paste0("./data/", sample_id, "/single_cell/call_variants/autosomal_variants/", stringency, "_autosome_variants.txt"))
    colnames(duplex_snv_candidate_sample) <- 
      c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "VALUE", "CELL_ID")
    depth <- data.frame(str_split(duplex_snv_candidate_sample$VALUE, ":", simplify = T)[, 3:4])
    mean_depth <- apply(depth, 2, function(x){mean(as.numeric(x))})
    per_5_depth <- apply(depth, 2, function(x){quantile(as.numeric(x), c(0.05))})
    per_95_depth <- apply(depth, 2, function(x){quantile(as.numeric(x), c(0.95))})
    
    duplex_snv_candidate_depth <- rbind(duplex_snv_candidate_depth, c(sample_id, stringency, mean_depth, per_5_depth, per_95_depth))
  }
}

duplex_snv_candidate_depth <- data.frame(duplex_snv_candidate_depth)
colnames(duplex_snv_candidate_depth) <- c("case_id", "stringency", "strand1_depth", "strand2_depth", 
                                          "strand1_depth_per_5", "strand2_depth_per_5",
                                          "strand1_depth_per_95", "strand2_depth_per_95")
duplex_snv_candidate_depth$strand1_depth <- as.numeric(duplex_snv_candidate_depth$strand1_depth)
duplex_snv_candidate_depth$strand2_depth <- as.numeric(duplex_snv_candidate_depth$strand2_depth)
duplex_snv_candidate_depth$strand1_depth_per_5 <- as.numeric(duplex_snv_candidate_depth$strand1_depth_per_5)
duplex_snv_candidate_depth$strand2_depth_per_5 <- as.numeric(duplex_snv_candidate_depth$strand2_depth_per_5)
duplex_snv_candidate_depth$strand1_depth_per_95 <- as.numeric(duplex_snv_candidate_depth$strand1_depth_per_95)
duplex_snv_candidate_depth$strand2_depth_per_95 <- as.numeric(duplex_snv_candidate_depth$strand2_depth_per_95)

# To plot strand 1 and 2 in one bar
duplex_snv_candidate_depth$strand1_depth_per_5 <- duplex_snv_candidate_depth$strand1_depth_per_5 + duplex_snv_candidate_depth$strand2_depth
duplex_snv_candidate_depth$strand1_depth_per_95 <- duplex_snv_candidate_depth$strand1_depth_per_95 + duplex_snv_candidate_depth$strand2_depth

reshape2::melt(duplex_snv_candidate_depth, id.vars = c("case_id", "stringency", 
                                                       "strand1_depth_per_5", "strand2_depth_per_5",
                                                       "strand1_depth_per_95", "strand2_depth_per_95")) -> duplex_snv_candidate_depth
duplex_snv_candidate_depth$case_id <- factor(duplex_snv_candidate_depth$case_id, 
                                             levels = c(brain_COLO829BLT_case_id_list),
                                             labels = c("UMB4428", "UMB1278", "UMB1864", "UMB4638", "UMB1465", 
                                                        "UMB5451", "Brain1901106", "UMB5657", "UMB5823", 
                                                        "COLO829BLT50_rep1", "COLO829BLT50_rep2"))

ggplot(duplex_snv_candidate_depth, aes(x = case_id, y = value, color = variable)) +
  geom_col(position = "stack", fill = "white") +
  geom_errorbar(aes(ymin = strand1_depth_per_5, ymax = strand1_depth_per_95),
                width = 0.2, position = position_nudge(x = 0.05), color = "#F8766D") +
  geom_errorbar(aes(ymin = strand2_depth_per_5, ymax = strand2_depth_per_95),
                width = 0.2, position = position_nudge(x = -0.05), color = "#00BFC4") +
  facet_wrap(.~stringency) +
  labs(x = "Sample", y = "Average depth of detected sSNVs", color = "Strand") +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
  )
ggsave("./figures/QC_metrics/detected_SNV_depth.pdf", width = 14, height = 8, units = "cm")


# ------------------------------------------------------------------------------
# coverage of cells
coverage_cell <- c()
for (sample_id in brain_COLO829BLT_case_id_list){
  for (a in c(2, 4, 6)){
    # with duplex
    s <- a %/% 2
    stringency <- str_c("a", a, "s", s)
    
    summary_file <- paste0("./data/", sample_id, "/single_cell/burden_estimation/Dan/",
                           stringency, "_burden_results.txt")
    burden_result <- read.table(summary_file, sep = "\t", header = F)
    cell_id <- str_extract(burden_result$V5, pattern = "[ATCG]{16}-1")
    burden_result$V6 <- cell_id
    
    coverage <- cbind(sample_id, stringency, burden_result$V3, burden_result$V6)
    
    coverage_cell <- rbind(coverage_cell, coverage)
    
    # without duplex
    s <- 0
    stringency <- str_c("a", a, "s", s)
    
    summary_file <- paste0("./data/", sample_id, "/single_cell/burden_estimation/Dan/",
                           stringency, "_burden_results.txt")
    burden_result <- read.table(summary_file, sep = "\t", header = F)
    cell_id <- str_extract(burden_result$V5, pattern = "[ATCG]{16}-1")
    burden_result$V6 <- cell_id
    
    coverage <- cbind(sample_id, stringency, burden_result$V3, burden_result$V6)
    
    coverage_cell <- rbind(coverage_cell, coverage)
  }
}

coverage_cell <- data.frame(coverage_cell)
colnames(coverage_cell) <- c("case_id", "stringency", "coverage", "cell_id")
coverage_cell$coverage <- as.integer(coverage_cell$coverage)
coverage_cell$case_id <- factor(coverage_cell$case_id, 
                                             levels = c(brain_COLO829BLT_case_id_list),
                                             labels = c("UMB4428", "UMB1278", "UMB1864", "UMB4638", "UMB1465", 
                                                        "UMB5451", "Brain1901106", "UMB5657", "UMB5823", 
                                                        "COLO829BLT50_rep1", "COLO829BLT50_rep2"))
# coverage_cell$stringency <- factor(coverage_cell$stringency, 
#                                    levels = c("a2s0", "a2s1"),
#                                    labels = c("w/o duplex", "w/ duplex"))

ggplot(coverage_cell, aes(x = case_id, y = coverage, fill = stringency)) +
  geom_boxplot() +
  facet_wrap(.~stringency, scales = "free_y") +
  labs(x = "Sample", y = "Callable regions per cell (bp)", fill = "Strand") +
  theme_classic() +
  scale_y_continuous(labels = scientific_10) +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
  )


coverage_cell %>% group_by(case_id, stringency) %>% 
  summarise("median" = median(coverage),
            "per_5" = quantile(coverage, 0.05),
            "per_95" = quantile(coverage, 0.95)) -> coverage_cell_summary
ggplot(coverage_cell_summary, aes(x = case_id, color = stringency)) +
  geom_pointrange(aes(y = median, ymin = per_5, ymax = per_95)) +
  facet_wrap(.~stringency, scales = "free_y") +
  labs(x = "Sample", y = "Callable regions per cell (bp)", fill = "Strand") +
  theme_classic() +
  scale_y_continuous(labels = scientific_10) +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
  )
ggsave("./figures/QC_metrics/callable_region_size.pdf", width = 14, height = 10, units = "cm")

genome_size <- 5.75e9
ggplot(coverage_cell_summary, aes(x = case_id, color = stringency)) +
  geom_pointrange(aes(y = median/genome_size, ymin = per_5/genome_size, ymax = per_95/genome_size)) +
  facet_wrap(.~stringency, scales = "free_y") +
  labs(x = "Sample", y = "Fraction of callable regions per cell", fill = "Strand") +
  theme_classic() +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
  )
ggsave("./figures/QC_metrics/callable_region_fraction.pdf", width = 14, height = 10, units = "cm")

# By cell type
brain_cell_type <- read.table("./results/celltype_annotation/annotation_table_20240512.txt", sep = "\t")
colnames(brain_cell_type) <- c("case_id", "cell_id", "cell_type")
brain_cell_type[!str_detect(brain_cell_type$cell_type, "doublet"),] -> brain_cell_type
brain_cell_type$cell_type[str_detect(brain_cell_type$cell_type, "EN")] <- "Excitatory Neuron"
brain_cell_type$cell_type[str_detect(brain_cell_type$cell_type, "IN")] <- "Inhibitory Neuron"
brain_cell_type$case_id <- factor(brain_cell_type$case_id, 
                                  levels = case_id_list,
                                  labels = c("UMB4428", "UMB1278", "UMB1864", "UMB4638", "UMB1465", 
                                             "UMB5451", "Brain1901106", "UMB5657", "UMB5823"))

coverage_cell_brain <- coverage_cell[coverage_cell$case_id %in% c("UMB4428", "UMB1278", "UMB1864", "UMB4638", "UMB1465", 
                                                                  "UMB5451", "Brain1901106", "UMB5657", "UMB5823"),]
coverage_cell_brain_cell_type <- inner_join(coverage_cell_brain, brain_cell_type)

colors_cell_type <- c("#2E2585", "#337538", "#7CAE00", "#E69F00", "#DCCD7D", "#00A9FF", "#CC79A7", "#7E2954")  
ggplot(coverage_cell_brain_cell_type, aes(x = cell_type, y = coverage, fill = cell_type)) +
  geom_violin() +
  stat_summary(fun = "median", 
               geom="point", color="black", size = 0.5,
               position = position_dodge(width = 0.8)) +
  facet_wrap(.~stringency, scales = "free_y") +
  labs(x = "Cell type", y = "Callable regions per cell", fill = "Cell type") +
  theme_classic() +
  scale_y_continuous(labels = scientific_10) +
  scale_fill_manual(values = colors_cell_type) +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
  )
ggsave("./figures/QC_metrics/callable_region_size_brain_celltype.pdf", width = 15, height = 10, units = "cm")


# COLO829BLT50
ggplot(coverage_cell_summary[coverage_cell_summary$case_id %in% COLO829BLT_case_id_list,], aes(x = case_id, color = stringency)) +
  geom_pointrange(aes(y = median, ymin = per_5, ymax = per_95)) +
  facet_wrap(.~stringency, scales = "free_y") +
  labs(x = "Sample", y = "Callable regions per cell (bp)", fill = "Strand") +
  theme_classic() +
  scale_y_continuous(labels = scientific_10) +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
  )
ggsave("./figures/QC_metrics/callable_region_size_COLO829BLT50.pdf", width = 8, height = 10, units = "cm")


coverage_cell_COLO829BLT50 <- coverage_cell[coverage_cell$case_id %in% COLO829BLT_case_id_list,]
COLO829BLT50_cell_type_annotation <- read.table("./results/celltype_annotation/annotation_COLO829BLT50_20240523.txt") 
colnames(COLO829BLT50_cell_type_annotation) <-  c("case_id", "cell_id", "cell_type")
coverage_cell_COLO829BLT50_cell_type <- inner_join(coverage_cell_COLO829BLT50, COLO829BLT50_cell_type_annotation)

ggplot(coverage_cell_COLO829BLT50_cell_type, aes(x = case_id, y = coverage, fill = cell_type)) +
  geom_violin(position = position_dodge(width = 0.8)) +
  stat_summary(fun = "median", 
               geom="point", color="black", size = 0.5,
               position = position_dodge(width = 0.8)) +
  facet_wrap(.~stringency, scales = "free_y") +
  labs(x = "Sample", y = "Callable regions per cell", fill = "Cell type") +
  theme_classic() +
  scale_y_continuous(labels = scientific_10) +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    #legend.position = "none",
  )
ggsave("./figures/QC_metrics/callable_region_size_COLO829BLT50_celltype.pdf", width = 14, height = 10, units = "cm")

# ------------------------------------------------------------------------------
# sensitivity
sensitivity_cell <- c()
for (sample_id in brain_COLO829BLT_case_id_list){
  for (a in c(6)){
    # with duplex
    s <- a %/% 2
    stringency <- str_c("a", a, "s", s)
    
    summary_file <- paste0("./data/", sample_id, "/single_cell/call_variants/germline_mutations/",
                           stringency, "_germline_mutation_results.txt")
    called_germline_mutation <- read.table(summary_file, sep = "\t", header = F)
    cell_id <- str_split(called_germline_mutation$V2, pattern = "[/]", simplify = T)
    cell_id <- cell_id[,dim(cell_id)[2]]
    cell_id <- str_split(cell_id, pattern = "[_]", simplify = T)[,1]
    rownames(called_germline_mutation) <- cell_id
    
    summary_file <- paste0("./data/", sample_id, "/single_cell/call_variants/all_germline_mutations_callable_region/",
                           "all_germline_mutations_callable_region_summary_", stringency, ".txt")
    all_germline_mutation <- read.table(summary_file, sep = "\t", header = F)
    cell_id <- str_split(all_germline_mutation$V2, pattern = "[/]", simplify = T)
    cell_id <- cell_id[,dim(cell_id)[2]]
    cell_id <- str_split(cell_id, pattern = "[_]", simplify = T)[,1]
    rownames(all_germline_mutation) <- cell_id
    
    sensitivity <- cbind(sample_id, stringency, called_germline_mutation[rownames(called_germline_mutation), "V1"]/all_germline_mutation[rownames(called_germline_mutation), "V1"])
    sensitivity <- na.omit(sensitivity)
    
    sensitivity_cell <- rbind(sensitivity_cell, sensitivity)
  }
}

sensitivity_cell <- data.frame(sensitivity_cell)
colnames(sensitivity_cell) <- c("case_id", "stringency", "sensitivity")
sensitivity_cell$sensitivity <- as.numeric(sensitivity_cell$sensitivity)
sensitivity_cell[is.finite(sensitivity_cell$sensitivity), ] -> sensitivity_cell
sensitivity_cell[sensitivity_cell$sensitivity <= 1,] -> sensitivity_cell
sensitivity_cell$case_id <- factor(sensitivity_cell$case_id, 
                                levels = c(brain_COLO829BLT_case_id_list),
                                labels = c("UMB4428", "UMB1278", "UMB1864", "UMB4638", "UMB1465", 
                                           "UMB5451", "Brain1901106", "UMB5657", "UMB5823", 
                                           "COLO829BLT50_rep1", "COLO829BLT50_rep2"))
sensitivity_cell %>% group_by(case_id, stringency) %>% 
  summarise("median" = median(sensitivity, na.rm = T),
            "per_5" = quantile(sensitivity, 0.05, na.rm = T),
            "per_95" = quantile(sensitivity, 0.95, na.rm = T)) -> sensitivity_cell_summary

ggplot(sensitivity_cell_summary[sensitivity_cell_summary$case_id %in% COLO829BLT_case_id_list,], aes(x = case_id, color = stringency)) +
  geom_pointrange(aes(y = median, ymin = per_5, ymax = per_95)) +
  labs(x = "Sample", 
       y = "Sensitivity\n(detected germline mutations/all germline mutations\nin callable regions)", 
       title = "a6s3") +  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
  )
ggsave("./figures/QC_metrics/germline_sensitivity_COLO829BLT50.pdf", width = 3, height = 8, units = "cm")


ggplot(sensitivity_cell_summary, aes(x = case_id, color = stringency)) +
  geom_pointrange(aes(y = median, ymin = per_5, ymax = per_95)) +
  labs(x = "Sample", 
       y = "Sensitivity\n(detected germline mutations/all germline mutations\nin callable regions)", 
       title = "a6s3") +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
  )


# ------------------------------------------------------------------------------
# actual number of mutation calls (by cell type)
duplex_snv_candidate_count <- c()
for (sample_id in case_id_list){
  for (a in c(6)){
    # with duplex
    s <- a %/% 2
    stringency <- str_c("a", a, "s", s)
    file_path <- paste0("./data/", sample_id, "/single_cell/call_variants/autosomal_variants/", stringency, "_autosome_variants.txt")
    if(file.info(file_path)$size==0){
      duplex_snv_candidate_sample_count <- matrix(c(sample_id, stringency, "Undetermined", 0), nrow = 1)
      colnames(duplex_snv_candidate_sample_count) <- c("sample_id", "stringency", "cell_type", "count")
      duplex_snv_candidate_count <- rbind(duplex_snv_candidate_count, duplex_snv_candidate_sample_count)
      next
    }
    duplex_snv_candidate_sample <- read.table(paste0("./data/", sample_id, "/single_cell/call_variants/autosomal_variants/", stringency, "_autosome_variants.txt"))
    colnames(duplex_snv_candidate_sample) <- 
      c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "VALUE", "CELL_ID")

    cell_type_info_sample <- brain_cell_type[brain_cell_type$case_id == str_remove(sample_id, "_deeper"),]
    duplex_snv_candidate_sample <- left_join(duplex_snv_candidate_sample, cell_type_info_sample,
              by = c("CELL_ID" = "cell_id"))
    duplex_snv_candidate_sample$cell_type[is.na(duplex_snv_candidate_sample$cell_type)] <- "Undetermined"
    duplex_snv_candidate_sample %>% group_by(cell_type) %>% summarise("count" = n()) -> duplex_snv_candidate_sample_count
    duplex_snv_candidate_sample_count <- cbind(sample_id, stringency, duplex_snv_candidate_sample_count)
    duplex_snv_candidate_count <- rbind(duplex_snv_candidate_count, duplex_snv_candidate_sample_count)
  }
}
duplex_snv_candidate_count <- data.frame(duplex_snv_candidate_count)
colnames(duplex_snv_candidate_count) <- c("case_id", "stringency", "cell_type", "snv")
duplex_snv_candidate_count$snv <- as.integer(duplex_snv_candidate_count$snv)
duplex_snv_candidate_count$case_id <- factor(duplex_snv_candidate_count$case_id, 
                                             levels = c(case_id_list),
                                             labels = c("UMB4428", "UMB1278", "UMB1864", "UMB4638", "UMB1465", 
                                                        "UMB5451", "Brain1901106", "UMB5657", "UMB5823"))

ggplot(duplex_snv_candidate_count, aes(x = case_id, y = snv)) +
  geom_col(position = "dodge") +
  geom_text(
    aes(label = snv),
    position = position_dodge(width = 0.9),
    vjust = -0.1, size = 2) +
  facet_wrap(.~cell_type) +
  theme_classic() +
  labs(x = "Sample", y = "Detected sSNV number", fill = "Stringency") +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
  )
ggsave("./figures/QC_metrics/absolute_SNV_count_per_sample.pdf", width = 14, height = 16, units = "cm")


# ------------------------------------------------------------------------------
# clonal mutation with duplex
# COLO829BLT50
# BLT50 bulk
clonal_mutation_list <- read.csv("./results/clonal_mutation/clonal_mutation_list_COLO829BLT50_BLT50_bulk_realigned_20241007.csv")

duplex_clonal_mutation <- c()
for (mut in unique(clonal_mutation_list$key)){
  clonal_mutation <- clonal_mutation_list[clonal_mutation_list$key == mut,]
  str_split(clonal_mutation$value, pattern = ":", simplify = T) -> strand_count
  duplex_clonal_mutation[mut] <- sum(as.integer(strand_count[,3]) > 0 & as.integer(strand_count[,4]) > 0)
}

# plot pie chart
df <- data.frame(
  group = c("w/ duplex support", "w/o duplex support"),
  value = c(sum(duplex_clonal_mutation > 0), length(duplex_clonal_mutation) - sum(duplex_clonal_mutation > 0))
)
df <- df %>%
  mutate(
    percent = value / sum(value) * 100,
    label = paste0(value, " (", round(percent, 1), "%)")
  )

# Create pie chart
p1 <- ggplot(df, aes(x = "", y = value, fill = group)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5)) +
  labs(fill = NULL, title = "COLO829BLT50\n(BLT50 bulk WGS)") +
  theme_void() +
  theme(legend.position = "right")

# BL bulk
clonal_mutation_list <- read.csv("./results/clonal_mutation/clonal_mutation_list_COLO829BLT50_BL_bulk_realigned_20241007.csv")

duplex_clonal_mutation <- c()
for (mut in unique(clonal_mutation_list$key)){
  clonal_mutation <- clonal_mutation_list[clonal_mutation_list$key == mut,]
  str_split(clonal_mutation$value, pattern = ":", simplify = T) -> strand_count
  duplex_clonal_mutation[mut] <- sum(as.integer(strand_count[,3]) > 0 & as.integer(strand_count[,4]) > 0)
}

# plot pie chart
df <- data.frame(
  group = c("w/ duplex support", "w/o duplex support"),
  value = c(sum(duplex_clonal_mutation > 0), length(duplex_clonal_mutation) - sum(duplex_clonal_mutation > 0))
)
df <- df %>%
  mutate(
    percent = value / sum(value) * 100,
    label = paste0(value, " (", round(percent, 1), "%)")
  )

# Create pie chart
p2 <- ggplot(df, aes(x = "", y = value, fill = group)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5)) +
  labs(fill = NULL, title = "COLO829BLT50\n(BL bulk WGS)") +
  theme_void() +
  theme(legend.position = "right")

ggarrange(p1, p2, ncol = 2, common.legend = T, legend = "right")
ggsave("./figures/QC_metrics/clonal_SNV_duplex_support_COLO829BLT50.pdf", width = 12, height = 4, units = "cm")

# Brain
clonal_mutation_list <- read.csv("./results/clonal_mutation/clonal_mutation_list_healthy_brain_realigned_20241028.csv")

duplex_clonal_mutation <- c()
for (mut in unique(clonal_mutation_list$key)){
  clonal_mutation <- clonal_mutation_list[clonal_mutation_list$key == mut,]
  str_split(clonal_mutation$value, pattern = ":", simplify = T) -> strand_count
  duplex_clonal_mutation[mut] <- sum(as.integer(strand_count[,3]) > 0 & as.integer(strand_count[,4]) > 0)
}

# plot pie chart
df <- data.frame(
  group = c("w/ duplex support", "w/o duplex support"),
  value = c(sum(duplex_clonal_mutation > 0), length(duplex_clonal_mutation) - sum(duplex_clonal_mutation > 0))
)
df <- df %>%
  mutate(
    percent = value / sum(value) * 100,
    label = paste0(value, " (", round(percent, 1), "%)")
  )

# Create pie chart
ggplot(df, aes(x = "", y = value, fill = group)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5)) +
  labs(fill = NULL, title = "Brain") +
  theme_void() +
  theme(legend.position = "right")

ggsave("./figures/QC_metrics/clonal_SNV_duplex_support_brain.pdf", width = 8, height = 4, units = "cm")


# ------------------------------------------------------------------------------
# Proportion of real duplex
duplex_rate <- c()
for (case_id in brain_COLO829BLT_case_id_list){
  dir_path <- paste0("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/", case_id, "/single_cell/read_families")
  dir(dir_path, pattern = "*families.bed$", full.names = T) -> bed_file_list

  for (bed_file in bed_file_list){
    read_family <- read.table(bed_file, sep = "\t")
    read_family <- read_family[!duplicated(read_family$V4),]
    duplex_rf <- read_family[read_family$V7 > 0 & read_family$V8 > 0,]
    duplex_rf_count <- sum(duplex_rf$V7 + duplex_rf$V8)%/%2
    rf_count <- sum(read_family$V7 + read_family$V8)%/%2
    if(rf_count != 0){
      duplex_rf_proportion <- duplex_rf_count/rf_count
    }else{
      duplex_rf_proportion <- NA
    }
    duplex_rate <- rbind(duplex_rate, c(case_id, duplex_rf_count, rf_count, duplex_rf_proportion))
  }
}
#saveRDS(duplex_rate, "./results/QC_metrics/duplex_rf_rate_all_reads.rds")
duplex_rate <- readRDS("./results/QC_metrics/duplex_rf_rate_all_reads.rds")

duplex_rate <- data.frame(duplex_rate)
colnames(duplex_rate) <- c("case_id", "duplex_rf_count", "rf_count", "duplex_rf_proportion")
duplex_rate$duplex_rf_count <- as.integer(duplex_rate$duplex_rf_count)
duplex_rate$rf_count <- as.integer(duplex_rate$rf_count)
duplex_rate$duplex_rf_proportion <- as.numeric(duplex_rate$duplex_rf_proportion)
duplex_rate$case_id <- factor(duplex_rate$case_id, 
                                 levels = c(brain_COLO829BLT_case_id_list),
                                 labels = c("UMB4428", "UMB1278", "UMB1864", "UMB4638", "UMB1465", 
                                            "UMB5451", "Brain1901106", "UMB5657", "UMB5823", 
                                            "COLO829BLT50_rep1", "COLO829BLT50_rep2"))

ggplot(duplex_rate, aes(x = case_id, y = duplex_rf_proportion)) +
  geom_violin() +
  labs(x = "Sample", y = "% of duplex reads per cell") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
  )
ggsave("./figures/QC_metrics/proportion_duplex_read.pdf", width = 8, height = 8, units = "cm")

ggplot(duplex_rate[duplex_rate$case_id %in% COLO829BLT_case_id_list,], aes(x = case_id, y = duplex_rf_proportion)) +
  geom_violin() +
  labs(x = "Sample", y = "% of duplex reads per cell") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
  )

ggsave("./figures/QC_metrics/proportion_duplex_read_COLO829BLT50.pdf", width = 3, height = 8, units = "cm")


# ------------------------------------------------------------------------------
# Downsampled PTA spectrum v.s. original spectrum
# PTA
spectrum_PTA_OL_EN <- readRDS("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/figures/manuscript_figures/Figure4/spectrum_PTA_OL_EN.rds")
# Duples Multiome
spectrum_duplex <- readRDS("./results/spectrum/corrected_mutation_spectrum_a6s3_20240627.rds")

mut_mat_EN <- cbind(spectrum_PTA_OL_EN[, "EN_PTA"], spectrum_duplex[["EN"]])
mut_mat_OL <- cbind(spectrum_PTA_OL_EN[, "OL_PTA"], spectrum_duplex[["Oligodendrocyte"]])

colnames(mut_mat_EN) <- c("PTA", "Duplex-Multiome")
colnames(mut_mat_OL) <- c("PTA", "Duplex-Multiome")

# calculate cosine similarity
cos_sim_OL_downsampled <- c()
cos_sim_EN_downsampled <- c()
cos_sim_EN_vs_OL_downsampled <- c()
for (i in 1:1000){
  # OL
  snv_type_vector_OL <- rep(names(mut_mat_OL[,"PTA"]), mut_mat_OL[,"PTA"])
  snv_type_vector_OL <- factor(snv_type_vector_OL, levels = row.names(mut_mat_OL))
  snv_type_vector_OL_downsampled <- sample(snv_type_vector_OL, size = colSums(mut_mat_OL)["Duplex-Multiome"], replace = F)
  cos_sim_OL_downsampled <- c(cos_sim_OL_downsampled, cos_sim(table(snv_type_vector_OL_downsampled), mut_mat_OL[, "PTA"]))
  
  # EN
  snv_type_vector_EN <- rep(names(mut_mat_EN[,"PTA"]), mut_mat_EN[,"PTA"])
  snv_type_vector_EN <- factor(snv_type_vector_EN, levels = row.names(mut_mat_EN))
  snv_type_vector_EN_downsampled <- sample(snv_type_vector_EN, size = round(colSums(mut_mat_EN)["Duplex-Multiome"]), replace = F)
  cos_sim_EN_downsampled <- c(cos_sim_EN_downsampled, cos_sim(table(snv_type_vector_EN_downsampled), mut_mat_EN[, "PTA"]))
  
  # EN vs OL
  cos_sim_EN_vs_OL_downsampled <- c(cos_sim_EN_vs_OL_downsampled, cos_sim(table(snv_type_vector_EN_downsampled), table(snv_type_vector_OL_downsampled)))
}

cos_sim_OL_downsampled <- cbind(cos_sim_OL_downsampled, "Oligodendrocyte")
cos_sim_EN_downsampled <- cbind(cos_sim_EN_downsampled, "Excitatory Neuron")
cos_sim_EN_vs_OL_downsampled <- cbind(cos_sim_EN_vs_OL_downsampled, "EN v.s. OL")
cos_sim_downsampled <- rbind(cos_sim_EN_downsampled, cos_sim_OL_downsampled, cos_sim_EN_vs_OL_downsampled)
colnames(cos_sim_downsampled) <- c("cos_sim", "cell_type")
cos_sim_downsampled <- data.frame(cos_sim_downsampled)
cos_sim_downsampled$cos_sim <- as.numeric(cos_sim_downsampled$cos_sim)

p1 <- ggplot(cos_sim_downsampled[cos_sim_downsampled$cell_type != "EN v.s. OL",], 
       aes(x = cell_type, y = cos_sim, fill = cell_type)) + 
  geom_boxplot() +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("#7CAE00", "#00A9FF")) +
  labs(x = "Cell type", y = "Cosine Similarity", title = "Original PTA spectrum v.s.\nDownsampled PTA spectrum\n(1000 times downsampling)")

cos_sim(mut_mat_OL[, "PTA"], mut_mat_EN[, "PTA"])
p2 <-ggplot(cos_sim_downsampled[cos_sim_downsampled$cell_type == "EN v.s. OL",], 
       aes(x = cell_type, y = cos_sim)) + 
  geom_boxplot(width = 0.4) +
  geom_hline(yintercept = cos_sim(mut_mat_OL[, "PTA"], mut_mat_EN[, "PTA"]), 
             linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = cos_sim(mut_mat_OL[, "PTA"], mut_mat_EN[, "PTA"])*1.02, 
           label = "cos sim. between the original spectra", color = "red", size = 2) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none"
  ) +
  labs(x = NULL, y = "Cosine Similarity", title = "Downsampled PTA spectra\nof OL and EN\n(1000 times downsampling)")
ggarrange(p1, p2, ncol = 2)
ggsave("./figures/QC_metrics/PTA_downsampled_spectrum_cos_sim.pdf", width = 10, height = 6, units = "cm")

# ------------------------------------------------------------------------------
# absolute number of clonal calls in each cell
# clonal mutation
brain_clonal_mutations <- read.csv("./results/clonal_mutation/clonal_mutation_list_healthy_brain_realigned_20241028.csv")
COLO829BLT50_clonal_mutations <- read.csv("./results/clonal_mutation/clonal_mutation_list_COLO829BLT50_BL_bulk_realigned_20241007.csv")
clonal_mutations <- rbind(COLO829BLT50_clonal_mutations, brain_clonal_mutations[,colnames(COLO829BLT50_clonal_mutations)])

clonal_mutation_count_cell <- c()
for (sample_id in brain_COLO829BLT_case_id_list){
  stringency <- "a2s0"
  
  summary_file <- paste0("./data/", sample_id, "/single_cell/burden_estimation/Dan/",
                         stringency, "_burden_results.txt")
  burden_result <- read.table(summary_file, sep = "\t", header = F)
  cell_id <- str_extract(burden_result$V5, pattern = "[ATCG]{16}-1")
  clonal_mutations_sample <- clonal_mutations[clonal_mutations$case_id == sample_id,]
  clonal_mutations_sample$cell <- factor(clonal_mutations_sample$cell, levels = cell_id)
  
  clonal_mutation_count <- cbind(sample_id, stringency, table(clonal_mutations_sample$cell))
  
  clonal_mutation_count_cell <- rbind(clonal_mutation_count_cell, clonal_mutation_count)
}

clonal_mutation_count_cell <- data.frame(clonal_mutation_count_cell)
colnames(clonal_mutation_count_cell) <- c("case_id", "stringency", "clonal_mutation")
clonal_mutation_count_cell$clonal_mutation <- as.integer(clonal_mutation_count_cell$clonal_mutation)
clonal_mutation_count_cell$case_id <- factor(clonal_mutation_count_cell$case_id, 
                                 levels = c(brain_COLO829BLT_case_id_list),
                                 labels = c("UMB4428", "UMB1278", "UMB1864", "UMB4638", "UMB1465", 
                                            "UMB5451", "Brain1901106", "UMB5657", "UMB5823", 
                                            "COLO829BLT50_rep1", "COLO829BLT50_rep2"))

ggplot(clonal_mutation_count_cell, aes(x = clonal_mutation)) +
  geom_histogram(binwidth = 1, boundary = 0.5, fill = "skyblue", color = "black") +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    vjust = -0.1, size = 2
  ) +
  facet_wrap(.~case_id) +
  scale_x_continuous(breaks = 0:10) +
  labs(x = "Absolute clonal sSNV number (in a single cell)", y = "Cell number") +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
  )
ggsave("./figures/QC_metrics/clonal_SNV_count_per_cell.pdf", width = 16, height = 16, units = "cm")

ggplot(clonal_mutation_count_cell[clonal_mutation_count_cell$clonal_mutation > 0,], aes(x = clonal_mutation)) +
  geom_histogram(binwidth = 1, boundary = 0.5, fill = "skyblue", color = "black") +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    vjust = -0.1, size = 2
  ) +
  facet_wrap(.~case_id) +
  scale_x_continuous(breaks = 0:10) +
  labs(x = "Absolute clonal sSNV number (in a single cell)", y = "Cell number") +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
  )

# By sample
clonal_mutation_count_sample <- clonal_mutations %>% group_by(key, case_id) %>% summarise() %>% group_by(case_id) %>% summarise("clonal_mutation" = n()) 
clonal_mutation_count_sample$case_id <- factor(clonal_mutation_count_sample$case_id, 
                                               levels = c(brain_COLO829BLT_case_id_list),
                                               labels = c("UMB4428", "UMB1278", "UMB1864", "UMB4638", "UMB1465", 
                                                          "UMB5451", "Brain1901106", "UMB5657", "UMB5823", 
                                                          "COLO829BLT50_rep1", "COLO829BLT50_rep2"))
ggplot(clonal_mutation_count_sample, aes(x = case_id, y = clonal_mutation)) +
  geom_col() +
  geom_text(
    aes(label = clonal_mutation),
    vjust = -0.1, size = 2
  ) +
  labs(y = "Absolute clonal sSNV number", x = "Sample") +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
  )
ggsave("./figures/QC_metrics/clonal_SNV_count_per_sample.pdf", width = 6, height = 6, units = "cm")


# ------------------------------------------------------------------------------
# coverage of cells in accessible regions
ENCODE_regulatory_elements <- read.table("../ref/ENCODE_hg38_regulatory_elements.bed", sep = "\t")
ENCODE_promoter <- ENCODE_regulatory_elements[ENCODE_regulatory_elements$V13 == "prom",]
ENCODE_promoter <- GRanges(seqnames = ENCODE_promoter$V1,
                           ranges = IRanges(start = ENCODE_promoter$V2,
                                            end = ENCODE_promoter$V3))
ENCODE_promoter_size <- sum(width(ENCODE_promoter))

ENCODE_regulatory_elements <- GRanges(seqnames = ENCODE_regulatory_elements$V1,
                           ranges = IRanges(start = ENCODE_regulatory_elements$V2,
                                            end = ENCODE_regulatory_elements$V3))
ENCODE_regulatory_element_size <- sum(width(ENCODE_regulatory_elements))

coverage_cell_accssible_regions <- c()
for (sample_id in brain_COLO829BLT_case_id_list){
  print(sample_id)
  gc()
  read_family_sample <- read.table(paste0("./results/QC_metrics/combined_covered_regions/", sample_id, "_open_region_a2s0_sorted.bed"), sep = "\t")
  gc()
  read_family_sample <- GRanges(seqnames = read_family_sample$V1, 
                         ranges = IRanges(start = read_family_sample$V2,
                                          end = read_family_sample$V3),
                         strand1 = read_family_sample$V7, 
                         strand2 = read_family_sample$V8)
  gc()
  
  for (a in c(2, 4, 6)){
    # with duplex
    s <- a %/% 2
    stringency <- str_c("a", a, "s", s)
    
    read_family_sample_stringency <- read_family_sample[read_family_sample$strand1 >= s*2 & read_family_sample$strand2 >= s*2, ]
    read_family_sample_stringency <- reduce(read_family_sample_stringency)
    
    #intersection <- findOverlaps(ENCODE_promoter, read_family_sample_stringency)
    intersected_ranges <- intersect(ENCODE_promoter, read_family_sample_stringency)
    #ENCODE_promoter[queryHits(intersection)]
    covered_promoter_size <- sum(width(reduce(intersected_ranges)))
    
    #intersection <- findOverlaps(ENCODE_regulatory_elements, read_family_sample_stringency)
    intersected_ranges <- intersect(ENCODE_regulatory_elements, read_family_sample_stringency)
    #ENCODE_regulatory_elements[queryHits(intersection)]
    covered_regulatory_element_size <- sum(width(reduce(intersected_ranges)))
    
    coverage_cell_accssible_regions <- rbind(coverage_cell_accssible_regions, c(sample_id, stringency, covered_promoter_size, covered_regulatory_element_size))
    
    # without duplex
    s <- 0
    stringency <- str_c("a", a, "s", s)
    
    read_family_sample_stringency <- read_family_sample[read_family_sample$strand1 + read_family_sample$strand2 >= a*2, ]
    read_family_sample_stringency <- reduce(read_family_sample_stringency)
    
    #intersection <- findOverlaps(ENCODE_promoter, read_family_sample_stringency)
    intersected_ranges <- intersect(ENCODE_promoter, read_family_sample_stringency)
    #ENCODE_promoter[queryHits(intersection)]
    covered_promoter_size <- sum(width(reduce(intersected_ranges)))
    
    #intersection <- findOverlaps(ENCODE_regulatory_elements, read_family_sample_stringency)
    intersected_ranges <- intersect(ENCODE_regulatory_elements, read_family_sample_stringency)
    #ENCODE_regulatory_elements[queryHits(intersection)]
    covered_regulatory_element_size <- sum(width(reduce(intersected_ranges)))
    
    coverage_cell_accssible_regions <- rbind(coverage_cell_accssible_regions, c(sample_id, stringency, covered_promoter_size, covered_regulatory_element_size))
  }
}  
saveRDS(coverage_cell_accssible_regions, "./results/QC_metrics/coverage_cell_accssible_regions.rds")

coverage_cell_accssible_regions <- readRDS("./results/QC_metrics/coverage_cell_accssible_regions.rds")
coverage_cell_accssible_regions <- data.frame(coverage_cell_accssible_regions)
colnames(coverage_cell_accssible_regions) <- c("case_id", "stringency", "covered_promoter_size", "covered_regulatory_element_size")
coverage_cell_accssible_regions$covered_promoter_size <- as.integer(coverage_cell_accssible_regions$covered_promoter_size)
coverage_cell_accssible_regions$covered_regulatory_element_size <- as.integer(coverage_cell_accssible_regions$covered_regulatory_element_size)
coverage_cell_accssible_regions[,"covered_promoter_proportion"] <- coverage_cell_accssible_regions$covered_promoter_size/ENCODE_promoter_size
coverage_cell_accssible_regions[,"covered_regulatory_element_proportion"] <- coverage_cell_accssible_regions$covered_regulatory_element_size/ENCODE_regulatory_element_size

coverage_cell_accssible_regions$case_id <- factor(coverage_cell_accssible_regions$case_id, 
                                levels = c(brain_COLO829BLT_case_id_list),
                                labels = c("UMB4428", "UMB1278", "UMB1864", "UMB4638", "UMB1465", 
                                           "UMB5451", "Brain1901106", "UMB5657", "UMB5823", 
                                           "COLO829BLT50_rep1", "COLO829BLT50_rep2"))

ggplot(coverage_cell_accssible_regions, aes(x = case_id, y = covered_promoter_proportion, fill = stringency)) +
  geom_col() +
  facet_wrap(.~stringency, scales = "free_y") +
  labs(x = "Sample", y = "Coverage of promoters (% of bp)") +
  theme_classic() +
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75, 1), limits = c(0, 1),
                     labels = scales::percent(c(0.25, 0.5, 0.75, 1))) +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
  )
ggsave("./figures/QC_metrics/coverage_promoter.pdf", width = 15, height = 10, units = "cm")

ggplot(coverage_cell_accssible_regions, aes(x = case_id, y = covered_regulatory_element_proportion, fill = stringency)) +
  geom_col() +
  facet_wrap(.~stringency, scales = "free_y") +
  labs(x = "Sample", y = "Coverage of regulatory elements (% of bp)") +
  theme_classic() +
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75, 1), limits = c(0, 1),
                     labels = scales::percent(c(0.25, 0.5, 0.75, 1))) +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
  )
ggsave("./figures/QC_metrics/coverage_regulatory_elements.pdf", width = 15, height = 10, units = "cm")

# ------------------------------------------------------------------------------
# coverage overlap across cell types
cell_type_intersection <- c()
for (case_id in case_id_list){
  print(case_id)
  
  cell_type <- c("EN", "IN", "Astrocyte", "Microglia", "Oligodendrocyte", "OPC")
  
  for (cell_type_1 in c("EN", "IN", "Astrocyte", "Microglia", "Oligodendrocyte", "OPC")){
    print(cell_type_1)
    #file_name <- paste0("./results/clonal_mutation/multicell_covered_region_celltype_ATAC_recovered/", case_id, "_", cell_type_1, "_open_region_a2s0_coverage_count.bed")
    file_name <- paste0("./results/QC_metrics/combined_covered_regions_a6s3/", case_id, "_", cell_type_1, "_open_region_a6s3_sorted.bed")
    if(!file.exists(file_name)){next}
    genomecov_cell_type_1 <- read.table(file_name, sep = "\t")
    genomecov_cell_type_1 <- GRanges(seqnames = genomecov_cell_type_1$V1,
                                     ranges = IRanges(start = genomecov_cell_type_1$V2,
                                                      end = genomecov_cell_type_1$V3))
    genomecov_cell_type_1 <- reduce(genomecov_cell_type_1)
    
    other_cell_type <- cell_type[cell_type != cell_type_1]
    for (cell_type_2 in other_cell_type){
      print(cell_type_2)
      #file_name <- paste0("./results/clonal_mutation/multicell_covered_region_celltype_ATAC_recovered/", case_id, "_", cell_type_2, "_open_region_a2s0_coverage_count.bed")
      file_name <- paste0("./results/QC_metrics/combined_covered_regions_a6s3/", case_id, "_", cell_type_2, "_open_region_a6s3_sorted.bed")
      if(!file.exists(file_name)){next}
      genomecov_cell_type_2 <- read.table(file_name, sep = "\t")
      genomecov_cell_type_2 <- GRanges(seqnames = genomecov_cell_type_2$V1,
                                       ranges = IRanges(start = genomecov_cell_type_2$V2,
                                                        end = genomecov_cell_type_2$V3))
      genomecov_cell_type_2 <- reduce(genomecov_cell_type_2)
      
      intersection <- intersect(genomecov_cell_type_1, genomecov_cell_type_2)
      cell_type_intersection[[paste0(case_id, "_", cell_type_1, "_", cell_type_2)]] <- reduce(intersection)
    }
  }
}
#save.image("./results/QC_metrics/coverage_cell_type_a2s0.RData")
#save.image("./results/QC_metrics/coverage_cell_type_a6s3.RData")

overlap_region_size <- c()
for (case_id in case_id_list){
  print(case_id)
  for (cell_type_1 in c("EN", "IN", "Astrocyte", "Microglia", "Oligodendrocyte", "OPC")){
    print(cell_type_1)
    #file_name <- paste0("./results/clonal_mutation/multicell_covered_region_celltype_ATAC_recovered/", case_id, "_", cell_type_1, "_open_region_a2s0_coverage_count.bed")
    file_name <- paste0("./results/QC_metrics/combined_covered_regions_a6s3/", case_id, "_", cell_type_1, "_open_region_a6s3_sorted.bed")
    if(!file.exists(file_name)){next}
    genomecov_cell_type_1 <- read.table(file_name, sep = "\t")
    genomecov_cell_type_1 <- GRanges(seqnames = genomecov_cell_type_1$V1,
                                     ranges = IRanges(start = genomecov_cell_type_1$V2,
                                                      end = genomecov_cell_type_1$V3))
    genomecov_cell_type_1 <- reduce(genomecov_cell_type_1)
    covered_size <- sum(width(genomecov_cell_type_1))
    
    other_cell_type <- cell_type[cell_type != cell_type_1]
    
    overlap_regions <- GRanges()
    for (cell_type_2 in other_cell_type){
      if(paste0(case_id, "_", cell_type_1, "_", cell_type_2) %in% names(cell_type_intersection)){
        overlap_regions <- GenomicRanges::union(overlap_regions, cell_type_intersection[[paste0(case_id, "_", cell_type_1, "_", cell_type_2)]])
      }
    }
    overlap_regions <- reduce(overlap_regions)
    overlap_size <- sum(width(overlap_regions))
    
    overlap_region_size <- rbind(overlap_region_size, c(case_id, cell_type_1, covered_size, overlap_size))
  }
}
#saveRDS(overlap_region_size, "./results/QC_metrics/overlap_coverage_size_cell_type_a2s0.rds")
#saveRDS(overlap_region_size, "./results/QC_metrics/overlap_coverage_size_cell_type_a6s3.rds")

overlap_region_size <- readRDS("./results/QC_metrics/overlap_coverage_size_cell_type_a2s0.rds")
overlap_region_size <- data.frame(overlap_region_size)
colnames(overlap_region_size) <- c("case_id", "cell_type", "covered_size", "overlap_size")
overlap_region_size$covered_size <- as.integer(overlap_region_size$covered_size)
overlap_region_size$overlap_size <- as.integer(overlap_region_size$overlap_size)
overlap_region_size[,"overlap_proportion"] <- overlap_region_size$overlap_size/overlap_region_size$covered_size

overlap_region_size$case_id <- factor(overlap_region_size$case_id, 
                                                  levels = c(case_id_list),
                                                  labels = c("UMB4428", "UMB1278", "UMB1864", "UMB4638", "UMB1465", 
                                                             "UMB5451", "Brain1901106", "UMB5657", "UMB5823"))

ggplot(overlap_region_size, aes(x = cell_type, y = overlap_proportion)) +
  geom_col() +
  facet_wrap(.~case_id, scales = "free_y") +
  labs(x = "Cell type", y = "Overlapping coverage of other cell types\n(% of covered regions)") +
  theme_classic() +
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75, 1), limits = c(0, 1),
                     labels = scales::percent(c(0.25, 0.5, 0.75, 1))) +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
  )
#ggsave("./figures/QC_metrics/overlap_coverage_size_cell_type_a2s0.pdf", width = 15, height = 10, units = "cm")
ggsave("./figures/QC_metrics/overlap_coverage_size_cell_type_a6s3.pdf", width = 15, height = 10, units = "cm")

overlap_region_size_summary <- overlap_region_size %>% group_by(cell_type) %>% 
  summarise("mean" = mean(overlap_proportion),
            "median" = median(overlap_proportion),
            "5_percentile" = quantile(overlap_proportion, 0.05),
            "95_percentile" = quantile(overlap_proportion, 0.95))

ggplot(overlap_region_size_summary, aes(x = cell_type, y = median)) +
  geom_col() +
  geom_errorbar(aes(ymin = `5_percentile`, ymax = `95_percentile`), width = 0.3) + 
  labs(x = "Cell type", y = "Overlapping coverage of other cell types\n(% of covered regions)") +
  theme_classic() +
  scale_y_continuous(breaks = c(0.25, 0.5), limits = c(0, 0.6),
                     labels = scales::percent(c(0.25, 0.5))) +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
  )
ggsave("./figures/QC_metrics/overlap_coverage_size_cell_type_a6s3_summary.pdf", width = 8, height = 6.5, units = "cm")
