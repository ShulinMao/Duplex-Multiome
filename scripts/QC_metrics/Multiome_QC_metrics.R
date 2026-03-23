# Libraries
.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(scales)
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

case_id_list <- c("UMB4428", "UMB1278_deeper", "UMB1864_deeper", "UMB4638", "UMB1465", 
                  "UMB5451", "Brain1901106_deeper", "UMB5657", "UMB5823_deeper")
COLO829BLT_case_id_list <- c("COLO829BLT50_rep1", "COLO829BLT50_rep2")
brain_COLO829BLT_case_id_list <- c(case_id_list, COLO829BLT_case_id_list)

# ------------------------------------------------------------------------------
# actual number of cells
brain_cell_type <- read.table("./results/celltype_annotation/annotation_table_20240512.txt", sep = "\t")
colnames(brain_cell_type) <- c("case_id", "cell_id", "cell_type")
brain_cell_type[!str_detect(brain_cell_type$cell_type, "doublet"),] -> brain_cell_type
brain_cell_type$cell_type[str_detect(brain_cell_type$cell_type, "EN")] <- "Excitatory\nNeuron"
brain_cell_type$cell_type[str_detect(brain_cell_type$cell_type, "IN")] <- "Inhibitory\nNeuron"
brain_cell_type$case_id <- factor(brain_cell_type$case_id, 
                                  levels = case_id_list,
                                  labels = c("UMB4428", "UMB1278", "UMB1864", "UMB4638", "UMB1465", 
                                             "UMB5451", "Brain1901106", "UMB5657", "UMB5823"))

COLO829BLT50_cell_type <- read.table("./results/celltype_annotation/annotation_COLO829BLT50_20240523.txt", sep = "\t")
colnames(COLO829BLT50_cell_type) <- c("case_id", "cell_id", "cell_type")

colors_cell_type <- c("#2E2585", "#337538", "#7CAE00", "#E69F00", "#DCCD7D", "#00A9FF", "#CC79A7", "#7E2954")  
ggplot(brain_cell_type, aes(x = case_id)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 2.5) +
  geom_bar(aes(fill = cell_type)) +
  theme_classic() +
  labs(x = "Sample", y = "Cell number", fill = "Cell type") +
  scale_fill_manual(values = colors_cell_type) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
  )
ggsave("./figures/QC_metrics/cell_type_count_brain.pdf", width = 10, height = 8, units = "cm")

ggplot(brain_cell_type, aes(x = case_id)) +
  geom_bar(aes(fill = cell_type), position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Sample", y = "Cell type composition (%)", fill = "Cell type") +
  scale_fill_manual(values = colors_cell_type) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
  )
ggsave("./figures/QC_metrics/cell_type_percentage_brain.pdf", width = 10, height = 8, units = "cm")

ggplot(COLO829BLT50_cell_type, aes(x = case_id)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 2.5) +
  geom_bar(aes(fill = cell_type)) +
  theme_classic() +
  labs(x = "Sample", y = "Cell number", fill = "Cell type") +
  scale_fill_manual(values = c("#A34C49", "lightskyblue1")) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
  )
ggsave("./figures/QC_metrics/cell_type_count_COLO829BLT50.pdf", width = 7, height = 8, units = "cm")

ggplot(COLO829BLT50_cell_type, aes(x = case_id)) +
  geom_bar(aes(fill = cell_type), position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_classic() +
  labs(x = "Sample", y = "Cell type composition (%)", fill = "Cell type") +
  scale_fill_manual(values = c("#A34C49", "lightskyblue1")) +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.text = element_text(color = "black", size = axis_font_size, family = font),
  )
ggsave("./figures/QC_metrics/cell_type_percentage_COLO829BLT50.pdf", width = 7, height = 8, units = "cm")

# ------------------------------------------------------------------------------
# reads in peak regions
data_path <- "/lab-share/Gene-Walsh-e2/Public/indData/AK/Multiome-strand-02212022/"
folder_names <- c(
  "UMB1465" = "Novogene_11022023/U1465_Mult_addl-RNA/outs",
  "UMB5451" = "Novogene_11022023/U5451_Mult_addl-RNA/outs",
  "UMB5657" = "Novogene_11022023/U5657_Mult_addl-RNA/outs",
  "UMB4638" = "Novogene_01312024/U4638_Mult/outs",
  "UMB4428" = "Novogene_04052024/4428_B/outs",
  "UMB1278_deeper" = "Novogene_04292024/U1278_deeper_05012024/outs",
  "UMB1864_deeper" = "Novogene_04292024/U1864_deeper_05012024/outs",
  "UMB5823_deeper" = "Novogene_04292024/U5823_deeper_05012024/outs",
  "Brain1901106_deeper" = "Novogene_04292024/B1901106_deeper_05012024/outs",
  "COLO829BLT50_rep1" = "Novogene_04292024/BLT50_A/outs", 
  "COLO829BLT50_rep2" = "Novogene_04292024/BLT50_B/outs" 
)

read_in_peak <- c()

for (sample_id in brain_COLO829BLT_case_id_list){
  print(sample_id)
  peak_region <- read.table(paste0(data_path, folder_names[sample_id], "/atac_peaks.bed"), sep="\t")
  peak_region <- peak_region[peak_region$V1 %in% paste0("chr", seq(1,22)),]
  peak_region <- GRanges(seqnames = peak_region$V1, 
                         ranges = IRanges(start = peak_region$V2,
                                          end = peak_region$V3))
  peak_region <- reduce(peak_region)
  peak_region_size <- sum(width(peak_region))
  
  dir_path <- paste0("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/", sample_id, "/single_cell/read_families/")
  dir(dir_path, pattern = "*families.bed$", full.names = T) -> bed_file_list
  
  for (bed_file in bed_file_list){
    read_family <- read.table(bed_file, sep = "\t")
    read_family <- GRanges(seqnames = read_family$V1, 
                           ranges = IRanges(start = read_family$V2,
                                            end = read_family$V3),
                           strand1 = read_family$V7, 
                           strand2 = read_family$V8,
                           rf_id = read_family$V4)
    
    intersection <- findOverlaps(read_family, peak_region)
    intersected_ranges <- read_family[queryHits(intersection)]
    
    for (a in c(2, 4, 6)){
      # with duplex
      s <- a %/% 2
      stringency <- str_c("a", a, "s", s)
      rf_in_peak <- intersected_ranges[intersected_ranges$strand1 >= s*2 & intersected_ranges$strand2 >= s*2, ]
      rf_total <- read_family[read_family$strand1 >= s*2 & read_family$strand2 >= s*2,]
      
      rf_count_in_peak <- length(unique(rf_in_peak$rf_id))
      rf_count_total <- length(unique(rf_total$rf_id))
      
      read_in_peak <- rbind(read_in_peak, c(sample_id, stringency, rf_count_in_peak, rf_count_total))
      
      # without duplex
      s <- 0
      stringency <- str_c("a", a, "s", s)
      rf_in_peak <- intersected_ranges[intersected_ranges$strand1 + intersected_ranges$strand2 >= a*2, ]
      rf_total <- read_family[read_family$strand1 + read_family$strand2 >= a*2,]
      
      rf_count_in_peak <- length(unique(rf_in_peak$rf_id))
      rf_count_total <- length(unique(rf_total$rf_id))
      
      read_in_peak <- rbind(read_in_peak, c(sample_id, stringency, rf_count_in_peak, rf_count_total))
    }
  }
}
#saveRDS(read_in_peak, "./results/QC_metrics/read_in_peak.rds")

read_in_peak <- readRDS("./results/QC_metrics/read_in_peak.rds")
read_in_peak <- data.frame(read_in_peak)
colnames(read_in_peak) <- c("case_id", "stringency", "rf_count_in_peak", "rf_count_total")
read_in_peak$rf_count_in_peak <- as.integer(read_in_peak$rf_count_in_peak)
read_in_peak$rf_count_total <- as.integer(read_in_peak$rf_count_total)
read_in_peak$rf_proportion_in_peak <- ifelse(read_in_peak$rf_count_total > 0,
                                             read_in_peak$rf_count_in_peak/read_in_peak$rf_count_total, 0)

ggplot(read_in_peak, aes(x = rf_proportion_in_peak)) +
  geom_histogram(fill = "skyblue", color = "black") +
  facet_wrap(.~stringency) +
  scale_x_continuous(labels = scales::percent_format()) +
  labs(x = "% of read families overlapping peak", y = "Cell number") +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
  )
ggsave("./figures/QC_metrics/proportion_rf_in_peak.pdf", width = 12, height = 10, units = "cm")

ggplot(read_in_peak[read_in_peak$stringency == "a6s3",], aes(x = rf_proportion_in_peak)) +
  geom_histogram(fill = "skyblue", color = "black") +
  facet_wrap(.~case_id) +
  scale_x_continuous(labels = scales::percent_format()) +
  labs(x = "% of read families overlapping peak", y = "Cell number", title = "a6s3") +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "none",
  )

# ------------------------------------------------------------------------------
# TSS enrichment
library(rvest)
library(xml2) 
library(jsonlite)

data_path <- "/lab-share/Gene-Walsh-e2/Public/indData/AK/Multiome-strand-02212022/"
folder_names <- c(
  "UMB1465" = "Novogene_11022023/U1465_Mult_addl-RNA/outs",
  "UMB5451" = "Novogene_11022023/U5451_Mult_addl-RNA/outs",
  "UMB5657" = "Novogene_11022023/U5657_Mult_addl-RNA/outs",
  "UMB4638" = "Novogene_01312024/U4638_Mult/outs",
  "UMB4428" = "Novogene_04052024/4428_B/outs",
  "UMB1278_deeper" = "Novogene_04292024/U1278_deeper_05012024/outs",
  "UMB1864_deeper" = "Novogene_04292024/U1864_deeper_05012024/outs",
  "UMB5823_deeper" = "Novogene_04292024/U5823_deeper_05012024/outs",
  "Brain1901106_deeper" = "Novogene_04292024/B1901106_deeper_05012024/outs",
  "COLO829BLT50_rep1" = "Novogene_04292024/BLT50_A/outs", 
  "COLO829BLT50_rep2" = "Novogene_04292024/BLT50_B/outs" 
)

TSS_enrichment_curve <- c()
for (sample_id in brain_COLO829BLT_case_id_list){
  print(sample_id)
  
  # Read the HTML content from the file
  web_summary <- read_html(paste0(data_path, folder_names[sample_id], "/web_summary.html"))
  # Use CSS selector to find all <script> tags with type="text/javascript"
  javascript_nodes <- html_nodes(web_summary, 'script[type="text/javascript"]')
  # Extract the text content from these nodes
  javascript_code <- html_text(javascript_nodes)
  # Find the starting position of the JSON
  json_start_pos <- str_locate(javascript_code, "const data = ")[[1, 2]] + 1
  json_string_partial <- str_sub(javascript_code, json_start_pos)
  json_object <- fromJSON(json_string_partial)
  x <- json_object$atac_tss_enrichment_plot$data$x[[1]]
  y <- json_object$atac_tss_enrichment_plot$data$y[[1]]
  
  TSS_enrichment_curve <- rbind(TSS_enrichment_curve, cbind(sample_id, x, y))
}

TSS_enrichment_curve <- data.frame(TSS_enrichment_curve)
colnames(TSS_enrichment_curve) <- c("case_id", "position", "enrichment")
TSS_enrichment_curve$position <- as.numeric(TSS_enrichment_curve$position)
TSS_enrichment_curve$enrichment <- as.numeric(TSS_enrichment_curve$enrichment)
TSS_enrichment_curve$case_id <- factor(TSS_enrichment_curve$case_id, 
                                       levels = c(brain_COLO829BLT_case_id_list),
                                       labels = c("UMB4428", "UMB1278", "UMB1864", "UMB4638", "UMB1465", 
                                                  "UMB5451", "Brain1901106", "UMB5657", "UMB5823", 
                                                  "COLO829BLT50_rep1", "COLO829BLT50_rep2"))

ggplot(TSS_enrichment_curve, aes(x = position, y = enrichment, color = case_id)) +
  geom_line(size = 0.05) +
  labs(y = "Relative enrichment", x = "Position (bp from TSS)",
       color = "Sample") +
  theme_classic() +
  theme(
    text = element_text(color = "black", size = font_size, family = font),
    axis.text.x = element_text(color = "black", size = axis_font_size, family = font, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = axis_font_size, family = font),
    legend.position = "right",
    legend.key.size = unit(0.3, "cm")
  )
ggsave("./figures/QC_metrics/TSS_enrichment.pdf", width = 8, height = 5, units = "cm")
