.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(stringr)
library(tidyr)
library(rngtools)
library(stringr)
library(ggpubr)
library(VennDiagram)

ref_genome="BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = T)
chr_orders <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages_sig/")
library(MutationalPatterns)

spectrum_mat <- function(clonal_mutation_list){
  clonal_mutation_list_grange <- GRanges(seqnames = clonal_mutation_list$chr, 
                                         ranges = IRanges(start = clonal_mutation_list$pos,
                                                          end = clonal_mutation_list$pos),
                                         ref = clonal_mutation_list$ref, alt = clonal_mutation_list$alt)
  chr_length <- seqlengths(Hsapiens)
  seqlengths(clonal_mutation_list_grange) <- chr_length[names(seqlengths(clonal_mutation_list_grange))]
  seqlevels(clonal_mutation_list_grange) <- seqlevels(clonal_mutation_list_grange)[order(factor(seqlevels(clonal_mutation_list_grange), levels = chr_orders))]
  genome(clonal_mutation_list_grange) <- "hg38"
  
  mut_mat <- mut_matrix(vcf_list = clonal_mutation_list_grange, ref_genome = ref_genome)
}

# generate mutation context
mutation_context <- function(clonal_mutation_list){
  clonal_mutation_list_grange <- GRanges(seqnames = clonal_mutation_list$chr, 
                                         ranges = IRanges(start = clonal_mutation_list$pos,
                                                          end = clonal_mutation_list$pos),
                                         ref = clonal_mutation_list$ref, alt = clonal_mutation_list$alt)
  chr_length <- seqlengths(Hsapiens)
  seqlengths(clonal_mutation_list_grange) <- chr_length[names(seqlengths(clonal_mutation_list_grange))]
  seqlevels(clonal_mutation_list_grange) <- seqlevels(clonal_mutation_list_grange)[order(factor(seqlevels(clonal_mutation_list_grange), levels = chr_orders))]
  genome(clonal_mutation_list_grange) <- "hg38"
  
  type <- type_context(clonal_mutation_list_grange, ref_genome = ref_genome)
}

# ------------------------------------------------------------------------------
# read clonal mutation list
clonal_mutation_list <- read.csv("results/clonal_mutation/clonal_mutation_list_COLO829BLT50_BLT50_bulk_realigned_20241007.csv")
clonal_mutation_list$key <- str_replace_all(clonal_mutation_list$key, "-", ":")
clonal_mutation_list <- na.omit(clonal_mutation_list)
clonal_mutation_list %>% group_by(key, cell_type) %>% summarise("count" = n()) -> clonal_mutation_count
clonal_mutation_count %>% group_by(key) %>% filter(count == max(count)) %>% select(cell_type) -> clonal_mutation_list
clonal_mutation_list <- separate(clonal_mutation_list, "key", c("chr", "pos", "ref", "alt"), sep = ":", remove = F)
clonal_mutation_list$pos <- as.integer(clonal_mutation_list$pos)

clonal_mutation_list_BL_bulk <- read.csv("results/clonal_mutation/clonal_mutation_list_COLO829BLT50_BL_bulk_realigned_20241007.csv")

clonal_mutation_list_BL_bulk <- na.omit(clonal_mutation_list_BL_bulk)
clonal_mutation_list_BL_bulk %>% group_by(key, cell_type) %>% summarise("count" = n()) -> clonal_mutation_count_BL_bulk
clonal_mutation_count_BL_bulk %>% group_by(key) %>% filter(count == max(count)) %>% select(cell_type) -> clonal_mutation_list_BL_bulk
clonal_mutation_list_BL_bulk <- separate(clonal_mutation_list_BL_bulk, "key", c("chr", "pos", "ref", "alt"), sep = ":", remove = F)
clonal_mutation_list_BL_bulk$pos <- as.integer(clonal_mutation_list_BL_bulk$pos)

clonal_mutation_list_tumor_bulk <- read.csv("results/clonal_mutation/clonal_mutation_list_COLO829BLT50_tumor_bulk_realigned_20241007.csv")

clonal_mutation_list_tumor_bulk <- na.omit(clonal_mutation_list_tumor_bulk)
clonal_mutation_list_tumor_bulk %>% group_by(key, cell_type) %>% summarise("count" = n()) -> clonal_mutation_count_tumor_bulk
clonal_mutation_count_tumor_bulk %>% group_by(key) %>% filter(count == max(count)) %>% select(cell_type) -> clonal_mutation_list_tumor_bulk
clonal_mutation_list_tumor_bulk <- separate(clonal_mutation_list_tumor_bulk, "key", c("chr", "pos", "ref", "alt"), sep = ":", remove = F)
clonal_mutation_list_tumor_bulk$pos <- as.integer(clonal_mutation_list_tumor_bulk$pos)

# ------------------------------------------------------------------------------
# melanoma_fibroblast
melanoma_fibroblast_clonal_mutation_list <- clonal_mutation_list[clonal_mutation_list$cell_type == "melanoma_fibroblast",]
melanoma_fibroblast_clonal_mutation_list_BL_bulk <- clonal_mutation_list_BL_bulk[clonal_mutation_list_BL_bulk$cell_type == "melanoma_fibroblast",]
melanoma_fibroblast_clonal_mutation_list_tumor_bulk <- clonal_mutation_list_tumor_bulk[clonal_mutation_list_tumor_bulk$cell_type == "melanoma_fibroblast",]

spectrum_mat(melanoma_fibroblast_clonal_mutation_list) -> melanoma_fibroblast_clonal_mutation_list_mat
colnames(melanoma_fibroblast_clonal_mutation_list_mat) <- "Melanoma"
plot_96_profile(melanoma_fibroblast_clonal_mutation_list_mat, ymax = 0.35)
ggsave("./figures/mutation_spectrum/COLO829BLT50_clonal_mutation/melanoma_fibroblast_BLT50_bulk_20241007.pdf", height = 3, width = 14)

plot_96_profile(spectrum_mat(melanoma_fibroblast_clonal_mutation_list_BL_bulk), ymax = 0.35)
ggsave("./figures/mutation_spectrum/COLO829BLT50_clonal_mutation/melanoma_fibroblast_BL_bulk_20241007.pdf", height = 3, width = 14)
melanoma_fibroblast_clonal_mutation_list_BL_bulk_spectrum_mat <- spectrum_mat(melanoma_fibroblast_clonal_mutation_list_BL_bulk)

plot_96_profile(spectrum_mat(melanoma_fibroblast_clonal_mutation_list_tumor_bulk))
ggsave("./figures/mutation_spectrum/COLO829BLT50_clonal_mutation/melanoma_fibroblast_tumor_bulk_20241007.pdf", height = 3, width = 14)

melanoma_fibroblast_clonal_mutation_list_all <- rbind(melanoma_fibroblast_clonal_mutation_list, melanoma_fibroblast_clonal_mutation_list_BL_bulk)
melanoma_fibroblast_clonal_mutation_list_all <- melanoma_fibroblast_clonal_mutation_list_all[!duplicated(melanoma_fibroblast_clonal_mutation_list_all),]

plot_96_profile(spectrum_mat(melanoma_fibroblast_clonal_mutation_list_all))
ggsave("./figures/mutation_spectrum/COLO829BLT50_clonal_mutation/melanoma_fibroblast_BLT50_bulk_BL_bulk_20241007.pdf", height = 3, width = 14)

# ------------------------------------------------------------------------------
# B_lymphoblast
B_lymphoblast_clonal_mutation_list <- clonal_mutation_list[clonal_mutation_list$cell_type == "B_lymphoblast",]
B_lymphoblast_clonal_mutation_list_tumor_bulk <- clonal_mutation_list_tumor_bulk[clonal_mutation_list_tumor_bulk$cell_type == "B_lymphoblast",]
B_lymphoblast_clonal_mutation_list_BL_bulk <- clonal_mutation_list_BL_bulk[clonal_mutation_list_BL_bulk$cell_type == "B_lymphoblast",]

spectrum_mat(B_lymphoblast_clonal_mutation_list) -> B_lymphoblast_clonal_mutation_list_mat
colnames(B_lymphoblast_clonal_mutation_list_mat) <- "B-cell"
plot_96_profile(B_lymphoblast_clonal_mutation_list_mat, ymax = 0.08) +
  scale_y_continuous(breaks=seq(0, 0.08, 0.02))
ggsave("./figures/mutation_spectrum/COLO829BLT50_clonal_mutation/B_lymphoblast_BLT50_bulk_20241007.pdf", height = 3, width = 14)

plot_96_profile(spectrum_mat(B_lymphoblast_clonal_mutation_list_tumor_bulk), ymax = 0.08) +
  scale_y_continuous(breaks=seq(0, 0.08, 0.02))
ggsave("./figures/mutation_spectrum/COLO829BLT50_clonal_mutation/B_lymphoblast_tumor_bulk_20241007.pdf", height = 3, width = 14)
B_lymphoblast_clonal_mutation_list_tumor_bulk_spectrum_mat <- spectrum_mat(B_lymphoblast_clonal_mutation_list_tumor_bulk)

plot_96_profile(spectrum_mat(B_lymphoblast_clonal_mutation_list_BL_bulk), ymax = 0.08) +
  scale_y_continuous(breaks=seq(0, 0.08, 0.02))
ggsave("./figures/mutation_spectrum/COLO829BLT50_clonal_mutation/B_lymphoblast_BL_bulk_20241007.pdf", height = 3, width = 14)

# Overlap between BLT50 bulk and tumor bulk
sum(!B_lymphoblast_clonal_mutation_list$key %in% B_lymphoblast_clonal_mutation_list_tumor_bulk$key)
B_lymphoblast_clonal_mutation_list[(B_lymphoblast_clonal_mutation_list$key %in% B_lymphoblast_clonal_mutation_list_tumor_bulk$key),] -> B_lymphoblast_clonal_mutation_list_overlap
B_lymphoblast_clonal_mutation_list[!(B_lymphoblast_clonal_mutation_list$key %in% B_lymphoblast_clonal_mutation_list_tumor_bulk$key),] -> B_lymphoblast_clonal_mutation_list_nonoverlap
left_join(B_lymphoblast_clonal_mutation_list_nonoverlap, clonal_mutation_count) -> B_lymphoblast_clonal_mutation_list_nonoverlap
left_join(B_lymphoblast_clonal_mutation_list_overlap, clonal_mutation_count) -> B_lymphoblast_clonal_mutation_list_overlap

summary(B_lymphoblast_clonal_mutation_list_nonoverlap$count)
summary(B_lymphoblast_clonal_mutation_list_overlap$count)

venn.diagram(x = list(B_lymphoblast_clonal_mutation_list$key, B_lymphoblast_clonal_mutation_list_tumor_bulk$key),
             category.names = c("BLT50 bulk WGS", "Tumor bulk WGS"),
             imagetype = "png",
             filename = "./figures/mutation_spectrum/COLO829BLT50_clonal_mutation/Venn_plot_B_lymphoblast_BLT50_bulk_overlap_tumor_bulk_20241007.png",
             output = F,
             height = 480 , 
             width = 480 , 
             resolution = 300,
             compression = "lzw",
             lwd = 1,
             col=c("#440154ff", '#21908dff'),
             fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
             cex = 0.5,
             fontfamily = "sans",
             cat.cex = 0.3,
             cat.pos = c(-27, 27),
             cat.dist = c(0.055, 0.055),
             cat.default.pos = "outer",
             cat.fontfamily = "sans",
             cat.col = c("#440154ff", '#21908dff'),
             margin = 0.1)

plot_96_profile(spectrum_mat(B_lymphoblast_clonal_mutation_list_nonoverlap))
ggsave("./figures/mutation_spectrum/COLO829BLT50_clonal_mutation/B_lymphoblast_BLT50_bulk_nonoverlap_tumor_bulk_20241007.pdf", height = 3, width = 14)

plot_96_profile(spectrum_mat(B_lymphoblast_clonal_mutation_list_overlap))
ggsave("./figures/mutation_spectrum/COLO829BLT50_clonal_mutation/B_lymphoblast_BLT50_bulk_overlap_tumor_bulk_20241007.pdf", height = 3, width = 14)

clonal_mutation_list_spectrum_mat <- list()
clonal_mutation_list_spectrum_mat[["B_lymphoblast_clonal_mutation_list_tumor_bulk"]] <- B_lymphoblast_clonal_mutation_list_tumor_bulk_spectrum_mat
clonal_mutation_list_spectrum_mat[["melanoma_fibroblast_clonal_mutation_list_BL_bulk"]] <- melanoma_fibroblast_clonal_mutation_list_BL_bulk_spectrum_mat
saveRDS(clonal_mutation_list_spectrum_mat, "./figures/manuscript_figures/Figure2/spectrum_clonal_mutation_COLO829BLT50.rds")

# ------------------------------------------------------------------------------
# Signature contribution
sbs96_selected <- get_known_signatures()

spectrum_mat(melanoma_fibroblast_clonal_mutation_list_BL_bulk) -> melanoma_fibroblast_clonal_mutation_list_BL_bulk_mat
colnames(melanoma_fibroblast_clonal_mutation_list_BL_bulk_mat) <- "Melanoma"

fit <- fit_to_signatures(melanoma_fibroblast_clonal_mutation_list_BL_bulk_mat, sbs96_selected)
plot_contribution(fit$contribution,
                  coord_flip = FALSE,
                  mode = "relative")+
  labs(title = "Vairants surviving BL cell\nbulk WGS filter")
ggsave("./figures/mutation_spectrum/COLO829BLT50_clonal_mutation/melanoma_fibroblast_BL_bulk_signature_analysis_20241007.pdf", height = 5, width = 3)

fit <- fit_to_signatures(spectrum_mat(B_lymphoblast_clonal_mutation_list), sbs96_selected)
plot_contribution(fit$contribution,
                  coord_flip = FALSE,
                  mode = "relative")

spectrum_mat(B_lymphoblast_clonal_mutation_list_tumor_bulk) -> B_lymphoblast_clonal_mutation_list_tumor_bulk_mat
colnames(B_lymphoblast_clonal_mutation_list_tumor_bulk_mat) <- "B-cell"
fit <- fit_to_signatures(B_lymphoblast_clonal_mutation_list_tumor_bulk_mat, sbs96_selected)
plot_contribution(fit$contribution,
                  coord_flip = FALSE,
                  mode = "relative")+
  labs(title = "Vairants surviving Tumor cell\nbulk WGS filter")
ggsave("./figures/mutation_spectrum/COLO829BLT50_clonal_mutation/B_lymphoblast_tumor_bulk_signature_analysis_20241007.pdf", height = 5, width = 3)

