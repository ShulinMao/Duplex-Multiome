# Libraries
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

ref_genome="BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = T)
chr_orders <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages_sig/")
library(MutationalPatterns)


# ------------------------------------------------------------------------------
# Healthy Human samples
case_id_list <- c("UMB1278_deeper", "UMB1864_deeper", "UMB4638", "UMB1465", "UMB5451", "Brain1901106_deeper", "UMB5657",  "UMB5823_deeper")
ssDNA_mut_grangelist <- GRangesList()

# convert vcf files to grange objects
for (case_id in case_id_list){
  print(case_id)
  ssDNA_mut <- readRDS(paste0("./results/ssDNA_mut/", case_id, "_spectrum.rds"))
  
  # cell type
  file_list <- dir(paste0("./data/", case_id, "/cell_type/"), pattern = "*.txt")
  cell_type_list <- str_remove(file_list, "[.]txt")
  
  # raw counts
  for(cell_type in cell_type_list){
    ssDNA_mut_cell_type <- ssDNA_mut[[cell_type]]
    ssDNA_mut_cell_type_grange <- GRanges(seqnames = ssDNA_mut_cell_type$V1, 
                                          ranges = IRanges(start = ssDNA_mut_cell_type$V2,
                                                           end = ssDNA_mut_cell_type$V2),
                                          ref = ssDNA_mut_cell_type$V4, 
                                          alt = ssDNA_mut_cell_type$V5)
    ssDNA_mut_grangelist[[paste0(case_id, "_", cell_type, "_Original")]] <- ssDNA_mut_cell_type_grange
  }
}


# signature analysis
chr_length <- seqlengths(Hsapiens)
seqlengths(ssDNA_mut_grangelist) <- chr_length[names(seqlengths(ssDNA_mut_grangelist))]
seqlevels(ssDNA_mut_grangelist) <- seqlevels(ssDNA_mut_grangelist)[order(factor(seqlevels(ssDNA_mut_grangelist), levels = chr_orders))]
genome(ssDNA_mut_grangelist) = 'hg38'
# calculate 96 type of snvs
mut_mat_summary <- mut_matrix(ssDNA_mut_grangelist, ref_genome = ref_genome)

saveRDS(mut_mat_summary, "./results/ssDNA_mut/spectrum_summary.rds")

# spectrum figures
case_id_list <- c("UMB1278_deeper", "UMB1864_deeper", "UMB4638", "UMB1465", "UMB5451", "Brain1901106_deeper", "UMB5657",  "UMB5823_deeper") # increase with ages
cell_type_list <- c("Astrocyte", "EN", "IN", "Microglia", "Oligodendrocyte", "OPC")
stringency <- "a2s0"
spectrum_type <- "Original"
for (cell_type in cell_type_list){
  mut_mat <- mut_mat_summary[,str_detect(colnames(mut_mat_summary), paste0("_", cell_type, "_", spectrum_type, "$"))]
  colnames(mut_mat) <- str_remove(colnames(mut_mat), paste0("_", cell_type, "_", spectrum_type, "$"))
  
  plot_96_profile(mut_mat[, case_id_list[case_id_list %in% colnames(mut_mat)]])
  ggsave(paste0("./figures/single_strand_damage_spectrum/", stringency, "/by_celltype/", cell_type, "_", spectrum_type, ".pdf"), height = 16, width = 14)
}

spectrum_type <- "Original"
for (case_id in case_id_list){
  mut_mat <- mut_mat_summary[,str_detect(colnames(mut_mat_summary), paste0("^", case_id, "_.*_", spectrum_type, "$"))]
  colnames(mut_mat) <- str_remove(colnames(mut_mat), paste0("_", spectrum_type, "$"))
  colnames(mut_mat) <- str_remove(colnames(mut_mat), paste0("^", case_id, "_"))
  
  plot_96_profile(mut_mat[,cell_type_list[cell_type_list %in% colnames(mut_mat)]])
  ggsave(paste0("./figures/single_strand_damage_spectrum/", stringency, "/by_case/", case_id, "_", spectrum_type, ".pdf"), height = 12, width = 14)
}

single_strand_damage_spectrum_df <- c()
for (cell_type in cell_type_list){
  mut_mat <- mut_mat_summary[,str_detect(colnames(mut_mat_summary), paste0("_", cell_type, "_", spectrum_type, "$"))]
  mut_mat <- rowSums(mut_mat)
  
  single_strand_damage_spectrum_df <- cbind(single_strand_damage_spectrum_df, mut_mat)
}
colnames(single_strand_damage_spectrum_df) <- c(cell_type_list)

plot_96_profile(single_strand_damage_spectrum_df, ymax = 0.1) +
  scale_y_continuous(breaks = seq(0, 0.1, 0.05)) 
ggsave(paste0("./figures/single_strand_damage_spectrum/", stringency, "/by_celltype/single_strand_damage_spectrum_all_combined_orginal.pdf"), height = 12, width = 14)

sbs30 <- get_known_signatures()[,c("SBS30"), drop=F]
colnames(sbs30) <- "SBS30"
rownames(sbs30) <- rownames(single_strand_damage_spectrum_df)

plot_96_profile(sbs30, ymax = 0.15) +
  scale_y_continuous(breaks = seq(0, 0.15, 0.05)) 
ggsave(paste0("./figures/single_strand_damage_spectrum/SBS30.pdf"), height = 2.7, width = 14)

df_cos_sim <- apply(single_strand_damage_spectrum_df, 2, function(X){cos_sim(X, sbs30[,1])})
df_cos_sim <- data.frame(df_cos_sim)
df_cos_sim[,"cell_type"] <- row.names(df_cos_sim)

ggplot(df_cos_sim, aes(x = cell_type, y = df_cos_sim)) + 
  geom_col(width = 0.5, fill = "skyblue") +
  geom_text(aes(label = round(df_cos_sim, 2), y = df_cos_sim + 0.05)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = "Cell type", y = "Cosine Similarity")
ggsave(paste0("./figures/single_strand_damage_spectrum/", stringency, "/by_celltype/single_strand_damage_spectrum_cos_sim_SBS30.pdf"), height = 3, width = 5)

# ------------------------------------------------------------------------------
# COLO829BLT50
case_id_list <- c("COLO829BLT50_rep1", "COLO829BLT50_rep2")
ssDNA_mut_grangelist <- GRangesList()

# convert vcf files to grange objects
for (case_id in case_id_list){
  print(case_id)
  ssDNA_mut <- readRDS(paste0("./results/ssDNA_mut/", case_id, "_spectrum.rds"))
  
  # cell type
  file_list <- dir(paste0("./data/", case_id, "/cell_type/"))
  cell_type_list <- str_remove(file_list, "[.]txt")
  
  # raw counts
  for(cell_type in cell_type_list){
    ssDNA_mut_cell_type <- ssDNA_mut[[cell_type]]
    ssDNA_mut_cell_type_grange <- GRanges(seqnames = ssDNA_mut_cell_type$V1, 
                                          ranges = IRanges(start = ssDNA_mut_cell_type$V2,
                                                           end = ssDNA_mut_cell_type$V2),
                                          ref = ssDNA_mut_cell_type$V4, 
                                          alt = ssDNA_mut_cell_type$V5)
    ssDNA_mut_grangelist[[paste0(case_id, "_", cell_type, "_Original")]] <- ssDNA_mut_cell_type_grange
  }
}


# signature analysis
chr_length <- seqlengths(Hsapiens)
seqlengths(ssDNA_mut_grangelist) <- chr_length[names(seqlengths(ssDNA_mut_grangelist))]
seqlevels(ssDNA_mut_grangelist) <- seqlevels(ssDNA_mut_grangelist)[order(factor(seqlevels(ssDNA_mut_grangelist), levels = chr_orders))]
genome(ssDNA_mut_grangelist) = 'hg38'
# calculate 96 type of snvs
mut_mat_summary <- mut_matrix(ssDNA_mut_grangelist, ref_genome = ref_genome)

saveRDS(mut_mat_summary, "./results/ssDNA_mut/COLO829BLT50_spectrum_summary.rds")

