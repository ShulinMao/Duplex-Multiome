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
library(pheatmap)

ref_genome="BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = T)
chr_orders <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages_sig/")
library(MutationalPatterns)

case_id_list <- c("Brain1901106_deeper", "UMB1278_deeper", "UMB1864_deeper", "UMB5823_deeper", "UMB4638", "UMB1465", "UMB5451", "UMB5657")

a <- 6
s <- a%/%2
stringency <- paste0("a", a, "s", s)

# ------------------------------------------------------------------------------
# collect spectrum from all cells
mut_mat_summary <- c()
for(case_id in case_id_list){
  # cell type
  file_list <- dir(paste0("./data/", case_id, "/cell_type/"), pattern = "*.txt")
  cell_type_list <- str_remove(file_list, "[.]txt")

  for (cell_type in cell_type_list){
    read.table(paste0("./data/", case_id, "/cell_type/", cell_type, ".txt"), sep = "\t") -> selected_cells
    
    normalized_spectrum <- rep(0, 96)
    original_spectrum <- rep(0, 96)
    # read in normalized spectrum
    for (cell_id in selected_cells$V2){
      specturm_cell <-
        read.table(paste0("./data/", case_id, "/single_cell/burden_estimation/Dan/", stringency, "/", stringency, "_filterGermline_burden.", cell_id, ".txt"), 
                   skip = 43, comment.char = "", header = T)
      specturm_cell$OriginalCount[is.na(specturm_cell$OriginalCount)] <- 0
      specturm_cell$Adjusted[is.na(specturm_cell$AdjustedCount)] <- 0
      original_spectrum <- original_spectrum + specturm_cell$OriginalCount
      normalized_spectrum <- normalized_spectrum + specturm_cell$AdjustedCount
    }
    
    mut_mat <- matrix(c(normalized_spectrum, original_spectrum), ncol = 2, byrow = F)
    colnames(mut_mat) <- c(paste(case_id, cell_type, "Normalized", sep="_"), paste(case_id, cell_type, "Original", sep="_"))
    row.names(mut_mat) <- 
      paste0(substring(specturm_cell$Context, 1, 1), 
             "[", specturm_cell$Variant, "]", 
             substring(specturm_cell$Context, 3, 3))
    
    mut_mat_summary <- cbind(mut_mat_summary, mut_mat)
  }
}

saveRDS(mut_mat_summary, paste0("./results/spectrum/mutation_spectrum", "_", stringency, "_20240627.rds"))
mut_mat_summary <- readRDS(paste0("./results/spectrum/mutation_spectrum", "_", stringency, "_20240627.rds"))
sum(is.na(mut_mat_summary))
# spectrum figures
case_id_list <- c("UMB1278_deeper", "UMB1864_deeper", "UMB4638", "UMB1465", "UMB5451", "Brain1901106_deeper", "UMB5657",  "UMB5823_deeper") # increase with ages
cell_type_list <- c("Astrocyte", "EN", "IN", "Microglia", "Oligodendrocyte", "OPC")

spectrum_type <- "Normalized"
for (cell_type in cell_type_list){
  mut_mat <- mut_mat_summary[,str_detect(colnames(mut_mat_summary), paste0("_", cell_type, "_", spectrum_type, "$"))]
  colnames(mut_mat) <- str_remove(colnames(mut_mat), paste0("_", cell_type, "_", spectrum_type, "$"))

  plot_96_profile(mut_mat[, case_id_list[case_id_list %in% colnames(mut_mat)]])
  ggsave(paste0("./figures/mutation_spectrum/", stringency, "/by_celltype/", cell_type, "_", spectrum_type, ".pdf"), height = 16, width = 14)
}

spectrum_type <- "Normalized"
for (case_id in case_id_list){
  mut_mat <- mut_mat_summary[,str_detect(colnames(mut_mat_summary), paste0("^", case_id, "_.*_", spectrum_type, "$"))]
  colnames(mut_mat) <- str_remove(colnames(mut_mat), paste0("_", spectrum_type, "$"))
  colnames(mut_mat) <- str_remove(colnames(mut_mat), paste0("^", case_id, "_"))
  colnames(mut_mat) <- str_remove(colnames(mut_mat), paste0("^deeper_"))
  
  plot_96_profile(mut_mat[,cell_type_list[cell_type_list %in% colnames(mut_mat)]])
  ggsave(paste0("./figures/mutation_spectrum/", stringency, "/by_case/", case_id, "_", spectrum_type, ".pdf"), height = 12, width = 14)
}

spectrum_type <- "Original"
uncorrected_spectrum_df <- c()
for (cell_type in cell_type_list){
  mut_mat <- mut_mat_summary[,str_detect(colnames(mut_mat_summary), paste0("_", cell_type, "_", spectrum_type, "$"))]
  mut_mat <- rowSums(mut_mat)
  
  uncorrected_spectrum_df <- cbind(uncorrected_spectrum_df, mut_mat)
}
colnames(uncorrected_spectrum_df) <- cell_type_list
plot_96_profile(uncorrected_spectrum_df, ymax = 0.08) +
  scale_y_continuous(breaks = seq(0, 0.08, 0.04)) 
ggsave(paste0("./figures/mutation_spectrum/", stringency, "/by_celltype/uncorrected_combined_spectrum.pdf"), height = 12, width = 14)


# ------------------------------------------------------------------------------
# spectrum analysis
mut_mat_summary <- readRDS(paste0("./results/spectrum/mutation_spectrum", "_", stringency, "_20240627.rds"))

mut_spectrum_corrected_celltype <- list()

# sbs96 <- read.table("../ref/COSMIC_v3.4_SBS_GRCh38.txt", header = T)
# sbs96$Type -> row.names(sbs96)
# sbs96 <- sbs96[,-1]
# sbs96 <- sbs96[row.names(mut_spectrum_corrected),]
# sbs96 <- as.matrix(sbs96)
# sbs96_selected <- sbs96[,c("SBS5", "SBS16", "SBS1", "SBS32", "SBS19")]
sbs96_selected <- get_known_signatures()[,c("SBS5", "SBS16", "SBS1", "SBS32", "SBS19")]

# 3 mer count in covered region / 3 mer count in the whole genome
context_3mer_ratio <- readRDS(paste0("./results/context_3mer_ratio/", "a", a, "s", s, "_celltype.rds"))
context_3mer_ratio <- context_3mer_ratio[c(1:16, 1:16, 1:16, 17:32, 17:32, 17:32),]

# ------------------------------------------------------------------------------
# remove ss damages
# EN
spectrum_type <- "Original"
cell_type <- "EN"
mut_mat <- mut_mat_summary[,str_detect(colnames(mut_mat_summary), paste0("_", cell_type, "_", spectrum_type, "$"))]
colnames(mut_mat) <- str_remove(colnames(mut_mat), paste0("_", cell_type, "_", spectrum_type, "$"))
mut_mat <- mut_mat[, case_id_list[case_id_list %in% colnames(mut_mat)]]
mut_spectrum <- matrix(rowSums(mut_mat))
rownames(mut_spectrum) <- row.names(mut_mat)
colnames(mut_spectrum) <- "EN_ds_mutation"
p1 <- plot_96_profile(mut_spectrum, ymax = 0.08) + 
  scale_y_continuous(breaks=seq(0, 0.08, 0.02))

# ssDNA damage spectrum
ssDNA_mut_mat_summary <- readRDS("./results/ssDNA_mut/spectrum_summary.rds")
spectrum_type <- "Original"
ssDNA_mut_mat <- ssDNA_mut_mat_summary[,str_detect(colnames(ssDNA_mut_mat_summary), paste0("_", cell_type, "_", spectrum_type, "$"))]
ssDNA_mut_spectrum <- matrix(rowSums(ssDNA_mut_mat))
rownames(ssDNA_mut_spectrum) <- row.names(ssDNA_mut_mat)
colnames(ssDNA_mut_spectrum) <- "EN_ss_damage"
p2 <- plot_96_profile(ssDNA_mut_spectrum, ymax = 0.08) + 
  scale_y_continuous(breaks=seq(0, 0.08, 0.02))

ssDNA_mut_spectrum / sum(ssDNA_mut_spectrum) -> ssDNA_mut_spectrum

sbs96 <- get_known_signatures()
fit <- fit_to_signatures(ssDNA_mut_spectrum, sbs96)
ssDNA_mut_spectrum <- #(ssDNA_mut_spectrum - fit$reconstructed) +
  rowSums(sweep(sbs96[,c("SBS2", "SBS7a", "SBS7b", "SBS11", "SBS19", "SBS23", "SBS30", "SBS32", "SBS40", "SBS42", "SBS84", "SBS44")], 2,
                fit$contribution[c("SBS2", "SBS7a", "SBS7b", "SBS11", "SBS19", "SBS23", "SBS30", "SBS32", "SBS40", "SBS42", "SBS84", "SBS44"),], "*"))
ssDNA_mut_spectrum <- matrix(ssDNA_mut_spectrum)
rownames(ssDNA_mut_spectrum) <- rownames(mut_spectrum)
colnames(ssDNA_mut_spectrum) <- "EN_ss_damage"
plot_96_profile(ssDNA_mut_spectrum, ymax = 0.08) + 
  scale_y_continuous(breaks=seq(0, 0.08, 0.02))

mut_spectrum_corrected <- mut_spectrum - sum(mut_spectrum[,1]) * 0.6 * ssDNA_mut_spectrum[,1]

mut_spectrum_corrected[mut_spectrum_corrected <= 0] <- 0
colnames(mut_spectrum_corrected) <- "EN_ds_mutation_corrected"

# normalized by 3-mer ratio
mut_spectrum_corrected <- mut_spectrum_corrected / context_3mer_ratio[,cell_type]

p3 <- plot_96_profile(mut_spectrum_corrected, ymax = 0.12) + 
  scale_y_continuous(breaks=seq(0, 0.08, 0.02))

ggarrange(p1, p2, p3, ncol=1)
ggsave(paste0("./figures/mutation_spectrum/", stringency, "_", cell_type, ".pdf"), height = 9, width = 14)

mut_spectrum_corrected_celltype[[cell_type]] <- mut_spectrum_corrected

# ------------------------------------------------------------------------------
# OL
spectrum_type <- "Original"
cell_type <- "Oligodendrocyte"
mut_mat <- mut_mat_summary[,str_detect(colnames(mut_mat_summary), paste0("_", cell_type, "_", spectrum_type, "$"))]
colnames(mut_mat) <- str_remove(colnames(mut_mat), paste0("_", cell_type, "_", spectrum_type, "$"))
mut_mat <- mut_mat[, case_id_list[case_id_list %in% colnames(mut_mat)]]
mut_spectrum <- matrix(rowSums(mut_mat))
rownames(mut_spectrum) <- row.names(mut_mat)
colnames(mut_spectrum) <- "OL_ds_mutation"
p1 <- plot_96_profile(mut_spectrum, ymax = 0.08) + 
  scale_y_continuous(breaks=seq(0, 0.08, 0.02))

# ssDNA damage spectrum
ssDNA_mut_mat_summary <- readRDS("./results/ssDNA_mut/spectrum_summary.rds")
spectrum_type <- "Original"
ssDNA_mut_mat <- ssDNA_mut_mat_summary[,str_detect(colnames(ssDNA_mut_mat_summary), paste0("_", cell_type, "_", spectrum_type, "$"))]
ssDNA_mut_spectrum <- matrix(rowSums(ssDNA_mut_mat))
rownames(ssDNA_mut_spectrum) <- row.names(ssDNA_mut_mat)
colnames(ssDNA_mut_spectrum) <- "OL_ss_damage"
p2 <- plot_96_profile(ssDNA_mut_spectrum, ymax = 0.08) + 
  scale_y_continuous(breaks=seq(0, 0.08, 0.02))

ssDNA_mut_spectrum / sum(ssDNA_mut_spectrum) -> ssDNA_mut_spectrum

sbs96 <- get_known_signatures()
fit <- fit_to_signatures(ssDNA_mut_spectrum, sbs96)
ssDNA_mut_spectrum <- #(ssDNA_mut_spectrum - fit$reconstructed) +
  rowSums(sweep(sbs96[,c("SBS2", "SBS7a", "SBS7b", "SBS11", "SBS19", "SBS23", "SBS30", "SBS31", "SBS32", "SBS40", "SBS42", "SBS84")], 2,
                fit$contribution[c("SBS2", "SBS7a", "SBS7b", "SBS11", "SBS19", "SBS23", "SBS30", "SBS31", "SBS32", "SBS40", "SBS42", "SBS84"),], "*"))
ssDNA_mut_spectrum <- matrix(ssDNA_mut_spectrum)
rownames(ssDNA_mut_spectrum) <- rownames(mut_spectrum)
colnames(ssDNA_mut_spectrum) <- "OL_ss_damage"

mut_spectrum_corrected <- mut_spectrum - sum(mut_spectrum[,1]) * 0.7 * ssDNA_mut_spectrum

mut_spectrum_corrected[mut_spectrum_corrected <= 0] <- 0
colnames(mut_spectrum_corrected) <- "OL_ds_mutation_corrected"

# normalized by 3-mer ratio
mut_spectrum_corrected <- mut_spectrum_corrected / context_3mer_ratio[,cell_type]

p3 <- plot_96_profile(mut_spectrum_corrected, ymax = 0.08) + 
  scale_y_continuous(breaks=seq(0, 0.08, 0.02))
ggarrange(p1, p2, p3, ncol=1)
ggsave(paste0("./figures/mutation_spectrum/", stringency, "_", cell_type, ".pdf"), height = 9, width = 14)

mut_spectrum_corrected_celltype[[cell_type]] <- mut_spectrum_corrected

# ------------------------------------------------------------------------------
# IN
spectrum_type <- "Original"
cell_type <- "IN"
mut_mat <- mut_mat_summary[,str_detect(colnames(mut_mat_summary), paste0("_", cell_type, "_", spectrum_type, "$"))]
colnames(mut_mat) <- str_remove(colnames(mut_mat), paste0("_", cell_type, "_", spectrum_type, "$"))
mut_mat <- mut_mat[, case_id_list[case_id_list %in% colnames(mut_mat)]]
mut_spectrum <- matrix(rowSums(mut_mat))
rownames(mut_spectrum) <- row.names(mut_mat)
colnames(mut_spectrum) <- "IN_ds_mutation"
p1 <- plot_96_profile(mut_spectrum, ymax = 0.08) +
  scale_y_continuous(breaks=seq(0, 0.08, 0.02))

# ssDNA damage spectrum
ssDNA_mut_mat_summary <- readRDS("./results/ssDNA_mut/spectrum_summary.rds")
spectrum_type <- "Original"
ssDNA_mut_mat <- ssDNA_mut_mat_summary[,str_detect(colnames(ssDNA_mut_mat_summary), paste0("_", cell_type, "_", spectrum_type, "$"))]
ssDNA_mut_spectrum <- matrix(rowSums(ssDNA_mut_mat))
rownames(ssDNA_mut_spectrum) <- row.names(ssDNA_mut_mat)
colnames(ssDNA_mut_spectrum) <- "IN_ss_damage"
p2 <- plot_96_profile(ssDNA_mut_spectrum, ymax = 0.08) +
  scale_y_continuous(breaks=seq(0, 0.08, 0.02))

ssDNA_mut_spectrum / sum(ssDNA_mut_spectrum) -> ssDNA_mut_spectrum

sbs96 <- get_known_signatures()
fit <- fit_to_signatures(ssDNA_mut_spectrum, sbs96)
ssDNA_mut_spectrum <- #(ssDNA_mut_spectrum - fit$reconstructed) +
  rowSums(sweep(sbs96[,c("SBS2", "SBS7a", "SBS7b", "SBS11", "SBS19", "SBS23", "SBS30", "SBS31", "SBS32", "SBS40", "SBS42", "SBS84")], 2,
                fit$contribution[c("SBS2", "SBS7a", "SBS7b", "SBS11", "SBS19", "SBS23", "SBS30", "SBS31", "SBS32", "SBS40", "SBS42", "SBS84"),], "*"))
ssDNA_mut_spectrum <- matrix(ssDNA_mut_spectrum)
rownames(ssDNA_mut_spectrum) <- rownames(mut_spectrum)
colnames(ssDNA_mut_spectrum) <- "IN_ss_damage"

mut_spectrum_corrected <- mut_spectrum - sum(mut_spectrum[,1]) * 0.7 * ssDNA_mut_spectrum

cor(sbs96, sbs96[,"SBS30"]) > 0.5
mut_spectrum_corrected[mut_spectrum_corrected <= 0] <- 0
colnames(mut_spectrum_corrected) <- "IN_ds_mutation_corrected"

# normalized by 3-mer ratio
mut_spectrum_corrected <- mut_spectrum_corrected / context_3mer_ratio[,cell_type]

p3 <- plot_96_profile(mut_spectrum_corrected, ymax = 0.08) +
  scale_y_continuous(breaks=seq(0, 0.08, 0.02))
ggarrange(p1, p2, p3, ncol=1)
ggsave(paste0("./figures/mutation_spectrum/", stringency, "_", cell_type, ".pdf"), height = 9, width = 14)

mut_spectrum_corrected_celltype[[cell_type]] <- mut_spectrum_corrected

# ------------------------------------------------------------------------------
# Astrocyte
spectrum_type <- "Original"
cell_type <- "Astrocyte"
mut_mat <- mut_mat_summary[,str_detect(colnames(mut_mat_summary), paste0("_", cell_type, "_", spectrum_type, "$"))]
colnames(mut_mat) <- str_remove(colnames(mut_mat), paste0("_", cell_type, "_", spectrum_type, "$"))
mut_mat <- mut_mat[, case_id_list[case_id_list %in% colnames(mut_mat)]]
mut_spectrum <- matrix(rowSums(mut_mat))
rownames(mut_spectrum) <- row.names(mut_mat)
colnames(mut_spectrum) <- "Astrocyte_ds_mutation"
p1 <- plot_96_profile(mut_spectrum, ymax = 0.08) +
  scale_y_continuous(breaks=seq(0, 0.08, 0.02))

# ssDNA damage spectrum
ssDNA_mut_mat_summary <- readRDS("./results/ssDNA_mut/spectrum_summary.rds")
spectrum_type <- "Original"
ssDNA_mut_mat <- ssDNA_mut_mat_summary[,str_detect(colnames(ssDNA_mut_mat_summary), paste0("_", cell_type, "_", spectrum_type, "$"))]
ssDNA_mut_spectrum <- matrix(rowSums(ssDNA_mut_mat))
rownames(ssDNA_mut_spectrum) <- row.names(ssDNA_mut_mat)
colnames(ssDNA_mut_spectrum) <- "Astrocyte_ss_damage"
p2 <- plot_96_profile(ssDNA_mut_spectrum, ymax = 0.08) +
  scale_y_continuous(breaks=seq(0, 0.08, 0.02))

ssDNA_mut_spectrum / sum(ssDNA_mut_spectrum) -> ssDNA_mut_spectrum

sbs96 <- get_known_signatures()
fit <- fit_to_signatures(ssDNA_mut_spectrum, sbs96)
ssDNA_mut_spectrum <- #(ssDNA_mut_spectrum - fit$reconstructed) +
  rowSums(sweep(sbs96[,c("SBS2", "SBS7a", "SBS7b", "SBS11", "SBS19", "SBS23", "SBS30", "SBS31", "SBS32", "SBS40", "SBS42", "SBS84")], 2,
                fit$contribution[c("SBS2", "SBS7a", "SBS7b", "SBS11", "SBS19", "SBS23", "SBS30", "SBS31", "SBS32", "SBS40", "SBS42", "SBS84"),], "*"))
ssDNA_mut_spectrum <- matrix(ssDNA_mut_spectrum)
rownames(ssDNA_mut_spectrum) <- rownames(mut_spectrum)
colnames(ssDNA_mut_spectrum) <- "Astrocyte_ss_damage"

mut_spectrum_corrected <- mut_spectrum - sum(mut_spectrum[,1]) * 0.6 * ssDNA_mut_spectrum

cor(sbs96, sbs96[,"SBS30"]) > 0.5
mut_spectrum_corrected[mut_spectrum_corrected <= 0] <- 0
colnames(mut_spectrum_corrected) <- "Astrocyte_ds_mutation_corrected"

# normalized by 3-mer ratio
mut_spectrum_corrected <- mut_spectrum_corrected / context_3mer_ratio[,cell_type]

p3 <- plot_96_profile(mut_spectrum_corrected, ymax = 0.14) +
  scale_y_continuous(breaks=seq(0, 0.14, 0.02))
ggarrange(p1, p2, p3, ncol=1)
ggsave(paste0("./figures/mutation_spectrum/", stringency, "_", cell_type, ".pdf"), height = 9, width = 14)

mut_spectrum_corrected_celltype[[cell_type]] <- mut_spectrum_corrected

# ------------------------------------------------------------------------------
# Microglia
spectrum_type <- "Original"
cell_type <- "Microglia"
mut_mat <- mut_mat_summary[,str_detect(colnames(mut_mat_summary), paste0("_", cell_type, "_", spectrum_type, "$"))]
colnames(mut_mat) <- str_remove(colnames(mut_mat), paste0("_", cell_type, "_", spectrum_type, "$"))
mut_mat <- mut_mat[, case_id_list[case_id_list %in% colnames(mut_mat)]]
mut_spectrum <- matrix(rowSums(mut_mat))
rownames(mut_spectrum) <- row.names(mut_mat)
colnames(mut_spectrum) <- "Microglia_ds_mutation"
p1 <- plot_96_profile(mut_spectrum, ymax = 0.08) +
  scale_y_continuous(breaks=seq(0, 0.08, 0.02))

# ssDNA damage spectrum
ssDNA_mut_mat_summary <- readRDS("./results/ssDNA_mut/spectrum_summary.rds")
spectrum_type <- "Original"
ssDNA_mut_mat <- ssDNA_mut_mat_summary[,str_detect(colnames(ssDNA_mut_mat_summary), paste0("_", cell_type, "_", spectrum_type, "$"))]
ssDNA_mut_spectrum <- matrix(rowSums(ssDNA_mut_mat))
rownames(ssDNA_mut_spectrum) <- row.names(ssDNA_mut_mat)
colnames(ssDNA_mut_spectrum) <- "Microglia_ss_damage"
p2 <- plot_96_profile(ssDNA_mut_spectrum, ymax = 0.08) +
  scale_y_continuous(breaks=seq(0, 0.08, 0.02))

ssDNA_mut_spectrum / sum(ssDNA_mut_spectrum) -> ssDNA_mut_spectrum

sbs96 <- get_known_signatures()
fit <- fit_to_signatures(ssDNA_mut_spectrum, sbs96)
ssDNA_mut_spectrum <- #(ssDNA_mut_spectrum - fit$reconstructed) +
  rowSums(sweep(sbs96[,c("SBS2", "SBS7a", "SBS7b", "SBS11", "SBS19", "SBS23", "SBS30", "SBS31", "SBS32", "SBS40", "SBS42", "SBS84")], 2,
                fit$contribution[c("SBS2", "SBS7a", "SBS7b", "SBS11", "SBS19", "SBS23", "SBS30", "SBS31", "SBS32", "SBS40", "SBS42", "SBS84"),], "*"))
ssDNA_mut_spectrum <- matrix(ssDNA_mut_spectrum)
rownames(ssDNA_mut_spectrum) <- rownames(mut_spectrum)
colnames(ssDNA_mut_spectrum) <- "Microglia_ss_damage"

mut_spectrum_corrected <- mut_spectrum - sum(mut_spectrum[,1]) * 0.7 * ssDNA_mut_spectrum

mut_spectrum_corrected[mut_spectrum_corrected <= 0] <- 0
colnames(mut_spectrum_corrected) <- "Microglia_ds_mutation_corrected"

# normalized by 3-mer ratio
mut_spectrum_corrected <- mut_spectrum_corrected / context_3mer_ratio[,cell_type]

p3 <- plot_96_profile(mut_spectrum_corrected, ymax = 0.08) +
  scale_y_continuous(breaks=seq(0, 0.08, 0.02))
ggarrange(p1, p2, p3, ncol=1)
ggsave(paste0("./figures/mutation_spectrum/", stringency, "_", cell_type, ".pdf"), height = 9, width = 14)

mut_spectrum_corrected_celltype[[cell_type]] <- mut_spectrum_corrected


# ------------------------------------------------------------------------------
# OPC
spectrum_type <- "Original"
cell_type <- "OPC"
mut_mat <- mut_mat_summary[,str_detect(colnames(mut_mat_summary), paste0("_", cell_type, "_", spectrum_type, "$"))]
colnames(mut_mat) <- str_remove(colnames(mut_mat), paste0("_", cell_type, "_", spectrum_type, "$"))
mut_mat <- mut_mat[, case_id_list[case_id_list %in% colnames(mut_mat)]]
mut_spectrum <- matrix(rowSums(mut_mat))
rownames(mut_spectrum) <- row.names(mut_mat)
colnames(mut_spectrum) <- "OPC_ds_mutation"
p1 <- plot_96_profile(mut_spectrum, ymax = 0.08) +
  scale_y_continuous(breaks=seq(0, 0.08, 0.02))

# ssDNA damage spectrum
ssDNA_mut_mat_summary <- readRDS("./results/ssDNA_mut/spectrum_summary.rds")
spectrum_type <- "Original"
ssDNA_mut_mat <- ssDNA_mut_mat_summary[,str_detect(colnames(ssDNA_mut_mat_summary), paste0("_", cell_type, "_", spectrum_type, "$"))]
ssDNA_mut_spectrum <- matrix(rowSums(ssDNA_mut_mat))
rownames(ssDNA_mut_spectrum) <- row.names(ssDNA_mut_mat)
colnames(ssDNA_mut_spectrum) <- "OPC_ss_damage"
p2 <- plot_96_profile(ssDNA_mut_spectrum, ymax = 0.08) +
  scale_y_continuous(breaks=seq(0, 0.08, 0.02))

ssDNA_mut_spectrum / sum(ssDNA_mut_spectrum) -> ssDNA_mut_spectrum

sbs96 <- get_known_signatures()
fit <- fit_to_signatures(ssDNA_mut_spectrum, sbs96)
ssDNA_mut_spectrum <- #(ssDNA_mut_spectrum - fit$reconstructed) +
  rowSums(sweep(sbs96[,c("SBS2", "SBS7a", "SBS7b", "SBS11", "SBS19", "SBS23", "SBS30", "SBS31", "SBS32", "SBS40", "SBS42", "SBS84")], 2,
                fit$contribution[c("SBS2", "SBS7a", "SBS7b", "SBS11", "SBS19", "SBS23", "SBS30", "SBS31", "SBS32", "SBS40", "SBS42", "SBS84"),], "*"))
ssDNA_mut_spectrum <- matrix(ssDNA_mut_spectrum)
rownames(ssDNA_mut_spectrum) <- rownames(mut_spectrum)
colnames(ssDNA_mut_spectrum) <- "OPC_ss_damage"

mut_spectrum_corrected <- mut_spectrum - sum(mut_spectrum[,1]) * 0.7 * ssDNA_mut_spectrum

mut_spectrum_corrected[mut_spectrum_corrected <= 0] <- 0
colnames(mut_spectrum_corrected) <- "OPC_ds_mutation_corrected"

# normalized by 3-mer ratio
mut_spectrum_corrected <- mut_spectrum_corrected / context_3mer_ratio[,cell_type]

p3 <- plot_96_profile(mut_spectrum_corrected, ymax = 0.08) +
  scale_y_continuous(breaks=seq(0, 0.08, 0.02))
ggarrange(p1, p2, p3, ncol=1)
ggsave(paste0("./figures/mutation_spectrum/", stringency, "_", cell_type, ".pdf"), height = 9, width = 14)

mut_spectrum_corrected_celltype[[cell_type]] <- mut_spectrum_corrected

saveRDS(mut_spectrum_corrected_celltype, paste0("./results/spectrum/corrected_mutation_spectrum_", stringency, "_20240627.rds"))

# ------------------------------------------------------------------------------
# cosine similarity between spectrums
mut_spectrum_corrected_celltype <- readRDS(paste0("./results/spectrum/corrected_mutation_spectrum_", stringency, "_20240627.rds"))
data.frame(mut_spectrum_corrected_celltype) -> mut_spectrum_corrected_celltype_df
colnames(mut_spectrum_corrected_celltype_df) <- c("Excitatory Neuron", 
                                                  "Oligodendrocyte",
                                                  "Inhibitory Neuron", 
                                                  "Astrocyte",
                                                  "Microglia",
                                                  "OPC")

pdf(paste0("./figures/mutation_spectrum/", stringency, "_cos_sim_between_celltype.pdf"), width = 5, height = 4)
pheatmap(cos_sim_matrix(mut_spectrum_corrected_celltype_df,
                        mut_spectrum_corrected_celltype_df),
         scale = "none", cluster_rows = F, cluster_cols = F, 
         display_numbers = T, fontsize_number = 12, cellheight = 30, cellwidth = 30)
dev.off()

# ------------------------------------------------------------------------------
# compare with PTA spectrum
EN_PTA_snv <- read.table("./data/others/PTA/neuron_snv_pass_hg38.bed")
EN_PTA_snv_grange <- GRanges(seqnames = EN_PTA_snv$V1,
                              ranges = IRanges(start = EN_PTA_snv$V2, end = EN_PTA_snv$V3),
                              ref = EN_PTA_snv$V4, alt = EN_PTA_snv$V5)
chr_length <- seqlengths(Hsapiens)
seqlengths(EN_PTA_snv_grange) <- chr_length[names(seqlengths(EN_PTA_snv_grange))]
seqlevels(EN_PTA_snv_grange) <- seqlevels(EN_PTA_snv_grange)[order(factor(seqlevels(EN_PTA_snv_grange), levels = chr_orders))]
genome(EN_PTA_snv_grange) <- "hg38"
EN_PTA_snv_mut_mat <- mut_matrix(EN_PTA_snv_grange, ref_genome = ref_genome)
plot_96_profile(EN_PTA_snv_mut_mat)

OL_PTA_snv <- read.table("./data/others/PTA/oligo_snv_pass_hg38.bed")
OL_PTA_snv_grange <-  GRanges(seqnames = OL_PTA_snv$V1,
                              ranges = IRanges(start = OL_PTA_snv$V2, end = OL_PTA_snv$V3),
                              ref = OL_PTA_snv$V4, alt = OL_PTA_snv$V5)
chr_length <- seqlengths(Hsapiens)
seqlengths(OL_PTA_snv_grange) <- chr_length[names(seqlengths(OL_PTA_snv_grange))]
seqlevels(OL_PTA_snv_grange) <- seqlevels(OL_PTA_snv_grange)[order(factor(seqlevels(OL_PTA_snv_grange), levels = chr_orders))]
genome(OL_PTA_snv_grange) <- "hg38"
OL_PTA_snv_mut_mat <- mut_matrix(OL_PTA_snv_grange, ref_genome = ref_genome)
plot_96_profile(OL_PTA_snv_mut_mat)

#EN_PTA_snv_mut_mat <- readRDS("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/results/spectrum/PTA_spectrum_EN_a4s2.rds")
#OL_PTA_snv_mut_mat <- readRDS("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/results/spectrum/PTA_spectrum_Oligodendrocyte_a4s2.rds")

data.frame(mut_spectrum_corrected_celltype) -> mut_spectrum_corrected_celltype_df
cbind(EN_PTA_snv_mut_mat, OL_PTA_snv_mut_mat) -> PTA_snv_mut_mat
colnames(PTA_snv_mut_mat) <- c("EN_PTA", "OL_PTA")
saveRDS(PTA_snv_mut_mat, "./figures/manuscript_figures/Figure4/spectrum_PTA_OL_EN.rds")

pdf(paste0("./figures/mutation_spectrum/", stringency, "_cos_sim_PTA.pdf"), width = 6, height = 4)
pheatmap(cos_sim_matrix(mut_spectrum_corrected_celltype_df, PTA_snv_mut_mat),
         scale = "none", cluster_rows = F, cluster_cols = F, 
         display_numbers = T, fontsize_number = 12, cellheight = 30, cellwidth = 30)
dev.off()

pdf(paste0("./figures/mutation_spectrum/", stringency, "_cos_sim_PTA_remove_CtoT.pdf"), width = 6, height = 4)
pheatmap(cos_sim_matrix(mut_spectrum_corrected_celltype_df[c(1:32, 49:96),], PTA_snv_mut_mat[c(1:32, 49:96),]),
         scale = "none", cluster_rows = F, cluster_cols = F, 
         display_numbers = T, fontsize_number = 12, cellheight = 30, cellwidth = 30)
dev.off()

cos_sim_matrix(PTA_snv_mut_mat[c(1:32, 49:96),], PTA_snv_mut_mat[c(1:32, 49:96),])

# ------------------------------------------------------------------------------
# Signature contribution
mut_spectrum_corrected_celltype <- readRDS(paste0("./results/spectrum/corrected_mutation_spectrum_a6s3_20240627.rds"))
mut_spectrum_corrected_celltype <- mut_spectrum_corrected_celltype[c("Oligodendrocyte", "EN")]

sbs96 <- get_known_signatures()
sbs96_selected <- get_known_signatures()[,c("SBS5", "SBS16", "SBS1", "SBS32")]

fit <- fit_to_signatures(data.frame(mut_spectrum_corrected_celltype), sbs96_selected)
fit$contribution[,1]/colSums(fit$contribution)[1]
fit$contribution[,2]/colSums(fit$contribution)[2]
plot_contribution(fit$contribution,
                  coord_flip = FALSE,
                  mode = "relative")
ggsave(paste0("./figures/mutation_spectrum/", stringency, "/signature/fitting_selected_signature.pdf"), width = 6, height = 6)

mut_spectrum_corrected_celltype[["Oligodendrocyte_OPC"]] <- mut_spectrum_corrected_celltype$Oligodendrocyte + mut_spectrum_corrected_celltype$OPC
fit <- fit_to_signatures(data.frame(mut_spectrum_corrected_celltype), sbs96)
fit$contribution["SBS1",]/colSums(fit$contribution)

plot_contribution(fit$contribution,
                  coord_flip = FALSE,
                  mode = "relative")
ggsave(paste0("./figures/mutation_spectrum/", stringency, "/signature/fitting_all_signature.pdf"), width = 6, height = 6)

fit <- fit_to_signatures_strict(data.frame(mut_spectrum_corrected_celltype), sbs96)
plot_contribution(fit$fit_res$contribution,
                  coord_flip = FALSE,
                  mode = "relative")
ggsave(paste0("./figures/mutation_spectrum/", stringency, "/signature/fitting_all_signature_strict.pdf"), width = 6, height = 6)


sbs16 <- sbs96[,c("SBS16"), drop = F]
for (mut_spectrum_corrected in mut_spectrum_corrected_celltype){
  fit <- fit_to_signatures(mut_spectrum_corrected, sbs16)
  print(fit$contribution/sum(mut_spectrum_corrected))
}

fit <- fit_to_signatures(PTA_snv_mut_mat, sbs16)
print(fit$contribution/colSums(PTA_snv_mut_mat))

sbs16 <- sbs96[,c("SBS16"), drop = F]
for (mut_spectrum_corrected in mut_spectrum_corrected_celltype){
  fit <- fit_to_signatures(mut_spectrum_corrected[65:80, , drop=F], sbs16[65:80, , drop=F])
  print(fit$contribution/sum(mut_spectrum_corrected))
}
