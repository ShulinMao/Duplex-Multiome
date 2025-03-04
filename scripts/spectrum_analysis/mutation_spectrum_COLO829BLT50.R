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

case_id_list <- c("COLO829BLT50_rep1", "COLO829BLT50_rep2")

a <- 6
s <- a%/%2
stringency <- paste0("a", a, "s", s)

# ------------------------------------------------------------------------------
# collect spectrum from all cells
mut_mat_summary <- c()
for(case_id in case_id_list){
  # cell type
  file_list <- dir(paste0("./data/", case_id, "/cell_type/"))
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

saveRDS(mut_mat_summary, paste0("./results/spectrum/COLO829BLT50_mutation_spectrum", "_", stringency, ".rds"))

# ------------------------------------------------------------------------------
# spectrum analysis
mut_mat_summary <- readRDS(paste0("./results/spectrum/COLO829BLT50_mutation_spectrum", "_", stringency, ".rds"))

mut_spectrum_corrected_celltype <- list()

# 3 mer count in covered region / 3 mer count in the whole genome
context_3mer_ratio <- readRDS(paste0("./results/context_3mer_ratio/COLO829BLT50_", "a", a, "s", s, "_celltype.rds"))
context_3mer_ratio <- context_3mer_ratio[c(1:16, 1:16, 1:16, 17:32, 17:32, 17:32),]

# ------------------------------------------------------------------------------
# remove ss damages
# melanoma_fibroblast
spectrum_type <- "Original"
cell_type <- "melanoma_fibroblast"
mut_mat <- mut_mat_summary[,str_detect(colnames(mut_mat_summary), paste0("_", cell_type, "_", spectrum_type, "$")), drop=F]
colnames(mut_mat) <- str_remove(colnames(mut_mat), paste0("_", cell_type, "_", spectrum_type, "$"))
mut_mat <- mut_mat[, case_id_list[case_id_list %in% colnames(mut_mat)]]
mut_spectrum <- matrix(rowSums(mut_mat))
rownames(mut_spectrum) <- row.names(mut_mat)
colnames(mut_spectrum) <- "melanoma_fibroblast_ds_mutation"
p1 <- plot_96_profile(mut_spectrum, ymax = 0.2)

# ssDNA damage spectrum
ssDNA_mut_mat_summary <- readRDS("./results/ssDNA_mut/COLO829BLT50_spectrum_summary.rds")
spectrum_type <- "Original"
ssDNA_mut_mat <- ssDNA_mut_mat_summary[,str_detect(colnames(ssDNA_mut_mat_summary), paste0("_", cell_type, "_", spectrum_type, "$")), drop=F]
ssDNA_mut_spectrum <- matrix(rowSums(ssDNA_mut_mat))
rownames(ssDNA_mut_spectrum) <- row.names(ssDNA_mut_mat)
colnames(ssDNA_mut_spectrum) <- "melanoma_fibroblast_ss_damage"
p2 <- plot_96_profile(ssDNA_mut_spectrum, ymax = 0.2)

ssDNA_mut_spectrum / sum(ssDNA_mut_spectrum) -> ssDNA_mut_spectrum

sbs96 <- get_known_signatures()
fit <- fit_to_signatures(ssDNA_mut_spectrum, sbs96)
ssDNA_mut_spectrum <- #(ssDNA_mut_spectrum - fit$reconstructed) +
  rowSums(sweep(sbs96[,c("SBS2", "SBS7a", "SBS7b", "SBS11", "SBS19", "SBS23", "SBS30", "SBS32", "SBS40", "SBS42", "SBS84", "SBS44")], 2,
                fit$contribution[c("SBS2", "SBS7a", "SBS7b", "SBS11", "SBS19", "SBS23", "SBS30", "SBS32", "SBS40", "SBS42", "SBS84", "SBS44"),], "*"))
ssDNA_mut_spectrum <- matrix(ssDNA_mut_spectrum)
rownames(ssDNA_mut_spectrum) <- rownames(mut_spectrum)
colnames(ssDNA_mut_spectrum) <- "melanoma_fibroblast_ss_damage"
plot_96_profile(ssDNA_mut_spectrum, ymax = 0.08) + 
  scale_y_continuous(breaks=seq(0, 0.08, 0.02))

mut_spectrum_corrected <- mut_spectrum - sum(mut_spectrum[,1]) * 0.5 * ssDNA_mut_spectrum[,1] # 0.52

mut_spectrum_corrected[mut_spectrum_corrected <= 0] <- 0
colnames(mut_spectrum_corrected) <- "melanoma_fibroblast_ds_mutation_corrected"

# normalized by 3-mer ratio
mut_spectrum_corrected <- mut_spectrum_corrected / context_3mer_ratio[,cell_type]

p3 <- plot_96_profile(mut_spectrum_corrected, ymax = 0.2)

ggarrange(p1, p2, p3, ncol=1)
ggsave(paste0("./figures/mutation_spectrum/COLO829BLT50_", stringency, "_", cell_type, ".pdf"), height = 9, width = 14)

mut_spectrum_corrected_celltype[[cell_type]] <- mut_spectrum_corrected

# ------------------------------------------------------------------------------
# B_lymphoblast
spectrum_type <- "Original"
cell_type <- "B_lymphoblast"
mut_mat <- mut_mat_summary[,str_detect(colnames(mut_mat_summary), paste0("_", cell_type, "_", spectrum_type, "$")), drop=F]
colnames(mut_mat) <- str_remove(colnames(mut_mat), paste0("_", cell_type, "_", spectrum_type, "$"))
mut_mat <- mut_mat[, case_id_list[case_id_list %in% colnames(mut_mat)]]
mut_spectrum <- matrix(rowSums(mut_mat))
rownames(mut_spectrum) <- row.names(mut_mat)
colnames(mut_spectrum) <- "B_lymphoblast_ds_mutation"
p1 <- plot_96_profile(mut_spectrum, ymax = 0.08) + 
  scale_y_continuous(breaks=seq(0, 0.08, 0.02))

# ssDNA damage spectrum
ssDNA_mut_mat_summary <- readRDS("./results/ssDNA_mut/COLO829BLT50_spectrum_summary.rds")
spectrum_type <- "Original"
ssDNA_mut_mat <- ssDNA_mut_mat_summary[,str_detect(colnames(ssDNA_mut_mat_summary), paste0("_", cell_type, "_", spectrum_type, "$")), drop=F]
ssDNA_mut_spectrum <- matrix(rowSums(ssDNA_mut_mat))
rownames(ssDNA_mut_spectrum) <- row.names(ssDNA_mut_mat)
colnames(ssDNA_mut_spectrum) <- "B_lymphoblast_ss_damage"
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
colnames(ssDNA_mut_spectrum) <- "B_lymphoblast_ss_damage"
#plot_96_profile(ssDNA_mut_spectrum, ymax = 0.08) + 
#  scale_y_continuous(breaks=seq(0, 0.08, 0.02))

mut_spectrum_corrected <- mut_spectrum - sum(mut_spectrum[,1]) * 0.7 * ssDNA_mut_spectrum[,1]

mut_spectrum_corrected[mut_spectrum_corrected <= 0] <- 0
colnames(mut_spectrum_corrected) <- "B_lymphoblast_ds_mutation_corrected"

# normalized by 3-mer ratio
mut_spectrum_corrected <- mut_spectrum_corrected / context_3mer_ratio[,cell_type]

p3 <- plot_96_profile(mut_spectrum_corrected, ymax = 0.08) + 
  scale_y_continuous(breaks=seq(0, 0.08, 0.02))

ggarrange(p1, p2, p3, ncol=1)
ggsave(paste0("./figures/mutation_spectrum/COLO829BLT50_", stringency, "_", cell_type, ".pdf"), height = 9, width = 14)

mut_spectrum_corrected_celltype[[cell_type]] <- mut_spectrum_corrected

saveRDS(mut_spectrum_corrected_celltype, "./figures/manuscript_figures/Figure2/spectrum_private_a6s3_COLO829BLT50.rds")
# ------------------------------------------------------------------------------
# Signature contribution
sbs96_selected <- get_known_signatures()
as.data.frame(mut_spectrum_corrected_celltype) -> mut_spectrum_corrected_celltype_df
colnames(mut_spectrum_corrected_celltype_df) <- c("Melanoma", "B-cell")
fit <- fit_to_signatures_strict(mut_spectrum_corrected_celltype_df, sbs96_selected)
plot_contribution(fit$fit_res$contribution,
                  coord_flip = FALSE,
                  mode = "relative") +
  labs(title = "Vairants surviving BLT-50 bulk WGS filter")
ggsave(paste0("figures/mutation_spectrum/COLO829BLT50_signature_analysis_", stringency, "_06112024.pdf"), width = 5, height = 5)

