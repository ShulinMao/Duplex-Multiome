# Libraries
.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(dplyr)
ref_genome="BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = T)
chr_orders <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
library(GenomicRanges)
library(ggplot2)
library(stringr)

.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages_sig/")
library(MutationalPatterns)

cell_type <- "EN"

# mutation
mutation_grange <- readRDS(paste0("./data/others/PTA/", cell_type, "_PTA_mutation_hg37.rds"))

# duplex-multiome covered regions
duplex_covered_region <- read.table(paste0("./data/others/covered_region_summary/", cell_type, "_open_region_a6s3_combined_merged_hg19.bed"), header = T)
duplex_covered_region_grange <- GRanges(seqnames = duplex_covered_region$seqnames,
                                        ranges = IRanges(start = duplex_covered_region$start,
                                                         end = duplex_covered_region$end))
rm(duplex_covered_region)

# mutations in or not in duplex-multiome covered regions
intersection <- findOverlaps(mutation_grange, duplex_covered_region_grange)
mutation_duplex_covered_region <- mutation_grange[queryHits(intersection)]
mutation_not_duplex_covered_region <- subsetByOverlaps(mutation_grange, mutation_duplex_covered_region, invert = T)

# settings for MutationPatterns
chr_length <- seqlengths(Hsapiens)
seqlengths(mutation_duplex_covered_region) <- chr_length[names(seqlengths(mutation_duplex_covered_region))]
seqlevels(mutation_duplex_covered_region) <- seqlevels(mutation_duplex_covered_region)[order(factor(seqlevels(mutation_duplex_covered_region), levels = chr_orders))]
genome(mutation_duplex_covered_region) = 'hg19'
seqlengths(mutation_not_duplex_covered_region) <- chr_length[names(seqlengths(mutation_not_duplex_covered_region))]
seqlevels(mutation_not_duplex_covered_region) <- seqlevels(mutation_not_duplex_covered_region)[order(factor(seqlevels(mutation_not_duplex_covered_region), levels = chr_orders))]
genome(mutation_not_duplex_covered_region) = 'hg19'
seqlengths(mutation_grange) <- chr_length[names(seqlengths(mutation_grange))]
seqlevels(mutation_grange) <- seqlevels(mutation_grange)[order(factor(seqlevels(mutation_grange), levels = chr_orders))]
genome(mutation_grange) = 'hg19'

# calculate 96 type of snvs
mutation_duplex_covered_region_mut_mat <- mut_matrix(mutation_duplex_covered_region, ref_genome = ref_genome)
mutation_not_duplex_covered_region_mut_mat <- mut_matrix(mutation_not_duplex_covered_region, ref_genome = ref_genome)
mutation_mut_mat <- mut_matrix(mutation_grange, ref_genome = ref_genome)

# 3 mer count in covered region / 3 mer count in the whole genome
context_3mer_ratio <- readRDS(paste0("./results/context_3mer_ratio/a6s3_celltype.rds"))
context_3mer_ratio <- context_3mer_ratio[c(1:16, 1:16, 1:16, 17:32, 17:32, 17:32),]

# normalized by 3-mer ratio
mutation_duplex_covered_region_mut_mat_corrected <- mutation_duplex_covered_region_mut_mat / context_3mer_ratio[, cell_type]

sbs96 <- get_known_signatures()[, c("SBS1", "SBS5", "SBS16", "SBS19", "SBS32") , drop = F]
#sbs96 <- get_known_signatures()[, c("SBS19") , drop = F]
mutation_duplex_covered_region_fit <- fit_to_signatures(mutation_duplex_covered_region_mut_mat, sbs96)
mutation_duplex_covered_region_corrected_fit <- fit_to_signatures(mutation_duplex_covered_region_mut_mat_corrected, sbs96)
mutation_not_duplex_covered_region_fit <- fit_to_signatures(mutation_not_duplex_covered_region_mut_mat, sbs96)
mutation_fit <- fit_to_signatures(mutation_mut_mat, sbs96)

selected_signature <- "SBS1"
mutation_duplex_covered_region_fit$contribution[selected_signature,]/sum(mutation_duplex_covered_region_mut_mat)
mutation_duplex_covered_region_corrected_fit$contribution[selected_signature,]/sum(mutation_duplex_covered_region_mut_mat_corrected)
mutation_fit$contribution[selected_signature,]/sum(mutation_mut_mat)

SBS_signature_enrichment <- data.frame("group" = c("Duplex-Multiome covered regions\n(Normalized)", "Genome"),
                                       "enrichment" = c((mutation_duplex_covered_region_corrected_fit$contribution[selected_signature,]/sum(mutation_duplex_covered_region_mut_mat_corrected))/
                                                          (mutation_fit$contribution[selected_signature,]/sum(mutation_mut_mat)), 1))

ggplot(SBS_signature_enrichment, aes(x = group, y = enrichment)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 0.5) +
  theme_classic() +
  geom_hline(yintercept=1, linetype="dashed") +
  labs(title = selected_signature,
       subtitle = "(EN, PTA sSNVs, Duplex a4s2 regions)",
       x = "Group",
       y = "Relative SBS % contribution") 

#saveRDS(list("covered_region" = mutation_duplex_covered_region_mut_mat,
#             "all" = mutation_mut_mat),
#        paste0("./figures/manuscript_figures/Figure4/signature_enrichment_analysis_", cell_type, "_a6s3.rds"))



