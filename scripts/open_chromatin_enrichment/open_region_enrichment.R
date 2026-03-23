# Libraries
.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(stringr)

# EN
# metadata
pta_metadata <- read.csv("./data/others/PTA/PTA_neuron_burden.csv")
duplex_metadata <- read.csv("./data/others/metadata_07012024.csv")

enrichment <- c()
for (sample in duplex_metadata$Sample){
  print(sample)
  if(!file.exists(paste0("./data/", sample, "/cell_type/subtype/EN.txt"))){
    next
  }
  
  for (a in seq(2, 10, 2)){
    stringency <- paste0("a", a, "s", a%/%2)
    print(stringency)
    # read in callable regions
    duplex_callable_region <- read.table(paste0("./data/", sample, "/PTA_results/EN_open_region_", stringency, "_sorted_merged.bed"))
    duplex_callable_region_grange <- GRanges(seqnames = duplex_callable_region$V1, 
                                             ranges = IRanges(start = duplex_callable_region$V2, 
                                                              end = duplex_callable_region$V3))
    duplex_callable_region_size <- sum(duplex_callable_region$V3 - duplex_callable_region$V2)
    
    # select age matched PTA cells
    duplex_sample_age <- duplex_metadata$Age[duplex_metadata$Sample == sample]
    selected_pta_cell <- pta_metadata$Cell_ID[pta_metadata$Age <= duplex_sample_age + 15 | pta_metadata$Age >= duplex_sample_age - 15]
    
    # read in PTA calls
    for (cell_id in selected_pta_cell){
      pta_snv <- read.table(paste0("./data/others/PTA/EN_snv/", cell_id, "-snv-hg38.bed"), sep = "\t")
      pta_snv_grange <- GRanges(seqnames = pta_snv$V1, 
                                ranges = IRanges(start = pta_snv$V2, 
                                                 end = pta_snv$V3))
      
      # extract overlaps
      intersection <- findOverlaps(pta_snv_grange, duplex_callable_region_grange)
      intersected_ranges <- pta_snv_grange[queryHits(intersection)]
      pta_mutation_num <- length(intersected_ranges)
      
      enrichment_ratio <- (pta_mutation_num/duplex_callable_region_size)/(dim(pta_snv)[1]/(5.75e9/2))
      
      enrichment <- rbind(enrichment, c(sample, cell_id, stringency, enrichment_ratio, pta_mutation_num, dim(pta_snv)[1], duplex_callable_region_size))
    }
  }
}

enrichment <- data.frame(enrichment)
colnames(enrichment) <- c("sample", "cell", "stringency", "enrichment_ratio", "overlapping_pta_mutation_num", "all_pta_mutation_num", "duplex_callable_region_size")
enrichment$overlapping_pta_mutation_num <- as.integer(enrichment$overlapping_pta_mutation_num)
enrichment$all_pta_mutation_num <- as.integer(enrichment$all_pta_mutation_num)
enrichment$duplex_callable_region_size <- as.integer(enrichment$duplex_callable_region_size)
enrichment_summary <- enrichment %>% 
  group_by(sample, stringency, duplex_callable_region_size) %>% 
  summarise("overlapping_pta_mutation_sum" = sum(overlapping_pta_mutation_num),
            "all_pta_mutation_sum" = sum(all_pta_mutation_num))
enrichment_summary$enrichment_ratio <- (enrichment_summary$overlapping_pta_mutation_sum/enrichment_summary$duplex_callable_region_size)/(enrichment_summary$all_pta_mutation_sum/2.875e9)

enrichment_summary$stringency <- factor(enrichment_summary$stringency,
                                        levels = c("a2s1", "a4s2", "a6s3", "a8s4", "a10s5"))
ggplot(enrichment_summary, aes(x = stringency, y = enrichment_ratio, fill = sample)) +
  geom_col(position = "dodge2") +
  theme_classic() +
  labs(x = "Strigency", y = "Mutation enrichment ratio \n (callable region/whole genome)", title = "EN")
ggsave("./figures/open_region_enrichment/EN.pdf", width = 8, height = 5)

saveRDS(enrichment_summary, "./results/open_region_enrichment/EN.rds")
write.csv(enrichment_summary, "./results/open_region_enrichment/EN.csv", quote = F, row.names = F)


# Oligodendrocyte
# metadata
pta_metadata <- read.csv("./data/others/PTA/PTA_oligo_burden.csv")
duplex_metadata <- read.csv("./data/others/metadata_07012024.csv")

enrichment <- c()
for (sample in duplex_metadata$Sample){
  print(sample)
  if(sample == "UMB4428"){next}
  for (a in seq(2, 10, 2)){
    stringency <- paste0("a", a, "s", a%/%2)
    print(stringency)
    # read in callable regions
    duplex_callable_region <- read.table(paste0("./data/", sample, "/PTA_results/Oligodendrocyte_open_region_", stringency, "_sorted_merged.bed"))
    duplex_callable_region_grange <- GRanges(seqnames = duplex_callable_region$V1, 
                                             ranges = IRanges(start = duplex_callable_region$V2, 
                                                              end = duplex_callable_region$V3))
    duplex_callable_region_size <- sum(duplex_callable_region$V3 - duplex_callable_region$V2)
    
    # select age matched PTA cells
    duplex_sample_age <- duplex_metadata$Age[duplex_metadata$Sample == sample]
    selected_pta_cell <- pta_metadata$Cell_ID[pta_metadata$Age <= duplex_sample_age + 15 | pta_metadata$Age >= duplex_sample_age - 15]
    
    # read in PTA calls
    for (cell_id in selected_pta_cell){
      pta_snv <- read.table(paste0("./data/others/PTA/Oligodendrocyte_snv/", cell_id, "-snv-hg38.bed"), sep = "\t")
      pta_snv_grange <- GRanges(seqnames = pta_snv$V1, 
                                ranges = IRanges(start = pta_snv$V2, 
                                                 end = pta_snv$V3))
      
      # extract overlaps
      intersection <- findOverlaps(pta_snv_grange, duplex_callable_region_grange)
      intersected_ranges <- pta_snv_grange[queryHits(intersection)]
      pta_mutation_num <- length(intersected_ranges)
      
      enrichment_ratio <- (pta_mutation_num/duplex_callable_region_size)/(dim(pta_snv)[1]/(5.75e9/2))
      
      enrichment <- rbind(enrichment, c(sample, cell_id, stringency, enrichment_ratio, pta_mutation_num, dim(pta_snv)[1], duplex_callable_region_size))
    }
  }
}

enrichment <- data.frame(enrichment)
colnames(enrichment) <- c("sample", "cell", "stringency", "enrichment_ratio", "overlapping_pta_mutation_num", "all_pta_mutation_num", "duplex_callable_region_size")
enrichment$overlapping_pta_mutation_num <- as.integer(enrichment$overlapping_pta_mutation_num)
enrichment$all_pta_mutation_num <- as.integer(enrichment$all_pta_mutation_num)
enrichment$duplex_callable_region_size <- as.integer(enrichment$duplex_callable_region_size)
enrichment_summary <- enrichment %>% 
  group_by(sample, stringency, duplex_callable_region_size) %>% 
  summarise("overlapping_pta_mutation_sum" = sum(overlapping_pta_mutation_num),
            "all_pta_mutation_sum" = sum(all_pta_mutation_num))
enrichment_summary$enrichment_ratio <- (enrichment_summary$overlapping_pta_mutation_sum/enrichment_summary$duplex_callable_region_size)/(enrichment_summary$all_pta_mutation_sum/2.875e9)

enrichment_summary$stringency <- factor(enrichment_summary$stringency,
                                        levels = c("a2s1", "a4s2", "a6s3", "a8s4", "a10s5"))
ggplot(enrichment_summary, aes(x = stringency, y = enrichment_ratio, fill = sample)) +
  geom_col(position = "dodge2") +
  theme_classic() +
  labs(x = "Strigency", y = "Mutation enrichment ratio \n (callable region/whole genome)", title = "Oligodendrocyte")
ggsave("./figures/open_region_enrichment/Oligodendrocyte.pdf", width = 8, height = 5)

saveRDS(enrichment_summary, "./results/open_region_enrichment/Oligodendrocyte.rds")
write.csv(enrichment_summary, "./results/open_region_enrichment/Oligodendrocyte.csv", quote = F, row.names = F)
