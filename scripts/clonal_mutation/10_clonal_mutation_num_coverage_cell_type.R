.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(reshape2)
library(rtracklayer)
library(patchwork)
library(GenomicRanges)
library(ggpubr)

neuronal <- c("EN", "IN")
glial <- c("Oligodendrocyte", "Microglia", "Astrocyte", "OPC")
#OL <- c("Oligodendrocyte", "OPC")

case_id_list <- c("Brain1901106_deeper", "UMB1278_deeper", "UMB1864_deeper", 
                  "UMB5823_deeper", "UMB4638", "UMB1465", "UMB5451", "UMB5657", "UMB4428")

metadata <- read.csv("./data/others/metadata_07012024.csv")

# clonal mutation
clonal_mutation <- read.csv("./results/clonal_mutation/clonal_mutation_list_healthy_brain_realigned_ATAC_celltype_20241028.csv")
clonal_mutation$vaf <- clonal_mutation$detected/clonal_mutation$covered
#clonal_mutation <- clonal_mutation[clonal_mutation$vaf<0.01,]
# validated_mutation <- read.csv("./data/amplicon-seq/healthy_brain_ampliconseq_valideated_01142025.csv", header = T)
# validated_mutation[,"key"] <- paste(validated_mutation$Chrom,
#                                     validated_mutation$Pos,
#                                     validated_mutation$Ref,
#                                     validated_mutation$Alt,
#                                     sep = "-")
# clonal_mutation <- clonal_mutation[clonal_mutation$key %in% validated_mutation$key,]
# clonal_mutation <- na.omit(clonal_mutation)
clonal_mutation$cell_type[is.na(clonal_mutation$cell_type)] <-
  clonal_mutation$celltype_ATAC[is.na(clonal_mutation$cell_type)]
clonal_mutation <- clonal_mutation[!is.na(clonal_mutation$cell_type),]
# clonal_mutation <- clonal_mutation[!is.na(clonal_mutation$cell_type),]
clonal_mutation$cell_type[str_detect(clonal_mutation$cell_type, "EN")] <- "EN"
clonal_mutation$cell_type[str_detect(clonal_mutation$cell_type, "IN")] <- "IN"
clonal_mutation$cell_type[str_detect(clonal_mutation$cell_type, "Astrocyte")] <- "Astrocyte"
clonal_mutation$cell_type[str_detect(clonal_mutation$cell_type, "Microglia")] <- "Microglia"

clonal_mutation$case_id <- factor(clonal_mutation$case_id, levels = case_id_list)
clonal_mutation_count <- 
  #clonal_mutation %>% group_by(cell_type) %>% summarise("count" = n())
  clonal_mutation %>% group_by(case_id, cell_type) %>% summarise("count" = n())
  #clonal_mutation %>% group_by(key, cell_type) %>% summarise() %>% group_by(cell_type, .drop=F) %>% summarise("count" = n())

clonal_mutation_covered_lower_bound <- 
  clonal_mutation %>% group_by(case_id) %>% summarise("min_covered" = min(covered))
#min(clonal_mutation$covered)
clonal_mutation_max_vaf <- clonal_mutation %>% group_by(key, case_id) %>% summarise("vaf" = vaf) %>% group_by(case_id) %>% summarise("max_vaf" = max(vaf))

# summary the size of covered regions that are covered >= 42 (clonal_mutation_covered_lower_bound) cells
genomecov_clone <- data.frame()
for (case_id in case_id_list){
  print(case_id)
  file_name <- paste0("./results/clonal_mutation/multicell_covered_region/", case_id, "_open_region_a2s0_coverage_count.bed")
  if(!file.exists(file_name)){next}
  genomecov_case <- import(file_name, format = "bed")
  
  if(sum(clonal_mutation_covered_lower_bound$case_id == case_id) == 0){
    clonal_mutation_covered_lower_bound_sample <- floor(mean(clonal_mutation_covered_lower_bound$min_covered))
  }else{
    clonal_mutation_covered_lower_bound_sample <- clonal_mutation_covered_lower_bound$min_covered[clonal_mutation_covered_lower_bound$case_id == case_id]
  }
  genomecov_case <- genomecov_case[genomecov_case$name > clonal_mutation_covered_lower_bound_sample,]
  
  for (cell_type in unique(clonal_mutation_count$cell_type)){
    gc()
    print(cell_type)
    file_name <- paste0("./results/clonal_mutation/multicell_covered_region_celltype_ATAC_recovered/", case_id, "_", cell_type, "_open_region_a2s0_coverage_count.bed")
    if(!file.exists(file_name)){next}
    genomecov <- import(file_name, format = "bed")
    intersection <- findOverlaps(genomecov, genomecov_case)
    intersected_ranges <- genomecov[queryHits(intersection)]
    
    genomecov_clone[case_id, cell_type] <- sum(intersected_ranges@ranges@width * as.integer(intersected_ranges$name))
  }
}
#save.image("./results/clonal_mutation/clonal_mutation_rate_celltype_ATAC_recovered_01242025.RData")
genomecov_clone[is.na(genomecov_clone)] <- 0
genomecov_clone$case_id <- row.names(genomecov_clone)
genomecov_clone$case_id <- factor(genomecov_clone$case_id, levels = case_id_list)
genomecov_clone <- melt(genomecov_clone, variable.name = "cell_type")
genomecov_clone %>% group_by(case_id, cell_type) %>% summarise("genomecov_clone" = sum(value)) -> genomecov_clone
genomecov_clone$cell_type <- as.character(genomecov_clone$cell_type)
# genomecov_clone$cell_type <- ifelse(genomecov_clone$cell_type %in% neuronal,
#                                     "Neuron", "Non-neuronal")
genomecov_clone$cell_type[genomecov_clone$cell_type %in% neuronal] <- "neuronal"
genomecov_clone$cell_type[genomecov_clone$cell_type %in% glial] <- "glial"
genomecov_clone <- genomecov_clone[genomecov_clone$cell_type %in% c("neuronal", "glial"),]
genomecov_clone$cell_type <- as.factor(genomecov_clone$cell_type)
genomecov_clone %>% group_by(case_id, cell_type) %>% summarise("genomecov_clone" = sum(genomecov_clone)) -> genomecov_clone

clonal_mutation_count$case_id <- factor(clonal_mutation_count$case_id, levels = case_id_list)
clonal_mutation_count %>% group_by(case_id, cell_type, .drop=F) %>% summarise("count" = sum(count)) -> clonal_mutation_count
clonal_mutation_count$cell_type <- as.character(clonal_mutation_count$cell_type)
# clonal_mutation_count$cell_type <- ifelse(clonal_mutation_count$cell_type %in% neuronal,
#                                           "Neuron", "Non-neuronal")
clonal_mutation_count$cell_type[clonal_mutation_count$cell_type %in% neuronal] <- "neuronal"
clonal_mutation_count$cell_type[clonal_mutation_count$cell_type %in% glial] <- "glial"
clonal_mutation_count <- clonal_mutation_count[clonal_mutation_count$cell_type %in% c("neuronal", "glial"),]
clonal_mutation_count$cell_type <- as.factor(clonal_mutation_count$cell_type)
clonal_mutation_count %>% group_by(case_id, cell_type) %>% summarise("count" = sum(count)) -> clonal_mutation_count

clonal_mutation_count_coverage <- inner_join(clonal_mutation_count, genomecov_clone, by = c("case_id", "cell_type"))
clonal_mutation_count_coverage$rate <- clonal_mutation_count_coverage$count/clonal_mutation_count_coverage$genomecov_clone

ggplot(clonal_mutation_count_coverage, aes(x = cell_type, y = rate)) +
  geom_col() +
  labs(x = "Cell type", y = "Detected clonal sSNVs rate per bp") +
  #scale_y_continuous(labels = scientific_10) +
  theme_classic()
#ggsave("./figures/burden_simulation/clonal_mutation_cell_type.pdf", width = 3)

#clonal_mutation_count_coverage <- clonal_mutation_count_coverage[clonal_mutation_count_coverage$case_id != "UMB5823_deeper",]
clonal_mutation_count_coverage %>% 
  group_by(cell_type) %>%
  summarise("mean_rate" = mean(rate),
            "sd_rate" = sd(rate)) -> clonal_mutation_count_coverage_rate


wilcox.test(x = clonal_mutation_count_coverage$rate[clonal_mutation_count_coverage$cell_type == "neuronal"],
            y = clonal_mutation_count_coverage$rate[clonal_mutation_count_coverage$cell_type == "glial"], 
       paired = T, alternative = "less") -> wilcox_result
wilcox_result

wilcox.test(x = c(4.010195e-10, 1.400655e-10, 2.511217e-10, 0.000000e+00, 0.000000e+00, 9.601248e-11, 2.889015e-10, 0.000000e+00),
      y = c(1.314647e-10, 7.502546e-11, 2.087542e-10, 0.000000e+00, 0.000000e+00, 5.298137e-11, 0.000000e+00, 0.000000e+00), 
      paired = T, alternative = "less")

clonal_mutation_count_coverage$rate[clonal_mutation_count_coverage$cell_type == "glial"]/
  clonal_mutation_count_coverage$rate[clonal_mutation_count_coverage$cell_type == "neuronal"]
inner_join(metadata, clonal_mutation_count_coverage, by = c("Sample" = "case_id")) -> metadata
regression_results <- lm(rate ~ Sample + cell_type, data = metadata)
summary(regression_results)
summary(lm(rate ~ cell_type + Age + Sex, data = metadata))
summary(lm(count ~ genomecov_clone + Age, data = metadata))
plot(metadata$genomecov_clone, metadata$count)

summary(lm(rate ~ Age, data = metadata[metadata$cell_type == "glial",]))
plot(metadata[metadata$cell_type == "glial", "Age"],
     metadata[metadata$cell_type == "glial", "rate"])

ggplot(metadata, aes(x=genomecov_clone, y=count, shape=cell_type, color=Sample))+
  geom_point(size=3) +
  labs(x = "Callable region size", y = "Number of clonal sSNVs detected") +
  theme_classic() 
#ggsave("./figures/burden_simulation/clonal_mutation_neuronal_glial.pdf")

for (sample in unique(metadata$Sample)){
  print(sample)
  print(metadata$rate[metadata$Sample==sample & metadata$cell_type=="glial"]/
          metadata$rate[metadata$Sample==sample & metadata$cell_type=="neuronal"])
}

ggplot(clonal_mutation_count_coverage, aes(x = cell_type, y = rate, color = cell_type)) +
  # geom_col() +
  # geom_errorbar(aes(ymin = mean_rate - sd_rate,
  #                   ymax = mean_rate + sd_rate)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  labs(x = "Callable region size", y = "Number of clonal sSNVs detected") +
  theme_classic()
  
ggplot(metadata, aes(x = cell_type, y = rate, fill = Sample)) +
    geom_col(position = "dodge") +
    labs(x = "Callable region size", y = "Number of clonal sSNVs detected") +
    theme_classic()

metadata$cell_type <- factor(metadata$cell_type, levels = c("neuronal", "glial"))
ggplot(metadata, aes(x = cell_type, y = rate, group = Sample)) +
  geom_point() +
  geom_line() +
  labs(x = "Callable region size", y = "Number of clonal sSNVs detected") +
  theme_classic() 
  

saveRDS(metadata, "./figures/manuscript_figures/Figure5/clonal_mutation_count_coverage_cell_type.rds")

metadata$Sample <- as.factor(metadata$Sample)
metadata %>% group_by(Sample, Age, Sex) %>% summarise("count" = sum(count),
                                            "genomecov_clone" = sum(genomecov_clone)) -> metadata_sample
plot(metadata_sample$genomecov_clone, metadata_sample$count)
summary(lm(count ~ genomecov_clone, data = metadata_sample))
saveRDS(metadata_sample, "./figures/manuscript_figures/Figure5/clonal_mutation_count_coverage.rds")
