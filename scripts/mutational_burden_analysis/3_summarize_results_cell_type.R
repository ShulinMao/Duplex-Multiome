.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/")

library(stringr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)
library(truncnorm)
set.seed(42)

library(httpgd)

source("./scripts/simulation_burdens/simulation_regression_func.R")


case_id_list <- #c("ABNCR6D", "AN06365")
  c("Brain1901106_deeper", "UMB1278_deeper", "UMB1864_deeper", "UMB5823_deeper", "UMB4638", "UMB1465", "UMB5451", "UMB5657")
  #c("Brain1901106", "UMB1278", "UMB1864", "UMB4638", "UMB1465", "UMB5451", "UMB5657", "UMB5823")

metadata <- read.csv("./data/others/metadata_07012024.csv")

# ------------------------------------------------------------------------------
# summarize simulation results
simulation_results_sample <- c()

for(case_id in case_id_list){
  print(case_id)
  
  # cell type
  file_list <- dir(paste0("./data/", case_id, "/cell_type/"), pattern = ".*txt")
  cell_type_list <- str_remove(file_list, "[.]txt")
  
  for (cell_type in cell_type_list){
    print(cell_type)
    
    if(!file.exists(paste0("./results/burden_simulation/", case_id, "_", cell_type, ".rds"))){next}

    # read in simulation results
    simulation_results <- readRDS(paste0("./results/burden_simulation/", case_id, "_", cell_type, ".rds"))
    
    simulation_results_summary <- c()
    for (a in seq(6, 10, 2)){
      s <- a%/%2
      stringency <- paste0("a", a, "s", s)
      # select the results with the lowest error
      simulation_results_summary <-
        rbind(simulation_results_summary, 
            c(stringency, simulation_results[[stringency]][which.max(simulation_results[[stringency]][,6]),]))
            #c(stringency, simulation_results[[stringency]][which.min(simulation_results[[stringency]][,7]),]))
    }
    colnames(simulation_results_summary) <- c("stringency", "count", "sd", "burden", "raw_coverage", "coverage", "likelihood", "error")
    simulation_results_summary <- data.frame(simulation_results_summary)
    simulation_results_summary["cell_type"] <- cell_type
    simulation_results_summary["case_id"] <- case_id
    simulation_results_summary$coverage <- as.integer(simulation_results_summary$coverage)
    
    # select the fittest situations
    # a6s3 or more stringent
    # covered regions > 20000000 bps
    covered_region_size <- 20000000
    if(sum(simulation_results_summary$coverage >= covered_region_size) <= 0){
      #simulation_results_summary <- c(NA, NA, NA, NA, NA, NA, NA, cell_type, case_id)
    } else {
      simulation_results_summary <- simulation_results_summary[simulation_results_summary$coverage >= covered_region_size,]
      simulation_results_summary <- simulation_results_summary[which.min(simulation_results_summary$burden),]
      simulation_results_sample <- rbind(simulation_results_sample, simulation_results_summary)
    }
  }
}

#simulation_results_sample[simulation_results_sample$case_id == "UMB1864_deeper" & simulation_results_sample$cell_type == "EN",] <-
#  c("a8s4", 0, 0, 0, 36861264, 32550349, NA, NA, "EN", "UMB1864_deeper")

# UMB4428: we detected 0 sSNVs from this sample
simulation_results_summary <- c(NA, 0, 0, NA, NA, 21577839, NA, NA, "OPC", "UMB4428")
simulation_results_sample <- rbind(simulation_results_sample, simulation_results_summary)
simulation_results_summary <- c(NA, 0, 0, NA, NA, 38896790, NA, NA, "IN", "UMB4428")
simulation_results_sample <- rbind(simulation_results_sample, simulation_results_summary)

simulation_results_sample <- data.frame(simulation_results_sample)
simulation_results_sample$count <- as.numeric(simulation_results_sample$count)
simulation_results_sample$sd <- as.numeric(simulation_results_sample$sd)
simulation_results_sample$burden <- as.numeric(simulation_results_sample$burden)
simulation_results_sample$coverage <- as.integer(simulation_results_sample$coverage)
simulation_results_sample[,"rate"] <- simulation_results_sample$count/simulation_results_sample$coverage
simulation_results_sample[,"rate_sd"] <- simulation_results_sample$sd/simulation_results_sample$coverage

simulation_results <- left_join(simulation_results_sample, metadata, c("case_id" = "Sample"))

saveRDS(simulation_results, "./results/burden_simulation/burden_summary_all_cell_type_20241018.rds")
readRDS("./results/burden_simulation/burden_summary_all_cell_type_20241018.rds") -> simulation_results

summary(lm(rate~Age:cell_type + Age, data=simulation_results))

ggplot(simulation_results, aes(x = Age, y = rate, color = cell_type)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  stat_cor() +
  theme_classic() +  
  geom_errorbar(aes(ymax = rate + rate_sd, ymin = rate - rate_sd), width = 1) +
  labs(y = "Mutation Burden (sSNVs/bp)")
ggsave("./figures/burden_simulation/simulation_results_20241018.pdf", width = 10, height = 5)

# ------------------------------------------------------------------------------
# EN & Oligodendrocyte Burden
simulation_results_EN_Oligo <- simulation_results[simulation_results$cell_type %in% c("EN", "Oligodendrocyte"),]

# Enrichment
EN_enrichment_summary <- readRDS("./results/open_region_enrichment/EN.rds")
EN_enrichment_summary[, "cell_type"] <- "EN"
Oligo_enrichment_summary <- readRDS("./results/open_region_enrichment/Oligodendrocyte.rds")
Oligo_enrichment_summary[, "cell_type"] <- "Oligodendrocyte"
enrichment_summary <- rbind(EN_enrichment_summary, Oligo_enrichment_summary)
colnames(enrichment_summary)[1] <- "case_id"

simulation_results_EN_Oligo <- left_join(simulation_results_EN_Oligo, enrichment_summary, by = c("case_id", "stringency", "cell_type"))
simulation_results_EN_Oligo[, "burden_corrected"] <-
  simulation_results_EN_Oligo$burden/simulation_results_EN_Oligo$enrichment_ratio
simulation_results_EN_Oligo$burden_sd <- 
(simulation_results_EN_Oligo$sd/simulation_results_EN_Oligo$coverage)*5.75e9
simulation_results_EN_Oligo[, "burden_sd_corrected"] <-
  simulation_results_EN_Oligo$burden_sd/simulation_results_EN_Oligo$enrichment_ratio


summary(lm(burden_corrected~Age, data=simulation_results_EN_Oligo[simulation_results_EN_Oligo$cell_type == "Oligodendrocyte",]))
summary(lm(burden_corrected~Age, data=simulation_results_EN_Oligo[simulation_results_EN_Oligo$cell_type == "EN",]))
summary(lm(burden_corrected~Age*cell_type, data=simulation_results_EN_Oligo))

ggplot(simulation_results_EN_Oligo, aes(x = Age, y = burden, color = cell_type)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  stat_cor() +
  geom_errorbar(aes(ymax = burden + burden_sd, ymin = burden - burden_sd)) +
  theme_classic()

ggplot(simulation_results_EN_Oligo, aes(x = Age, y = burden_corrected, color = cell_type)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  stat_cor() +
  geom_errorbar(aes(ymax = burden_corrected + burden_sd_corrected, ymin = burden_corrected - burden_sd_corrected), width = 3) +
  theme_classic() +
  labs(y = "Autosomal Mutation Burden")
ggsave("./figures/burden_simulation/EN_OL_20241018.pdf", width = 6, height = 5)

saveRDS(simulation_results_EN_Oligo, "./results/burden_simulation/burden_summary_EN_OL_20241018.rds")
write.csv(simulation_results_sample, "./results/burden_simulation/burden_summary_20241018.csv", quote = F, row.names = F)

