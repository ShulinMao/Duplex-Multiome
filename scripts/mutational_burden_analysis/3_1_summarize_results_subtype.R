.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")


library(stringr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)
library(truncnorm)
set.seed(42)

library(httpgd)
#hgd()

source("./scripts/simulation_burdens/simulation_regression_func.R")

case_id_list <- c("Brain1901106_deeper", "UMB1278_deeper", "UMB1864_deeper", "UMB5823_deeper", 
                  "UMB4638", "UMB1465", "UMB5451", "UMB5657")
  #c("Brain1901106", "UMB1278", "UMB1864", "UMB4638", "UMB1465", "UMB5451", "UMB5657", "UMB5823")

metadata <- read.csv("./data/others/metadata_07012024.csv")

# ------------------------------------------------------------------------------
# summarize simulation results
simulation_results_sample <- c()

for(case_id in case_id_list){
  print(case_id)
  
  # sub cell type
  file_list <- dir(paste0("./data/", case_id, "/cell_type/subtype"), pattern = ".*txt")
  cell_type_list <- str_remove(file_list, "[.]txt")
  
  for (cell_type in cell_type_list){
    print(cell_type)

    # read in simulation results
    simulation_results <- readRDS(paste0("./results/burden_simulation/subtype/", case_id, "_", cell_type, ".rds"))
    
    simulation_results_summary <- c()
    for (a in seq(6, 10, 2)){
      s <- a%/%2
      stringency <- paste0("a", a, "s", s)
      # select the results with the lowest error
      simulation_results_summary <-
        rbind(simulation_results_summary, 
            #c(stringency, simulation_results[[stringency]][which.max(simulation_results[[stringency]][,6]),]))
            c(stringency, simulation_results[[stringency]][which.min(simulation_results[[stringency]][,7]),]))
    }
    colnames(simulation_results_summary) <- c("stringency", "count", "sd", "burden", "raw_coverage", "coverage", "likelihood", "error")
    simulation_results_summary <- data.frame(simulation_results_summary)
    simulation_results_summary["cell_type"] <- cell_type
    simulation_results_summary["case_id"] <- case_id
    simulation_results_summary$coverage <- as.integer(simulation_results_summary$coverage)
    
    # select the fittest situations
    # a6s3 or more stringent
    # simulation_results_summary <- simulation_results_summary[3:5,]
    # covered regions > 20000000 bps
    covered_region_size <- 20000000
    if(sum(simulation_results_summary$coverage >= covered_region_size) == 0){
      #simulation_results_summary <- c(NA, NA, NA, NA, NA, NA, NA, cell_type, case_id)
    } else {
      simulation_results_summary <- simulation_results_summary[simulation_results_summary$coverage >= covered_region_size,]
      simulation_results_summary <- simulation_results_summary[which.min(simulation_results_summary$burden),]
      simulation_results_sample <- rbind(simulation_results_sample, simulation_results_summary)
    }
  }
}

simulation_results_sample <- data.frame(simulation_results_sample)
simulation_results_sample$count <- as.numeric(simulation_results_sample$count)
simulation_results_sample$sd <- as.numeric(simulation_results_sample$sd)
simulation_results_sample$burden <- as.numeric(simulation_results_sample$burden)
simulation_results_sample$coverage <- as.integer(simulation_results_sample$coverage)
simulation_results_sample[,"rate"] <- simulation_results_sample$count/simulation_results_sample$coverage
simulation_results_sample[,"rate_sd"] <- simulation_results_sample$sd/simulation_results_sample$coverage

simulation_results <- left_join(simulation_results_sample, metadata, c("case_id" = "Sample"))

saveRDS(simulation_results, "./results/burden_simulation/burden_summary_subtype_20241018.rds")
simulation_results <- readRDS("./results/burden_simulation/burden_summary_subtype_20241018.rds")

ggplot(simulation_results, aes(x = Age, y = rate, color = cell_type)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  stat_cor() +
  theme_classic() +
  geom_errorbar(aes(ymax = rate + rate_sd, ymin = rate - rate_sd), width = 1) +
  labs(y = "Mutation Burden (sSNVs/bp)")
ggsave("./figures/burden_simulation/simulation_results_subtype_20241018.pdf", width = 10, height = 5)

ggplot(simulation_results[str_detect(simulation_results$cell_type, "IN"),], aes(x = Age, y = rate, color = cell_type)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  stat_cor() +
  theme_classic() +
  geom_errorbar(aes(ymax = rate + rate_sd, ymin = rate - rate_sd), width = 1) +
  labs(y = "Mutation Burden (sSNVs/bp)")
