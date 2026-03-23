.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(stringr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(MASS)
library(truncnorm)
set.seed(42)

# read in sample & cell type from command lines
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

case_id <- args[1]
cell_type <- args[2]
subtype <- args[3]
rm(args)

print(case_id)
print(cell_type)
print(subtype)

# ------------------------------------------------------------------------------
# read in metadata
# cell type
read.table(paste0("./data/", case_id, "/cell_type/subtype/", cell_type, "_", subtype, ".txt"), sep = "\t") -> selected_cells

# the # of incorrect read families
rf_num <- read.table(paste0("./data/", case_id, "/single_cell/read_family_strand_barcode/read_family_strand_barcode_summary.txt"), header = T)

# error rate for incorrect read families
error_rate_incorrect_rf_summary_raw <- readRDS(paste0("./results/error_rate_incorrect_rf/", case_id, ".rds"))[[cell_type]]
error_rate_incorrect_rf_summary_raw <- as.data.frame(error_rate_incorrect_rf_summary_raw)
error_rate_incorrect_rf_summary_raw$levels <- c(1:5)
error_rate_incorrect_rf_summary <- 
  cbind(lm(ssnv_rate_strand2 ~ levels, data = error_rate_incorrect_rf_summary_raw)$fitted.values, 
        lm(ssnv_rate_strand12 ~ levels, data = error_rate_incorrect_rf_summary_raw)$fitted.values)

# ------------------------------------------------------------------------------
# Collect mutation data (require deplex consensus)
path <- paste0("./data/", case_id, "/single_cell/burden_estimation/Dan/")
cell_id <- read.table(paste0("./data/", case_id, "/single_cell/per_barcode_metrics.cells"), header = F)
burden_summary <- data.frame(row.names = cell_id$V1)
sSNV_number_summary <- data.frame(row.names = cell_id$V1)
coverage_summary <- data.frame(row.names = cell_id$V1)
for (summary_file in paste0(path, "/", 
                            c("a2s1", "a4s2", "a6s3", "a8s4", "a10s5"),
                            "_burden_results.txt")){
  burden_result <- read.table(summary_file, sep = "\t", header = F)
  cell_id <- str_split(burden_result$V5, pattern = "[.]", simplify = T)[,2]
  burden_summary[cell_id, str_extract(summary_file, pattern = "a[0-9]*s[0-9]*")] <- burden_result$V4 * 5.75e9
  sSNV_number_summary[cell_id, str_extract(summary_file, pattern = "a[0-9]*s[0-9]*")] <- burden_result$V1
  coverage_summary[cell_id, str_extract(summary_file, pattern = "a[0-9]*s[0-9]*")] <- burden_result$V3
}

# germline sensitivity
path <- paste0("./data/", case_id, "/single_cell/call_variants/germline_mutations") 
cell_id <- read.table(paste0("./data/", case_id, "/single_cell/per_barcode_metrics.cells"), header = F)
germline_mutation_summary <- data.frame(row.names = cell_id$V1)
for (summary_file in paste0(path, "/", 
                            c("a2s1", "a4s2", "a6s3", "a8s4", "a10s5"),
                            "_germline_mutation_results.txt")){
  germline_mutation_result <- read.table(summary_file, sep = "\t", header = F)
  cell_id <- str_split(germline_mutation_result$V2, pattern = "[/]", simplify = T)
  cell_id <- cell_id[,dim(cell_id)[2]]
  cell_id <- str_split(cell_id, pattern = "[_]", simplify = T)[,1]
  germline_mutation_summary[cell_id, str_extract(summary_file, pattern = "a[0-9]*s[0-9]*")] <- germline_mutation_result$V1
}

path <- paste0("./data/", case_id, "/single_cell/call_variants/all_germline_mutations_callable_region") 
cell_id <- read.table(paste0("./data/", case_id, "/single_cell/per_barcode_metrics.cells"), header = F)
all_possible_germline_mutation_callabel_region = data.frame(row.names = cell_id$V1)
for (summary_file in paste0(path, "/all_germline_mutations_callable_region_summary_", 
                            c("a2s1", "a4s2", "a6s3", "a8s4", "a10s5"), ".txt")){
  germline_mutation_result <- read.table(summary_file, sep = "\t", header = F)
  cell_id <- str_split(germline_mutation_result$V2, pattern = "[/]", simplify = T)
  cell_id <- cell_id[,dim(cell_id)[2]]
  cell_id <- str_split(cell_id, pattern = "[_]", simplify = T)[,1]
  all_possible_germline_mutation_callabel_region[cell_id, str_extract(summary_file, pattern = "a[0-9]*s[0-9]*")] <- germline_mutation_result$V1
}

# ------------------------------------------------------------------------------
# extract cells from the selected cell type
detection_efficiency_all_cell <- germline_mutation_summary/all_possible_germline_mutation_callabel_region
detection_efficiency_all_cell["cell_id"] <- row.names(detection_efficiency_all_cell)
coverage_summary["cell_id"] <- row.names(coverage_summary)
sSNV_number_summary["cell_id"] <- row.names(sSNV_number_summary)
coverage_detection_efficiency_all_cell <- inner_join(detection_efficiency_all_cell, coverage_summary, by="cell_id", suffix=c(".efficiency", ".coverage"))
coverage_detection_efficiency_all_cell <- inner_join(coverage_detection_efficiency_all_cell, sSNV_number_summary, by="cell_id")

rf_num_celltype <- rf_num[rf_num$cell %in% selected_cells$V2,]
colSums(rf_num_celltype[,1:20])

prop_rf_strand2 <- data.frame(rf_num_celltype$a2s0/rf_num_celltype$a2s1,
                              rf_num_celltype$a4s0/rf_num_celltype$a4s2, 
                              rf_num_celltype$a6s0/rf_num_celltype$a6s3, 
                              rf_num_celltype$a8s0/rf_num_celltype$a8s4, 
                              rf_num_celltype$a10s0/rf_num_celltype$a10s5)
prop_rf_strand2 <- reshape2::melt(prop_rf_strand2)
prop_rf_strand2$value[is.infinite(prop_rf_strand2$value)] <- NA
prop_rf_strand2 <- na.omit(prop_rf_strand2)
prop_rf_strand2$strigency <- as.numeric(as.character(factor(prop_rf_strand2$variable, levels(prop_rf_strand2$variable), c(1:5))))
prop_rf_strand2[prop_rf_strand2$value == 0, "value"] <- NA
prop_rf_strand2_lm <- lm(log(value) ~ strigency, data = prop_rf_strand2)
summary(prop_rf_strand2_lm)

prop_rf_strand12 <- data.frame(rf_num_celltype[, paste0("a4", "s", 1:1)]/rf_num_celltype$a4s2,
                               rowSums(rf_num_celltype[, paste0("a6", "s", 1:2)])/rf_num_celltype$a6s3,
                               rowSums(rf_num_celltype[, paste0("a8", "s", 1:3)])/rf_num_celltype$a8s4,
                               rowSums(rf_num_celltype[, paste0("a10", "s", 1:4)])/rf_num_celltype$a10s5)
prop_rf_strand12 <- reshape2::melt(prop_rf_strand12)
prop_rf_strand12$value[is.infinite(prop_rf_strand12$value)] <- NA
prop_rf_strand12 <- na.omit(prop_rf_strand12)

prop_rf_strand12$strigency <- as.numeric(as.character(factor(prop_rf_strand12$variable, levels(prop_rf_strand12$variable), c(2:5))))
prop_rf_strand12[prop_rf_strand12$value == 0, "value"] <- NA
prop_rf_strand12_lm <- lm(log(value) ~ strigency, data = prop_rf_strand12)
summary(prop_rf_strand12_lm)

rf_num_celltype <- cbind(exp(predict(prop_rf_strand2_lm, data.frame("strigency" = 1:5))), c(0, exp(predict(prop_rf_strand12_lm, data.frame("strigency" = 2:5)))))

simulation_results <- list()
for(a in seq(2, 10, 2)){
  print(paste0("a", a))
  s <- a%/%2
  
  # ------------------------------------------------------------------------------
  # calculate burden
  stringency_sim <- paste0("a", a, "s", s)
  if (a == 2){
    simulation_time <- 20000
  }else{
    simulation_time <- 50000
  }
  
  
  coverage_detection_efficiency <- coverage_detection_efficiency_all_cell[coverage_detection_efficiency_all_cell$cell_id %in% selected_cells$V2,]
  colSums(coverage_detection_efficiency[, c("a2s1", "a4s2", "a6s3", "a8s4", "a10s5")])
  open_chromatin <- colSums(coverage_detection_efficiency[, paste0(c("a2s1", "a4s2", "a6s3", "a8s4", "a10s5"), ".coverage")])
  open_chromatin <- open_chromatin[paste0(stringency_sim, ".coverage")]
  
  error_rate_incorrect_rf <- error_rate_incorrect_rf_summary[paste0("a", a),]
  
  if (a <= 8){
    proportion_rf <- rf_num_celltype[s,]/0.75
  }else{
    proportion_rf <- rf_num_celltype[s,]/0.75/0.87 
  }
  sensitivity <- 2*median(coverage_detection_efficiency[, paste0(stringency_sim, ".efficiency")], na.rm = T)
  
  init <- (sum(coverage_detection_efficiency[, stringency_sim]) - sum(proportion_rf*open_chromatin*error_rate_incorrect_rf))/sensitivity
  print(init)
  #5.75e9*(sum(coverage_detection_efficiency[, stringency_sim]) - sum(proportion_rf*open_chromatin*error_rate_incorrect_rf))/sensitivity/(open_chromatin*(1-sum(proportion_rf)))
  
  results <- c()
  # skip if 0 snv were detected under a stringency
  if (sum(coverage_detection_efficiency[, stringency_sim]) == 0 | init < 0){
    simulation_results[[stringency_sim]] <- rbind(results, c(0, 0, 0, as.numeric(open_chromatin), (open_chromatin*(1-sum(proportion_rf))), -Inf, Inf))
  }else{
    for (lambda in init * seq(0.1, 2, 0.1)){
      for (sd in c(0.5, 1, 1.5, 2)){
        detected_snv_num_simulation <- c()
        for (i in 1:simulation_time){
          
          cell_selected <- sample(1:dim(coverage_detection_efficiency)[1], 1) 
          covered_region <- coverage_detection_efficiency[cell_selected, paste0(stringency_sim, ".coverage")]
          snv_num <- rtruncnorm(1, a = 0, b = Inf, mean = lambda, sd = sd)
          
          rate <- error_rate_incorrect_rf
          if(length(rate) == 2 & rate[2] < snv_num/(open_chromatin*(1-sum(proportion_rf))) & (1-sum(proportion_rf))){rate[2] <- snv_num/(open_chromatin*(1-sum(proportion_rf)))}
          
          covered_snv_num <- sum(runif(round(covered_region*(1-sum(proportion_rf)))) <= (snv_num/(open_chromatin*(1-sum(proportion_rf)))))
          
          detected_snv_num <- 0
          for (rf in 1:length(proportion_rf)){
            detected_snv_num <- detected_snv_num + sum(runif(round(covered_region*proportion_rf[rf])) <= rate[rf])
          }
          
          if(!is.na(sensitivity)){
            if (covered_snv_num > 0){
              detected_snv_num <- sum(runif(covered_snv_num) <= as.numeric(sensitivity)) + detected_snv_num
            }
          }
          detected_snv_num_simulation <- c(detected_snv_num_simulation, detected_snv_num)
        }
        
        freq <- table(coverage_detection_efficiency[, stringency_sim])
        poisson_model <- fitdistr(detected_snv_num_simulation, densfun = "poisson")
        prob <- dpois(as.integer(names(freq)), poisson_model$estimate)
        likelihood <- sum(log(prob)*freq)
        
        
        sampling_results <- data.frame(table(detected_snv_num_simulation))
        sampling_results$detected_snv_num_simulation <- as.integer(sampling_results$detected_snv_num_simulation)
        sampling_results$Freq <- (sampling_results$Freq / simulation_time) * dim(coverage_detection_efficiency)[1]
        sampling_results_Freq <- sampling_results$Freq
        if(length(sampling_results_Freq) > length(freq)){
          freq <- c(freq, rep(0, length(sampling_results$Freq) - length(freq)))
        }else if(length(sampling_results_Freq) < length(freq)){
          sampling_results_Freq <- c(sampling_results_Freq, rep(0, length(freq) - length(sampling_results$Freq)))
        }
        print(sampling_results_Freq)
        print(freq)
        difference <- sum((sampling_results_Freq - freq)**2)
        #print(paste(snv_num, difference))
        print(paste(lambda, sd, lambda/(open_chromatin*(1-sum(proportion_rf)))*5.75e9, likelihood, difference))
        results <- rbind(results, c(lambda, sd, lambda/(open_chromatin*(1-sum(proportion_rf)))*5.75e9, open_chromatin, (open_chromatin*(1-sum(proportion_rf))), likelihood, difference))
      }
    }
    simulation_results[[stringency_sim]] <- results
    results <- c()
  }
}

saveRDS(simulation_results, paste0("./results/burden_simulation/subtype/", case_id, "_", cell_type, "_", subtype, ".rds"))
