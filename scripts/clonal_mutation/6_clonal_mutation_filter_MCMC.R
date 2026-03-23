.libPaths(c("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackage_rstan/"))
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(stringr)
library(dplyr)
library(tidyr)


# ------------------------------------------------------------------------------
# All brain samples (control + ASD)
case_id_list_ASD <- c("ABNCR6D", "AN06365")
case_id_list_healthy_brain <- c("Brain1901106_deeper", "UMB1278_deeper", "UMB1864_deeper", 
                                "UMB5823_deeper", "UMB4638", "UMB1465", "UMB5451", "UMB5657")
case_id_list <- c(case_id_list_ASD, case_id_list_healthy_brain)

clonal_mutation_list_healthy_brain <- read.csv("results/clonal_mutation/clonal_mutation_list_20240812.csv")
clonal_mutation_list_ASD <- read.csv("results/clonal_mutation/clonal_mutation_list_ASD_20241028.csv")
clonal_mutation_list <- rbind(clonal_mutation_list_healthy_brain, clonal_mutation_list_ASD)
clonal_mutation_list$key <- str_replace_all(clonal_mutation_list$key, pattern = "[:]", replacement = "-")

# ------------------------------------------------------------------------------
# Germline filter
cell_cover_clonal_mutation_site_summary <- c()
for(case_id in case_id_list){
  print(case_id)
  
  # Read in clonal mutations
  cell_cover_clonal_mutation_site <- 
    readRDS(paste0("./results/clonal_mutation/clonal_mutation_VAF_", case_id, ".rds"))
  if(dim(cell_cover_clonal_mutation_site)[1] == 0){next}
  cell_cover_clonal_mutation_site["case_id"] <- case_id
  cell_cover_clonal_mutation_site["key"] <- row.names(cell_cover_clonal_mutation_site)
  
  cell_cover_clonal_mutation_site <- cell_cover_clonal_mutation_site[cell_cover_clonal_mutation_site[,"key"] %in% clonal_mutation_list[clonal_mutation_list$case_id == case_id, "key"],]
  
  cell_cover_clonal_mutation_site_summary <- rbind(cell_cover_clonal_mutation_site_summary, cell_cover_clonal_mutation_site)
}
cell_cover_clonal_mutation_site_summary <- data.frame(cell_cover_clonal_mutation_site_summary)

# remove mutations shared in multiple samples
selected_site <- (cell_cover_clonal_mutation_site_summary %>% group_by(key) %>% summarise("count" = n()))
selected_site <- selected_site[selected_site$count == 1, "key", drop = T]
cell_cover_clonal_mutation_site_summary <- 
  cell_cover_clonal_mutation_site_summary[cell_cover_clonal_mutation_site_summary$key %in% selected_site,]

# remove mutations in < 3 cells
cell_cover_clonal_mutation_site_summary <- 
  cell_cover_clonal_mutation_site_summary[cell_cover_clonal_mutation_site_summary$detected > 2,]

cell_cover_clonal_mutation_site_filtered_out_summary <- cell_cover_clonal_mutation_site_summary

# remove potential germline mutations
cell_cover_clonal_mutation_site_summary <- 
  cell_cover_clonal_mutation_site_summary[cell_cover_clonal_mutation_site_summary$detected < floor(cell_cover_clonal_mutation_site_summary$bound_0.001.),]

# only keep mutations VAF < 
# plot(density(cell_cover_clonal_mutation_site_summary$detected / cell_cover_clonal_mutation_site_summary$covered))
# cell_cover_clonal_mutation_site_summary <- 
#   cell_cover_clonal_mutation_site_summary[(cell_cover_clonal_mutation_site_summary$detected / cell_cover_clonal_mutation_site_summary$covered < 0.125*(0.35/0.5)),]
# plot(density(cell_cover_clonal_mutation_site_summary$detected / cell_cover_clonal_mutation_site_summary$covered))

cell_cover_clonal_mutation_site_filtered_out_summary <- 
  cell_cover_clonal_mutation_site_filtered_out_summary %>% anti_join(cell_cover_clonal_mutation_site_summary)

clonal_mutation_list <- rbind(clonal_mutation_list_healthy_brain, clonal_mutation_list_ASD)
clonal_mutation_list$key <- str_replace_all(clonal_mutation_list$key, pattern = "[:]", replacement = "-")
clonal_mutation_list_filtered <- inner_join(clonal_mutation_list, cell_cover_clonal_mutation_site_summary, by = c("key", "case_id"))
clonal_mutation_list_filtered[,"recovered"] <- FALSE
clonal_mutation_list_filtered <- 
  clonal_mutation_list_filtered[!(duplicated(str_c(clonal_mutation_list_filtered$cell, clonal_mutation_list_filtered$key, clonal_mutation_list_filtered$case_id))),]

# ------------------------------------------------------------------------------
# Recover clonal mutations that is found under standard stringency but removed by our germline filter
# stringency: a6s3, VAF==0 reads
clone <- c()
for (case_id in case_id_list){
  # stringency: a6s3, VAF==0 reads
  somatic_mutation <- read.table(paste0("./data/", case_id, "/single_cell/call_variants/autosomal_variants/a6s3_autosome_variants.txt"), header = F, sep = "\t")
  colnames(somatic_mutation) <- c("chr", "pos", "id", "ref", "alt", "qual", "filter", "info", "fromat", "value_cell")
  somatic_mutation <- separate(somatic_mutation, col = "value_cell", into = c("value", "cell"), sep = " ")
  somatic_mutation <- somatic_mutation[,c("chr", "pos", "cell", "ref", "alt", "qual", "filter", "info", "fromat", "value")]
  somatic_mutation[,"key"] <- unite(somatic_mutation, col = "key", c("chr", "pos", "ref", "alt"), sep = "-")[,"key"]
  somatic_mutation[,"case_id"] <- case_id
  
  filtered_out_clonal_mutation <-
    cell_cover_clonal_mutation_site_filtered_out_summary[cell_cover_clonal_mutation_site_filtered_out_summary$case_id == case_id,]
  somatic_clonal_mutation_list <- filtered_out_clonal_mutation[filtered_out_clonal_mutation$key %in% unique(somatic_mutation$key),]
  
  # for (large_clonal_snv in somatic_clonal_mutation_list){
  #   clone <- rbind(clone, somatic_mutation[somatic_mutation$key == large_clonal_snv, ])
  # }
  clone <- rbind(clone, somatic_clonal_mutation_list)
}

cell_type_healthy_brain <- read.table("./results/celltype_annotation/annotation_table_20240512.txt",
                                      sep = "\t")
cell_type_ASD <- read.table("./results/celltype_annotation/annotation_table_20240813_ASD.txt",
                            sep = "\t")
cell_type <- rbind(cell_type_healthy_brain, cell_type_ASD)
colnames(cell_type) <- c("case_id", "cell", "cell_type")
clonal_mutation_list_recovered <- inner_join(clone, clonal_mutation_list, by = c("key", "case_id"))
clonal_mutation_list_recovered[,"recovered"] <- TRUE
clonal_mutation_list_recovered <- clonal_mutation_list_recovered[!duplicated(str_c(clonal_mutation_list_recovered$cell, clonal_mutation_list_recovered$key, clonal_mutation_list_recovered$case_id)),]

# Combine
cell_cover_clonal_mutation_site_summary_combined <- rbind(cell_cover_clonal_mutation_site_summary, clone)
clonal_mutation_list_combined <- rbind(clonal_mutation_list_recovered, clonal_mutation_list_filtered)

# remove mutations in multiple samples
clonal_mutation_list_combined_count <- 
  clonal_mutation_list_combined %>% group_by(case_id, key) %>% 
  summarise() %>% group_by(key) %>% summarise("count" = n())
unique_clonal_mutation <- clonal_mutation_list_combined_count[clonal_mutation_list_combined_count$count == 1, "key", drop=T]
clonal_mutation_list_combined <- clonal_mutation_list_combined[clonal_mutation_list_combined$key %in% unique_clonal_mutation,]
clonal_mutation_list_combined_count <- clonal_mutation_list_combined %>% group_by(key) %>% summarise("count" = n())
selected_clonal_mutation <- clonal_mutation_list_combined_count[clonal_mutation_list_combined_count$count > 2, "key", drop = T]
clonal_mutation_list_combined <- clonal_mutation_list_combined[clonal_mutation_list_combined$key %in% selected_clonal_mutation,]

# save control and ASD clonal mutations separately
clonal_mutation_list_combined_ASD <- clonal_mutation_list_combined[clonal_mutation_list_combined$case_id %in% case_id_list_ASD,]
clonal_mutation_list_combined_healthy_brain <- clonal_mutation_list_combined[clonal_mutation_list_combined$case_id %in% case_id_list_healthy_brain,]

write.csv(clonal_mutation_list_combined_ASD,
          "./results/clonal_mutation/clonal_mutation_list_ASD_filtered_20241028.csv",
          quote = F, row.names = F)

write.csv(clonal_mutation_list_combined_healthy_brain,
          "./results/clonal_mutation/clonal_mutation_list_filtered_20241028.csv",
          quote = F, row.names = F)

# ------------------------------------------------------------------------------
# COLO829BLT50
case_id_list <- c("COLO829BLT50_rep1", "COLO829BLT50_rep2")

# ------------------------------------------------------------------------------
# Germline filter
cell_cover_clonal_mutation_site_summary <- c()
for(case_id in case_id_list){
  print(case_id)
  
  # Read in clonal mutations
  cell_cover_clonal_mutation_site <- 
    readRDS(paste0("./results/clonal_mutation/clonal_mutation_VAF_", case_id, ".rds"))
  if(dim(cell_cover_clonal_mutation_site)[1] == 0){next}
  cell_cover_clonal_mutation_site["case_id"] <- case_id
  cell_cover_clonal_mutation_site["key"] <- row.names(cell_cover_clonal_mutation_site)
  
  cell_cover_clonal_mutation_site_summary <- rbind(cell_cover_clonal_mutation_site_summary, cell_cover_clonal_mutation_site)
}
cell_cover_clonal_mutation_site_summary <- data.frame(cell_cover_clonal_mutation_site_summary)

cell_cover_clonal_mutation_site_filtered_out_summary <- 
  cell_cover_clonal_mutation_site_summary[cell_cover_clonal_mutation_site_summary$detected >= cell_cover_clonal_mutation_site_summary$bound_5.,]
cell_cover_clonal_mutation_site_summary <- 
  cell_cover_clonal_mutation_site_summary[cell_cover_clonal_mutation_site_summary$detected < cell_cover_clonal_mutation_site_summary$bound_5.,]

cell_cover_clonal_mutation_site_summary <- 
  cell_cover_clonal_mutation_site_summary[cell_cover_clonal_mutation_site_summary$detected > 2,]

clonal_mutation_list <- read.csv("results/clonal_mutation/clonal_mutation_list_COLO829BLT50_20240523.csv")
clonal_mutation_list$key <- str_replace_all(clonal_mutation_list$key, pattern = "[:]", replacement = "-")
clonal_mutation_list_filtered <- inner_join(clonal_mutation_list, cell_cover_clonal_mutation_site_summary, by = c("key", "case_id"))
clonal_mutation_list_filtered[,"recovered"] <- FALSE
clonal_mutation_list_filtered <- 
  clonal_mutation_list_filtered[!(duplicated(str_c(clonal_mutation_list_filtered$cell, clonal_mutation_list_filtered$key, clonal_mutation_list_filtered$case_id))),]

# ------------------------------------------------------------------------------
# Recover clonal mutations that is found under standard stringency but removed by our germline filter
# stringency: a6s3, VAF==0 reads
clone <- c()
for (case_id in case_id_list){
  # stringency: a6s3, VAF==0 reads
  somatic_mutation <- read.table(paste0("./data/", case_id, "/single_cell/call_variants/autosomal_variants/a6s3_autosome_variants.txt"), header = F, sep = "\t")
  colnames(somatic_mutation) <- c("chr", "pos", "id", "ref", "alt", "qual", "filter", "info", "fromat", "value_cell")
  somatic_mutation <- separate(somatic_mutation, col = "value_cell", into = c("value", "cell"), sep = " ")
  somatic_mutation <- somatic_mutation[,c("chr", "pos", "cell", "ref", "alt", "qual", "filter", "info", "fromat", "value")]
  somatic_mutation[,"key"] <- unite(somatic_mutation, col = "key", c("chr", "pos", "ref", "alt"), sep = "-")[,"key"]
  somatic_mutation[,"case_id"] <- case_id
  
  filtered_out_clonal_mutation <-
    cell_cover_clonal_mutation_site_filtered_out_summary[cell_cover_clonal_mutation_site_filtered_out_summary$case_id == case_id,]
  somatic_clonal_mutation_list <- filtered_out_clonal_mutation$key[filtered_out_clonal_mutation$key %in% somatic_mutation$key]
  
  for (large_clonal_snv in somatic_clonal_mutation_list){
    clone <- rbind(clone, somatic_mutation[somatic_mutation$key == large_clonal_snv, ])
  }
}

cell_type <- read.table("./results/celltype_annotation/annotation_COLO829BLT50_20240523.txt",
                        sep = "\t")
colnames(cell_type) <- c("case_id", "cell", "cell_type")
clone <- left_join(clone, cell_type, by = c("case_id", "cell"))
clone[,"recovered"] <- TRUE
clone <- clone[!duplicated(str_c(clone$cell, clone$key, clone$case_id)),]

# Combine
clonal_mutation_list_combined <- rbind(clone, clonal_mutation_list_filtered[, c("chr", "pos", "cell", "ref", "alt", "qual", "filter", "info", "fromat", "value", "key", "case_id", "cell_type", "recovered")])

write.csv(clonal_mutation_list_combined,
          "./results/clonal_mutation/clonal_mutation_list_filtered_COLO829BLT50_05292024.csv",
          quote = F, row.names = F)


