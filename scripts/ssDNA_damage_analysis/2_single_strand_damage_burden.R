# Libraries
.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(stringr)
library(ggplot2)
library(dplyr)
library(ggpubr)

duplex_metadata <- read.csv("./data/others/metadata_07012024.csv")

case_id_list <- c("Brain1901106_deeper", "UMB1278_deeper", "UMB1864_deeper", "UMB1465", "UMB5451", "UMB5657", "UMB4638", "UMB5823_deeper")

ssDNA_mut_summary <- c()
# convert vcf files to grange objects
for (case_id in case_id_list){
  print(case_id)
  ssDNA_mut <- readRDS(paste0("./results/ssDNA_mut/", case_id, ".rds"))
  
  # cell type
  file_list <- dir(paste0("./data/", case_id, "/cell_type/"), pattern = "*.txt")
  cell_type_list <- str_remove(file_list, "[.]txt")
  
  # raw counts
  for(cell_type in cell_type_list){
    ssDNA_mut_cell_type <- ssDNA_mut[[cell_type]]
    ssDNA_mut_summary <- rbind(ssDNA_mut_summary, c(case_id, cell_type, ssDNA_mut_cell_type[1, 1]))
  }
}
ssDNA_mut_summary <- data.frame(ssDNA_mut_summary)
colnames(ssDNA_mut_summary) <- c("case_id", "cell_type", "rate")

ssDNA_mut_summary <- left_join(ssDNA_mut_summary, duplex_metadata, by = c("case_id" = "Sample"))
ssDNA_mut_summary$rate <- as.numeric(ssDNA_mut_summary$rate)

ggplot(ssDNA_mut_summary, aes(x=Age, y=rate, group=cell_type, color=cell_type)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  stat_cor() +
  labs(x = "Age", y = "ssDNA damage burden (per base)", color = "Cell type") +
  theme_classic()
ggsave("./figures/single_strand_damage_spectrum/single_strand_damage_rate.pdf", width = 6, height = 4)
