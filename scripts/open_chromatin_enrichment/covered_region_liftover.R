library(rtracklayer)

chain <- import("./ref/hg38ToHg19.over.chain")

for (cell_type in cell_type_list){
  for(stringency in c("a2s1", "a4s2", "a6s3", "a8s4", "a10s5")){
    regions <- read.table(paste0("./Duplex_multiome/data/others/covered_region_summary/", cell_type, "_open_region_", stringency, "_combined_merged.bed"))
    region_grange <- GRanges(seqnames = regions$V1, 
                             ranges = IRanges(start = regions$V2,
                                              end = regions$V3))
    
    liftover_results <- liftOver(region_grange, chain)
    df_liftover_results <- data.frame(liftover_results)
    write.table(df_liftover_results[,3:5], paste0("./Duplex_multiome/data/others/covered_region_summary/", cell_type, "_open_region_", stringency, "_combined_merged_hg19.bed"), quote = F, row.names = F)
  }
}
