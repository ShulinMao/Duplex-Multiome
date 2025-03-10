.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(stringr)

case_id_list <- 
  c("Brain1901106_deeper", "UMB1278_deeper", "UMB1864_deeper", "UMB4428", "UMB5823_deeper", "UMB4638", "UMB1465", "UMB5451", "UMB5657")
  #c("COLO829BLT50_rep1", "COLO829BLT50_rep2")

directory <- "/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/scripts/single_strand_burden/scripts/"
for(case_id in case_id_list){
  print(case_id)

  file_name <- paste0(directory, case_id, ".sh")
  sink(file_name)
  cat("#!/bin/bash \n")
  cat("#SBATCH --partition=bch-compute	 			# queue to be used \n")
  cat("#SBATCH --time=2-00:00:00	 			# Running time (in hours-minutes-seconds) \n")
  cat(paste0("#SBATCH --job-name=", case_id, "		# Job name \n"))
  cat("#SBATCH --mail-type=BEGIN,END,FAIL 		# send and email when the job begins, ends or fails \n")
  cat("#SBATCH --mail-user=shulin_mao@fas.harvard.edu	 	# Email address to send the job status \n")
  cat(paste0("#SBATCH --output=output_", case_id, ".txt 			# Name of the output file \n"))
  cat("#SBATCH --mem=12G \n")
  cat("#SBATCH --nodes=1				# Number of compute nodes \n")
  cat("#SBATCH --ntasks=1			# Number of threads/tasks on one node \n")
  
  cat("source /programs/biogrids.shrc \n")
  cat("export R_X=4.0.2 \n")
  
  cat(paste("Rscript", 
            "/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/scripts/single_strand_burden/template.R",
            case_id, "\n"))
  sink()

}

file_name <- paste0(directory, "run.sh")
sink(file_name)
for(case_id in case_id_list){
  cat(paste0("sbatch ", case_id, ".sh \n"))
}
sink()

