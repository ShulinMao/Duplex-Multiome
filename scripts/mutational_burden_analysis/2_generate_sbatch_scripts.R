.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(stringr)

case_id_list <- 
  c("UMB1932", "ABNCR6D", "AN06365")
  #c("Brain1901106_deeper", "UMB1278_deeper", "UMB1864_deeper", "UMB5823_deeper", "UMB4638", "UMB1465", "UMB5451", "UMB5657")
#c("Brain1901106", "UMB1278", "UMB1864", "UMB4638", "UMB1465", "UMB5451", "UMB5657", "UMB5823")

# ------------------------------------------------------------------------------
# cell type
directory <- "/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/scripts/simulation_burdens/cell_type/"
for(case_id in case_id_list){
  print(case_id)
  
  # cell type
  file_list <- dir(paste0("./data/", case_id, "/cell_type/"))
  cell_type_list <- str_remove(file_list, "[.]txt")
  
  for(cell_type in cell_type_list){
    file_name <- paste0(directory, case_id, "_", cell_type, ".sh")
    sink(file_name)
    cat("#!/bin/bash \n")
    cat("#SBATCH --partition=bch-compute	 			# queue to be used \n")
    cat("#SBATCH --time=3-00:00:00	 			# Running time (in hours-minutes-seconds) \n")
    cat(paste0("#SBATCH --job-name=", case_id, "_", cell_type, "		# Job name \n"))
    cat("#SBATCH --mail-type=BEGIN,END,FAIL 		# send and email when the job begins, ends or fails \n")
    cat("#SBATCH --mail-user=shulin_mao@fas.harvard.edu	 	# Email address to send the job status \n")
    cat(paste0("#SBATCH --output=output_", case_id, "_", cell_type, ".txt 			# Name of the output file \n"))
    cat("#SBATCH --mem=12G \n")
    cat("#SBATCH --nodes=1				# Number of compute nodes \n")
    cat("#SBATCH --ntasks=1			# Number of threads/tasks on one node \n")
    
    cat("source /programs/biogrids.shrc \n")
    cat("export R_X=4.0.2 \n")
    
    cat(paste("Rscript", 
              "/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/scripts/simulation_burdens/template.R",
              case_id, cell_type, "\n"))
    sink()
  }
}

file_name <- paste0(directory, "run.sh")
sink(file_name)
for(case_id in case_id_list){
  
  # cell type
  file_list <- dir(paste0("./data/", case_id, "/cell_type/"))
  cell_type_list <- str_remove(file_list, "[.]txt")
  
  for(cell_type in cell_type_list){
    cat(paste0("sbatch ", case_id, "_", cell_type, ".sh \n"))
  }
}
sink()


# ------------------------------------------------------------------------------
# sub-cell-type
directory <- "/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/scripts/simulation_burdens/subtype/"
for(case_id in case_id_list){
  print(case_id)
  
  # cell type
  file_list <- dir(paste0("./data/", case_id, "/cell_type/subtype/"))
  cell_type_list <- str_remove(file_list, "[.]txt")
  
  for(cell_type in cell_type_list){
    file_name <- paste0(directory, case_id, "_", cell_type, ".sh")
    sink(file_name)
    cat("#!/bin/bash \n")
    cat("#SBATCH --partition=bch-compute	 			# queue to be used \n")
    cat("#SBATCH --time=3-00:00:00	 			# Running time (in hours-minutes-seconds) \n")
    cat(paste0("#SBATCH --job-name=", case_id, "_", cell_type, "		# Job name \n"))
    cat("#SBATCH --mail-type=BEGIN,END,FAIL 		# send and email when the job begins, ends or fails \n")
    cat("#SBATCH --mail-user=shulin_mao@fas.harvard.edu	 	# Email address to send the job status \n")
    cat(paste0("#SBATCH --output=output_", case_id, "_", cell_type, ".txt 			# Name of the output file \n"))
    cat("#SBATCH --mem=12G \n")
    cat("#SBATCH --nodes=1				# Number of compute nodes \n")
    cat("#SBATCH --ntasks=1			# Number of threads/tasks on one node \n")
    
    cat("source /programs/biogrids.shrc \n")
    cat("export R_X=4.0.2 \n")
    
    cat(paste("Rscript", 
              "/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/scripts/simulation_burdens/template_subtype.R",
              case_id, str_split_fixed(cell_type, "_", n = 2)[1,1], str_split_fixed(cell_type, "_", n = 2)[1,2], "\n"))
    sink()
  }
}

file_name <- paste0(directory, "run_all_sbatch_scripts.sh")
sink(file_name)
for(case_id in case_id_list){
  
  # cell type
  file_list <- dir(paste0("./data/", case_id, "/cell_type/subtype/"))
  cell_type_list <- str_remove(file_list, "[.]txt")
  
  for(cell_type in cell_type_list){
    cat(paste0("sbatch ", case_id, "_", cell_type, ".sh \n"))
  }
}
sink()

