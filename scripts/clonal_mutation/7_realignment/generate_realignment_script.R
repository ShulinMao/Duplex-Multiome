.libPaths("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/Rpackages/")
setwd("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome")

library(stringr)

clonal_mutation <- read.csv("./results/clonal_mutation/clonal_mutation_list_ASD_filtered_20241028.csv")

# ------------------------------------------------------------------------------
directory <- "/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/scripts/clonal_mutation/realignment/ASD/"
output_dir <- "/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/results/clonal_mutation/realigned/ASD/"

for (i in seq(dim(clonal_mutation)[1])) {
  case_id <- clonal_mutation$case_id[i]
  variant <- clonal_mutation$key[i]
  cell <- clonal_mutation$cell[i]
  wd <- paste0("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/", case_id, "/single_cell/")
  RF <- str_split(clonal_mutation$value[i], ":", simplify = T)[1, 5]
  
  file_name <- paste0(directory, case_id, "_", cell, "_", variant, ".sh")
  sink(file_name)
  cat("#!/bin/bash \n")
  cat("#SBATCH --partition=bch-compute	 			# queue to be used \n")
  cat("#SBATCH --time=00:30:00	 			# Running time (in hours-minutes-seconds) \n")
  cat(paste0("#SBATCH --job-name=", case_id, "_", cell, "_", variant, "		# Job name \n"))
  cat("#SBATCH --mail-type=FAIL 		# send and email when the job begins, ends or fails \n")
  cat("#SBATCH --mail-user=shulin_mao@fas.harvard.edu	 	# Email address to send the job status \n")
  cat(paste0("#SBATCH --output=output_", case_id, "_", cell, "_", variant, ".txt 			# Name of the output file \n"))
  cat("#SBATCH --mem=8G \n")
  cat("#SBATCH --nodes=1				# Number of compute nodes \n")
  cat("#SBATCH --ntasks=2			# Number of threads/tasks on one node \n")
  
  cat("source /programs/biogrids.shrc\n\n")
  
  cat(paste0("(samtools view -h ", 
             wd, "read_families/", cell, ".families.bam | \n",
             "awk '{if(match($0, \"^@\") || match($0, \"RF:Z:", RF, "$\")){print $0}}' | \n",
             "samtools view -b - | \n",
             "samtools collate -Oun128  - | \n",
             "samtools fastq -OT \"*\" - | \n",
             "bwa mem -p -O 4,4 -C /lab-share/Gene-Lee-e2/Public/home/shulin/ref/hg38/genome.fa - | \n",
             "samtools sort -o ", output_dir, case_id, ".", cell, ".", variant, ".families.realigned.bam -) &&\n",
             "(samtools index -b ", output_dir, case_id, ".", cell, ".", variant, ".families.realigned.bam)\n\n"))
  
  cat(paste0(
    "/home/ch241780/go/bin/mcsCallVariants ",
		"-r /lab-share/Gene-Lee-e2/Public/home/shulin/ref/hg38/genome.fa ",
		"-i ", output_dir, case_id, ".", cell, ".", variant, ".families.realigned.bam ",
		"-o ", output_dir, case_id, ".", cell, ".", variant, ".realigned.vcf ",
		"-b ", wd, "read_families/", cell, ".families.bed ",
		"-e /lab-share/Gene-Lee-e2/Public/home/shulin/ref/mask/hg38_mask.bed ",
		"-a 2 -s 0 -threads 1 -minReadFamilyLength 0\n\n"
  ))

  cat(paste0("(samtools view -h ",
            wd, "read_families/", cell, ".families.bam | \n",
            "awk '{if(match($0, \"^@\") || match($0, \"RF:Z:", RF, "$\")){print $0}}' | \n",
            "samtools view -b - | \n",
            "samtools reheader -P /lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/results/clonal_mutation/realigned/bam_header/header_", case_id, " - > ", output_dir, case_id, "_", cell, "_", variant, ".bam) && \n",
            "(samtools index ", output_dir, case_id, "_", cell, "_", variant, ".bam)\n\n"))

  cat(paste0("gatk --java-options \"-Xmx4g\" HaplotypeCaller ",
             "-R /lab-share/Gene-Lee-e2/Public/home/shulin/ref/hg38/genome.fa ",
             "-I ", output_dir, case_id, "_", cell, "_", variant, ".bam ",
             "-O ", output_dir, case_id, "_", cell, "_", variant, ".vcf ",
             "--disable-read-filter NotDuplicateReadFilter\n"))
  sink()
}


file_name <- paste0(directory, "run.sh")
sink(file_name)
dir(output_dir, pattern = ".realigned.vcf$") -> file_list
for (i in seq(dim(clonal_mutation)[1])) {
  case_id <- clonal_mutation$case_id[i]
  variant <- clonal_mutation$key[i]
  cell <- clonal_mutation$cell[i]
  wd <- paste0("/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/", case_id, "/single_cell/")
  RF <- str_split(clonal_mutation$value[i], ":", simplify = T)[1, 5]
  
  file_name <- paste0(directory, case_id, "_", cell, "_", variant, ".sh")
  cat(paste0("sbatch ", file_name, "\n"))
}
sink()

