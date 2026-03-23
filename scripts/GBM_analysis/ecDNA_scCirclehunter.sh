#!/bin/bash
#SBATCH --partition=bch-compute
#SBATCH --time=5-00:00:00	 			# Running time (in hours-minutes-seconds)
#SBATCH --job-name=ecDNA			# Job name
#SBATCH --mail-type=BEGIN,END,FAIL 		# send and email when the job begins, ends or fails
#SBATCH --mail-user=shulin_mao@fas.harvard.edu	 	# Email address to send the job status
#SBATCH --output=scCirclehunter_ecDNA_output.txt 			# Name of the output file
#SBATCH --mem=64G
#SBATCH --nodes=8				# Number of compute nodes

source /lab-share/Gene-Lee-e2/Public/home/shulin/softwares/anaconda3/etc/profile.d/conda.sh
conda activate /lab-share/Gene-Lee-e2/Public/home/shulin/softwares/mambaforge/envs/scCirclehunter

circlehunter2 -p 8 --blacklist /lab-share/Gene-Lee-e2/Public/home/shulin/ref/mask/hg38_mask.bed /lab-share/Gene-Walsh-e2/Public/indData/AK/Multiome-strand-02212022/Novogene_11242025/Mult_3828T_addl_seq/outs/atac_possorted_bam.bam /lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/results/ecDNA/3828T_addl_scCirclehunter_ecDNA.bed > scCirclehunter_ecDNA_output.txt 2>&1
circlehunter2 -p 8 --blacklist /lab-share/Gene-Lee-e2/Public/home/shulin/ref/mask/hg38_mask.bed /lab-share/Gene-Walsh-e2/Public/indData/AK/Multiome-strand-02212022/Novogene_11242025/Mult_3828N_addl_seq/outs/atac_possorted_bam.bam /lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/results/ecDNA/3828N_addl_scCirclehunter_ecDNA.bed > scCirclehunter_ecDNA_output.txt 2>&1

circlehunter2 -p 8 --blacklist /lab-share/Gene-Lee-e2/Public/home/shulin/ref/mask/hg38_mask.bed /lab-share/Gene-Walsh-e2/Public/indData/AK/Multiome-strand-02212022/Novogene_05152025/Mult_4397N/outs/atac_possorted_bam.bam /lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/results/ecDNA/UMB4397_normal_scCirclehunter_ecDNA.bed > scCirclehunter_ecDNA_output.txt 2>&1
circlehunter2 -p 8 --blacklist /lab-share/Gene-Lee-e2/Public/home/shulin/ref/mask/hg38_mask.bed /lab-share/Gene-Walsh-e2/Public/indData/AK/Multiome-strand-02212022/Novogene_05152025/Mult_4397T/outs/atac_possorted_bam.bam /lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/results/ecDNA/UMB4397_tumor_scCirclehunter_ecDNA.bed > scCirclehunter_ecDNA_output.txt 2>&1
