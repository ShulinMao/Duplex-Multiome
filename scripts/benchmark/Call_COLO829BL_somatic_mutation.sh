#!/bin/bash
#SBATCH --partition=bch-compute	 			# queue to be used
#SBATCH --time=10-00:00:00	 			# Running time (in hours-minutes-seconds)
#SBATCH --job-name=Mutect2_COLO829BL		# Job name
#SBATCH --mail-type=BEGIN,END,FAIL 		# send and email when the job begins, ends or fails
#SBATCH --mail-user=shulin_mao@fas.harvard.edu	 	# Email address to send the job status
#SBATCH --output=output_Mutect2_COLO829BL.txt 			# Name of the output file
#SBATCH --mem=128G
#SBATCH --nodes=2				# Number of compute nodes

source /programs/biogrids.shrc

gatk Mutect2 \
  -R /lab-share/Gene-Lee-e2/Public/home/shulin/ref/hg38/genome.fa \
  -I /lab-share/Gene-Walsh-e2/Public/indData/AK/bulk_bam/COLO829BL/SMHTCOLO829BL-X-X-M45-A001-uwsc-SMAFIAN1924G-sentieon_bwamem_202308.01_GRCh38.aligned.sorted.bam \
  -I /lab-share/Gene-Walsh-e2/Public/indData/AK/bulk_bam/COLO829T/SMHTCOLO829T-X-X-M45-A001-uwsc-SMAFIZO17Q9I-sentieon_bwamem_202308.01_GRCh38.aligned.sorted.bam \
  -normal SMACU1SULTKP \
  --germline-resource /lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/others/COLO829BLT50_truth_set/colo829bl_germline.vcf \
  --panel-of-normals /lab-share/Gene-Lee-e2/Public/home/shulin/ref/gatk_pon/1000g_pon.hg38.vcf.gz \
  -O /lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/others/COLO829BLT50_truth_set/COLO829BL_somatic_07252024.vcf.gz
