#!/bin/bash
#SBATCH --partition=bch-compute	 			# queue to be used
#SBATCH --time=01:00:00	 			# Running time (in hours-minutes-seconds)
#SBATCH --job-name=amplicon_seq			# Job name
#SBATCH --mail-type=FAIL 		# send and email when the job begins, ends or fails
#SBATCH --mail-user=shulin_mao@fas.harvard.edu	 	# Email address to send the job status
#SBATCH --output=/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/scripts/amplicon-seq/slurm-%j.out 			# Name of the output file
#SBATCH --mem=32G
#SBATCH --ntasks=1				# Number of compute nodes
#SBATCH --cpus-per-task=2
#SBATCH -A bch


# software version control
source /programs/biogrids.shrc
export PICARD_X=2.26.10;SAMTOOLS_X=1.13;JAVA_X=jdk1.8.0_144

gatk_3_6=/lab-share/Gene-Lee-e2/Public/home/shulin/softwares/gatk-3.6/gatk
NaiveCaller=/lab-share/Gene-Lee-e2/Public/home/shulin/softwares/NaiveCaller

fastq_dir=$1
BARCODE_ID=$2
fastq_R1_file=${fastq_dir}/${BARCODE_ID}_R1_001.fastq.gz
fastq_R2_file=${fastq_dir}/${BARCODE_ID}_R2_001.fastq.gz
ref_genome=/lab-share/Gene-Lee-e2/Public/home/shulin/ref/hg38/genome.fa
output_dir=$3
candidate_pos=${output_dir}/candidate_pos.bed
# surrounding_pos=${output_dir}/surrounding_pos.bed


mkdir -p ${output_dir}/tmp 
tmp_dir=${output_dir}/tmp

# Read alignment
echo "Read alignment started"

mkdir -p ${output_dir}/bam
bwa mem -M -t 1 ${ref_genome} ${fastq_R1_file} ${fastq_R2_file} | samtools view -bSho ${output_dir}/bam/${BARCODE_ID}_bwa.bam -

picard SortSam \
	-I ${output_dir}/bam/${BARCODE_ID}_bwa.bam \
	-O ${output_dir}/bam/${BARCODE_ID}_bwa_sorted.bam \
	--TMP_DIR ${tmp_dir} \
	--SORT_ORDER coordinate \
	--VALIDATION_STRINGENCY LENIENT \
	--CREATE_INDEX false
samtools index ${output_dir}/bam/${BARCODE_ID}_bwa_sorted.bam

picard AddOrReplaceReadGroups \
	-I ${output_dir}/bam/${BARCODE_ID}_bwa_sorted.bam \
	-O ${output_dir}/bam/${BARCODE_ID}_bwa_sorted_addrg.bam \
	--RGID ${BARCODE_ID} \
	--RGLB ${BARCODE_ID} \
	--RGPL illumina \
	--RGPU unit1 \
	--RGSM ${BARCODE_ID}
samtools index ${output_dir}/bam/${BARCODE_ID}_bwa_sorted_addrg.bam

java -XX:+UseSerialGC -Xmx64G -Djava.io.tmpdir=${tmp_dir} -jar ${gatk_3_6}/GenomeAnalysisTK.jar \
	-T RealignerTargetCreator \
	-nt 2 \
	-R ${ref_genome} \
	-I ${output_dir}/bam/${BARCODE_ID}_bwa_sorted_addrg.bam \
	-o ${output_dir}/bam/${BARCODE_ID}.realigned.intervals

java -XX:+UseSerialGC -Xmx64G -Djava.io.tmpdir=${tmp_dir} -jar ${gatk_3_6}/GenomeAnalysisTK.jar \
	-T IndelRealigner \
	-R ${ref_genome} \
	-I ${output_dir}/bam/${BARCODE_ID}_bwa_sorted_addrg.bam \
	-targetIntervals ${output_dir}/bam/${BARCODE_ID}.realigned.intervals \
	-o ${output_dir}/bam/${BARCODE_ID}.final.bam
samtools index ${output_dir}/bam/${BARCODE_ID}.final.bam

rm -f ${output_dir}/bam/${BARCODE_ID}_bwa*

echo "Read alignment finished"

# Calling
rm -rf ${output_dir}/NaiveCaller
mkdir -p ${output_dir}/NaiveCaller

java -Xmx64g -jar ${NaiveCaller}/NaiveCaller.jar \
	-r ${ref_genome} \
	-b ${output_dir}/bam/${BARCODE_ID}.final.bam \
	-e ${candidate_pos} \
	-P 0 -q 30 -Q 20 -o ${output_dir}/NaiveCaller/${BARCODE_ID}.NC.tsv
	
# java -Xmx64g -jar ${NaiveCaller}/NaiveCaller.jar \
# 	-r ${ref_genome} \
# 	-b ${output_dir}/bam/${BARCODE_ID}.final.bam \
# 	-e ${surrounding_pos} \
# 	-P 0 -q 30 -Q 20 -o ${output_dir}/NaiveCaller/${BARCODE_ID}.NC.surrounding.tsv
