#!/bin/bash

# software version control
source /programs/biogrids.shrc
export PICARD_X=2.26.10;GATK_X=4.3.0.0;SAMTOOLS_X=1.13

# path
input_dir="/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/Brain1901106/bulk"
sample_name="Brain1901106"
fastq_R1_file="/lab-share/Gene-Walsh-e2/Public/indData/AK/for_Shulin/fastq/CON1901106_bulk_R1.fastq.gz"
fastq_R2_file="/lab-share/Gene-Walsh-e2/Public/indData/AK/for_Shulin/fastq/CON1901106_bulk_R2.fastq.gz"
mkdir ${input_dir}/bam
output_dir=${input_dir}/bam

ref_genome="/lab-share/Gene-Lee-e2/Public/home/shulin/ref/hg38/genome.fa"
known_sites="/lab-share/Gene-Lee-e2/Public/home/shulin/ref/dbsnp_147_b38_common_all_20160601.vcf"


# (1) Preprocessing
echo "Preprocessing start"

gzip -c ${input_dir}/${fastq_R1_file} > ${input_dir}/${fastq_R1_file}.gz
gzip -c ${input_dir}/${fastq_R2_file} > ${input_dir}/${fastq_R2_file}.gz

fastq_R1_file=${fastq_R1_file}.gz
fastq_R2_file=${fastq_R2_file}.gz

echo "Preprocessing done"

# (2) Alignment

echo "Alignment start"

bwa mem -M -t 8 ${ref_genome} ${fastq_R1_file} ${fastq_R2_file} | samtools view -bSho ${output_dir}/${sample_name}_bwa.bam -

mkdir ${output_dir}/tmp 
picard SortSam \
	-I ${output_dir}/${sample_name}_bwa.bam \
	-O ${output_dir}/${sample_name}_bwa_sorted.bam \
	--TMP_DIR ${output_dir}/tmp \
	--SORT_ORDER coordinate \
	--VALIDATION_STRINGENCY LENIENT \
	--CREATE_INDEX false

samtools index ${output_dir}/${sample_name}_bwa_sorted.bam
samtools idxstats ${output_dir}/${sample_name}_bwa_sorted.bam > ${output_dir}/${sample_name}_bwa_sorted.bam.stat

echo "Alignment done"

# (3) Add read group
echo "Add read group start"

picard AddOrReplaceReadGroups \
	-I ${output_dir}/${sample_name}_bwa_sorted.bam \
	-O ${output_dir}/${sample_name}_bwa_sorted_addrg.bam \
	--RGID ${sample_name} \
	--RGLB ${sample_name} \
	--RGPL illumina \
	--RGPU unit1 \
	--RGSM ${sample_name}

samtools index ${output_dir}/${sample_name}_bwa_sorted_addrg.bam
samtools idxstats ${output_dir}/${sample_name}_bwa_sorted_addrg.bam > ${output_dir}/${sample_name}_bwa_sorted_addrg.bam.stat

echo "Add read group done"

# (4) Mark duplicates reads
echo "Mark duplicates start"

picard MarkDuplicates \
	-I ${output_dir}/${sample_name}_bwa_sorted_addrg.bam \
	-O ${output_dir}/${sample_name}_bwa_sorted_addrg_mkdup.bam \
	-M ${output_dir}/${sample_name}_dup.matrix \
	--COMPRESSION_LEVEL 1 \
	--ASSUME_SORTED true \
	--VALIDATION_STRINGENCY LENIENT \
	--REMOVE_DUPLICATES false

samtools index ${output_dir}/${sample_name}_bwa_sorted_addrg_mkdup.bam
samtools idxstats ${output_dir}/${sample_name}_bwa_sorted_addrg_mkdup.bam > ${output_dir}/${sample_name}_bwa_sorted_addrg_mkdup.bam.stat

echo "Mark duplicates done"

# (5) BQSR
echo "BQSR start"

gatk BaseRecalibrator \
	-I ${output_dir}/${sample_name}_bwa_sorted_addrg_mkdup.bam \
	-R ${ref_genome} \
	--known-sites ${known_sites} \
	-O ${output_dir}/${sample_name}_recal_data.table

gatk ApplyBQSR \
   -R ${ref_genome} \
   -I ${output_dir}/${sample_name}_bwa_sorted_addrg_mkdup.bam \
   --bqsr-recal-file ${output_dir}/${sample_name}_recal_data.table \
   -O ${output_dir}/${sample_name}_bwa_sorted_addrg_mkdup_bqsr.bam

samtools index ${output_dir}/${sample_name}_bwa_sorted_addrg_mkdup_bqsr.bam
samtools idxstats ${output_dir}/${sample_name}_bwa_sorted_addrg_mkdup_bqsr.bam > ${output_dir}/${sample_name}_bwa_sorted_addrg_mkdup_bqsr.bam.stat

echo "BQSR done"

# (6) Remove intermediated files
rm -f ${output_dir}/${sample_name}_bwa_sorted.bam*
rm -f ${output_dir}/${sample_name}_bwa_sorted_addrg.bam*
rm -f ${output_dir}/${sample_name}_bwa_sorted_addrg_mkdup.bam*

