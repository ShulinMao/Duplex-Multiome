#!/bin/bash

# software version control
source /programs/biogrids.shrc
export PICARD_X=2.26.10;GATK_X=4.3.0.0;SAMTOOLS_X=1.13

# path
input_dir=/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/Brain1901106/bulk
sample_name=Brain1901106
mkdir ${input_dir}/germline_mutation
output_dir=${input_dir}/germline_mutation
bulk_bam=/lab-share/Gene-Walsh-e2/Public/indData/AK/bulk_bam/Brain1901106/bam/Brain1901106_bwa_sorted_addrg_mkdup_bqsr.bam

ref_genome=/lab-share/Gene-Lee-e2/Public/home/shulin/ref/hg38/genome.fa
known_sites=/lab-share/Gene-Lee-e2/Public/home/shulin/ref/dbsnp_147_b38_common_all_20160601.vcf.gz

gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R ${ref_genome} \
   -I ${bulk_bam} \
   -O ${output_dir}/${sample_name}.g.vcf.gz \
   #-ERC GVCF

bcftools isec ${output_dir}/${sample_name}.g.vcf.gz ${known_sites} -n =2 -w 1 -O z -o ${output_dir}/${sample_name}_intersection_with_dbsnp.vcf.gz
tabix ${output_dir}/${sample_name}_intersection_with_dbsnp.vcf.gz

bcftools view -R /lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/temp/hg38_autosomes.bed -i 'FORMAT/GT=="0/1"' ${output_dir}/${sample_name}_intersection_with_dbsnp.vcf.gz -Oz -o ${output_dir}/${sample_name}_germline_het_variants.vcf.gz
tabix ${output_dir}/${sample_name}_germline_het_variants.vcf.gz

rm -f ${output_dir}/${sample_name}_intersection_with_dbsnp.vcf.gz*

