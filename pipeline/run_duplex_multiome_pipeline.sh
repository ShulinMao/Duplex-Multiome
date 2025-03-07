#!/bin/bash
#SBATCH --partition=bch-compute	 			# queue to be used
#SBATCH --time=6-00:00:00	 			# Running time (in hours-minutes-seconds)
#SBATCH --job-name=UMB1117_Duplex_Multiome_Heart			# Job name
#SBATCH --mail-type=BEGIN,END,FAIL 		# send and email when the job begins, ends or fails
#SBATCH --mail-user=shulin_mao@fas.harvard.edu	 	# Email address to send the job status
#SBATCH --output=output.txt 			# Name of the output file
#SBATCH --mem=8G
#SBATCH --nodes=1				# Number of compute nodes

source /lab-share/Gene-Lee-e2/Public/home/shulin/softwares/anaconda3/etc/profile.d/conda.sh
conda activate /lab-share/Gene-Lee-e2/Public/home/shulin/softwares/mambaforge/envs/tn5
source /programs/biogrids.shrc

snakefile=/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/scripts/pipeline_preprocessing/Snakefile.Duplex_Multiome
profile=/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/scripts/slurm_config

# basic info of samples
dir=/lab-share/Gene-Choudhury-e2/Public/Shulin/Duplex_multiome/data

sample_name=UMB1117
base_dir=${dir}/${sample_name}/single_cell && mkdir -p ${base_dir}

cellranger_dir=/lab-share/Gene-Walsh-e2/Public/indData/AK/Multiome-strand-02212022/Novogene_04052024/1117_H/outs
bulk_bam=/lab-share/Gene-Walsh-e2/Public/indData/AK/bulk_bam/UMB1117/UMB1117_bwa_sorted_addrg_mkdup_bqsr.bam
het_germline_mutations=/lab-share/Gene-Choudhury-e2/Public/Shulin/Duplex_multiome/data/UMB1117/bulk/germline_mutation/UMB1117_germline_het_variants.vcf.gz
bulk_bam_list=${base_dir}/bulk_bam_list.txt

go_dir=/home/ch241780/go/bin
reference_genome=/lab-share/Gene-Lee-e2/Public/home/shulin/ref/hg38/genome.fa
known_snp=/lab-share/Gene-Lee-e2/Public/home/shulin/ref/hg38_gnomad312_genome_snp_AF0.01_sorted.vcf.gz
masked_region=/lab-share/Gene-Lee-e2/Public/home/shulin/ref/mask/hg38_mask.bed
max_VAF_clonal_mutation=0.1
autosomes_bed=/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/temp/hg38_autosomes.bed

# print cell barcode list
sed "s/,/\t/g" ${cellranger_dir}/per_barcode_metrics.csv | awk '{if ($4==1) print}' | awk '{print $1}' > ${base_dir}/per_barcode_metrics.cells

# generate strand barcode
echo -e "ATGTCCAG\nCGACGTCA\nGCTATAGC\nTACGAGTT" > ${base_dir}/strand1
echo -e "ACACCTAA\nCGTTTGGG\nGACAAACC\nTTGGGCTT" > ${base_dir}/strand2

# generate config file for Snakefile
echo "sample_name: \"${sample_name}\"" > ${base_dir}/config.split_bam.yaml
echo "cellranger_dir: \"${cellranger_dir}\"" >> ${base_dir}/config.split_bam.yaml
echo "base_dir: \"${base_dir}\"" >> ${base_dir}/config.split_bam.yaml
echo "cell_barcode_list: \"${base_dir}/per_barcode_metrics.cells\"" >> ${base_dir}/config.split_bam.yaml
echo "go_dir: \"${go_dir}\"" >> ${base_dir}/config.split_bam.yaml
echo "reference_genome: \"${reference_genome}\"" >> ${base_dir}/config.split_bam.yaml
echo "bulk_bam: \"${bulk_bam}\"" >> ${base_dir}/config.split_bam.yaml
echo "known_snp: \"${known_snp}\"" >> ${base_dir}/config.split_bam.yaml
echo "strand1_bc: \"${base_dir}/strand1\"" >> ${base_dir}/config.split_bam.yaml
echo "strand2_bc: \"${base_dir}/strand2\"" >> ${base_dir}/config.split_bam.yaml
echo "bulk_bam_list: \"${bulk_bam_list}\"" >> ${base_dir}/config.split_bam.yaml
echo "masked_region: \"${masked_region}\"" >> ${base_dir}/config.split_bam.yaml
echo "max_VAF_clonal_mutation: \"${max_VAF_clonal_mutation}\"" >> ${base_dir}/config.split_bam.yaml
echo "het_germline_mutations:  \"${het_germline_mutations}\"" >> ${base_dir}/config.split_bam.yaml
echo "autosomes_bed: \"${autosomes_bed}\"" >> ${base_dir}/config.split_bam.yaml


mkdir ${dir}/${sample_name}/tmp
configfile=${base_dir}/config.split_bam.yaml
directory=${base_dir}

snakemake -s ${snakefile} --profile ${profile} --configfile ${configfile} --directory ${directory} --unlock # unlock a directory if there was a kill signal
#snakemake -s ${snakefile} --profile ${profile} --configfile ${configfile} --directory ${directory} -n # dry run
snakemake -s ${snakefile} --profile ${profile} --configfile ${configfile} --directory ${directory} --jobs 300 --default-resources "tmpdir='${dir}/${sample_name}/tmp'"

