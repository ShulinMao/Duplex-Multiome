#!/bin/bash
#SBATCH --partition=bch-compute	 			# queue to be used
#SBATCH --time=1-00:00:00	 			# Running time (in hours-minutes-seconds)
#SBATCH --job-name=demo_Duplex_Multiome			# Job name
#SBATCH --mail-type=BEGIN,END,FAIL 		# send and email when the job begins, ends or fails
#SBATCH --mail-user=shulin_mao@fas.harvard.edu	 	# Email address to send the job status
#SBATCH --output=output.txt 			# Name of the output file
#SBATCH --mem=8G
#SBATCH --nodes=1				# Number of compute nodes

source /lab-share/Gene-Lee-e2/Public/home/shulin/softwares/anaconda3/etc/profile.d/conda.sh
conda activate /lab-share/Gene-Lee-e2/Public/home/shulin/softwares/mambaforge/envs/Duplex-Multiome_2

# -----------------------------------------------------------------------------
# basic info of samples
dir=$(pwd)
sample_name=Duplex_multiome_demo
base_dir=${dir}/${sample_name}/single_cell && mkdir -p ${base_dir}
tmp=${dir}/${sample_name}/tmp && mkdir -p ${tmp}

# cellranger output directory for the sample
cellranger_dir=${dir}/CellRanger_output

# bulk WGS bam and germline variants called from the WGS data
bulk_bam=${dir}/bulk_data/demo_bulk.bam
het_germline_mutations=${dir}/bulk_data/demo_germline_variants.vcf.gz
bulk_bam_list=${base_dir}/bulk_bam_list.txt # bulk WGS bam for other samples in the batch

# reference
reference_genome=/lab-share/Gene-Lee-e2/Public/home/shulin/ref/hg38/genome.fa
known_snp=/lab-share/Gene-Lee-e2/Public/home/shulin/ref/hg38_gnomad312_genome_snp_AF0.01_sorted.vcf.gz
masked_region=/lab-share/Gene-Lee-e2/Public/home/shulin/ref/mask/hg38_mask.bed
autosomes_bed=/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/temp/hg38_autosomes.bed

# Sankemake pipeline and slurm settings
snakefile=${dir}/../pipeline/pipeline_preprocessing/Snakefile.Duplex_Multiome
# please update your slurm system job submission settings in "cluster-config.yml" in the following folder
profile=${dir}/../pipeline/slurm_config

# Go dir for duplexTools
go_dir=/home/ch241780/go/bin

# cut-off for filtering out variants present in the bulk WGS bam
max_VAF_clonal_mutation=0.1 # only keep variants with <= 10% VAF shown in the bulk WGS bam

# cell barcode list
sed "s/,/\t/g" ${cellranger_dir}/per_barcode_metrics.csv | awk '{if ($4==1) print}' | awk '{print $1}' > ${base_dir}/per_barcode_metrics.cells

# strand barcode
echo -e "AGACTTTC\nCCGAGGCA\nGATGCAGT\nTTCTACAG" > ${base_dir}/strand1
echo -e "AGCTGCGT\nCAACCATC\nGTGGAGCA\nTCTATTAG" > ${base_dir}/strand2


# -----------------------------------------------------------------------------
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

# set up snakemake pipeline
configfile=${base_dir}/config.split_bam.yaml
directory=${base_dir}


# snakemake -s ${snakefile} --profile ${profile} --configfile ${configfile} --directory ${directory} -n # dry run

# if running on a machine with a slurm system
snakemake -s ${snakefile} --profile ${profile} --configfile ${configfile} --directory ${directory} --unlock # unlock a directory if there was a kill signal
snakemake -s ${snakefile} --profile ${profile} --configfile ${configfile} --directory ${directory} --jobs 10 --default-resources "tmpdir='${tmp}'"

# if running on a machine without a slurm system
# snakemake -s ${snakefile} --configfile ${configfile} --directory ${directory} --unlock # unlock a directory if there was a kill signal
# snakemake -s ${snakefile} --configfile ${configfile} --directory ${directory} --jobs 1 --default-resources "tmpdir='${tmp}'"
