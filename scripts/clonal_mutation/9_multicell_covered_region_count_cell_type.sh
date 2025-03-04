#!/bin/bash
#SBATCH --partition=bch-compute	 			# queue to be used
#SBATCH --time=2-00:00:00	 			# Running time (in hours-minutes-seconds)
#SBATCH --job-name=covered_region_count			# Job name
#SBATCH --mail-type=END,FAIL 		# send and email when the job begins, ends or fails
#SBATCH --mail-user=shulin_mao@fas.harvard.edu	 	# Email address to send the job status
#SBATCH --output=output.txt 			# Name of the output file
#SBATCH --mem=16G

source /programs/biogrids.shrc
sample_list=("UMB4428" "Brain1901106_deeper" "UMB1278_deeper" "UMB5823_deeper" "UMB1864_deeper" "UMB4638" "UMB1465" "UMB5451" "UMB5657")
stringency="a2s0"

for sample in "${sample_list[@]}"
do
	echo $sample

	bed_file_dir=/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/${sample}/single_cell/read_families/a2s0_0minReadFamilyLength
	output_dir=/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/results/clonal_mutation/multicell_covered_region_celltype_ATAC_recovered
	cell_type_dir=/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/${sample}/cell_type/celltype_ATAC
	
	mkdir -p ${output_dir}

	cd ${cell_type_dir}
	for cell_type in $(ls *.txt | cut -d. -f1)
	do
		echo $cell_type
		cell_barcode_list=${cell_type_dir}/${cell_type}.txt

		rm -rf ${output_dir}/${sample}_${cell_type}_open_region_${stringency}.bed

		while read -r f1 cell f3; do
			cat ${bed_file_dir}/${cell}.families.analysis.calledSites.bed >> ${output_dir}/${sample}_${cell_type}_open_region_${stringency}.bed
		done < ${cell_barcode_list}

		sort -k 1,1 ${output_dir}/${sample}_${cell_type}_open_region_${stringency}.bed > ${output_dir}/${sample}_${cell_type}_open_region_${stringency}_sorted.bed

		bedtools genomecov -i ${output_dir}/${sample}_${cell_type}_open_region_${stringency}_sorted.bed -g /lab-share/Gene-Lee-e2/Public/home/shulin/ref/hg38/hg38_main_chr_length.txt > ${output_dir}/${sample}_${cell_type}_open_region_${stringency}_genomecov.txt

		bedtools genomecov -bg -i ${output_dir}/${sample}_${cell_type}_open_region_${stringency}_sorted.bed -g /lab-share/Gene-Lee-e2/Public/home/shulin/ref/hg38/hg38_main_chr_length.txt > ${output_dir}/${sample}_${cell_type}_open_region_${stringency}_coverage_count.bed

		rm -f ${output_dir}/${sample}_${cell_type}_open_region_${stringency}.bed ${output_dir}/${sample}_${cell_type}_open_region_${stringency}_sorted.bed
	
	done
done
