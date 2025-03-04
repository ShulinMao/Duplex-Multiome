#!/bin/bash
#SBATCH --partition=bch-compute	 			# queue to be used
#SBATCH --time=06:00:00	 			# Running time (in hours-minutes-seconds)
#SBATCH --job-name=run_R_script		# Job name
#SBATCH --mail-type=BEGIN,END,FAIL 		# send and email when the job begins, ends or fails
#SBATCH --mail-user=shulin_mao@fas.harvard.edu	 	# Email address to send the job status
#SBATCH --output=output_callable_region_merge.txt 			# Name of the output file
#SBATCH --mem=2G
#SBATCH --nodes=1				# Number of compute nodes
#SBATCH --ntasks=1			# Number of threads/tasks on one node

source /programs/biogrids.shrc

sample_list=("Brain1901106_deeper" "UMB1278_deeper" "UMB1864_deeper" "UMB5823_deeper" "UMB4638" "UMB1465" "UMB5451" "UMB5657" "UMB4428")

# for sample in "${sample_list[@]}"
# do
# 	bed_file_dir=/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/${sample}/single_cell/call_variants/autosomal_callable_region
# 	cell_barcode_list_dir=/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/${sample}/cell_type/
# 	output_dir=/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/${sample}/PTA_results
# 
# 	mkdir -p ${output_dir}
# 
# 	#cell_type_list=("EN" "Oligodendrocyte")
# 	cell_type_list=("IN")
# 	for cell_type in "${cell_type_list[@]}"
# 	do
# 		cell_barcode_list=${cell_barcode_list_dir}/${cell_type}.txt
# 		
#     if [[ ! -f "$cell_barcode_list" ]]; then
#         continue  # Skip to the next cell type
#     fi
# 
# 		stringency_list=("a2s1" "a4s2" "a6s3" "a8s4" "a10s5")
# 		for stringency in "${stringency_list[@]}"
# 		do
# 			rm -rf ${output_dir}/${cell_type}_open_region_${stringency}.bed
# 			while read f1 f2 f3; do
# 				cat ${bed_file_dir}/${stringency}/${f2}.families.analysis.calledSites.bed >> ${output_dir}/${cell_type}_open_region_${stringency}.bed
# 			done < ${cell_barcode_list}
# 			wait
# 			sort -k1,1 -k2,2n ${output_dir}/${cell_type}_open_region_${stringency}.bed | bedtools merge -c 1 -o count -i - > ${output_dir}/${cell_type}_open_region_${stringency}_sorted_merged.bed
# 		done
# 	done
# done


# subtypes of neurons
for sample in "${sample_list[@]}"
do
	bed_file_dir=/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/${sample}/single_cell/call_variants/autosomal_callable_region
	cell_barcode_list_dir=/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/${sample}/cell_type/subtype
	output_dir=/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/${sample}/PTA_results

	mkdir -p ${output_dir}

	cell_type_list=("EN_deep_layer")
	for cell_type in "${cell_type_list[@]}"
	do
		cell_barcode_list=${cell_barcode_list_dir}/${cell_type}.txt
		
    if [[ ! -f "$cell_barcode_list" ]]; then
        continue  # Skip to the next cell type
    fi

		stringency_list=("a2s1" "a4s2" "a6s3" "a8s4" "a10s5")
		for stringency in "${stringency_list[@]}"
		do
			rm -rf ${output_dir}/${cell_type}_open_region_${stringency}.bed
			while read f1 f2 f3; do
				cat ${bed_file_dir}/${stringency}/${f2}.families.analysis.calledSites.bed >> ${output_dir}/${cell_type}_open_region_${stringency}.bed
			done < ${cell_barcode_list}
			wait
			sort -k1,1 -k2,2n ${output_dir}/${cell_type}_open_region_${stringency}.bed | bedtools merge -c 1 -o count -i - > ${output_dir}/${cell_type}_open_region_${stringency}_sorted_merged.bed
		done
	done
done
