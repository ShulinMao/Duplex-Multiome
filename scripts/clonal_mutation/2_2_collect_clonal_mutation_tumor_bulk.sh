sample_list=("COLO829BLT50_rep1" "COLO829BLT50_rep2")

for case_id in "${sample_list[@]}"
do
	echo $case_id
	dir=/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/${case_id}

	while read cell; do  
		awk -v cell_id="$cell" '{if($0!~/^#/){printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, cell_id, $4, $5, $6, $7, $8, $9, $10)}}' ${dir}/single_cell/call_variants/filter_germline_mutation/a2s0_tumor_bulk/a2s0_filterGermline.${cell}.vcf
	done < ${dir}/single_cell/per_barcode_metrics.cells > ${dir}/single_cell/call_variants/filter_germline_mutation/a2s0_tumor_bulk.vcf
done
