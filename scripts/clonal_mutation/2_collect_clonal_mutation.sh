sample_list=("ABNCR6D")
#("ABNCR6D" "AN06365")
#("Brain1901106_deeper" "UMB1278_deeper" "UMB5823_deeper" "UMB1864_deeper" "UMB4638" "UMB1465" "UMB5451" "UMB5657")
#("COLO829BLT50_rep1" "COLO829BLT50_rep2")

germline_VAF=0.1

for case_id in "${sample_list[@]}"
do
	echo $case_id
	dir=/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/${case_id}

	while read cell; do  
		awk -v cell_id="$cell" '{if($0!~/^#/){printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, cell_id, $4, $5, $6, $7, $8, $9, $10)}}' ${dir}/single_cell/call_variants/filter_germline_mutation/a2s0_AF${germline_VAF}/a2s0_filterGermline.${cell}.vcf
	done < ${dir}/single_cell/per_barcode_metrics.cells > ${dir}/single_cell/call_variants/filter_germline_mutation/a2s0_AF${germline_VAF}.vcf
done
