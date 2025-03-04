source /programs/biogrids.shrc

sample_list=("ABNCR6D" "UMB1932")
#("Brain1901106" "UMB1278" "UMB1465" "UMB1864" "UMB5451" "UMB5657" "UMB4638" "UMB5823")

for case_id in "${sample_list[@]}"
do

	#case_id=UMB5657
	dir=/lab-share/Gene-Lee-e2/Public/home/shulin/Duplex_multiome/data/${case_id}
	het_germline_mutation=${dir}/bulk/germline_mutation/${case_id}_germline_het_variants.vcf.gz

	mkdir ${dir}/single_cell/call_variants/single_strand_damage

	# summarize all a2s0 calls
	while read cell; do 
		awk -v cell_id="$cell" '{if($0!~/^#/){printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, cell_id, $4, $5, $6, $7, $8, $9, $10)}}' ${dir}/single_cell/call_variants/Dan/a2s0/${cell}.vcf
	done < ${dir}/single_cell/per_barcode_metrics.cells > ${dir}/single_cell/call_variants/single_strand_damage/a2s0.txt

	# sort
	# we need a header of a vcf file, so I randomly copy it from a file, $(head -n 1 ${dir}/single_cell/per_barcode_metrics.cells).vcf 
	bcftools view -h ${dir}/single_cell/call_variants/Dan/a2s0/$(head -n 1 ${dir}/single_cell/per_barcode_metrics.cells).vcf | cat - ${dir}/single_cell/call_variants/single_strand_damage/a2s0.txt | bcftools sort - -Oz -o ${dir}/single_cell/call_variants/single_strand_damage/a2s0.vcf.gz
	tabix ${dir}/single_cell/call_variants/single_strand_damage/a2s0.vcf.gz

	# intersect with high quality het germline mutations
	bcftools isec -n=2 -w1 ${dir}/single_cell/call_variants/single_strand_damage/a2s0.vcf.gz ${het_germline_mutation} -Oz -o ${dir}/single_cell/call_variants/single_strand_damage/a2s0_germline_mutation.vcf.gz
	tabix ${dir}/single_cell/call_variants/single_strand_damage/a2s0_germline_mutation.vcf.gz

	#rm -rf ${dir}/single_cell/call_variants/single_strand_damage/a2s0.txt
	#rm -rf ${dir}/single_cell/call_variants/single_strand_damage/a2s0.vcf.gz*

done

