

while read sample project tech tissue path; do
	echo $sample
	printf "${sample}\t" >> depths_IDT.tsv
	samtools depth -@ 9 -a -b /data/bbg/datasets/prominent/metadata/regions/data/pipeline_bedfile/pancancer/PanCancerPanel.original.hg38.chr.bed4.bed ${path}/sortbamduplexcons/${sample}.*.bam | awk '{sum+=$3; n++} END{print sum/n}' >> depths_IDT.tsv
	
done < IDT_samples.tsv
