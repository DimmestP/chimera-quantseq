rm /homes/sjhaynes/chimera-quantseq/data/output/sorted_bam/polyA_reads.bed

for f in /homes/sjhaynes/chimera-quantseq/data/output/sorted_bam/*/*.bam;
do
	echo $f;
	sample_name=$(expr match $f '.*\([ABE][0-9]\+\)_')
	bedtools intersect -a $f -b\
	../data/input/Scer_ref_genome/construct_integrated_genome/construct_genome_gffs/${sample_name}_sample_longest_full_ORF_with_constructs.gff\
	-bed | grep -w - | grep  mCherry >> /homes/sjhaynes/chimera-quantseq/data/output/sorted_bam/${sample_name}/${sample_name}_polyA_reads.bed
done;
