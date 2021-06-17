1. run generate_genome.py to create modified fasta and gtf files from S.cerevisiae R64 genome
2. S.cerevisiae genome fasta and GTF: https://www.ncbi.nlm.nih.gov/assembly/GCF_000146045.2/ (manually removed last line from original GTF file containing "###" only)
3. packages needed:
	gff3sort https://anaconda.org/bioconda/gff3sort
	samtools
4. plasmid_features.csv - plasmid features extracted from original snapgene files
5. the script generates two sets of GTF files: not sorted (*_unsorted.gtf) and sorted (*.gtf). Sorted is probably needed for some tools
