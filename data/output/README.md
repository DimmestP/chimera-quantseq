# Output folder
This folder holds the output from the nextflow quantseq rev pipeline (available at https://github.com/DimmestP/chimera-quantseq) ran on the 3'UTR motif construct data set (raw reads available on datastore bigdata/fastq/EdWallace-030521-data/). 


## multiqc.html

Output from the MultiQC program for the entire chimera project QuantSeq dataset. 
It combines alignment and initial FastQC results.

## counts
Output from the FeatureCounts program counting all reads mapped to the full transcript (5'UTR- ORF-3'UTR) of all genes in chimera project QuantSeq dataset.

## bedgraph
Output from running bedtools to make bedgraphs showing coverage of aligned reads

## alignments
Output BAM files from running the HISAT2 aligner.

## FastQC
Output from running FastQC on the raw fastq files (before adapter trimming).
