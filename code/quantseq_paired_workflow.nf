/*
QuantSeq paired read analysis pipeline for the Chimera Project S.cer motif insertion experiment.
*/

/*
Define dataset-specific input parameters.
These still need documentation and flexibility.
To find where they are used, search the document for the name, 
e.g. "params.featurename" is used in the featureCounts call.
*/
params.read_1_adapter = 'TTTTTTTTTTTTTTTTTT'
params.read_adapters_1 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
params.read_adapters_2 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
params.read_adapters_3 = 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
params.read_2_adapter = 'AAAAAAAAAAAAAAAAAA'
params.index_dir = '../data/input/Scer_ref_genome/construct_integrated_genome/construct_genome_fastas/indexed_genome/'
params.index_prefix = '_sample_with_saccharomyces_cerevisiae_R64'
params.mRNAgff_dir = '../data/input/Scer_ref_genome/construct_integrated_genome/construct_genome_gffs/'
params.input_fq_dir = '/homes/wallacelab/datastore/wallace_rna/bigdata/fastq/EdWallace-030521-data/' 
params.output_dir = '/homes/wallacelab/datastore/wallace_rna/data/2021/07-July/Sam/chimera_quantseq_pipeline_output/'
params.featuretype = 'primary_transcript'
params.featurename = 'ID'
params.num_processes = 4

/* Flatten nested list of file names (used after grouping by sample) */

flattenFileList = {
    list_of_paired_files = it[1]
    flattened_file_list = list_of_paired_files.flatten()
    it[1] = flattened_file_list
    it
}

/* Extract sample code from file name */
extractSampleCode = {
    filename = it[0]
    sample_name = (filename =~ /\w\d*(?=_\w\d+_L)/)[0]
    it[0] = sample_name
    it
}

/*
Define the input fastq.gz files, pairing forward and reverse reads and grouping across lanes by sample anme
*/

multi_lane_input_fq = Channel
    .fromFilePairs("${params.input_fq_dir}/*_R{1,2}_001.fastq.gz", size: 2)
    .map(extractSampleCode)
    .groupTuple(size: 4)
    .map(flattenFileList)

process combineLanesAcrossSamples {
    errorStrategy 'retry'
    maxRetries 3
    tag "${sample_id}"
    input:
    set sample_id, file(seq) from multi_lane_input_fq

    output:
    tuple val(sample_id), file("${sample_id}_R*.fastq.gz") into input_fq

    """
    cat ${seq.findAll{it =~/_R1_/}.asType(nextflow.util.BlankSeparatedList)} > ${sample_id + '_R1.fastq.gz'}
    cat ${seq.findAll{it =~/_R2_/}.asType(nextflow.util.BlankSeparatedList)} > ${sample_id + '_R2.fastq.gz'}
   """
}

/* split input_fq into two separate channels */
input_fq
    .tap{input_fq_qc}
    .tap{input_fq_cut}

/*
Run FastQC to produce a quality control report for the input data for every sample
*/

process runFastQC{
    conda 'bioconda::fastqc=0.11.9'
    errorStrategy 'retry'
    maxRetries 3
    tag "${sample_id}"
    publishDir "${params.output_dir}/FastQC/${sample_id}", mode: 'copy', overwrite: true
    input:
        set sample_id, file(paired_sample_fq) from input_fq_qc

    output:
        file("${sample_id}_fastqc/*.zip") into fastqc_files
   

    """
    mkdir ${sample_id}_fastqc
    fastqc --outdir ${sample_id}_fastqc \
    -t ${params.num_processes} \
    ${paired_sample_fq}
    """
}

/*
Cut sequencing adapters from 3' end of gene
*/


process cutAdapters {
    conda 'bioconda::cutadapt=1.18'
    errorStrategy 'retry'
    maxRetries 3
    tag "${sample_id}"
    input:
        set sample_id, file(sample_fq) from input_fq_cut
    output:
        tuple val(sample_id), file("trim_*.fq") into cut_fq
    shell:
        """
        cutadapt --trim-n -O 1 -m 20 -a ${params.read_1_adapter} -A ${params.read_2_adapter}\
            -A ${params.read_adapters_1} -A ${params.read_adapters_2} -A ${params.read_adapters_3}\
            -a ${params.read_adapters_1} -a ${params.read_adapters_2} -a ${params.read_adapters_3}\
            -o trim_1.fq -p trim_2.fq -j ${params.num_processes} ${sample_fq[0]} ${sample_fq[1]}
        """
}

/*
Define the aligner indexes for each construct
*/

extract_sample_name = {
    sample_name = it =~ /(?<=\/)\w\d*(?=_)/
    sample_file_tuple = [sample_name[0],it]
    sample_file_tuple
}

indexed_genomes = Channel.fromPath( "${params.index_dir}*.ht2" )
    .map(extract_sample_name)
    .groupTuple(size: 8)

reads_genome_tuple = indexed_genomes
    .join(cut_fq)

/*
Align trimmed reads to the genome with hisat2
*/

process alignHisat2 {
    conda 'bioconda::hisat2=2.1.0'
    errorStrategy 'retry'
    maxRetries 3
    tag "${sample_id}"
    publishDir "${params.output_dir}/alignment/${sample_id}", pattern: "*.hisat2_summary.txt", mode: 'copy', overwrite: true
    input:
        tuple val(sample_id), path(index_ht2_parts), path(sample_fq) from reads_genome_tuple
    output:
        file("unaligned.fq") into unaligned_fq
        file("${sample_id}.hisat2_summary.txt") into alignment_logs
        tuple val(sample_id), file("aligned.sam") into aligned_sam
    shell:
        """
        hisat2 --version
        hisat2 -p ${params.num_processes} -k 2 \
            --no-spliced-alignment \
            --no-unal \
            --un unaligned.fq -x ${sample_id}${params.index_prefix} \
            -S aligned.sam \
	    -1 ${sample_fq[0]} -2 ${sample_fq[1]} \
            --summary-file ${sample_id}.hisat2_summary.txt --maxins 1500
        """
}

/*
Turn unsorted aligned samfiles into sorted indexed compressed bamfiles
*/

process samViewSort {
    conda 'bioconda::samtools=1.11'
    errorStrategy 'retry'
    maxRetries 3
    tag "${sample_id}"
    input:
        set val(sample_id), file(sample_sam) from aligned_sam
    output:
        tuple val(sample_id), file("aligned_sorted.bam"), \
            file("aligned_sorted.bam.bai") into aligned_sorted_bam
    shell:
        """
        samtools --version
        samtools view -b ${sample_sam} | samtools sort \
            -@ ${params.num_processes} -O bam -o aligned_sorted.bam -
        samtools index aligned_sorted.bam
        """
}

// Split channel for use in multiple downstream processes.
aligned_sorted_bam.into { bedgraph_bam; htscount_bam }

/*
Make bedgraphs showing coverage of aligned reads
*/

process makeBedgraphs {
    conda 'bioconda::bedtools=2.30.0'
    errorStrategy 'retry'
    maxRetries 3
    tag "${sample_id}"
    publishDir "${params.output_dir}/bedgraph/${sample_id}", mode: 'copy', overwrite: true
    input:
        tuple val(sample_id), file(sample_bam), file(sample_bam_bai) \
            from bedgraph_bam
    output:
        tuple file("plus.bedgraph.gz"), \
            file("minus.bedgraph.gz") into bedgraph
    shell:
        """
        bedtools --version
        bedtools genomecov -ibam ${sample_bam} -trackline -bga \
            -strand + | gzip > plus.bedgraph.gz
        bedtools genomecov -ibam ${sample_bam} -trackline -bga \
            -strand - | gzip > minus.bedgraph.gz
        """
}

/*
Run rename Bam files by sample, for input into featureCounts.
*/

process renameBamSample {
    errorStrategy 'retry'
    maxRetries 3
    tag "${sample_id}"
    input:
        tuple val(sample_id), file(sample_bam), file(sample_bam_bai) \
            from htscount_bam
    output:
        tuple val(sample_id),  file("${sample_id}_aln.bam") into sampleid_aln_bam
    shell:
        """
        ln -s ${sample_bam} ${sample_id}_aln.bam
        """
}

/*
Define the gffs for each construct  
*/


mRNAgff = Channel.fromPath("${params.mRNAgff_dir}*.gff")
          .map(extract_sample_name)

gff_bam_tuple = mRNAgff
                .join(sampleid_aln_bam)

/*
Run featureCounts to count aligned reads to genes for all processed samples
*/

process countAllmRNA {
    conda 'bioconda::subread=2.0.0'
    errorStrategy 'retry'
    maxRetries 3
    tag "${sample_id}"
    publishDir "${params.output_dir}/counts/${sample_id}", pattern:"*.txt", mode: 'copy', overwrite: true
    publishDir "${params.output_dir}/counts/${sample_id}", pattern:"*.summary", mode: 'copy', overwrite: true
    publishDir "${params.output_dir}/sorted_bam/${sample_id}", pattern:"*.bam",  mode: 'copy', overwrite: true
    input:
        tuple val(sample_id), file(mRNAgff), file(sampleid_bams) from gff_bam_tuple
    output:
        file("${sample_id}_counts.txt") into counts
        file("${sample_id}_counts.txt.summary") into counts_summary
        file("${sample_id}_aln.bam.featureCounts.bam") into feature_assigned_bam
    shell:
        """
        featureCounts -p -T ${params.num_processes} -s 2 -t ${params.featuretype} -g ${params.featurename} -a ${mRNAgff} -o "${sample_id}_counts.txt" ${sampleid_bams.join(" ")} -R BAM
        """
}

/*
Run multiQC to collate single quality control report across all samples.
*/

process runMultiQC{
    errorStrategy 'retry'
    maxRetries 3
    tag { "multiQC" }
    publishDir "${params.output_dir}", mode: 'copy', overwrite: true
    input:
        file ('*') from fastqc_files.collect()
        file ('*') from alignment_logs.collect()
    output:
        file('multiqc_report.html')

    """
    multiqc .
    """
}

