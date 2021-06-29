/*
QuantSeq paired read analysis pipeline for the Chimera Project S.cer motif insertion experiment.
*/

/*
Define dataset-specific input parameters.
These still need documentation and flexibility.
To find where they are used, search the document for the name, 
e.g. "params.featurename" is used in the featureCounts call.
*/
params.read_1_adapters_1 = 'AAAAAAAAAAAAAAAAAA'
params.read_1_adapters_2 = 'TTTTTTTTTTTTTTTTTT'
params.read_2_adapters_1 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
params.read_2_adapters_2 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
params.read_2_adapters_3 = 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
params.index_dir = '../data/input/Scer_ref_genome/'
params.index_prefix = 'saccharomyces_cerevisiae_R64'
params.mRNAgff = '../data/input/Scer_ref_genome/longest_three_prime_UTRs.gff'
params.input_fq_dir = '../data/input/EdWallace-030521-data/'
params.output_dir = '../data/output/'
params.featuretype = 'three_prime_UTR'
params.featurename = 'ID'
params.num_processes = 4

/* Flatten nested list of file names (used after grouping by sample) */

flattenFileList = {
    list_of_paired_files = it[1]
    flattened_file_list = list_of_paired_files.flatten()
    it[1] = flattened_file_list
    it
}

/*
Define the input fastq.gz files, pairing forward and reverse reads and grouping across lanes by sample anme
*/

multi_lane_input_fq = Channel
    .fromFilePairs("${params.input_fq_dir}/*_R{1,2}_001.fastq.gz", size: 2) {file -> (file =~ /\w\d*_\w\d+(?=_L)/)[0]} /* closure extracts sample name from file name */
    .groupTuple(size: 4) /* group lanes */
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
Define the aligner index and feature file (gff)
*/

Channel
        .fromPath("${params.index_dir}/${params.index_prefix}.*.ht2",
                  checkIfExists: true)
        .collect().set { index_ht2_parts }

mRNAgff = Channel.fromPath(params.mRNAgff)


/*
Run FastQC to produce a quality control report for the input data for every sample
*/

process runFastQC{
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
    errorStrategy 'retry'
    maxRetries 3
    tag "${sample_id}"
    input:
        set sample_id, file(sample_fq) from input_fq_cut
    output:
        tuple val(sample_id), file("trim_*.fq") into cut_fq
    shell:
        """
        cutadapt --trim-n -O 1 -m 20 -a ${params.read_1_adapters_1} -a ${params.read_1_adapters_2}\
	    -A ${params.read_2_adapters_1} -A ${params.read_2_adapters_2} -A ${params.read_2_adapters_3}\
            -o trim_1.fq -p trim_2.fq -j ${params.num_processes} ${sample_fq[0]} ${sample_fq[1]}
        """
}

/*
Align trimmed reads to the genome with hisat2
*/

process alignHisat2 {
    errorStrategy 'retry'
    maxRetries 3
    tag "${sample_id}"
    publishDir "${params.output_dir}/alignment/${sample_id}", pattern: '*.hisat2_summary.txt', mode: 'copy', overwrite: true
    input:
        set sample_id, file(sample_fq) from cut_fq
        file(index_ht2_parts) from index_ht2_parts
    output:
        file("unaligned.fq") into unaligned_fq
        file("${sample_id}.hisat2_summary.txt") into alignment_logs
        tuple val(sample_id), file("aligned.sam") into aligned_sam
    shell:
        """
        hisat2 --version
        hisat2 -p ${params.num_processes} -k 2 \
            --pen-cansplice 4 --pen-noncansplice 12 --min-intronlen 40  --max-intronlen 200 \
            --no-unal \
            --un unaligned.fq -x ${params.index_prefix} \
            -S aligned.sam \
	    -1 ${sample_fq[0]} -2 ${sample_fq[1]} \
            --summary-file ${sample_id}.hisat2_summary.txt
        """
}

/*
Turn unsorted aligned samfiles into sorted indexed compressed bamfiles
*/

process samViewSort {
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
        file("${sample_id}_aln.bam") into sampleid_aln_bam
    shell:
        """
        ln -s ${sample_bam} ${sample_id}_aln.bam
        """
}

/*
Run featureCounts to count aligned reads to genes for all processed samples
*/

process countAllmRNA {
    errorStrategy 'retry'
    maxRetries 3
    publishDir "${params.output_dir}/counts/", mode: 'copy'
    input:
        file(sampleid_bams) from sampleid_aln_bam.collect()
        file(mRNAgff)
    output:
        file("counts.txt") into counts
    shell:
        """
        featureCounts -T ${params.num_processes} -s 1 -t ${params.featuretype} -g ${params.featurename} -a ${mRNAgff} -o counts.txt ${sampleid_bams.join(" ")} 
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

