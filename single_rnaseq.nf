nextflow.enable.dsl = 2

/*
 * pipeline input parameters
 */
params.study_id = "default_study"  // Study identifier, should be provided by user input
params.single_end_reads = "Data/${params.study_id}/*_pass.fastq"
params.transcriptome = "t_indxs/pao1_cnda.fa"
params.outdir = "Output/${params.study_id}"
params.kmer_size = 15  // Default value, can be overridden by user input

log.info """\
         R N A S E Q - N F   P I P E L I N E
         ===================================
         transcriptome: ${params.transcriptome}
         single-end reads: ${params.single_end_reads}
         outdir       : ${params.outdir}
         kmer size    : ${params.kmer_size}
         """
         .stripIndent()

/*
 * Define the `INDEX` process that creates a binary index
 * given the transcriptome file
 */
process INDEX {
    cpus 4
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path transcriptome

    output:
    path 'index'

    script:
    """
    salmon index --threads ${task.cpus} -t ${transcriptome} -i index -k ${params.kmer_size}
    """
}

/*
 * Define TrimGalore process for trimming single-end reads
 */
process TRIMGALORE_SINGLE {
    cpus 2
    tag "TrimGalore on $sample_id (single-end)"
    publishDir "${params.outdir}/${sample_id}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq"), path("*.txt")

    script:
    """
    trim_galore -o . ${reads}
    if [ "${reads.baseName}_trimmed.fq" != "${sample_id}_trimmed.fastq" ]; then
        mv ${reads.baseName}_trimmed.fq ${sample_id}_trimmed.fastq
    fi
    """
}

/*
 * Run Salmon to perform the quantification of expression using
 * the index and the matched read files
 */
process QUANT_SINGLE {
    cpus 4
    tag "Quantification on $sample_id (single-end)"
    publishDir "${params.outdir}/${sample_id}/quant", mode: 'copy'

    input:
    each index
    tuple val(sample_id), path(reads)

    output:
    path(sample_id)

    script:
    """
    salmon quant --threads ${task.cpus} -i ${index} --libType A -r ${reads} --validateMappings -o ${sample_id}
    """
}

/*
 * Run FastQC to check quality of reads files
 */
process FASTQC_SINGLE {
    cpus 2
    tag "FASTQC on $sample_id (single-end)"
    publishDir "${params.outdir}/${sample_id}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path("fastqc_${sample_id}_logs")

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

/*
 * Create a report using MultiQC for the quantification,
 * FastQC, and TrimGalore processes
 */
process MULTIQC {
    cpus 2
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path 'multiqc_inputs/*'

    output:
    path('multiqc_report.html')

    script:
    """
    multiqc multiqc_inputs
    """
}

workflow {
    single_end_reads_ch = Channel.fromPath("Data/${params.study_id}/fastq/*_pass.fastq", checkIfExists: true)
                            .map { file -> 
                                tuple(file.simpleName, file) 
                            }

    transcriptome_ch = Channel.fromPath(params.transcriptome, checkIfExists: true)

    index_ch = INDEX(transcriptome_ch)

    trimmed_single_end_reads_ch = TRIMGALORE_SINGLE(single_end_reads_ch)

    // Separate the outputs of TRIMGALORE_SINGLE
    trimmed_fq_ch = trimmed_single_end_reads_ch.map { tuple(it[0], it[1]) }
    trimmed_html_ch = trimmed_single_end_reads_ch.map { it -> it[2] }

    // Pass trimmed reads to FASTQC and QUANT
    fastqc_single_ch = FASTQC_SINGLE(trimmed_fq_ch)
    quant_single_ch = QUANT_SINGLE(index_ch, trimmed_fq_ch)

    // Collect the outputs for MultiQC
    multiqc_htmls_ch = trimmed_html_ch
    multiqc_quant_ch = quant_single_ch
    multiqc_fastqc_ch = fastqc_single_ch

    multiqc_ch = MULTIQC(quant_single_ch.mix(fastqc_single_ch, trimmed_html_ch).collect())
}


workflow.onComplete {

    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        Study ID    : ${params.study_id} 
        CPUs Used   : ${params.cpus} 
        """
        .stripIndent()

    // Save the summary message to a file
    def summaryFile = new File("${params.outdir}/pipeline_summary.txt")
    summaryFile.text = msg

    log.info(msg)

    log.info(workflow.success ? "\nDone! Open the following report in your browser --> ${params.outdir}/multiqc_report.html\n" : "Oops .. something went wrong")
}