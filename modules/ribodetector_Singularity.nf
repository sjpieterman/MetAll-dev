nextflow.enable.dsl=2

process RIBODETECTOR {
    tag "$sample_id"
    label 'process_high'
    container "${baseDir}/ribodetector.sif"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_ribofree_{1,2}.fq.gz"), emit: reads, optional: true
    tuple val(sample_id), path("${sample_id}_rrna_{1,2}.fq.gz"), emit: rrna_reads, optional: true
    path "${sample_id}_ribodetector.log"                           , emit: log

    publishDir "${params.outdir}/ribodetector", mode: 'copy'

    script:
    def output_args = ""
    if (params.rd_model == "both") {
        output_args = "-o ${sample_id}_ribofree_1.fq.gz ${sample_id}_ribofree_2.fq.gz --rrna ${sample_id}_rrna_1.fq.gz ${sample_id}_rrna_2.fq.gz"
    } else {
        output_args = "-o ${sample_id}_ribofree_1.fq.gz ${sample_id}_ribofree_2.fq.gz"
    }
    """
    # Running RiboDetector with parameterized configuration
    ribodetector \\
        -t ${params.rd_threads} \\
        -l ${params.rd_len} \\
        -i ${reads[0]} ${reads[1]} \\
        -e ${params.rd_model} \\
        -m ${params.rd_gpu_mem} \\
        ${output_args} \\
        2> ${sample_id}_ribodetector.log
    """
}