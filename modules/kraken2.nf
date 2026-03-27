process KRAKEN2 {
    tag "$sample_id"
    label 'process_high'
    container 'quay.io/biocontainers/kraken2:2.17.1--pl5321h077b44d_0'

    input:
    tuple val(sample_id), path(reads)
    path db

    output:
    tuple val(sample_id), path("${sample_id}.kraken2.report.txt"), emit: report
    path "${sample_id}.kraken2.output.txt"                       , emit: output

    publishDir "${params.outdir}/kraken2", mode: 'copy'

    script:
    """
    kraken2 \\
        --db ${db} \\
        --threads ${task.cpus} \\
        --confidence ${params.k2_confidence} \\
        --minimum-base-quality ${params.k2_min_hit_groups} \\
        --paired \\
        --gzip-compressed \\
        --output ${sample_id}.kraken2.output.txt \\
        --report ${sample_id}.kraken2.report.txt \\
        ${reads}
    """
}