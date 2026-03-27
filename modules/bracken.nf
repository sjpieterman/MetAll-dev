process BRACKEN {
    tag "$sample_id"
    container 'quay.io/biocontainers/bracken:3.1--h9948957_0'

    input:
    tuple val(sample_id), path(kraken_report)
    path db

    output:
    tuple val(sample_id), path("${sample_id}.bracken.tsv"), emit: report

    publishDir "${params.outdir}/bracken", mode: 'copy'

    script:
    """
    bracken \\
        -d ${db} \\
        -i ${kraken_report} \\
        -o ${sample_id}.bracken.tsv \\
        -r ${params.bracken_read_len} \\
        -l ${params.bracken_level} \\
        -t ${params.bracken_threshold}
    """
}