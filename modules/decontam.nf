process IDENTIFY_CONTAMINANTS {
    label 'process_low'
    container 'quay.io/biocontainers/bioconductor-decontam:1.30.0--r45hdfd78af_1'

    input:
    path(reports)
    path(metadata)
    val(neg_col)
    val(method)
    val(threshold)

    output:
    path "contaminants.txt", emit: contaminants

    script:
    """
    decontam_kraken.R ${metadata} "${neg_col}" "${method}" "${threshold}" ${reports}
    """
}

process APPLY_DECONTAM {
    tag "$sample_id"
    label 'process_low'
    container 'python:3.9-slim'

    input:
    tuple val(sample_id), path(report)
    path(contaminants)

    output:
    tuple val(sample_id), path("${sample_id}.decontam.report"), emit: report

    publishDir "${params.outdir}/kraken2_decontam", mode: 'copy'

    script:
    """
    fix_kraken_report.py ${report} ${contaminants} ${sample_id}.decontam.report
    """
}
