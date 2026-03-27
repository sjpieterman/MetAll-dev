nextflow.enable.dsl=2

process FASTQC {
    tag "$sample_id"
    label 'process_low'

    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(sample_id), path(reads)
    val suffix

    output:
    path("*.html"), emit: html
    path("*.zip") , emit: zip
    path("*_fastqc"), emit: qc_dir, optional: true

    publishDir "${params.outdir}/fastqc${suffix}", mode: 'copy', overwrite: true

    script:
    """
    fastqc -o . -t ${task.cpus} --extract ${reads}
    """
}