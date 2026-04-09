process MULTIQC {
    label 'process_low'
    container 'quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0'

    input:
    path 'qc_results/*', stageAs: 'results_??/*' 
    val suffix // e.g., "_raw" or ""

    output:
    path("multiqc_report${suffix}.html"), emit: report
    path("multiqc_data")                , emit: data, optional: true
    path("*_plots")                     , emit: plots, optional: true

    publishDir "${params.outdir}/multiqc${suffix}", mode: 'copy', overwrite: true

    script:
    """
    multiqc \\
        -f \\
        -n multiqc_report${suffix}.html \\
        -o . \\
        -p ${task.cpus} \\
        .
    """
}

