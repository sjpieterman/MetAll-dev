process FASTP {
    tag "$sample_id"
    label 'process_medium'
    container 'quay.io/biocontainers/fastp:0.24.0--h81e0983_0'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_fastp_1.fq.gz"), path("${sample_id}_fastp_2.fq.gz"), emit: reads
    path "${sample_id}.fastp.json"                                                            , emit: json
    path "${sample_id}.fastp.html"                                                            , emit: html

    publishDir "${params.outdir}/fastp", mode: 'copy'

    script:
    """
    fastp \\
        -i ${reads[0]} \\
        -I ${reads[1]} \\
        -o ${sample_id}_fastp_1.fq.gz \\
        -O ${sample_id}_fastp_2.fq.gz \\
        -j ${sample_id}.fastp.json \\
        -h ${sample_id}.fastp.html \\
        --thread ${task.cpus} \\
        --detect_adapter_for_pe \\
        --trim_poly_g \\
        --trim_poly_x
    """
}
