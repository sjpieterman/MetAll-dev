process TRIMGALORE {
    tag "$sample_id"
    label 'process_low'

    container "quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.fq.gz")      , emit: reads
    path "*trimming_report.txt"                , emit: report
    path "*_fastqc.{zip,html}"                 , emit: fastqc, optional: true

    publishDir "${params.outdir}/trimgalore", mode: 'copy'

    script:
    // Build the command piece by piece
    def args = []
    args << "--quality ${params.tg_quality}"
    args << "--length ${params.tg_min_length}"
    args << "--stringency ${params.tg_stringency}"
    args << "--e ${params.tg_error_rate}"
    
    if (params.tg_fastqc) args << "--fastqc"
    if (params.tg_adapter) args << "--adapter ${params.tg_adapter}"
    
    // Only add clipping if value is greater than 0
    if (params.tg_clip_r1 > 0) args << "--clip_R1 ${params.tg_clip_r1}"
    if (params.tg_clip_r2 > 0) args << "--clip_R2 ${params.tg_clip_r2}"
    if (params.tg_three_prime_clip_r1 > 0) args << "--three_prime_clip_R1 ${params.tg_three_prime_clip_r1}"
    if (params.tg_three_prime_clip_r2 > 0) args << "--three_prime_clip_R2 ${params.tg_three_prime_clip_r2}"

    """
    trim_galore \\
        --cores ${task.cpus} \\
        --paired \\
        --gzip \\
        --basename ${sample_id} \\
        ${args.join(' ')} \\
        ${reads}
    """
}