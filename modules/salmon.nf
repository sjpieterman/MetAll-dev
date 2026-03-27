nextflow.enable.dsl=2

process SALMON_INDEX {
    tag "$fasta"
    label 'process_medium'
    container 'quay.io/biocontainers/salmon:1.11.4--h7f96273_0'

    input:
    path fasta

    output:
    path "salmon_index", emit: index

    script:
    """
    salmon index -t ${fasta} -i salmon_index
    """
}

process SALMON_QUANT {
    tag "$sample_id"
    label 'process_medium'
    container 'quay.io/biocontainers/salmon:1.11.4--h7f96273_0'

    input:
    tuple val(sample_id), path(reads)
    path index
    path gtf

    output:
    tuple val(sample_id), path("${sample_id}_salmon"), emit: results
    path "${sample_id}_salmon/quant.sf", emit: quant

    publishDir "${params.outdir}/salmon", mode: 'copy'

    script:
    """
    salmon quant \\
        -i ${index} \\
        -l A \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -p ${task.cpus} \\
        -o ${sample_id}_salmon \\
        --validateMappings \\
        -g ${gtf}
    """
}
