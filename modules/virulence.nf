nextflow.enable.dsl=2

process SALMON_VIRULENCE_INDEX {
    tag "$fasta"
    label 'process_medium'
    container 'quay.io/biocontainers/salmon:1.10.3--h45fbf2d_5'

    input:
    path fasta

    output:
    path "salmon_virulence_index", emit: index

    script:
    """
    salmon index -t ${fasta} -i salmon_virulence_index
    """
}

process SALMON_VIRULENCE_QUANT {
    tag "$sample_id"
    label 'process_medium'
    container 'quay.io/biocontainers/salmon:1.10.3--h45fbf2d_5'

    input:
    tuple val(sample_id), path(reads)
    path index

    output:
    tuple val(sample_id), path("${sample_id}_salmon_virulence"), emit: results
    path "${sample_id}_salmon_virulence/quant.sf", emit: quant

    publishDir "${params.outdir}/virulence", mode: 'copy'

    script:
    """
    salmon quant \\
        -i ${index} \\
        -l A \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -p ${task.cpus} \\
        -o ${sample_id}_salmon_virulence \\
        --validateMappings
    """
}
