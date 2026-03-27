process RIBODETECTOR {
    tag "$sample_id"
    label 'process_high'

    // Use a direct string path. Nextflow handles the absolute path conversion.
    container "${baseDir}/ribodetector.sif"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_ribofree_{1,2}.fq.gz"), emit: reads
    path "${sample_id}_ribodetector.log"                           , emit: log

    publishDir "${params.outdir}/ribodetector", mode: 'copy'

    script:
        """
        # Running RiboDetector using your proven parameters
        ribodetector \\
            -t ${task.cpus} \\
            -l 151 \\
            -i ${reads[0]} ${reads[1]} \\
            -e rrna \\
            -m 20 \\
            -o ${sample_id}_ribofree_1.fq.gz ${sample_id}_ribofree_2.fq.gz \\
            2> ${sample_id}_ribodetector.log
        """
}