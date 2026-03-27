nextflow.enable.dsl=2


process STAR_ALIGN {
    tag "$sample_id"
    label 'process_high'
    // Updated to a robust biocontainers image
    container 'quay.io/biocontainers/star:2.7.3a--h5ca1c30_1'

    input:
    tuple val(sample_id), path(reads)
    path star_index
    path gtf

    output:
    tuple val(sample_id), path("${sample_id}.Aligned.sortedByCoord.out.bam"), emit: bam
    path "${sample_id}.Log.final.out"                                     , emit: log
    path "${sample_id}.SJ.out.tab"                                        , emit: sj
    path "${sample_id}.ReadsPerGene.out.tab"                              , emit: counts

    publishDir "${params.outdir}/star", mode: 'copy'

    script:
    """
    STAR \\
        --genomeDir ${star_index} \\
        --readFilesIn ${reads[0]} ${reads[1]} \\
        --readFilesCommand zcat \\
        --runThreadN ${task.cpus} \\
        --outFileNamePrefix ${sample_id}. \\
        --outSAMtype BAM SortedByCoordinate \\
        --sjdbGTFfile ${gtf} \\
        --quantMode GeneCounts
    """
}

process FEATURECOUNTS {
    tag "$sample_id"
    label 'process_medium'
    // Already using a good biocontainers image
    container 'quay.io/biocontainers/subread:2.0.6--he4a0461_0'

    input:
    tuple val(sample_id), path(bam)
    path gtf

    output:
    tuple val(sample_id), path("${sample_id}.featureCounts.txt")        , emit: counts
    tuple val(sample_id), path("${sample_id}.featureCounts.txt.summary"), emit: summary

    publishDir "${params.outdir}/featurecounts", mode: 'copy'

    script:
    """
    featureCounts \\
        -T ${task.cpus} \\
        -p \\
        -a ${gtf} \\
        -o ${sample_id}.featureCounts.txt \\
        ${bam}
    """
}

process MERGE_COUNTS {
    tag "merge_counts"
    label 'process_medium'
    container 'quay.io/biocontainers/pandas:1.5.2'

    input:
    path(fc_files)
    path(star_files)
    path(salmon_files)

    output:
    path "gene_counts_matrix.tsv", emit: matrix
    path "*_gene_counts_matrix.tsv", emit: all_matrices, optional: true

    publishDir "${params.outdir}/counts", mode: 'copy'

    script:
    def fc_arg = fc_files ? "--featurecounts ${fc_files}" : ""
    def star_arg = star_files ? "--star ${star_files}" : ""
    def salmon_arg = salmon_files ? "--salmon ${salmon_files}" : ""
    """
    python3 ${projectDir}/bin/merge_counts_enhanced.py \\
        ${fc_arg} \\
        ${star_arg} \\
        ${salmon_arg} \\
        --output gene_counts_matrix.tsv
    """
}