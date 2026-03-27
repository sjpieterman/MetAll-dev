nextflow.enable.dsl=2

process PARABRICKS_RNA {
    tag "$sample_id"
    label 'process_high'
    container "nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1"

    input:
    tuple val(sample_id), path(reads)
    path  index_dir
    path  ref_fasta
    path  gtf

    output:
    tuple val(sample_id), path("${sample_id}.bam")          , emit: bam
    tuple val(sample_id), path("${sample_id}.bam.bai")      , emit: bai
    path "${sample_id}_featurecounts.txt"                   , emit: counts
    path "${sample_id}_featurecounts.txt.summary"           , emit: summary

    publishDir "${params.outdir}/parabricks", mode: 'copy'

    script:
    def is_gz = reads[0].name.endsWith('.gz')
    def fq1   = is_gz ? "read_1.fq" : "${reads[0]}"
    def fq2   = is_gz ? "read_2.fq" : "${reads[1]}"

    """
    if [ "${is_gz}" = "true" ]; then
        gunzip -c ${reads[0]} > read_1.fq
        gunzip -c ${reads[1]} > read_2.fq
    fi

    # Accelerated Fastq to BAM for RNA (Parabricks 4.6.0-1)
    pbrun rna_fq2bam \\
        --ref ${ref_fasta} \\
        --genome-lib-dir ${index_dir} \\
        --in-fq ${fq1} ${fq2} \\
        --output-dir . \\
        --out-bam ${sample_id}.bam \\
        --num-gpus 1 \\
        --low-memory

    # Accelerated FeatureCounts
    pbrun featurecounts \\
        --bam ${sample_id}.bam \\
        --gtf ${gtf} \\
        --output ${sample_id}_featurecounts.txt \\
        --num-gpus 1
    """
}