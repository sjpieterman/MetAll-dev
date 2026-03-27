process SORTMERNA {
    tag "$sample_id"
    label 'process_high'
    container 'quay.io/biocontainers/sortmerna:4.3.6--h9ee0642_0'
    
    // If the process fails, don't stop the whole pipeline, just skip this sample
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(reads)
    path db_dir

    output:
    tuple val(sample_id), path("${sample_id}_non_rRNA_*.fq.gz"), emit: reads
    path "${sample_id}_sortmerna.log"                          , emit: log
    stdout emit: status // We'll use stdout to signal if it's a bad sample

    publishDir "${params.outdir}/sortmerna", mode: 'copy'

    script:
    // Construct dynamic arguments
    def args = []
    args << "--num_alignments ${params.smr_num_alignments}"
    args << "--e ${params.smr_evalue}"
    args << "--mismatch ${params.smr_mismatch}"
    args << "--id ${params.smr_coverage}"
    args << "--seed ${params.smr_rand_seed}"

    """
    # Find databases
    refs=\$(find ${db_dir} -name "*.fasta" | sed 's/^/--ref /' | tr '\\n' ' ')

    sortmerna \\
        \$refs \\
        --reads ${reads[0]} --reads ${reads[1]} \\
        --threads ${task.cpus} \\
        --workdir . \\
        ${args.join(' ')} \\
        --aligned ${sample_id}_rRNA \\
        --other ${sample_id}_non_rRNA \\
        --fastx \\
        --paired_in \\
        --out2 || true
    
    if [ -f ${sample_id}_non_rRNA_fwd.fq ] && [ -s ${sample_id}_non_rRNA_fwd.fq ]; then
        gzip ${sample_id}_non_rRNA*.fq
        echo "PASS" 
    else
        echo "${sample_id}"
        echo "" | gzip > ${sample_id}_non_rRNA_fwd.fq.gz
        echo "" | gzip > ${sample_id}_non_rRNA_rev.fq.gz
    fi

    [ -f sortmerna.log ] && mv sortmerna.log ${sample_id}_sortmerna.log || touch ${sample_id}_sortmerna.log
    """
}