nextflow.enable.dsl=2

process FASTQC {
    tag "$sample_id"
    label 'process_low'

    container { params.use_falco ? 'quay.io/biocontainers/falco:1.2.5--h077b44d_0' : 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0' }

    input:
    tuple val(sample_id), path(reads)
    val suffix

    output:
    path("*.html"), emit: html
    path("*.zip") , emit: zip, optional: true
    path("*.txt") , emit: txt, optional: true
    path("*_fastqc"), emit: qc_dir, optional: true

    publishDir "${params.outdir}/fastqc${suffix}", mode: 'copy', overwrite: true

    script:
    if (params.use_falco) {
        """
        # falco processes files and puts results in their own directories if -o is used
        # to match FastQC behavior closely, we'll run it on each file
        for f in ${reads}; do
            falco -t ${task.cpus} "\$f"
        done
        """
    } else {
        """
        fastqc -o . -t ${task.cpus} --extract ${reads}
        """
    }
}