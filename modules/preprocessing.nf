process SHEETMAKER {
    label 'process_low'

    input:
    val search_dir

    output:
    path "samplesheet.csv", emit: csv

    publishDir "${params.outdir}/preprocessing", mode: 'copy'

    script:
    """
    sheetmaker.sh ${search_dir} > samplesheet.csv
    """
}

process CHECK_GPU {
    label 'process_low'
    debug true

    script:
    """
    echo "--- NVIDIA HARDWARE CHECK ---"
    #nvidia-smi || echo "GPU NOT FOUND"
    """
}