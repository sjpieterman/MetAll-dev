nextflow.enable.dsl=2

process FASTQC {
    tag "$sample_id"
    label 'process_low'

    container { 
        params.use_rastqc ? params.rastqc_container :
        'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0' 
    }

    input:
    tuple val(sample_id), path(reads)
    val suffix

    output:
    path("*.html"), emit: html
    path("*.zip") , emit: zip, optional: true
    path("*.txt") , emit: txt, optional: true
    path("*.json"), emit: json, optional: true
    path("*_fastqc"), emit: qc_dir, optional: true

    publishDir "${params.outdir}/fastqc${suffix}", mode: 'copy', overwrite: true

    script:
    if (params.use_rastqc) {
        """
        # rastqc processes each file and outputs to the current directory
        for f in ${reads}; do
            # Check if file is empty (even if gzipped)
            if [ "\$(zcat "\$f" | head -n 4 | wc -l)" -eq 0 ]; then
                echo "Skipping empty file: \$f"
                continue
            fi

            ${params.rastqc_bin} -t ${task.cpus} --multiqc-json -o . "\$f"
            
            # RastQC uses filename minus extensions (e.g., .fastq.gz, .fq.gz, etc.)
            prefix=\$(basename "\$f" | sed 's/\\.fastq\\.gz\$//; s/\\.fq\\.gz\$//; s/\\.fastq\$//; s/\\.fq\$//')
            zip_file="\${prefix}_fastqc.zip"
            
            if [ -f "\$zip_file" ]; then
                # Create a temporary directory for fixing
                fix_dir="fix_\${prefix}"
                mkdir -p "\$fix_dir"
                unzip -q -o "\$zip_file" -d "\$fix_dir"
                
                # Standardize structure: all files directly in fix_dir first
                find "\$fix_dir" -type f -exec mv -t "\$fix_dir" {} + 2>/dev/null || true
                find "\$fix_dir" -mindepth 1 -type d -exec rm -rf {} + 2>/dev/null || true
                
                # Fix PASS/WARN/FAIL case to match what MultiQC expects (lowercase)
                # It crashes on 'WARN' with KeyError in MultiQC 1.19
                # We'll use a more robust sed that preserves the tab separator
                find "\$fix_dir" -type f -name "fastqc_data.txt" -exec sed -i 's/\\tPASS/\\tpass/g; s/\\tWARN/\\twarn/g; s/\\tFAIL/\\tfail/g' {} +
                find "\$fix_dir" -type f -name "summary.txt" -exec sed -i 's/^PASS/pass/; s/^WARN/warn/; s/^FAIL/fail/' {} +
                
                # Re-create the standard FastQC ZIP structure: a single subdirectory <prefix>_fastqc/
                mkdir -p "\$fix_dir/\${prefix}_fastqc"
                find "\$fix_dir" -maxdepth 1 -type f -exec mv -t "\$fix_dir/\${prefix}_fastqc" {} + 2>/dev/null || true
                
                # Re-zip with the standard structure
                rm "\$zip_file"
                (cd "\$fix_dir" && zip -q -r "../\$zip_file" *)
                rm -rf "\$fix_dir"
            fi
        done
        """
    } else {
        """
        fastqc -o . -t ${task.cpus} --extract ${reads}
        """
    }
}