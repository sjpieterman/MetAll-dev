process KRONA {
    tag "$sample_id"
    container 'quay.io/biocontainers/krona:2.8.1--pl5321hdfd78af_1'

    input:
    tuple val(sample_id), path(kraken_report)

    output:
    path "${sample_id}.krona.html", emit: html

    publishDir "${params.outdir}/kraken2", mode: 'copy'

    script:
    """
    # Convert hierarchical Kraken2 report to Krona text format
    # This preserves the full taxonomy hierarchy
    python3 -c "
import sys

def convert_kreport(infile, outfile):
    with open(infile, 'r') as f:
        lines = f.readlines()
    
    with open(outfile, 'w') as f:
        hierarchy = []
        for line in lines:
            line = line.rstrip('\\n')
            if not line: continue
            parts = line.split('\\t')
            if len(parts) < 6: continue
            
            # taxon_reads is column 3 (index 2)
            try:
                taxon_reads = int(parts[2])
            except ValueError:
                continue
                
            if taxon_reads == 0:
                continue
            
            name_col = parts[5]
            stripped_name = name_col.lstrip()
            indent = (len(name_col) - len(stripped_name))
            depth = indent // 2
            
            if depth < len(hierarchy):
                hierarchy = hierarchy[:depth]
            
            hierarchy.append(stripped_name)
            
            # Format: count\tlevel1\tlevel2...
            f.write(str(taxon_reads) + '\\t' + '\\t'.join(hierarchy) + '\\n')

if __name__ == '__main__':
    convert_kreport('${kraken_report}', '${sample_id}.krona.txt')
"

    ktImportText \\
        -o ${sample_id}.krona.html \\
        ${sample_id}.krona.txt
    """
}