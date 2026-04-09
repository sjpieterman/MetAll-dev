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
    # Using perl because it is guaranteed to be in the Krona container
    perl -ne '
        @p = split /\\t/;
        next if \$p[2] == 0;
        \$name = \$p[5];
        \$stripped = \$name;
        \$stripped =~ s/^\\s+//;
        \$depth = (length(\$name) - length(\$stripped)) / 2;
        @hierarchy = @hierarchy[0..\$depth-1];
        push @hierarchy, \$stripped;
        print \$p[2] . "\\t" . join("\\t", @hierarchy) . "\\n";
    ' ${kraken_report} > ${sample_id}.krona.txt

    ktImportText \\
        -o ${sample_id}.krona.html \\
        ${sample_id}.krona.txt
    """
}