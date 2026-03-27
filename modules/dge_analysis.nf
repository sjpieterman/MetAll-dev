nextflow.enable.dsl=2


process DGE_ANALYSIS {
    tag "${comparison_name}"
    label 'process_medium'
    publishDir "${params.outdir}/dge_analysis/${comparison_name}", mode: 'copy'

    container 'docker://python:3.10-slim'
    // For Singularity, use:
    // container "${params.singularity_cache}/dge_analysis.sif"

    input:
    path counts_file
    path metadata_file
    val control_group
    val treatment_group
    val comparison_name

    output:
    path "dge_results.csv", emit: results
    path "volcano_plot.html", emit: volcano
    path "ma_plot.html", emit: ma_plot
    path "pca_plot.html", emit: pca
    path "pca_scree_plot.html", emit: scree
    path "heatmap.html", emit: heatmap
    path "dge_summary.json", emit: summary
    path "top50_degs.csv", emit: top50
    path "versions.yml", emit: versions

    script:
    def p_thresh = params.dge_p_threshold ?: 0.05
    def fc_thresh = params.dge_fc_threshold ?: 1.0
    """
    # Install required packages if not in container
    pip3 install --quiet numpy pandas scipy scikit-learn plotly kaleido 2>/dev/null || true

    # Run DGE analysis
    python3 ${projectDir}/bin/dge_analysis.py \\
        --counts ${counts_file} \\
        --metadata ${metadata_file} \\
        --control_group ${control_group} \\
        --treatment_group ${treatment_group} \\
        --output_dir . \\
        --p_threshold ${p_thresh} \\
        --fc_threshold ${fc_thresh}

    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        scipy: \$(python3 -c "import scipy; print(scipy.__version__)")
        scikit-learn: \$(python3 -c "import sklearn; print(sklearn.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch dge_results.csv
    touch volcano_plot.html
    touch ma_plot.html
    touch pca_plot.html
    touch pca_scree_plot.html
    touch heatmap.html
    touch dge_summary.json
    touch top50_degs.csv
    touch versions.yml
    """
}


/*
 * Multi-comparison DGE Analysis
 * Allows multiple treatment vs control comparisons
 */
process DGE_MULTI_COMPARISON {
    tag "multi_comparison"
    label 'process_medium'
    publishDir "${params.outdir}/dge_analysis/multi_comparison", mode: 'copy'

    container 'docker://python:3.10-slim'

    input:
    path counts_file
    path metadata_file
    path comparison_table  // CSV with columns: comparison_name, control_group, treatment_group

    output:
    path "*/dge_results.csv", emit: all_results
    path "*/volcano_plot.html", emit: all_volcanos
    path "*/pca_plot.html", emit: all_pcas
    path "multi_comparison_summary.csv", emit: summary

    script:
    """
    # Install required packages
    pip3 install --quiet numpy pandas scipy scikit-learn plotly kaleido 2>/dev/null || true

    # Read comparison table and run each comparison
    python3 -c "
import pandas as pd
import subprocess
import os

comparisons = pd.read_csv('${comparison_table}')

summary_data = []
for idx, row in comparisons.iterrows():
    comp_name = row['comparison_name']
    ctrl = row['control_group']
    treat = row['treatment_group']

    os.makedirs(comp_name, exist_ok=True)

    cmd = [
        'python3', '${projectDir}/bin/dge_analysis.py',
        '--counts', '${counts_file}',
        '--metadata', '${metadata_file}',
        '--control_group', ctrl,
        '--treatment_group', treat,
        '--output_dir', comp_name,
        '--p_threshold', '${params.dge_p_threshold ?: 0.05}',
        '--fc_threshold', '${params.dge_fc_threshold ?: 1.0}'
    ]

    subprocess.run(cmd, check=True)

    # Read summary
    import json
    with open(f'{comp_name}/dge_summary.json') as f:
        summary = json.load(f)

    summary_data.append({
        'comparison': comp_name,
        'control': ctrl,
        'treatment': treat,
        'total_genes': summary['summary_statistics']['total_genes'],
        'significant_genes': summary['summary_statistics']['significant_genes'],
        'upregulated': summary['summary_statistics']['upregulated'],
        'downregulated': summary['summary_statistics']['downregulated']
    })

# Save multi-comparison summary
pd.DataFrame(summary_data).to_csv('multi_comparison_summary.csv', index=False)
"
    """
}
