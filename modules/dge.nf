nextflow.enable.dsl=2

process DGE_R_ANALYSIS {
    tag "${comparison_name}"
    label 'process_medium'
    publishDir "${params.outdir}/dge_analysis/${comparison_name}", mode: 'copy'

    container 'quay.io/biocontainers/mulled-v2-88ad01a0658428173499f668f12143003e670d8a:afaa2c0fcd83014cb964175390099898f82877a5-0'

    input:
    path counts_file
    path metadata_file
    val control_group
    val treatment_group
    val comparison_name
    val tool
    val batch_method
    val batch_col
    val group_col
    val covariates

    output:
    path "dge_results.csv", emit: results
    path "volcano_plot.png", emit: volcano_png
    path "volcano_plot.html", emit: volcano_html, optional: true
    path "ma_plot.html", emit: ma_plot, optional: true
    path "pca_plot.png", emit: pca_png
    path "pca_plot.html", emit: pca_html, optional: true
    path "heatmap.html", emit: heatmap, optional: true
    path "dge_summary.json", emit: summary
    path "top50_degs.csv", emit: top50
    path "normalized_counts.tsv", emit: normalized_counts, optional: true
    path "dataset_qc/*", emit: dataset_qc, optional: true

    script:
    def batch_arg = batch_method != 'none' && batch_col && batch_col != 'null' ? "--batch_method ${batch_method} --batch_col ${batch_col}" : ""
    def covar_arg = covariates && covariates != 'null' ? "--covariates ${covariates}" : ""
    def proxy_env = params.http_proxy ? "export http_proxy=${params.http_proxy}; export https_proxy=${params.http_proxy};" : ""
    def install_env = params.dge_autoinstall ? "" : "export METALL_NO_R_AUTOINSTALL=1;"
    
    def go_arg = params.dge_run_go ? "--run_go" : ""
    def gsea_arg = params.dge_run_gsea ? "--run_gsea" : ""
    def biomart_arg = params.dge_use_biomart ? "--use_biomart" : ""
    
    """
    ${proxy_env}
    ${install_env}
    Rscript ${projectDir}/bin/run_dge.R \\
        --counts ${counts_file} \\
        --metadata ${metadata_file} \\
        --control ${control_group} \\
        --treatment ${treatment_group} \\
        --outdir . \\
        --tool ${tool} \\
        --group_col ${group_col} \\
        ${batch_arg} \\
        ${covar_arg} \\
        ${go_arg} \\
        ${gsea_arg} \\
        ${biomart_arg} \\
        --organism_db ${params.dge_organism_db} \\
        --ontology ${params.dge_ontology} \\
        --keytype ${params.dge_keytype} \\
        --p_threshold ${params.dge_p_threshold ?: 0.05} \\
        --fc_threshold ${params.dge_fc_threshold ?: 1.0}
    """
}
