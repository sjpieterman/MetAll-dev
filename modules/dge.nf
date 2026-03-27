nextflow.enable.dsl=2

process DGE_R_ANALYSIS {
    tag "${comparison_name}"
    label 'process_medium'
    publishDir "${params.outdir}/dge_analysis/${comparison_name}", mode: 'copy'

    container 'quay.io/biocontainers/bioconductor-deseq2:1.40.0--r43hf17093f_0'

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
    path "pca_plot.png", emit: pca
    path "volcano_plot.png", emit: volcano

    script:
    def batch_arg = batch_method != 'none' && batch_col ? "--batch_method ${batch_method} --batch_col ${batch_col}" : ""
    def covar_arg = covariates != 'null' ? "--covariates ${covariates}" : ""
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
