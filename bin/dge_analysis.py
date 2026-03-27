#!/usr/bin/env python3
"""
Differential Gene Expression Analysis Script
Standalone script for DGE analysis with comprehensive visualizations
"""

import os
import sys
import argparse
import json
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')


class DGEAnalyzer:
    """Differential Gene Expression Analysis Class"""

    def __init__(self, counts_file, metadata_file, control_group, treatment_group,
                 output_dir, p_threshold=0.05, fc_threshold=1.0):
        """
        Initialize DGE Analyzer

        Parameters:
        -----------
        counts_file : str
            Path to featureCounts output file
        metadata_file : str
            Path to metadata CSV file (sample, group)
        control_group : str
            Name of control group in metadata
        treatment_group : str
            Name of treatment group in metadata
        output_dir : str
            Directory for output files
        p_threshold : float
            P-value threshold for significance
        fc_threshold : float
            Log2 fold change threshold
        """
        self.counts_file = counts_file
        self.metadata_file = metadata_file
        self.control_group = control_group
        self.treatment_group = treatment_group
        self.output_dir = output_dir
        self.p_threshold = p_threshold
        self.fc_threshold = fc_threshold

        # Create output directory
        os.makedirs(output_dir, exist_ok=True)

        # Data containers
        self.counts_df = None
        self.metadata_df = None
        self.results_df = None
        self.control_samples = []
        self.treatment_samples = []

    def load_data(self):
        """Load counts and metadata"""
        print(f"Loading counts from: {self.counts_file}")

        # Read featureCounts output
        self.counts_df = pd.read_csv(self.counts_file, sep='\t', comment='#', skiprows=1)
        print(f"  Loaded {len(self.counts_df)} genes")

        # Read metadata
        print(f"Loading metadata from: {self.metadata_file}")
        self.metadata_df = pd.read_csv(self.metadata_file)
        print(f"  Loaded {len(self.metadata_df)} samples")

        # Get sample columns (after first 6 annotation columns in featureCounts)
        available_samples = self.counts_df.columns[6:].tolist()

        # Match samples from metadata
        control_meta = self.metadata_df[self.metadata_df['group'] == self.control_group]['sample'].tolist()
        treatment_meta = self.metadata_df[self.metadata_df['group'] == self.treatment_group]['sample'].tolist()

        # Find matching samples in counts file
        self.control_samples = [s for s in control_meta if s in available_samples]
        self.treatment_samples = [s for s in treatment_meta if s in available_samples]

        print(f"  Control samples ({len(self.control_samples)}): {self.control_samples}")
        print(f"  Treatment samples ({len(self.treatment_samples)}): {self.treatment_samples}")

        if len(self.control_samples) == 0 or len(self.treatment_samples) == 0:
            raise ValueError("No matching samples found between metadata and counts file!")

    def perform_dge(self):
        """Perform differential gene expression analysis"""
        print("\nPerforming DGE analysis...")

        gene_ids = self.counts_df.iloc[:, 0]
        counts_control = self.counts_df[self.control_samples].values
        counts_treatment = self.counts_df[self.treatment_samples].values

        results = []
        for i in range(len(gene_ids)):
            ctrl_vals = counts_control[i, :]
            treat_vals = counts_treatment[i, :]

            # Calculate means
            mean_ctrl = np.mean(ctrl_vals)
            mean_treat = np.mean(treat_vals)

            # Add pseudocount for log transformation
            mean_ctrl_pc = mean_ctrl + 1
            mean_treat_pc = mean_treat + 1

            # Log2 fold change
            log2fc = np.log2(mean_treat_pc / mean_ctrl_pc)

            # T-test
            if np.std(ctrl_vals) > 0 or np.std(treat_vals) > 0:
                t_stat, p_val = stats.ttest_ind(ctrl_vals, treat_vals, equal_var=False)
            else:
                p_val = 1.0
                t_stat = 0.0

            results.append({
                'gene_id': gene_ids[i],
                'baseMean_control': mean_ctrl,
                'baseMean_treatment': mean_treat,
                'log2FoldChange': log2fc,
                'pvalue': p_val,
                'tstat': t_stat
            })

        self.results_df = pd.DataFrame(results)

        # Multiple testing correction (Bonferroni)
        self.results_df['padj'] = self.results_df['pvalue'] * len(self.results_df)
        self.results_df['padj'] = self.results_df['padj'].clip(upper=1.0)

        # Calculate -log10(pvalue) for volcano plot
        self.results_df['neg_log10_pval'] = -np.log10(self.results_df['pvalue'] + 1e-300)

        # Classify genes
        self.results_df['significance'] = 'Not Significant'
        up_mask = (self.results_df['padj'] < self.p_threshold) & (self.results_df['log2FoldChange'] > self.fc_threshold)
        down_mask = (self.results_df['padj'] < self.p_threshold) & (self.results_df['log2FoldChange'] < -self.fc_threshold)

        self.results_df.loc[up_mask, 'significance'] = 'Upregulated'
        self.results_df.loc[down_mask, 'significance'] = 'Downregulated'

        n_up = up_mask.sum()
        n_down = down_mask.sum()
        n_total = n_up + n_down

        print(f"  Total genes analyzed: {len(self.results_df)}")
        print(f"  Significant genes: {n_total}")
        print(f"  Upregulated: {n_up}")
        print(f"  Downregulated: {n_down}")

        # Save results
        output_file = os.path.join(self.output_dir, 'dge_results.csv')
        self.results_df.to_csv(output_file, index=False)
        print(f"  Results saved to: {output_file}")

    def create_volcano_plot(self):
        """Generate volcano plot"""
        print("\nGenerating volcano plot...")

        fig = px.scatter(
            self.results_df,
            x='log2FoldChange',
            y='neg_log10_pval',
            color='significance',
            color_discrete_map={
                'Upregulated': '#00ff88',
                'Downregulated': '#ff6b6b',
                'Not Significant': '#888888'
            },
            hover_data=['gene_id', 'pvalue', 'padj'],
            labels={
                'log2FoldChange': 'Log2 Fold Change',
                'neg_log10_pval': '-Log10(p-value)'
            },
            title=f'Volcano Plot: {self.treatment_group} vs {self.control_group}'
        )

        # Add threshold lines
        fig.add_hline(y=-np.log10(self.p_threshold), line_dash="dash",
                     line_color="gray", opacity=0.5, annotation_text=f"p={self.p_threshold}")
        fig.add_vline(x=self.fc_threshold, line_dash="dash",
                     line_color="gray", opacity=0.5)
        fig.add_vline(x=-self.fc_threshold, line_dash="dash",
                     line_color="gray", opacity=0.5)

        fig.update_layout(
            template="plotly_white",
            width=1000,
            height=700,
            font=dict(size=12)
        )

        output_file = os.path.join(self.output_dir, 'volcano_plot.html')
        fig.write_html(output_file)
        print(f"  Volcano plot saved to: {output_file}")

    def create_ma_plot(self):
        """Generate MA plot"""
        print("\nGenerating MA plot...")

        self.results_df['baseMean_log2'] = np.log2(
            (self.results_df['baseMean_control'] + self.results_df['baseMean_treatment']) / 2 + 1
        )

        fig = px.scatter(
            self.results_df,
            x='baseMean_log2',
            y='log2FoldChange',
            color='significance',
            color_discrete_map={
                'Upregulated': '#00ff88',
                'Downregulated': '#ff6b6b',
                'Not Significant': '#888888'
            },
            hover_data=['gene_id', 'pvalue', 'padj'],
            labels={
                'baseMean_log2': 'Log2 Mean Expression',
                'log2FoldChange': 'Log2 Fold Change'
            },
            title=f'MA Plot: {self.treatment_group} vs {self.control_group}'
        )

        fig.add_hline(y=0, line_dash="dash", line_color="black", opacity=0.5)

        fig.update_layout(
            template="plotly_white",
            width=1000,
            height=700,
            font=dict(size=12)
        )

        output_file = os.path.join(self.output_dir, 'ma_plot.html')
        fig.write_html(output_file)
        print(f"  MA plot saved to: {output_file}")

    def create_pca_plot(self):
        """Generate PCA plot"""
        print("\nGenerating PCA plot...")

        # Get all sample columns
        all_samples = self.control_samples + self.treatment_samples
        counts_matrix = self.counts_df[all_samples].values.T  # Transpose: samples x genes

        # Log transform and standardize
        counts_log = np.log2(counts_matrix + 1)
        scaler = StandardScaler()
        counts_scaled = scaler.fit_transform(counts_log)

        # PCA
        pca = PCA(n_components=min(len(all_samples), 10))
        pca_coords = pca.fit_transform(counts_scaled)

        # Create dataframe
        pca_df = pd.DataFrame({
            'PC1': pca_coords[:, 0],
            'PC2': pca_coords[:, 1],
            'Sample': all_samples,
            'Group': [self.control_group] * len(self.control_samples) +
                    [self.treatment_group] * len(self.treatment_samples)
        })

        fig = px.scatter(
            pca_df,
            x='PC1',
            y='PC2',
            color='Group',
            text='Sample',
            color_discrete_map={
                self.control_group: '#00d9ff',
                self.treatment_group: '#ff6b6b'
            },
            title=f'PCA Plot (PC1: {pca.explained_variance_ratio_[0]:.1%}, PC2: {pca.explained_variance_ratio_[1]:.1%})',
            labels={
                'PC1': f'PC1 ({pca.explained_variance_ratio_[0]:.1%})',
                'PC2': f'PC2 ({pca.explained_variance_ratio_[1]:.1%})'
            }
        )

        fig.update_traces(textposition='top center', marker=dict(size=12))
        fig.update_layout(
            template="plotly_white",
            width=1000,
            height=700,
            font=dict(size=12)
        )

        output_file = os.path.join(self.output_dir, 'pca_plot.html')
        fig.write_html(output_file)
        print(f"  PCA plot saved to: {output_file}")

        # Save scree plot
        fig_scree = go.Figure()
        fig_scree.add_trace(go.Bar(
            x=[f'PC{i+1}' for i in range(len(pca.explained_variance_ratio_))],
            y=pca.explained_variance_ratio_ * 100,
            marker_color='#00d9ff'
        ))
        fig_scree.update_layout(
            title='PCA Scree Plot',
            xaxis_title='Principal Component',
            yaxis_title='Variance Explained (%)',
            template="plotly_white",
            width=800,
            height=500
        )

        output_file = os.path.join(self.output_dir, 'pca_scree_plot.html')
        fig_scree.write_html(output_file)
        print(f"  Scree plot saved to: {output_file}")

    def create_heatmap(self, n_genes=50):
        """Generate heatmap of top DEGs"""
        print(f"\nGenerating heatmap (top {n_genes} genes)...")

        # Get top significant genes
        sig_genes = self.results_df[self.results_df['significance'] != 'Not Significant'].copy()
        sig_genes = sig_genes.sort_values('padj').head(n_genes)

        if len(sig_genes) == 0:
            print("  No significant genes found for heatmap")
            return

        # Get counts for these genes
        all_samples = self.control_samples + self.treatment_samples
        gene_indices = [self.counts_df.iloc[:, 0].tolist().index(g)
                       for g in sig_genes['gene_id']
                       if g in self.counts_df.iloc[:, 0].values]

        heatmap_data = self.counts_df.iloc[gene_indices][all_samples].values

        # Log transform
        heatmap_log = np.log2(heatmap_data + 1)

        # Z-score normalize
        heatmap_zscore = (heatmap_log - heatmap_log.mean(axis=1, keepdims=True)) / \
                        (heatmap_log.std(axis=1, keepdims=True) + 1e-10)

        fig = go.Figure(data=go.Heatmap(
            z=heatmap_zscore,
            x=all_samples,
            y=[sig_genes.iloc[i]['gene_id'] for i in range(len(gene_indices))],
            colorscale='RdBu_r',
            zmid=0,
            colorbar=dict(title="Z-score")
        ))

        fig.update_layout(
            title=f'Heatmap: Top {len(gene_indices)} DEGs',
            xaxis_title='Samples',
            yaxis_title='Genes',
            template="plotly_white",
            width=1000,
            height=max(700, len(gene_indices) * 15),
            font=dict(size=10)
        )

        output_file = os.path.join(self.output_dir, 'heatmap.html')
        fig.write_html(output_file)
        print(f"  Heatmap saved to: {output_file}")

    def create_summary_report(self):
        """Generate summary statistics and report"""
        print("\nGenerating summary report...")

        n_total = len(self.results_df)
        n_sig = len(self.results_df[self.results_df['significance'] != 'Not Significant'])
        n_up = len(self.results_df[self.results_df['significance'] == 'Upregulated'])
        n_down = len(self.results_df[self.results_df['significance'] == 'Downregulated'])

        # Get top genes
        top_up = self.results_df[self.results_df['significance'] == 'Upregulated'].sort_values('padj').head(10)
        top_down = self.results_df[self.results_df['significance'] == 'Downregulated'].sort_values('padj').head(10)

        summary = {
            'analysis_parameters': {
                'control_group': self.control_group,
                'treatment_group': self.treatment_group,
                'control_samples': self.control_samples,
                'treatment_samples': self.treatment_samples,
                'p_threshold': self.p_threshold,
                'fc_threshold': self.fc_threshold
            },
            'summary_statistics': {
                'total_genes': n_total,
                'significant_genes': n_sig,
                'upregulated': n_up,
                'downregulated': n_down,
                'percent_significant': round(100 * n_sig / n_total, 2)
            },
            'top_upregulated': top_up[['gene_id', 'log2FoldChange', 'padj']].to_dict('records'),
            'top_downregulated': top_down[['gene_id', 'log2FoldChange', 'padj']].to_dict('records')
        }

        # Save JSON report
        output_file = os.path.join(self.output_dir, 'dge_summary.json')
        with open(output_file, 'w') as f:
            json.dump(summary, f, indent=2)
        print(f"  Summary saved to: {output_file}")

        # Save top 50 DEGs table
        top50 = self.results_df[self.results_df['significance'] != 'Not Significant'].sort_values('padj').head(50)
        output_file = os.path.join(self.output_dir, 'top50_degs.csv')
        top50.to_csv(output_file, index=False)
        print(f"  Top 50 DEGs saved to: {output_file}")

    def run_complete_analysis(self):
        """Run complete DGE analysis pipeline"""
        print("="*60)
        print("Starting DGE Analysis Pipeline")
        print("="*60)

        self.load_data()
        self.perform_dge()
        self.create_volcano_plot()
        self.create_ma_plot()
        self.create_pca_plot()
        self.create_heatmap(n_genes=50)
        self.create_summary_report()

        print("\n" + "="*60)
        print("DGE Analysis Complete!")
        print(f"All outputs saved to: {self.output_dir}")
        print("="*60)


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description='Differential Gene Expression Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python dge_analysis.py \\
    --counts featureCounts_output.txt \\
    --metadata samples.csv \\
    --control_group Control \\
    --treatment_group Treatment \\
    --output_dir dge_results \\
    --p_threshold 0.05 \\
    --fc_threshold 1.0

Metadata file format (CSV):
  sample,group
  Sample1,Control
  Sample2,Control
  Sample3,Treatment
  Sample4,Treatment
        """
    )

    parser.add_argument('--counts', required=True,
                       help='Path to featureCounts output file')
    parser.add_argument('--metadata', required=True,
                       help='Path to metadata CSV file (sample, group columns)')
    parser.add_argument('--control_group', required=True,
                       help='Name of control group in metadata')
    parser.add_argument('--treatment_group', required=True,
                       help='Name of treatment group in metadata')
    parser.add_argument('--output_dir', required=True,
                       help='Output directory for results')
    parser.add_argument('--p_threshold', type=float, default=0.05,
                       help='P-value threshold (default: 0.05)')
    parser.add_argument('--fc_threshold', type=float, default=1.0,
                       help='Log2 fold change threshold (default: 1.0)')

    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.counts):
        print(f"Error: Counts file not found: {args.counts}")
        sys.exit(1)

    if not os.path.exists(args.metadata):
        print(f"Error: Metadata file not found: {args.metadata}")
        sys.exit(1)

    # Run analysis
    analyzer = DGEAnalyzer(
        counts_file=args.counts,
        metadata_file=args.metadata,
        control_group=args.control_group,
        treatment_group=args.treatment_group,
        output_dir=args.output_dir,
        p_threshold=args.p_threshold,
        fc_threshold=args.fc_threshold
    )

    analyzer.run_complete_analysis()


if __name__ == '__main__':
    main()
