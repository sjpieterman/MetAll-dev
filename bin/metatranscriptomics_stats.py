import os
import pandas as pd
import streamlit as st
import plotly.express as px


def run_metatranscriptomics_analysis(outdir):
    st.subheader("Metatranscriptomics: Taxonomic Composition")

    kraken_dir = os.path.join(outdir, "kraken2")

    if not os.path.exists(kraken_dir):
        st.info(f"Kraken2 results directory not found at: {os.path.abspath(kraken_dir)}")
        return

    # Look for common Kraken2 report extensions
    extensions = [".k2report", ".report.txt", ".kraken2.report.txt", ".report"]
    reports = [f for f in os.listdir(kraken_dir) if any(f.endswith(ext) for ext in extensions)]
    
    if not reports:
        st.info(f"No Kraken2 reports found in {kraken_dir}. Checked for extensions: {', '.join(extensions)}")
        # Debug: show what files ARE there
        all_files = os.listdir(kraken_dir)
        if all_files:
            st.write("Files found in directory:", all_files)
        return

    all_data = []
    for report in reports:
        # Clean up sample name by removing any of the known extensions
        sample_name = report
        for ext in extensions:
            sample_name = sample_name.replace(ext, "")
            
        report_path = os.path.join(kraken_dir, report)

        try:
            # Kraken2 Report Format:
            # 0: percentage, 1: clades_count, 2: taxon_count, 3: rank, 4: taxid, 5: name
            df_report = pd.read_csv(report_path, sep='\t', header=None)

            def get_stats(tax_id):
                row = df_report[df_report[4] == tax_id]
                if not row.empty:
                    return row.iloc[0][1], row.iloc[0][0]  # count, percentage
                return 0, 0.0

            bact_cnt, bact_pct = get_stats(2)  # Bacteria
            vir_cnt, vir_pct = get_stats(10239)  # Viruses
            human_cnt, human_pct = get_stats(9606)  # Homo sapiens

            all_data.append({
                "Sample": sample_name,
                "Bacteria (Abs)": bact_cnt, "Bacteria (%)": bact_pct,
                "Viral (Abs)": vir_cnt, "Viral (%)": vir_pct,
                "Human (Abs)": human_cnt, "Human (%)": human_pct
            })
        except Exception as e:
            st.warning(f"Could not parse {report}: {e}")

    summary_df = pd.DataFrame(all_data)

    if not summary_df.empty:
        col1, col2 = st.columns([1, 1])

        with col1:
            st.write("### Absolute Reads")
            fig_abs = px.bar(summary_df, x="Sample", y=["Bacteria (Abs)", "Viral (Abs)", "Human (Abs)"],
                             title="Absolute Read Distribution", barmode="group",
                             color_discrete_sequence=["#1f77b4", "#ff7f0e", "#2ca02c"])
            st.plotly_chart(fig_abs, use_container_width=True)

        with col2:
            st.write("### Percentages")
            fig_pct = px.bar(summary_df, x="Sample", y=["Bacteria (%)", "Viral (%)", "Human (%)"],
                             title="Percentage Composition", barmode="stack",
                             color_discrete_sequence=["#1f77b4", "#ff7f0e", "#2ca02c"])
            st.plotly_chart(fig_pct, use_container_width=True)

        st.write("### Detailed Statistics")
        st.dataframe(summary_df, use_container_width=True)