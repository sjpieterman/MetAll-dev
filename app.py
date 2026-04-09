import os
import subprocess
import re
import io
import base64
import importlib.util
import sys

#"""python -m pip install "cachetools<6.0,>=2.0.0" "keras<2.15,>=2.14.0" "ml-dtypes==0.2.0" "wrapt<1.15,>=1.11.0""""

def _ensure_packages(required):
    """
    Ensure required packages are installed into the SAME environment as this interpreter.
    Set environment variable METALL_NO_AUTOINSTALL=1 to disable.
    """
    if os.environ.get("METALL_NO_AUTOINSTALL", "").strip() in {"1", "true", "True", "yes", "YES"}:
        return

    missing = []
    for import_name, pip_name in required:
        if importlib.util.find_spec(import_name) is None:
            missing.append(pip_name)

    if not missing:
        return

    # Use the current Python interpreter to install into the current environment.
    cmd = [sys.executable, "-m", "pip", "install", *missing]

    try:
        print(f"[MetaLL] Missing Python packages detected: {', '.join(missing)}")
        print(f"[MetaLL] Installing with: {' '.join(cmd)}")
        subprocess.check_call(cmd)
        return
    except Exception as first_err:
        # Retry with --user (helps when site-packages is not writable)
        cmd_user = [sys.executable, "-m", "pip", "install", "--user", *missing]
        try:
            print(f"[MetaLL] Initial install failed ({first_err}). Retrying with: {' '.join(cmd_user)}")
            subprocess.check_call(cmd_user)
            return
        except Exception as second_err:
            raise SystemExit(
                "\n[MetaLL] Failed to auto-install required Python packages.\n"
                f"Missing: {', '.join(missing)}\n\n"
                "Tried:\n"
                f"  {' '.join(cmd)}\n"
                f"  {' '.join(cmd_user)}\n\n"
                "Common causes:\n"
                "  • No internet access\n"
                "  • Insufficient permissions to install packages\n"
                "  • pip not available for this interpreter\n\n"
                "You can disable auto-install by setting:\n"
                "  METALL_NO_AUTOINSTALL=1\n"
            ) from None


# Packages this app needs (import name -> pip name)
_ensure_packages([
    ("dash", "dash"),
    ("dash_bootstrap_components", "dash-bootstrap-components"),
    ("plotly", "plotly"),
    ("flask", "flask"),
    ("pandas", "pandas"),
    ("numpy", "numpy"),
    ("scipy", "scipy"),
    ("sklearn", "scikit-learn"),
    ("dash_ace", "dash-ace"),
    ("streamlit", "streamlit"),
])

import os
import subprocess
import pandas as pd
import numpy as np
import plotly.express as px
import re
import io
import base64
from dash import Dash, html, dcc, Input, Output, State, dash_table, callback_context, no_update
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import importlib.util
import sys
import flask
import urllib.parse
import glob
from dash_ace import DashAceEditor
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import time

# --- 1. Load Analysis Module ---
mt_stats = None
try:
    stats_path = os.path.join(os.path.dirname(__file__), "bin", "metatranscriptomics_stats.py")
    if os.path.exists(stats_path):
        spec = importlib.util.spec_from_file_location("metatranscriptomics_stats", stats_path)
        mt_stats = importlib.util.module_from_spec(spec)
        sys.modules["metatranscriptomics_stats"] = mt_stats
        spec.loader.exec_module(mt_stats)
except Exception as e:
    print(f"Warning: Failed to load analysis module: {e}")

# --- 2. Constants & Initial State ---
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def get_run_outdir(run_name, base_outdir=None):
    """
    Resolves the output directory for a given run name.
    If run_name is an absolute path, it returns it directly.
    Otherwise, it resolves it relative to the 'results' directory.
    """
    if not run_name:
        run_name = "run_001"
    run_name = run_name.strip()
    
    # If run_name is absolute, return it directly
    if os.path.isabs(run_name):
        return os.path.normpath(run_name)
        
    # Otherwise, resolve relative to base_outdir or default results dir
    if not base_outdir:
        base_outdir = os.path.join(BASE_DIR, "results")
    else:
        base_outdir = os.path.expanduser(base_outdir)
        if not os.path.isabs(base_outdir):
            base_outdir = os.path.join(BASE_DIR, base_outdir)
            
    return os.path.normpath(os.path.join(base_outdir, run_name))

def find_counts_file(run_outdir):
    """
    Tries to find the gene_counts_matrix.tsv file in expected locations and variations.
    """
    if not run_outdir or not os.path.exists(run_outdir):
        return None

    # Priority search subdirectories
    search_subdirs = ["counts", "featurecounts", "salmon", "star", ""]
    
    # Priority filenames (exact matches)
    target_names = [
        "gene_counts_matrix.tsv",
        "combined_gene_counts_matrix.tsv",
        "featurecounts_gene_counts_matrix.tsv",
        "salmon_gene_counts_matrix.tsv",
        "star_gene_counts_matrix.tsv",
        "gene_counts.tsv"
    ]
    
    for subdir in search_subdirs:
        d = os.path.join(run_outdir, subdir)
        if not os.path.exists(d):
            continue
            
        # Try exact matches first
        for name in target_names:
            p = os.path.join(d, name)
            if os.path.exists(p) and os.path.isfile(p):
                return p
        
        # Try any file ending with gene_counts_matrix.tsv
        try:
            for f in os.listdir(d):
                if f.endswith("gene_counts_matrix.tsv") and os.path.isfile(os.path.join(d, f)):
                    return os.path.join(d, f)
        except Exception:
            pass
            
    return None

INITIAL_PARAMS = {
    "reads": "/media/baseripper/volume_28TB_1/P11065_test/",
    "outdir": "/media/baseripper/volume_28TB_1/P11065_test/results",
    "profile": "singularity",
    "cpus": 40,
    "memory": "200.GB",
    "resume_pipeline": True,
    "run_trimgalore": True,
    "run_ribodetector": True,
    "run_sortmerna": False,
    "run_kraken2": True,
    "tg_quality": 20,
    "tg_min_length": 20,
    "tg_fastqc": True,
    "tg_clip_r1": 0,
    "tg_clip_r2": 0,
    "use_fastp": False,
    "use_falco": True,
    "rd_model": "norrna",
    "rd_len": 150,
    "rd_gpu_mem": 24,
    "rd_threads": 16,
    "rd_device": "gpu",
    "smr_coverage": 0.97,
    "smr_mismatch": 2,
    "smr_num_alignments": 1,
    "smr_db_dir": "/media/baseripper/volume_28TB_2/rRNA_databases",
    "kraken2_db": "/media/baseripper/2TB_Baseripper/Stef/k2_standard_20231009",
    "k2_confidence": 0.0,
    "bracken_level": "S",
    "bracken_read_len": 150,
    "run_parabricks": True,
    "parabricks_index": os.path.join(BASE_DIR, "assets/GRCh38_indexed/"),
    "genome_fasta": os.path.join(BASE_DIR, "assets/GRCh38_indexed/GRCh38.primary_assembly.genome.fa"),
    "genome_gtf": os.path.join(BASE_DIR, "assets/GRCh38_indexed/gencode.v44.annotation.gtf"),
    "run_star": False,
    "star_index": os.path.join(BASE_DIR, "assets/GRCh38_indexed/"),
    "run_featurecounts": False,
    "run_salmon": False,
    "salmon_index": "",
    "run_virulence": False,
    "virulence_db": "/home/baseripper/programs/vfdb/VFDB_setB_pro.fas",
    "virulence_index": "",
    "run_dge": False,
    "dge_tool": "deseq2",
    "dge_control": None,
    "dge_treatment": None,
    "dge_comparison_name": "treatment_vs_control",
    "batch_method": "none",
    "batch_col": "",
    "dge_p_threshold": 0.05,
    "dge_fc_threshold": 1.0,
    "norm_method": "voom"  # Default to voom; can be overridden
}

# --- 3. App Initialization ---
app = Dash(__name__, external_stylesheets=[dbc.themes.SLATE], suppress_callback_exceptions=True)
app.title = "MetaLL Pipeline Controller"

# --- 4. UI Helper Components ---

def build_sidebar():
    return dbc.Card([
        dbc.CardBody([
            html.H4("Settings", className="card-title mb-4", style={"fontWeight": "600", "color": "#00d9ff"}),
            html.Label("Input Reads Path", style={"fontWeight": "500", "fontSize": "0.9rem"}),
            dbc.Input(id="input-reads", value=INITIAL_PARAMS["reads"], type="text", className="mb-3"),

            html.Label("Max CPUs", className="mt-2", style={"fontWeight": "500", "fontSize": "0.9rem"}),
            dcc.Dropdown(
                options=[1, 4, 8, 16, 24, 32, 40, 48, 64],
                value=INITIAL_PARAMS["cpus"], id="slider-cpus",
                className="mb-3"
            ),

            html.Label("Max Memory", className="mt-2", style={"fontWeight": "500", "fontSize": "0.9rem"}),
            dcc.Dropdown(
                options=["16.GB", "64.GB", "128.GB", "200.GB", "248.GB"],
                value=INITIAL_PARAMS["memory"], id="drop-mem",
                className="mb-3"
            ),

            html.Label("Base Results Directory", className="mt-2", style={"fontWeight": "500", "fontSize": "0.9rem"}),
            dbc.Input(id="base-outdir", value=INITIAL_PARAMS["outdir"], type="text", debounce=True, className="mb-3"),

            html.Label("Output Run Name", className="mt-2", style={"fontWeight": "500", "fontSize": "0.9rem"}),
            dbc.Input(id="output-run-name", value="run_001", type="text", debounce=True, className="mb-3"),

            dbc.Checklist(
                options=[{"label": "Resume Pipeline", "value": 1}],
                value=[1] if INITIAL_PARAMS["resume_pipeline"] else [],
                id="check-resume", switch=True, className="mt-3"
            ),

            html.Hr(style={"borderColor": "#444"}),
            html.H4("Workflow", style={"fontWeight": "600", "color": "#00d9ff"}),
            dbc.Checklist(
                options=[
                    {"label": "TrimGalore", "value": "run_trimgalore"},
                    {"label": "RiboDetector", "value": "run_ribodetector"},
                    {"label": "SortMeRNA", "value": "run_sortmerna"},
                    {"label": "Kraken2 & Bracken", "value": "run_kraken2"},
                    {"label": "Parabricks (GPU)", "value": "run_parabricks"},
                    {"label": "STAR (CPU)", "value": "run_star"},
                    {"label": "featureCounts", "value": "run_featurecounts"},
                    {"label": "Salmon", "value": "run_salmon"},
                    {"label": "Virulence Factors", "value": "run_virulence"},
                    {"label": "DGE Analysis", "value": "run_dge"},
                ],
                value=["run_trimgalore", "run_ribodetector", "run_kraken2", "run_star", "run_featurecounts", "run_salmon", "run_virulence", "run_dge"],
                id="workflow-checks", switch=True
            ),

            dbc.Button("Run Pipeline", id="run-btn", color="info", className="mt-4 w-100",
                      style={"fontWeight": "600", "fontSize": "1.1rem", "padding": "12px"}),
        ])
    ], className="h-100", style={"boxShadow": "0 4px 6px rgba(0, 0, 0, 0.3)", "borderRadius": "8px"})


def tool_layout(title, settings_content, display_content):
    return dbc.Row([
        dbc.Col([html.H4(title), html.Hr(), settings_content], width=3, className="border-end"),
        dbc.Col(display_content, width=9)
    ], className="mt-3")


# --- 5. Tab Content Definitions ---

def build_main_controller_tab():
    return dbc.Row([
        dbc.Col(build_sidebar(), width=3),
        dbc.Col([
            dbc.Card([
                dbc.CardBody([
                    html.H4("Execution Log", className="mb-3", style={"fontWeight": "600", "color": "#00d9ff"}),
                    html.Div([
                        html.Pre(id="log-display", style={
                            "backgroundColor": "#1a1a1a", "color": "#00ff88",
                            "padding": "20px", "height": "500px", "overflowY": "scroll",
                            "borderRadius": "8px", "fontFamily": "monospace", "fontSize": "0.9rem",
                            "border": "1px solid #333", "boxShadow": "inset 0 2px 4px rgba(0, 0, 0, 0.4)"
                        }),
                    ], id="log-container"),
                ])
            ], style={"boxShadow": "0 4px 6px rgba(0, 0, 0, 0.3)", "borderRadius": "8px"}),
            dcc.Interval(id='log-update', interval=2000, n_intervals=0, disabled=True),
            dcc.Store(id='pipeline-running', data=False, storage_type='session')
        ], width=9)
    ], className="mt-3")


def build_preprocessing_tab():
    return html.Div([
        html.H4("Data Quality Control", className="mt-3 mb-4"),
        dbc.Tabs([
            dbc.Tab(label="Samplesheet Discovery", children=html.Div([
                html.H5("Discovered Samples", className="mt-3"),
                html.P("Listing all paired-end FASTQ files identified for processing.", className="text-muted small"),
                html.Div(id="samplesheet-table-container")
            ], className="p-3")),
            dbc.Tab(label="Raw Data QC", children=html.Div(id="raw-qc-container", className="p-3")),
            dbc.Tab(label="Post-Trimming QC", children=html.Div(id="trimmed-qc-container", className="p-3"))
        ])
    ])


def build_trimgalore_tab():
    return tool_layout(
        "TrimGalore Settings",
        html.Div([
            html.Label("Quality Threshold", style={"fontWeight": "500", "fontSize": "0.9rem"}),
            dbc.Input(id="tg-quality", type="number", value=INITIAL_PARAMS["tg_quality"], className="mb-2"),
            html.Label("Min Read Length", className="mt-2", style={"fontWeight": "500", "fontSize": "0.9rem"}),
            dbc.Input(id="tg-min-len", type="number", value=INITIAL_PARAMS["tg_min_length"], className="mb-2"),
            html.Label("5' Clip R1", className="mt-2", style={"fontWeight": "500", "fontSize": "0.9rem"}),
            dbc.Input(id="tg-clip-r1", type="number", value=INITIAL_PARAMS["tg_clip_r1"], className="mb-2"),
            html.Label("5' Clip R2", className="mt-2", style={"fontWeight": "500", "fontSize": "0.9rem"}),
            dbc.Input(id="tg-clip-r2", type="number", value=INITIAL_PARAMS["tg_clip_r2"], className="mb-2"),
            dbc.Checklist(
                options=[{"label": "Run Internal FastQC", "value": 1}],
                value=[1] if INITIAL_PARAMS["tg_fastqc"] else [],
                id="tg-fastqc", switch=True, className="mt-3"
            ),
            html.Hr(),
            html.H5("Performance Options", className="mt-3", style={"fontSize": "1rem"}),
            dbc.Checklist(
                options=[
                    {"label": "Use Fastp (Faster QC & Trimming)", "value": "use_fastp"},
                    {"label": "Use Falco (Faster FastQC)", "value": "use_falco"},
                ],
                value=(["use_fastp"] if INITIAL_PARAMS["use_fastp"] else []) + (["use_falco"] if INITIAL_PARAMS["use_falco"] else []),
                id="perf-checks", switch=True, className="mt-2"
            ),
        ]),
        html.Div([
            html.H5("Trimming Statistics", className="mt-3 mb-3", style={"fontWeight": "600", "color": "#00d9ff"}),
            html.Div(id="tg-stats-plots-container"),
            dcc.Interval(id='tg-stats-update', interval=5000, n_intervals=0, disabled=False)
        ])
    )


def build_ribodetector_tab():
    return tool_layout(
        "RiboDetector Settings",
        html.Div([
            html.Label("Model Read Length"),
            dcc.Dropdown(options=[50, 70, 100, 150], value=INITIAL_PARAMS["rd_len"], id="rd-len"),
            html.Label("Detection Mode", className="mt-2"),
            dbc.RadioItems(
                options=[
                    {"label": "No rRNA", "value": "norrna"},
                    {"label": "rRNA Only", "value": "rrna"},
                    {"label": "Both", "value": "both"}
                ],
                value=INITIAL_PARAMS["rd_model"],
                id="rd-mode"
            ),
            html.Label("Device", className="mt-2"),
            dbc.RadioItems(
                options=[{"label": "GPU", "value": "gpu"}, {"label": "CPU", "value": "cpu"}],
                value=INITIAL_PARAMS["rd_device"],
                id="rd-device"
            ),
            html.Label("GPU Memory (GB)", className="mt-2"),
            dbc.Input(id="rd-gpu-mem", type="number", value=INITIAL_PARAMS["rd_gpu_mem"]),
            html.Label("I/O Threads", className="mt-2"),
            dbc.Input(id="rd-threads", type="number", value=INITIAL_PARAMS["rd_threads"]),
        ]),
        html.Div([
            html.H5("Live Processing Status"),
            dbc.Tabs([
                dbc.Tab(label="Plots", children=html.Div(id="rd-plots-container", className="mt-3")),
                dbc.Tab(label="Statistics", children=html.Div(id="rd-stats-container", className="mt-3")),
            ]),
            dcc.Interval(id='rd-stats-update', interval=5000, n_intervals=0)
        ])
    )


def build_sortmerna_tab():
    return tool_layout(
        "SortMeRNA Settings",
        html.Div([
            html.Label("Min Coverage (ID)"),
            dcc.Slider(0.0, 1.0, 0.01, value=INITIAL_PARAMS["smr_coverage"],
                      id="smr-coverage", marks={0: '0', 0.5: '0.5', 1: '1'}),
            html.Label("Max Mismatches", className="mt-2"),
            dbc.Input(id="smr-mismatch", type="number", value=INITIAL_PARAMS["smr_mismatch"]),
            html.Label("Num Alignments", className="mt-2"),
            dcc.Dropdown(options=[1, 5, 10], value=INITIAL_PARAMS["smr_num_alignments"], id="smr-num-align"),
            html.Label("Database Path", className="mt-2"),
            dbc.Input(id="smr-db-dir", type="text", value=INITIAL_PARAMS["smr_db_dir"]),
        ]),
        html.Div([dbc.Alert("SortMeRNA configuration updated. Parameters will be used on next pipeline run.", color="info")])
    )


def build_virulence_tab():
    return tool_layout(
        "Virulence Factor Settings",
        html.Div([
            html.P("Align reads to a curated virulence gene database to confirm expression.", className="small text-muted"),
            html.Label("Virulence DB (FASTA)", style={"fontWeight": "500"}),
            dbc.Input(id="virulence-db", type="text", value=INITIAL_PARAMS["virulence_db"], placeholder="/path/to/vfdb.fasta", className="mb-2"),
            html.Label("Pre-computed Index", className="mt-2", style={"fontWeight": "500"}),
            dbc.Input(id="virulence-index", type="text", value=INITIAL_PARAMS["virulence_index"], placeholder="/path/to/salmon_index", className="mb-2"),
            html.Hr(),
            html.P("Common Databases:", className="small mb-1"),
            html.Ul([
                html.Li("VFDB (Virulence Factor Database)", className="small"),
                html.Li("PATRIC / Virulence", className="small"),
                html.Li("CARD (Resistance Factors)", className="small"),
                html.Li("BacMet (Metal Resistance)", className="small"),
            ], className="small text-muted")
        ]),
        html.Div([
            html.H5("Virulence Results", className="mt-3"),
            html.P("Abundance estimates (TPM) will be available in the results viewer after the run.", className="text-muted"),
            dbc.Alert("Ensure you have enabled 'Virulence Factors' in the sidebar workflow settings.", color="warning", className="mt-3")
        ])
    )


def build_results_viewer_tab():
    return tool_layout(
        "Pipeline Outputs",
        html.Div([
            html.H6("Select Tool Output"),
            dcc.Dropdown(
                id="results-folder-dropdown",
                options=[
                    {"label": "FastQC Raw", "value": "fastqc_raw"},
                    {"label": "FastQC Trimmed", "value": "fastqc_trimmed"},
                    {"label": "Trim Galore", "value": "trimgalore"},
                    {"label": "RiboDetector", "value": "ribodetector"},
                    {"label": "SortMeRNA", "value": "sortmerna"},
                    {"label": "Kraken2", "value": "kraken2"},
                    {"label": "Bracken", "value": "bracken"},
                    {"label": "Virulence", "value": "virulence"},
                    {"label": "STAR", "value": "star"},
                    {"label": "featureCounts", "value": "featurecounts"},
                    {"label": "Salmon", "value": "salmon"},
                    {"label": "Gene Counts", "value": "counts"},
                    {"label": "DGE Analysis", "value": "dge_analysis"},
                ],
                placeholder="Choose directory..."
            ),
            html.Hr(),
            html.H6("MultiQC Hub"),
            dbc.ButtonGroup([
                dbc.Button("Raw Data QC", id="btn-mqc-raw", color="secondary", outline=True, size="sm"),
                dbc.Button("Post-Trim QC", id="btn-mqc-trim", color="secondary", outline=True, size="sm"),
                dbc.Button("RiboDetector QC", id="btn-mqc-ribo", color="secondary", outline=True, size="sm"),
                dbc.Button("SortMeRNA QC", id="btn-mqc-smr", color="secondary", outline=True, size="sm"),
            ], vertical=True, className="w-100")
        ]),
        html.Div([
            html.Div(id="results-files-container", children=[
                dash_table.DataTable(
                    id='results-files-table',
                    columns=[{"name": "Filename", "id": "name"}, {"name": "Size (MB)", "id": "size"}],
                    data=[],
                    style_cell={'textAlign': 'left'},
                    page_size=20
                )
            ]),
            html.Iframe(
                id="mqc-frame-v2",
                src="",
                style={"width": "100%", "height": "800px", "border": "none", "display": "none"}
            )
        ])
    )


def build_transcriptomics_settings():
    return html.Div([
        html.H4("Quantification & Alignment Settings", className="mb-4"),
        dbc.Row([
            dbc.Col([
                html.H5("STAR Alignment (CPU)", style={"color": "#00d9ff"}),
                html.Label("STAR Index Path"),
                dbc.Input(id="star-index", type="text", value=INITIAL_PARAMS["star_index"], className="mb-2"),
                html.Label("Genome Fasta"),
                dbc.Input(id="genome-fasta", type="text", value=INITIAL_PARAMS["genome_fasta"], className="mb-2"),
                html.Label("Genome GTF"),
                dbc.Input(id="genome-gtf", type="text", value=INITIAL_PARAMS["genome_gtf"], className="mb-2"),
                html.Hr(),
                html.H5("Parabricks (GPU)", style={"color": "#00d9ff"}),
                html.Label("Parabricks Index Path"),
                dbc.Input(id="parabricks-index", type="text", value=INITIAL_PARAMS["parabricks_index"], className="mb-2"),
            ], width=6, className="border-end"),
            dbc.Col([
                html.H5("Salmon & featureCounts", style={"color": "#00d9ff"}),
                html.H6("Salmon"),
                html.Label("Salmon Index Path (Optional)"),
                dbc.Input(id="salmon-index", type="text", value=INITIAL_PARAMS["salmon_index"], className="mb-3", placeholder="Leave empty to build from fasta"),

                html.H6("featureCounts"),
                html.Label("Genome GTF"),
                dbc.Input(id="fc-gtf", type="text", value=INITIAL_PARAMS["genome_gtf"], className="mb-3"),

                dbc.Alert("These parameters configure the quantification tools. Enable them in the Workflow settings on the Main Controller tab.", color="info", className="mt-4")
            ], width=6),
        ])
    ], className="p-3")


def build_metatranscriptomics_tab():
    return html.Div([
        html.H4("Metatranscriptomics Analysis", className="mt-3 mb-4"),
        dbc.Row([
            dbc.Col([
                html.H6("Taxonomy Settings", style={"fontWeight": "600", "color": "#00d9ff"}),
                html.Label("Kraken2 DB Path", className="mt-2"),
                dbc.Input(id="k2-db-path", type="text", value=INITIAL_PARAMS["kraken2_db"], className="mb-2"),
                html.Label("Confidence Threshold"),
                dcc.Slider(0.0, 1.0, 0.05, value=INITIAL_PARAMS["k2_confidence"], id="k2-confidence", marks={0: '0', 0.5: '0.5', 1: '1'}),
                html.Label("Bracken Taxonomic Level", className="mt-2"),
                dcc.Dropdown(options=["D", "P", "C", "O", "F", "G", "S"], value=INITIAL_PARAMS["bracken_level"], id="bracken-level"),
                html.Hr(),
                html.Label("Select Sample for Detail"),
                dcc.Dropdown(id="taxonomy-sample-dropdown", placeholder="Select a sample...", persistence=True, persistence_type='session'),
                html.Div(id="mt-update-indicator", style={"fontSize": "0.7rem", "color": "#888", "marginTop": "10px"})
            ], width=3, className="border-end"),
            dbc.Col([
                dbc.Tabs([
                    dbc.Tab(label="Classification Overview", children=[
                        html.Div([
                            html.H5("Host vs Microbial Classification", className="mt-3 mb-4", style={"fontWeight": "600", "color": "#00d9ff"}),
                            html.Div(id="mt-pie-charts-container"),
                        ], className="p-3")
                    ]),
                    dbc.Tab(label="Species Abundance", children=[
                        html.Div([
                            dcc.Graph(id="bracken-bar-plot")
                        ], className="p-3")
                    ]),
                    dbc.Tab(label="Overall Abundance", children=[
                        html.Div([
                            dcc.Graph(id="overall-abundance-plot")
                        ], className="p-3")
                    ]),
                    dbc.Tab(label="Krona Plot", children=[
                        html.Div(id="krona-container", className="p-3")
                    ])
                ])
            ], width=9)
        ])
    ])


# --- 5. UI Helper Components ---

def build_dge_analysis_tab():
    return html.Div([
        html.H4("Differential Gene Expression Analysis", className="mt-3 mb-4"),
        dbc.Row([
            dbc.Col([
                html.H6("Run Configuration", style={"fontWeight": "600", "color": "#00d9ff"}),
                html.Label("DGE Tool", className="mt-2"),
                dcc.Dropdown(
                    id="dge-tool",
                    options=[
                        {"label": "DESeq2", "value": "deseq2"},
                        {"label": "edgeR", "value": "edger"},
                        {"label": "limma-voom", "value": "limma"}
                    ],
                    value=INITIAL_PARAMS["dge_tool"]
                ),
                html.Label("Group Column", className="mt-2"),
                dcc.Dropdown(id="dge-group-col", placeholder="Select metadata column..."),
                html.Label("Control Group Value", className="mt-2"),
                dcc.Dropdown(id="dge-control", placeholder="Select control value..."),
                html.Label("Treatment Group Value", className="mt-2"),
                dcc.Dropdown(id="dge-treatment", placeholder="Select treatment value..."),
                html.Label("Covariates", className="mt-2"),
                dcc.Dropdown(id="dge-covariates", multi=True, placeholder="Select covariates..."),
                html.Label("Comparison Name", className="mt-2"),
                dbc.Input(id="dge-comparison", type="text", value=INITIAL_PARAMS["dge_comparison_name"]),
                dbc.Button("Run DGE Analysis", id="btn-run-dge", color="primary", className="mt-3 w-100"),

                dbc.Accordion([
                    dbc.AccordionItem([
                        html.Label("Batch Correction Method"),
                        dcc.Dropdown(
                            id="batch-method",
                            options=[
                                {"label": "None", "value": "none"},
                                {"label": "ComBat", "value": "combat"},
                                {"label": "Limma-style", "value": "limma"}
                            ],
                            value=INITIAL_PARAMS["batch_method"]
                        ),
                        html.Label("Batch Column", className="mt-2"),
                        dcc.Dropdown(id="batch-col", placeholder="Select batch column..."),
                        html.Label("P-value Threshold", className="mt-2"),
                        dbc.Input(id="dge-p-thresh", type="number", value=INITIAL_PARAMS["dge_p_threshold"], step=0.01),
                        html.Label("Log2 FC Threshold", className="mt-2"),
                        dbc.Input(id="dge-fc-thresh", type="number", value=INITIAL_PARAMS["dge_fc_threshold"], step=0.1),
                        html.Hr(),
                        html.Label("Functional Enrichment"),
                        dbc.Checklist(
                            options=[
                                {"label": "Run GO (clusterProfiler)", "value": "go"},
                                {"label": "Run GSEA (clusterProfiler)", "value": "gsea"}
                            ],
                            value=[],
                            id="dge-enrichment-checks",
                            switch=True,
                            className="mt-2"
                        ),
                        html.Label("Ontology", className="mt-2"),
                        dcc.Dropdown(
                            id="dge-ontology",
                            options=[
                                {"label": "Biological Process (BP)", "value": "BP"},
                                {"label": "Cellular Component (CC)", "value": "CC"},
                                {"label": "Molecular Function (MF)", "value": "MF"},
                                {"label": "All (BP, CC, MF)", "value": "ALL"}
                            ],
                            value="BP"
                        ),
                        html.Label("Organism DB", className="mt-2"),
                        dcc.Dropdown(
                            id="dge-organism-db",
                            options=[{"label": "Human (org.Hs.eg.db)", "value": "org.Hs.eg.db"}],
                            value="org.Hs.eg.db"
                        ),
                        html.Label("Gene ID Type (for enrichment)", className="mt-2"),
                        dcc.Dropdown(
                            id="dge-keytype",
                            options=[
                                {"label": "Auto-detect", "value": "SYMBOL"},
                                {"label": "SYMBOL", "value": "SYMBOL"},
                                {"label": "ENSEMBL", "value": "ENSEMBL"},
                                {"label": "ENTREZID", "value": "ENTREZID"}
                            ],
                            value="SYMBOL"
                        ),
                        dbc.Checklist(
                            options=[{"label": "Use biomaRt to label plots with gene symbols (if counts use Ensembl)", "value": 1}],
                            value=[],
                            id="dge-use-biomart",
                            switch=True,
                            className="mt-2"
                        ),
                        html.Hr(),
                        html.H6("R Environment", style={"fontWeight": "600", "fontSize": "0.85rem"}),
                        dbc.Checklist(
                            options=[{"label": "Auto-install missing R packages", "value": 1}],
                            value=[1],
                            id="dge-auto-install",
                            switch=True,
                            className="mt-2"
                        ),
                        html.Label("HTTP Proxy (if needed)", className="mt-2"),
                        dbc.Input(id="dge-proxy", type="text", placeholder="e.g. http://proxy:8080", size="sm"),
                        html.P("Set this if your server is behind a proxy and lacks direct internet access.", 
                               className="text-muted", style={"fontSize": "0.75rem"}),
                    ], title="Advanced DGE Parameters")
                ], start_collapsed=True, className="mt-3", style={"fontSize": "0.9rem"}),

                html.Hr(),
                html.H6("Results Selection", style={"fontWeight": "600", "color": "#00d9ff"}),
                html.Label("Select Comparison", className="mt-2"),
                dcc.Dropdown(id="dge-comparison-dropdown", placeholder="Choose analysis..."),
                dbc.Button("Refresh Results", id="btn-refresh-dge", color="secondary", size="sm", className="mt-2 w-100"),
                html.Div(id="dge-analysis-info", className="mt-2"),

                html.Hr(),
                dbc.Alert(
                    [
                        html.Div("Markdown reports: place a rendered HTML report at:"),
                        html.Code("dge_analysis/<comparison>/dge_report.html"),
                        html.Div("Then open the “Markdown Analysis” tab."),
                    ],
                    color="info",
                    className="mt-2",
                    style={"fontSize": "0.85rem"}
                ),
            ], width=3, className="border-end"),
            dbc.Col([
                dbc.Tabs([
                    dbc.Tab(label="Summary", children=[
                        html.Div(id="dge-summary-stats", className="p-3")
                    ]),
                    dbc.Tab(label="Plots", children=[
                        html.Div([
                            dbc.Row([
                                dbc.Col([
                                    html.H5("Volcano Plot", className="text-center"),
                                    html.Iframe(id="dge-volcano-iframe", style={"width": "100%", "height": "700px", "border": "none", "overflow": "hidden"})
                                ], width=6),
                                dbc.Col([
                                    html.H5("MA Plot", className="text-center"),
                                    html.Iframe(id="dge-ma-iframe", style={"width": "100%", "height": "700px", "border": "none", "overflow": "hidden"})
                                ], width=6),
                            ], className="mb-4"),
                            dbc.Row([
                                dbc.Col([
                                    html.H5("PCA Plot", className="text-center"),
                                    html.Iframe(id="dge-pca-iframe", style={"width": "100%", "height": "700px", "border": "none", "overflow": "hidden"})
                                ], width=6),
                                dbc.Col([
                                    html.H5("Scree Plot", className="text-center"),
                                    html.Iframe(id="dge-scree-iframe", style={"width": "100%", "height": "500px", "border": "none", "overflow": "hidden"})
                                ], width=6),
                            ], className="mb-4"),
                            dbc.Row([
                                dbc.Col([
                                    html.H5("Heatmap", className="text-center"),
                                    html.Iframe(id="dge-heatmap-iframe", style={"width": "100%", "height": "800px", "border": "none"})
                                ], width=12),
                            ], className="mb-4"),
                        ], className="p-3")
                    ]),
                    dbc.Tab(label="Enrichment", children=[
                        html.Div([
                            dbc.Row([
                                dbc.Col([
                                    html.H5("GO Enrichment (UP/DOWN)", className="text-center"),
                                    html.Div(id="dge-go-figs", className="text-center")
                                ], width=6),
                                dbc.Col([
                                    html.H5("GSEA GO", className="text-center"),
                                    html.Div(id="dge-gsea-figs", className="text-center")
                                ], width=6),
                            ])
                        ], className="p-3")
                    ]),
                    dbc.Tab(label="Top 50 DEGs", children=[
                        html.Div([
                            html.H5("Top 50 Differentially Expressed Genes"),
                            dash_table.DataTable(
                                id='dge-top50-table',
                                columns=[],
                                data=[],
                                page_size=50,
                                style_cell={'textAlign': 'left', 'padding': '8px'},
                                style_header={'backgroundColor': '#2b3e50', 'color': 'white', 'fontWeight': 'bold'},
                                style_data_conditional=[
                                    {'if': {'column_id': 'log2FoldChange', 'filter_query': '{log2FoldChange} > 0'},
                                     'backgroundColor': '#00ff8820'},
                                    {'if': {'column_id': 'log2FoldChange', 'filter_query': '{log2FoldChange} < 0'},
                                     'backgroundColor': '#ff6b6b20'}
                                ],
                                sort_action='native',
                                filter_action='native'
                            )
                        ], className="p-3")
                    ]),
                    dbc.Tab(label="Markdown Analysis", children=[
                        html.Div([
                            html.P("You can edit the R Markdown code below and click 'Run Analysis' to re-render the report.", className="mb-2"),
                            dbc.Row([
                                dbc.Col([
                                    DashAceEditor(
                                        id='dge-markdown-editor',
                                        mode='r',
                                        theme='monokai',
                                        width='100%',
                                        height='600px',
                                        value='',
                                        tabSize=2,
                                        enableBasicAutocompletion=True,
                                        enableLiveAutocompletion=True,
                                    ),
                                    dbc.Button("Run Analysis", id="btn-run-markdown", color="primary", className="mt-2 w-100"),
                                    html.Div(id="markdown-run-status", className="mt-2")
                                ], width=6),
                                dbc.Col([
                                    html.Iframe(
                                        id="dge-markdown-iframe",
                                        style={"width": "100%", "height": "600px", "border": "none"},
                                    ),
                                    dbc.Alert(
                                        id="dge-markdown-missing",
                                        children="No markdown report found for this comparison yet.",
                                        color="warning",
                                        className="mt-2",
                                    ),
                                ], width=6)
                            ])
                        ], className="p-3")
                    ])
                ])
            ], width=9, style={"overflowY": "auto", "maxHeight": "90vh"})
        ]),
        dcc.Interval(id='dge-refresh-interval', interval=10000, n_intervals=0),
        dcc.Interval(id='dge-log-interval', interval=2000, n_intervals=0, disabled=True),
        html.Hr(),
        html.H6("Live DGE Log", style={"fontWeight": "600", "color": "#00d9ff"}),
        html.Pre(id="dge-live-log", style={"height": "240px", "overflowY": "auto", "backgroundColor": "#111", "color": "#ddd", "padding": "10px", "borderRadius": "4px"}),
        dcc.Store(id="dge-run-state", data={"running": False})
    ])


def build_host_transcriptomics_tab():
    return html.Div([
        html.H4("Host Transcriptomics: Transcriptome Coverage", className="mt-3 mb-4", style={"fontWeight": "600", "color": "#00d9ff"}),
        dbc.Row([
            dbc.Col([
                html.H6("Analysis Controls", style={"fontWeight": "500", "color": "#00d9ff"}),
                html.P("Calculates average and median transcriptome coverage using gene counts and lengths from featureCounts or Salmon.", style={"fontSize": "0.85rem"}),
                html.Hr(),
                dbc.Button("Refresh Coverage Data", id="btn-refresh-host", color="primary", className="mb-3 w-100"),
                html.Label("Min Reads Threshold", className="mt-2", style={"fontSize": "0.9rem"}),
                dcc.Slider(0, 10, 1, value=1, id="host-min-reads", marks={0: '0', 5: '5', 10: '10'}),
                html.Hr(),
                html.Div(id="host-analysis-info", style={"fontSize": "0.75rem", "color": "#888"})
            ], width=3, className="border-end"),
            dbc.Col([
                dbc.Tabs([
                    dbc.Tab(label="Coverage Overview", children=[
                        html.Div([
                            dcc.Graph(id="host-coverage-bar-plot")
                        ], className="p-3")
                    ]),
                    dbc.Tab(label="Per-Sample Distribution", children=[
                        html.Div([
                            dcc.Graph(id="host-coverage-dist-plot")
                        ], className="p-3")
                    ])
                ])
            ], width=9)
        ])
    ])


def build_analysis_tab():
    return html.Div([
        html.H3("Experimental Data Analysis", className="mt-3 mb-4"),
        dbc.Tabs([
            dbc.Tab(label="Metadata Manager", children=[
                tool_layout(
                    "Table Controls",
                    html.Div([
                        dcc.Upload(
                            id="upload-metadata",
                            children=html.Div(["Drag and Drop or ", html.A("Select Metadata")]),
                            style={
                                "borderWidth": "1px",
                                "borderStyle": "dashed",
                                "borderRadius": "5px",
                                "textAlign": "center",
                                "padding": "20px",
                                "cursor": "pointer",
                            },
                        ),
                        dbc.Button("Download CSV", id="btn-csv-download", color="success", className="mt-3 w-100"),
                        dcc.Download(id="download-metadata-csv"),
                        html.Hr(),
                        dbc.Checklist(
                            options=[{"label": "Show Column Statistics", "value": 1}],
                            value=[],
                            id="show-stats-check",
                            switch=True,
                            className="mt-2",
                        ),
                    ]),
                    html.Div([
                        dash_table.DataTable(
                            id="metadata-datatable",
                            columns=[],
                            data=[],
                            page_size=10,
                            editable=True,
                            style_cell={"textAlign": "left"},
                            style_table={"overflowX": "auto"},
                        ),
                        html.Div(id="metadata-stats-container", className="mt-3"),
                    ]),
                )
            ]),
            dbc.Tab(label="Metatranscriptomics", children=build_metatranscriptomics_tab()),
            dbc.Tab(label="Quantification Settings", children=build_transcriptomics_settings()),
            dbc.Tab(label="DGE Analysis", children=build_dge_analysis_tab()),
            dbc.Tab(label="Microbial transcriptomics", children=html.Div("Content for Microbial transcriptomics", className="p-3")),
            dbc.Tab(label="Host transcriptomics", children=build_host_transcriptomics_tab()),
        ])
    ])
# --- 6. Main Layout ---
app.layout = dbc.Container([
    dbc.Row(dbc.Col(html.H1("MetaLL Pipeline Controller", className="my-4",
                           style={"fontWeight": "700", "color": "#00d9ff", "textShadow": "0 2px 4px rgba(0, 217, 255, 0.3)"}))),
    dcc.Interval(id='status-update', interval=3000, n_intervals=0),
    dcc.Interval(id='mt-update', interval=5000, n_intervals=0),
    html.Div([
        html.Div(id="status-badges-container", children=[
            html.Span(id="status-last-update", children="", style={"marginRight": "15px", "fontSize": "0.7rem", "color": "#888"}),
            html.Span("Preprocessing: ", style={"marginRight": "5px", "fontWeight": "500"}),
            html.Span(id="status-preprocessing-proc", children="0",
                     style={"backgroundColor": "#ff9800", "color": "#000", "padding": "3px 8px",
                            "borderRadius": "4px", "fontWeight": "600", "fontSize": "0.75rem", "marginRight": "5px", "display": "inline-block"}),
            html.Span(id="status-preprocessing-comp", children="0",
                     style={"backgroundColor": "#00ff88", "color": "#000", "padding": "3px 8px",
                            "borderRadius": "4px", "fontWeight": "600", "fontSize": "0.75rem", "marginRight": "5px", "display": "inline-block"}),
            html.Span(id="status-preprocessing-fail", children="0",
                     style={"backgroundColor": "#ff6b6b", "color": "#fff", "padding": "3px 8px",
                            "borderRadius": "4px", "fontWeight": "600", "fontSize": "0.75rem", "marginRight": "20px", "display": "inline-block"}),

            html.Span("TrimGalore: ", style={"marginRight": "5px", "fontWeight": "500"}),
            html.Span(id="status-trimgalore-proc", children="0",
                     style={"backgroundColor": "#ff9800", "color": "#000", "padding": "3px 8px",
                            "borderRadius": "4px", "fontWeight": "600", "fontSize": "0.75rem", "marginRight": "5px", "display": "inline-block"}),
            html.Span(id="status-trimgalore-comp", children="0",
                     style={"backgroundColor": "#00ff88", "color": "#000", "padding": "3px 8px",
                            "borderRadius": "4px", "fontWeight": "600", "fontSize": "0.75rem", "marginRight": "5px", "display": "inline-block"}),
            html.Span(id="status-trimgalore-fail", children="0",
                     style={"backgroundColor": "#ff6b6b", "color": "#fff", "padding": "3px 8px",
                            "borderRadius": "4px", "fontWeight": "600", "fontSize": "0.75rem", "marginRight": "20px", "display": "inline-block"}),

            html.Span("RiboDetector: ", style={"marginRight": "5px", "fontWeight": "500"}),
            html.Span(id="status-ribodetector-proc", children="0",
                     style={"backgroundColor": "#ff9800", "color": "#000", "padding": "3px 8px",
                            "borderRadius": "4px", "fontWeight": "600", "fontSize": "0.75rem", "marginRight": "5px", "display": "inline-block"}),
            html.Span(id="status-ribodetector-comp", children="0",
                     style={"backgroundColor": "#00ff88", "color": "#000", "padding": "3px 8px",
                            "borderRadius": "4px", "fontWeight": "600", "fontSize": "0.75rem", "marginRight": "5px", "display": "inline-block"}),
            html.Span(id="status-ribodetector-fail", children="0",
                     style={"backgroundColor": "#ff6b6b", "color": "#fff", "padding": "3px 8px",
                            "borderRadius": "4px", "fontWeight": "600", "fontSize": "0.75rem", "marginRight": "20px", "display": "inline-block"}),

            html.Span("SortMeRNA: ", style={"marginRight": "5px", "fontWeight": "500"}),
            html.Span(id="status-sortmerna-proc", children="0",
                     style={"backgroundColor": "#ff9800", "color": "#000", "padding": "3px 8px",
                            "borderRadius": "4px", "fontWeight": "600", "fontSize": "0.75rem", "marginRight": "5px", "display": "inline-block"}),
            html.Span(id="status-sortmerna-comp", children="0",
                     style={"backgroundColor": "#00ff88", "color": "#000", "padding": "3px 8px",
                            "borderRadius": "4px", "fontWeight": "600", "fontSize": "0.75rem", "marginRight": "5px", "display": "inline-block"}),
            html.Span(id="status-sortmerna-fail", children="0",
                     style={"backgroundColor": "#ff6b6b", "color": "#fff", "padding": "3px 8px",
                            "borderRadius": "4px", "fontWeight": "600", "fontSize": "0.75rem", "marginRight": "20px", "display": "inline-block"}),

            html.Span("Kraken2: ", style={"marginRight": "5px", "fontWeight": "500"}),
            html.Span(id="status-kraken-proc", children="0",
                     style={"backgroundColor": "#ff9800", "color": "#000", "padding": "3px 8px",
                            "borderRadius": "4px", "fontWeight": "600", "fontSize": "0.75rem", "marginRight": "5px", "display": "inline-block"}),
            html.Span(id="status-kraken-comp", children="0",
                     style={"backgroundColor": "#00ff88", "color": "#000", "padding": "3px 8px",
                            "borderRadius": "4px", "fontWeight": "600", "fontSize": "0.75rem", "marginRight": "5px", "display": "inline-block"}),
            html.Span(id="status-kraken-fail", children="0",
                     style={"backgroundColor": "#ff6b6b", "color": "#fff", "padding": "3px 8px",
                            "borderRadius": "4px", "fontWeight": "600", "fontSize": "0.75rem", "display": "inline-block"})
        ], style={"textAlign": "left", "marginBottom": "15px", "fontSize": "0.85rem", "lineHeight": "2"})
    ]),
    dbc.Tabs([
        dbc.Tab(label="Main Controller", children=build_main_controller_tab(), tab_style={"fontWeight": "500"}),
        dbc.Tab(label="Preprocessing", children=build_preprocessing_tab(), tab_style={"fontWeight": "500"}),
        dbc.Tab(label="TrimGalore", children=build_trimgalore_tab(), tab_style={"fontWeight": "500"}),
        dbc.Tab(label="RiboDetector", children=build_ribodetector_tab(), tab_style={"fontWeight": "500"}),
        dbc.Tab(label="SortMeRNA", children=build_sortmerna_tab(), tab_style={"fontWeight": "500"}),
        dbc.Tab(label="Virulence", children=build_virulence_tab(), tab_style={"fontWeight": "500"}),
        dbc.Tab(label="Results Viewer", children=build_results_viewer_tab(), tab_style={"fontWeight": "500"}),
        dbc.Tab(label="Analysis", children=build_analysis_tab(), tab_style={"fontWeight": "500"}),
    ])
], fluid=True, style={"fontFamily": "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif"})


# --- 6.5. Flask route to serve reports safely ---
@app.server.route('/reports/serve')
def serve_reports():
    run_name = flask.request.args.get('run_name')
    path = flask.request.args.get('path')
    base_outdir = flask.request.args.get('base_outdir')
    
    if not run_name or not path:
        return "Missing run_name or path parameters", 400
        
    report_dir = get_run_outdir(run_name, base_outdir)
    return flask.send_from_directory(report_dir, path)


# --- 7. Callbacks ---

# Callback to populate DGE metadata dropdowns
@app.callback(
    [Output("dge-group-col", "options"),
     Output("batch-col", "options"),
     Output("dge-covariates", "options")],
    [Input("metadata-datatable", "columns")]
)
def update_dge_dropdowns(columns):
    if not columns:
        return [], [], []
    options = [{"label": col["name"], "value": col["id"]} for col in columns]
    return options, options, options

@app.callback(
    [Output("dge-control", "options"),
     Output("dge-treatment", "options")],
    [Input("dge-group-col", "value")],
    [State("metadata-datatable", "data")]
)
def update_dge_group_values(group_col, data):
    if not group_col or not data:
        return [], []
    df = pd.DataFrame(data)
    if group_col not in df.columns:
        return [], []
    unique_values = df[group_col].dropna().unique()
    options = [{"label": str(val), "value": str(val)} for val in unique_values]
    return options, options

# Global variable to track pipeline process
pipeline_process = None
pipeline_log_file = None

# Pipeline execution callback
@app.callback(
    [Output("log-display", "children"),
     Output("pipeline-running", "data"),
     Output("log-update", "disabled")],
    [Input("run-btn", "n_clicks")],
    [State("input-reads", "value"),
     State("slider-cpus", "value"),
     State("drop-mem", "value"),
     State("check-resume", "value"),
     State("workflow-checks", "value"),
     State("output-run-name", "value"),
     State("tg-quality", "value"),
     State("tg-min-len", "value"),
     State("tg-clip-r1", "value"),
     State("tg-clip-r2", "value"),
     State("tg-fastqc", "value"),
     State("rd-len", "value"),
     State("rd-mode", "value"),
     State("rd-device", "value"),
     State("rd-gpu-mem", "value"),
     State("rd-threads", "value"),
     State("smr-coverage", "value"),
     State("smr-mismatch", "value"),
     State("smr-num-align", "value"),
     State("smr-db-dir", "value"),
     State("k2-db-path", "value"),
     State("k2-confidence", "value"),
     State("bracken-level", "value"),
     State("salmon-index", "value"),
     State("virulence-db", "value"),
     State("virulence-index", "value"),
     State("star-index", "value"),
     State("parabricks-index", "value"),
     State("genome-fasta", "value"),
     State("genome-gtf", "value"),
     State("dge-tool", "value"),
     State("dge-control", "value"),
     State("dge-treatment", "value"),
     State("dge-comparison", "value"),
     State("batch-method", "value"),
     State("batch-col", "value"),
     State("dge-p-thresh", "value"),
     State("dge-fc-thresh", "value"),
     State("dge-group-col", "value"),
     State("dge-covariates", "value"),
     State("dge-enrichment-checks", "value"),
     State("dge-ontology", "value"),
     State("dge-organism-db", "value"),
     State("dge-keytype", "value"),
     State("dge-use-biomart", "value"),
     State("dge-auto-install", "value"),
     State("dge-proxy", "value"),
     State("base-outdir", "value"),
     State("perf-checks", "value")],
    prevent_initial_call=True
)
def run_pipeline(n_clicks, reads, cpus, memory, resume, workflow, run_name,
                 tg_quality, tg_min_len, tg_clip_r1, tg_clip_r2, tg_fastqc,
                 rd_len, rd_mode, rd_device, rd_gpu_mem, rd_threads,
                 smr_coverage, smr_mismatch, smr_num_align, smr_db_dir,
                 k2_db_path, k2_confidence, bracken_level,
                 salmon_index, virulence_db, virulence_index, star_index, parabricks_index, genome_fasta, genome_gtf,
                 dge_tool, dge_control, dge_treatment, dge_comparison,
                 batch_method, batch_col, dge_p_thresh, dge_fc_thresh,
                 dge_group_col, dge_covariates, enrichment_checks, ontology, organism_db, keytype,
                 use_biomart, auto_install, proxy, base_outdir, perf_checks):
    global pipeline_process, pipeline_log_file

    if n_clicks is None:
        return no_update, no_update, no_update

    def coalesce(value, default):
        return default if value is None else value

    # If a tab hasn't been visited yet, many State values can be None.
    reads = coalesce(reads, INITIAL_PARAMS["reads"])
    cpus = coalesce(cpus, INITIAL_PARAMS["cpus"])
    memory = coalesce(memory, INITIAL_PARAMS["memory"])
    workflow = coalesce(workflow, [])  # checklist values list

    tg_quality = coalesce(tg_quality, INITIAL_PARAMS["tg_quality"])
    tg_min_len = coalesce(tg_min_len, INITIAL_PARAMS["tg_min_length"])
    tg_clip_r1 = coalesce(tg_clip_r1, INITIAL_PARAMS["tg_clip_r1"])
    tg_clip_r2 = coalesce(tg_clip_r2, INITIAL_PARAMS["tg_clip_r2"])
    tg_fastqc = coalesce(tg_fastqc, [1] if INITIAL_PARAMS["tg_fastqc"] else [])

    rd_len = coalesce(rd_len, INITIAL_PARAMS["rd_len"])
    rd_mode = coalesce(rd_mode, INITIAL_PARAMS["rd_model"])
    rd_device = coalesce(rd_device, INITIAL_PARAMS["rd_device"])
    rd_gpu_mem = coalesce(rd_gpu_mem, INITIAL_PARAMS["rd_gpu_mem"])
    rd_threads = coalesce(rd_threads, INITIAL_PARAMS["rd_threads"])

    smr_coverage = coalesce(smr_coverage, INITIAL_PARAMS["smr_coverage"])
    smr_mismatch = coalesce(smr_mismatch, INITIAL_PARAMS["smr_mismatch"])
    smr_num_align = coalesce(smr_num_align, INITIAL_PARAMS["smr_num_alignments"])
    smr_db_dir = coalesce(smr_db_dir, INITIAL_PARAMS["smr_db_dir"])

    k2_db_path = coalesce(k2_db_path, INITIAL_PARAMS["kraken2_db"])
    k2_confidence = coalesce(k2_confidence, INITIAL_PARAMS["k2_confidence"])
    bracken_level = coalesce(bracken_level, INITIAL_PARAMS["bracken_level"])

    salmon_index = coalesce(salmon_index, INITIAL_PARAMS["salmon_index"])
    virulence_db = coalesce(virulence_db, INITIAL_PARAMS["virulence_db"])
    virulence_index = coalesce(virulence_index, INITIAL_PARAMS["virulence_index"])
    star_index = coalesce(star_index, INITIAL_PARAMS["star_index"])
    parabricks_index = coalesce(parabricks_index, INITIAL_PARAMS["parabricks_index"])
    genome_fasta = coalesce(genome_fasta, INITIAL_PARAMS["genome_fasta"])
    genome_gtf = coalesce(genome_gtf, INITIAL_PARAMS["genome_gtf"])

    dge_tool = coalesce(dge_tool, INITIAL_PARAMS["dge_tool"])
    dge_control = coalesce(dge_control, INITIAL_PARAMS["dge_control"])
    dge_treatment = coalesce(dge_treatment, INITIAL_PARAMS["dge_treatment"])
    dge_comparison = coalesce(dge_comparison, INITIAL_PARAMS["dge_comparison_name"])
    batch_method = coalesce(batch_method, INITIAL_PARAMS["batch_method"])
    dge_p_thresh = coalesce(dge_p_thresh, INITIAL_PARAMS["dge_p_threshold"])
    dge_fc_thresh = coalesce(dge_fc_thresh, INITIAL_PARAMS["dge_fc_threshold"])
    dge_group_col = coalesce(dge_group_col, "group")
    dge_covariates = coalesce(dge_covariates, [])
    perf_checks = coalesce(perf_checks, [])
    use_fastp = "use_fastp" in perf_checks
    use_falco = "use_falco" in perf_checks

    # --- Resolve output directory ---
    outdir = get_run_outdir(run_name, base_outdir)

    if not run_name.strip():
        return "Please provide an output run name (e.g. run_001).", False, True

    # If run_name is not absolute, it should be a simple folder name
    if not os.path.isabs(run_name.strip()):
        bad = (
            "/" in run_name
            or "\\" in run_name
            or run_name.startswith(".")
            or ".." in run_name
        )
        if bad:
            return "Invalid run name. Use a simple folder name like: run_001 (no slashes) or an absolute path.", False, True

    samplesheet_path = os.path.join(outdir, "samplesheet.csv")

    cmd = [
        "nextflow", "run", "./main.nf",
        "-profile", INITIAL_PARAMS["profile"],
        "--reads", str(reads),
        "--outdir", str(outdir),
        "--max_cpus", str(cpus),
        "--max_mem", str(memory),
        "--run_trimgalore", str("run_trimgalore" in workflow).lower(),
        "--run_ribodetector", str("run_ribodetector" in workflow).lower(),
        "--run_sortmerna", str("run_sortmerna" in workflow).lower(),
        "--run_kraken2", str("run_kraken2" in workflow).lower(),
        "--run_parabricks", str("run_parabricks" in workflow).lower(),
        "--run_star", str("run_star" in workflow).lower(),
        "--run_featurecounts", str("run_featurecounts" in workflow).lower(),
        "--run_salmon", str("run_salmon" in workflow).lower(),
        "--run_virulence", str("run_virulence" in workflow).lower(),
        "--run_dge", str("run_dge" in workflow).lower(),

        "--tg_quality", str(tg_quality),
        "--tg_min_length", str(tg_min_len),
        "--tg_clip_r1", str(tg_clip_r1),
        "--tg_clip_r2", str(tg_clip_r2),
        "--tg_fastqc", str(1 in tg_fastqc).lower() if tg_fastqc else "false",

        "--rd_model", str(rd_mode),
        "--rd_len", str(rd_len),
        "--rd_gpu_mem", str(rd_gpu_mem),
        "--rd_threads", str(rd_threads),
        "--rd_device", str(rd_device),

        "--smr_coverage", str(smr_coverage),
        "--smr_mismatch", str(smr_mismatch),
        "--smr_num_alignments", str(smr_num_align),
        "--smr_db_dir", str(smr_db_dir),

        "--kraken2_db", str(k2_db_path),
        "--k2_confidence", str(k2_confidence),
        "--bracken_level", str(bracken_level)
    ]

    # Optional parameters: only add if they are not empty/None
    if salmon_index:
        cmd.extend(["--salmon_index", str(salmon_index)])
    if virulence_db:
        cmd.extend(["--virulence_db", str(virulence_db)])
    if virulence_index:
        cmd.extend(["--virulence_index", str(virulence_index)])
    
    cmd.extend([
        "--star_index", str(star_index),
        "--parabricks_index", str(parabricks_index),
        "--genome_fasta", str(genome_fasta),
        "--genome_gtf", str(genome_gtf)
    ])

    if dge_tool:
        cmd.extend(["--dge_tool", str(dge_tool)])
    if dge_control:
        cmd.extend(["--dge_control", str(dge_control)])
    if dge_treatment:
        cmd.extend(["--dge_treatment", str(dge_treatment)])
    if dge_comparison:
        cmd.extend(["--dge_comparison_name", str(dge_comparison)])
    
    cmd.extend([
        "--group_col", str(dge_group_col) if dge_group_col else "group",
    ])

    if dge_covariates:
        cmd.extend(["--covariates", ",".join(dge_covariates)])
    
    cmd.extend([
        "--batch_method", str(batch_method or "none"),
    ])

    if batch_col:
        cmd.extend(["--batch_col", str(batch_col)])

    cmd.extend([
        "--dge_p_threshold", str(dge_p_thresh or 0.05),
        "--dge_fc_threshold", str(dge_fc_thresh or 1.0),

        "--dge_run_go", str(enrichment_checks and "go" in enrichment_checks).lower(),
        "--dge_run_gsea", str(enrichment_checks and "gsea" in enrichment_checks).lower(),
        "--dge_organism_db", str(organism_db or "org.Hs.eg.db"),
        "--dge_ontology", str(ontology or "BP"),
        "--dge_keytype", str(keytype or "SYMBOL"),
        "--dge_use_biomart", str(use_biomart and 1 in use_biomart).lower(),
        "--dge_autoinstall", str(auto_install and 1 in auto_install).lower()
    ])

    if proxy:
        cmd.extend(["--http_proxy", str(proxy)])

    cmd.extend([
        "--use_fastp", str(use_fastp).lower(),
        "--use_falco", str(use_falco).lower()
    ])

    if resume and 1 in resume:
        cmd.append("-resume")

    try:
        log_path = os.path.join(outdir, "pipeline_dash.log")  # CHANGED
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        pipeline_log_file = open(log_path, "w", buffering=1)

        pipeline_process = subprocess.Popen(
            cmd,
            stdout=pipeline_log_file,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            cwd=os.path.dirname(os.path.abspath(__file__))
        )

        return (
            f"Pipeline started!\n\nCommand: {' '.join(cmd)}\n\nLog file: {log_path}\n\nWaiting for output...",
            True,
            False
        )
    except Exception as e:
        return f"Error starting pipeline: {e}", False, True


# Log update callback - with error protection
@app.callback(
    [Output("log-display", "children", allow_duplicate=True),
     Output("pipeline-running", "data", allow_duplicate=True),
     Output("log-update", "disabled", allow_duplicate=True)],
    [Input("log-update", "n_intervals"),
     Input("pipeline-running", "data"),
     Input("output-run-name", "value")],
    [State("base-outdir", "value")],
    prevent_initial_call=True
)
def update_log(n, is_running, run_name, base_outdir):
    try:
        current_outdir = get_run_outdir(run_name, base_outdir)
        log_path = os.path.join(current_outdir, "pipeline_dash.log")

        # Check if log file exists
        if not os.path.exists(log_path):
            return "Waiting for pipeline output...", True, False

        # Try to read the log file
        try:
            with open(log_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
        except Exception:
            # If we can't read, just skip this update
            return no_update, no_update, no_update

        if not content:
            return "Waiting for pipeline output...", True, False

        # Get last 500 lines to avoid overwhelming the display
        lines = content.split('\n')
        display_lines = lines[-500:] if len(lines) > 500 else lines
        display_text = '\n'.join(display_lines)

        # Check if pipeline is still running
        global pipeline_process
        if pipeline_process is not None:
            poll_result = pipeline_process.poll()
            if poll_result is not None:
                # Pipeline finished - disable interval
                return display_text + "\n\n[Pipeline completed]", False, True

        # Pipeline still running - keep interval active
        return display_text, True, False

    except Exception as e:
        # Catch all exceptions to prevent page refresh
        print(f"Log update error: {e}")
        return no_update, no_update, no_update


# Status counter callback - updates all tab badges
@app.callback(
    [Output("status-last-update", "children"),
     Output("status-preprocessing-proc", "children"),
     Output("status-preprocessing-comp", "children"),
     Output("status-preprocessing-fail", "children"),
     Output("status-trimgalore-proc", "children"),
     Output("status-trimgalore-comp", "children"),
     Output("status-trimgalore-fail", "children"),
     Output("status-ribodetector-proc", "children"),
     Output("status-ribodetector-comp", "children"),
     Output("status-ribodetector-fail", "children"),
     Output("status-sortmerna-proc", "children"),
     Output("status-sortmerna-comp", "children"),
     Output("status-sortmerna-fail", "children"),
     Output("status-kraken-proc", "children"),
     Output("status-kraken-comp", "children"),
     Output("status-kraken-fail", "children")],
    [Input("status-update", "n_intervals"),
     Input("output-run-name", "value")],
    [State("base-outdir", "value")],
    prevent_initial_call=False
)
def update_status_counters(n, run_name, base_outdir):
    try:
        current_outdir = get_run_outdir(run_name, base_outdir)

        import datetime
        timestamp = datetime.datetime.now().strftime("%H:%M:%S")

        # Check samplesheet for total samples
        samplesheet_path = os.path.join(current_outdir, "preprocessing", "samplesheet.csv")
        if not os.path.exists(samplesheet_path):
            return tuple([f"Updated: {timestamp}"] + ["0"] * 15)

        df = pd.read_csv(samplesheet_path)
        total_samples = len(df)

        # Count completed for each stage
        # Preprocessing - samplesheet exists means completed
        prep_comp = total_samples if os.path.exists(samplesheet_path) else 0
        prep_proc = 0
        prep_fail = 0

        # TrimGalore - count completed trimming reports
        tg_results_dir = os.path.join(current_outdir, "trimgalore")
        tg_comp = 0
        tg_started = False
        if os.path.exists(tg_results_dir):
            try:
                reports = [f for f in os.listdir(tg_results_dir) if f.endswith('_trimming_report.txt')]
                tg_comp = len(reports) // 2  # Divide by 2 for R1/R2 pairs
                tg_started = len(reports) > 0  # Stage has started if any reports exist
            except Exception:
                tg_comp = 0
        tg_proc = max(0, total_samples - tg_comp) if tg_started else 0
        tg_fail = 0

        # RiboDetector - count completed log files (one per sample, not per R1/R2)
        rd_results_dir = os.path.join(current_outdir, "ribodetector")
        rd_comp = 0
        rd_started = False
        if os.path.exists(rd_results_dir):
            try:
                logs = [f for f in os.listdir(rd_results_dir) if f.endswith('.log')]
                rd_comp = len(logs)  # One log per sample
                rd_started = len(logs) > 0
            except Exception:
                rd_comp = 0
        rd_proc = max(0, total_samples - rd_comp) if rd_started else 0
        rd_fail = 0

        # SortMeRNA - count completed log files (one per sample, not per R1/R2)
        smr_results_dir = os.path.join(current_outdir, "sortmerna")
        smr_comp = 0
        smr_started = False
        if os.path.exists(smr_results_dir):
            try:
                logs = [f for f in os.listdir(smr_results_dir) if f.endswith('.log')]
                smr_comp = len(logs)  # One log per sample
                smr_started = len(logs) > 0
            except Exception:
                smr_comp = 0
        smr_proc = max(0, total_samples - smr_comp) if smr_started else 0
        smr_fail = 0

        # Kraken2 - count completed report files (one per sample, not per read pair)
        k2_results_dir = os.path.join(current_outdir, "kraken2")
        k2_comp = 0
        k2_started = False
        if os.path.exists(k2_results_dir):
            try:
                # Be flexible with extensions
                reports = [f for f in os.listdir(k2_results_dir) 
                           if f.endswith('.kraken2.report.txt') or f.endswith('.report.txt') or f.endswith('.report')]
                k2_comp = len(reports)
                k2_started = len(reports) > 0
            except Exception:
                k2_comp = 0
        k2_proc = max(0, total_samples - k2_comp) if k2_started else 0
        k2_fail = 0

        return (f"Updated: {timestamp}",
                str(prep_proc), str(prep_comp), str(prep_fail),
                str(tg_proc), str(tg_comp), str(tg_fail),
                str(rd_proc), str(rd_comp), str(rd_fail),
                str(smr_proc), str(smr_comp), str(smr_fail),
                str(k2_proc), str(k2_comp), str(k2_fail))

    except Exception as e:
        # Return zeros on error to prevent refresh
        import datetime
        timestamp = datetime.datetime.now().strftime("%H:%M:%S")
        print(f"Status counter error: {e}")
        return tuple([f"Error: {timestamp}"] + ["0"] * 15)


# Samplesheet display callback
@app.callback(
    Output("samplesheet-table-container", "children"),
    [Input("status-update", "n_intervals"),
     Input("output-run-name", "value")],
    [State("base-outdir", "value")]
)
def update_samplesheet(_, run_name, base_outdir):
    current_outdir = get_run_outdir(run_name, base_outdir)
    samplesheet_path = os.path.join(current_outdir, "preprocessing", "samplesheet.csv")
    if os.path.exists(samplesheet_path):
        try:
            df = pd.read_csv(samplesheet_path)
            return dash_table.DataTable(
                data=df.to_dict('records'),
                columns=[{"name": i, "id": i} for i in df.columns],
                style_cell={'textAlign': 'left', 'padding': '10px'},
                style_header={
                    'backgroundColor': 'rgb(230, 230, 230)',
                    'fontWeight': 'bold'
                },
                page_size=10,
                style_table={'overflowX': 'auto'}
            )
        except Exception as e:
            return dbc.Alert(f"Error reading samplesheet: {e}", color="danger")
    return dbc.Alert("Samplesheet not yet generated. Run the pipeline to discover files.", color="info")


# QC reports callbacks (Updated to use URL routing instead of srcDoc)
@app.callback(
    Output("raw-qc-container", "children"),
    [Input("status-update", "n_intervals"),
     Input("output-run-name", "value")],
    [State("base-outdir", "value")]
)
def update_raw_qc(_, run_name, base_outdir):
    rel_path = "multiqc_raw/multiqc_report_raw.html"
    current_outdir = get_run_outdir(run_name, base_outdir)
    full_path = os.path.join(current_outdir, rel_path)
    if os.path.exists(full_path):
        params = urllib.parse.urlencode({"run_name": run_name, "path": rel_path, "base_outdir": base_outdir})
        return html.Iframe(
            src=f"/reports/serve?{params}",
            style={"width": "100%", "height": "1000px", "border": "none"}
        )
    return dbc.Alert("Raw QC report not found.", color="info")


@app.callback(
    Output("trimmed-qc-container", "children"),
    [Input("status-update", "n_intervals"),
     Input("output-run-name", "value")],
    [State("base-outdir", "value")]
)
def update_trimmed_qc(_, run_name, base_outdir):
    rel_path = "multiqc_trimmed/multiqc_report_trimmed.html"
    current_outdir = get_run_outdir(run_name, base_outdir)
    full_path = os.path.join(current_outdir, rel_path)
    if os.path.exists(full_path):
        params = urllib.parse.urlencode({"run_name": run_name, "path": rel_path, "base_outdir": base_outdir})
        return html.Iframe(
            src=f"/reports/serve?{params}",
            style={"width": "100%", "height": "1000px", "border": "none"}
        )
    return dbc.Alert("Trimmed QC report not found.", color="info")


# RiboDetector stats callback
@app.callback(
    [Output("rd-plots-container", "children"),
     Output("rd-stats-container", "children")],
    [Input("rd-stats-update", "n_intervals"),
     Input("output-run-name", "value")],
    [State("base-outdir", "value")]
)
def update_ribodetector_stats(_, run_name, base_outdir):
    current_outdir = get_run_outdir(run_name, base_outdir)
    rd_results_dir = os.path.join(current_outdir, "ribodetector")
    if not os.path.exists(rd_results_dir):
        return None, dbc.Alert("RiboDetector directory not created yet.", color="info")

    log_files = [f for f in sorted(os.listdir(rd_results_dir)) if f.endswith('.log')]
    if not log_files:
        return None, dbc.Alert("Waiting for RiboDetector logs...", color="info")

    # Patterns for both current and potentially older RiboDetector versions
    total_re = re.compile(r"(?:Processed|Total reads processed[:]?)\s+([\d,]+)", re.IGNORECASE)
    rrna_re = re.compile(r"(?:Detected\s+([\d,]+)\s+rRNA|rRNA reads[:]?\s+([\d,]+))", re.IGNORECASE)
    non_rrna_re = re.compile(r"Detected\s+([\d,]+)\s+non-rRNA", re.IGNORECASE)
    ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
    
    stats = []
    plot_data = []
    
    for log in log_files:
        try:
            log_path = os.path.join(rd_results_dir, log)
            t_val, r_val, nr_val = None, None, None
            
            with open(log_path, 'r', encoding='utf-8', errors='replace') as f:
                content = ansi_escape.sub('', f.read())
            
            t_matches = total_re.findall(content)
            if t_matches:
                t_val = int(t_matches[-1].replace(',', ''))
            
            r_matches = rrna_re.findall(content)
            if r_matches:
                last_match = r_matches[-1]
                if isinstance(last_match, tuple):
                    r_str = last_match[0] or last_match[1]
                else:
                    r_str = last_match
                if r_str:
                    r_val = int(r_str.replace(',', ''))
            
            nr_matches = non_rrna_re.findall(content)
            if nr_matches:
                nr_val = int(nr_matches[-1].replace(',', ''))
            
            if t_val is not None:
                if nr_val is None and r_val is not None:
                    nr_val = t_val - r_val
                elif r_val is None and nr_val is not None:
                    r_val = t_val - nr_val
                
                if r_val is None: r_val = 0
                if nr_val is None: nr_val = 0
                
                sample_name = log.replace("_ribodetector.log", "").replace(".log", "")
                perc = (r_val / t_val) * 100 if t_val > 0 else 0
                
                stats.append({
                    "Sample": sample_name,
                    "Total Reads": f"{t_val:,}",
                    "non-rRNA Reads": f"{nr_val:,}",
                    "rRNA Reads": f"{r_val:,}",
                    "rRNA %": f"{perc:.2f}%"
                })
                
                plot_data.append({"Sample": sample_name, "Reads": nr_val, "Type": "non-rRNA"})
                plot_data.append({"Sample": sample_name, "Reads": r_val, "Type": "rRNA"})
                
        except Exception as e:
            print(f"Error parsing RiboDetector log {log}: {e}")
            continue

    if stats:
        # Table
        df = pd.DataFrame(stats)
        table = dash_table.DataTable(
            data=df.to_dict('records'),
            columns=[{"name": i, "id": i} for i in df.columns],
            sort_action="native",
            filter_action="native",
            style_table={'overflowX': 'auto'},
            style_cell={
                'textAlign': 'left',
                'backgroundColor': '#1e1e1e',
                'color': 'white',
                'padding': '10px',
                'fontSize': '12px'
            },
            style_header={
                'backgroundColor': '#333',
                'color': '#00d9ff',
                'fontWeight': 'bold',
                'border': '1px solid #444'
            },
            style_data={
                'border': '1px solid #444'
            }
        )
        
        # Plot
        df_plot = pd.DataFrame(plot_data)
        fig = px.bar(
            df_plot, x="Sample", y="Reads", color="Type",
            title="RiboDetector Read Distribution",
            barmode="stack",
            template="plotly_dark",
            color_discrete_map={"non-rRNA": "#00d9ff", "rRNA": "#ff4d4d"}
        )
        fig.update_layout(
            font_family="Inter, sans-serif",
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor='rgba(0,0,0,0)',
            title_font_size=20,
            xaxis_title="",
            yaxis_title="Number of Reads"
        )
        
        plot = dcc.Graph(figure=fig, config={'displayModeBar': True})
        
        return plot, table

    return None, dbc.Alert(f"Parsing {len(log_files)} logs in {rd_results_dir}... (No summary data found yet)", color="info")


# TrimGalore stats plot callback - creates multiple plots in grid
@app.callback(
    Output("tg-stats-plots-container", "children"),
    [Input("tg-stats-update", "n_intervals"),
     Input("output-run-name", "value")],
    [State("base-outdir", "value")]
)
def update_trimgalore_stats(_, run_name, base_outdir):
    try:
        current_outdir = get_run_outdir(run_name, base_outdir)
        tg_results_dir = os.path.join(current_outdir, "trimgalore")

        if not os.path.exists(tg_results_dir):
            return dbc.Alert("TrimGalore directory not created yet", color="info")

        # Look for trimming report files
        report_files = [f for f in os.listdir(tg_results_dir) if f.endswith('_trimming_report.txt')]

        if not report_files:
            return dbc.Alert("Waiting for TrimGalore reports...", color="info")

        # Regex patterns for summary parsing
        total_reads_re = re.compile(r"Total reads processed:\s+([\d,]+)")
        reads_adapters_re = re.compile(r"Reads with adapters:\s+[\d,]+\s+\(([\d.]+)%\)")
        reads_written_re = re.compile(r"Reads written \(passing filters\):\s+([\d,]+)\s+\(([\d.]+)%\)")
        quality_trimmed_re = re.compile(r"Quality-trimmed:\s+([\d,]+)\s+bp\s+\(([\d.]+)%\)")
        total_written_re = re.compile(r"Total written \(filtered\):\s+[\d,]+\s+bp\s+\(([\d.]+)%\)")

        # Parse trimming reports and calculate file sizes
        sample_stats = {}

        # Load samplesheet from CURRENT run to get actual input file paths
        samplesheet_path = os.path.join(current_outdir, "preprocessing", "samplesheet.csv")
        file_path_map = {}
        if os.path.exists(samplesheet_path):
            ss_df = pd.read_csv(samplesheet_path)
            for _, row in ss_df.iterrows():
                sample = row['sample']
                file_path_map[sample] = {
                    'R1': row['fastq_1'],
                    'R2': row['fastq_2']
                }

        for report in sorted(report_files):
            # Extract sample name (remove R1/R2 and file extension parts)
            sample_name = report.replace('_trimming_report.txt', '')
            # Remove R1/R2 identifiers
            base_name = re.sub(r'_R[12]_', '_R_', sample_name)

            with open(os.path.join(tg_results_dir, report), 'r') as f:
                content = f.read()

                # Determine if R1 or R2
                read_type = "R1" if "_R1_" in report else "R2"

                # Extract summary stats
                total_reads = 0
                adapters_pct = 0.0
                reads_written_pct = 0.0
                quality_trimmed_bp = 0
                total_written_pct = 0.0

                total_reads_match = total_reads_re.search(content)
                if total_reads_match:
                    total_reads = int(total_reads_match.group(1).replace(",", ""))
                
                reads_adapters_match = reads_adapters_re.search(content)
                if reads_adapters_match:
                    adapters_pct = float(reads_adapters_match.group(1))
                
                reads_written_match = reads_written_re.search(content)
                if reads_written_match:
                    reads_written_pct = float(reads_written_match.group(2))
                
                quality_trimmed_match = quality_trimmed_re.search(content)
                if quality_trimmed_match:
                    quality_trimmed_bp = int(quality_trimmed_match.group(1).replace(",", ""))
                
                total_written_match = total_written_re.search(content)
                if total_written_match:
                    total_written_pct = float(total_written_match.group(1))

                # Find input file size using samplesheet
                input_size_mb = 0

                # Try to match sample from report name to samplesheet
                sample_from_report = report.replace('_trimming_report.txt', '').replace('.fastq.gz', '').replace('.fq.gz', '')
                sample_from_report = re.sub(r'_R[12]_\d+$', '', sample_from_report)

                # Look for matching sample in samplesheet
                input_file_path = None
                for sample_key, paths in file_path_map.items():
                    if sample_key in sample_from_report or sample_from_report in sample_key:
                        input_file_path = paths.get(read_type)
                        break

                if input_file_path and os.path.exists(input_file_path):
                    input_size_mb = os.path.getsize(input_file_path) / (1024 * 1024)

                # Find output file size (trimmed file)
                output_size_mb = 0
                val_suffix = "_val_1.fq.gz" if read_type == "R1" else "_val_2.fq.gz"

                for f in os.listdir(tg_results_dir):
                    if f.endswith(val_suffix) and sample_from_report in f:
                        output_path = os.path.join(tg_results_dir, f)
                        if os.path.exists(output_path):
                            output_size_mb = os.path.getsize(output_path) / (1024 * 1024)
                            break

                size_reduced_mb = max(0, input_size_mb - output_size_mb)

                if base_name not in sample_stats:
                    sample_stats[base_name] = {}

                sample_stats[base_name][read_type] = {
                    "kept_size_mb": output_size_mb,
                    "reduced_size_mb": size_reduced_mb,
                    "input_size_mb": input_size_mb,
                    "total_reads": total_reads,
                    "adapters_pct": adapters_pct,
                    "reads_written_pct": reads_written_pct,
                    "quality_trimmed_bp": quality_trimmed_bp,
                    "total_written_pct": total_written_pct
                }

        if not sample_stats:
            return dbc.Alert("No valid trimming statistics found", color="warning")

        # Build summary data for table
        summary_data = []
        for sample_name, reads_data in sorted(sample_stats.items()):
            for read_type in ['R1', 'R2']:
                if read_type in reads_data:
                    d = reads_data[read_type]
                    summary_data.append({
                        "Sample": sample_name.split('/')[-1],
                        "Read": read_type,
                        "Total Reads": f"{d['total_reads']:,}",
                        "Adapters": f"{d['adapters_pct']:.1f}%",
                        "Reads Written": f"{d['reads_written_pct']:.1f}%",
                        "Quality Trimmed": f"{d['quality_trimmed_bp']:,} bp",
                        "Written Bases": f"{d['total_written_pct']:.1f}%"
                    })

        # Create individual plots for each sample
        plots = []
        for sample_name, reads_data in sorted(sample_stats.items()):
            # Create figure for this sample
            fig = go.Figure()

            read_types = []
            kept_size_vals = []
            reduced_size_vals = []
            input_size_vals = []
            hover_text = []

            for read_type in ['R1', 'R2']:
                if read_type in reads_data:
                    d = reads_data[read_type]
                    read_types.append(read_type)
                    kept_size_vals.append(d['kept_size_mb'])
                    reduced_size_vals.append(d['reduced_size_mb'])
                    input_size_vals.append(d['input_size_mb'])
                    
                    hover_info = (f"Reads: {d['total_reads']:,}<br>"
                                 f"Adapters: {d['adapters_pct']:.1f}%<br>"
                                 f"Written: {d['reads_written_pct']:.1f}%<br>"
                                 f"Qual Trim: {d['quality_trimmed_bp']:,} bp")
                    hover_text.append(hover_info)

            # Add bars - show size after trimming first (bottom), then reduction (top)
            fig.add_trace(go.Bar(
                name='Size After Trimming',
                x=read_types,
                y=kept_size_vals,
                marker_color='#00ff88',
                text=[f"{v:.1f}" for v in kept_size_vals],
                customdata=hover_text,
                textposition='inside',
                hovertemplate='<b>After Trimming</b><br>%{y:.2f} MB<br>%{customdata}<extra></extra>'
            ))

            fig.add_trace(go.Bar(
                name='Size Reduced',
                x=read_types,
                y=reduced_size_vals,
                marker_color='#ff6b6b',
                text=[f"{v:.1f}" for v in reduced_size_vals],
                customdata=hover_text,
                textposition='inside',
                hovertemplate='<b>Size Reduced</b><br>%{y:.2f} MB<br>%{customdata}<extra></extra>'
            ))

            # Create debug title showing input/output sizes
            debug_info = []
            for i, rt in enumerate(read_types):
                debug_info.append(f"{rt}: {input_size_vals[i]:.0f}→{kept_size_vals[i]:.0f}MB")
            
            # Use a cleaner sample name for the title
            display_name = sample_name.split('/')[-1]
            title_text = f"{display_name}<br><sub>({', '.join(debug_info)})</sub>"

            fig.update_layout(
                barmode='stack',
                template="plotly_dark",
                paper_bgcolor='rgba(0,0,0,0)',
                plot_bgcolor='rgba(0,0,0,0)',
                font=dict(color='#ffffff', size=11),
                title=dict(text=title_text, font=dict(size=12, color='#00d9ff')),
                xaxis_title="Read",
                yaxis_title="File Size (MB)",
                showlegend=True,
                legend=dict(
                    orientation="h",
                    yanchor="top",
                    y=-0.15,
                    xanchor="center",
                    x=0.5,
                    font=dict(size=10)
                ),
                height=350,
                margin=dict(l=50, r=20, t=50, b=60)
            )

            plots.append(dbc.Col(dcc.Graph(figure=fig, config={'displayModeBar': False}), width=4, className="mb-3"))

        # Return plots in rows of 3
        rows = []
        for i in range(0, len(plots), 3):
            rows.append(dbc.Row(plots[i:i+3]))

        # Final layout assembly
        return html.Div([
            dbc.Card([
                dbc.CardHeader(html.H6("Aggregated Trimming Summary", className="mb-0")),
                dbc.CardBody([
                    dash_table.DataTable(
                        data=summary_data,
                        columns=[{"name": i, "id": i} for i in summary_data[0].keys()],
                        sort_action="native",
                        filter_action="native",
                        page_size=10,
                        style_table={'overflowX': 'auto'},
                        style_cell={
                            'textAlign': 'left', 
                            'backgroundColor': '#1e1e1e', 
                            'color': 'white', 
                            'fontSize': '11px',
                            'padding': '5px'
                        },
                        style_header={
                            'backgroundColor': '#333',
                            'color': 'white',
                            'fontWeight': 'bold',
                            'border': '1px solid #444'
                        },
                        style_data={
                            'border': '1px solid #444'
                        }
                    )
                ])
            ], className="mb-4", style={"border": "1px solid #444", "backgroundColor": "#1e1e1e"}),
            
            html.H6("Per-Sample Trimming Details", className="mb-3", style={"color": "#00d9ff"}),
            html.Div(rows)
        ])

    except Exception as e:
        import traceback
        print(traceback.format_exc())
        return dbc.Alert(f"Error loading trimming stats: {str(e)[:200]}", color="danger")


# Taxonomy sample dropdown callback
@app.callback(
    [Output("taxonomy-sample-dropdown", "options"),
     Output("taxonomy-sample-dropdown", "value")],
    [Input("mt-update", "n_intervals"),
     Input("output-run-name", "value")],
    [State("taxonomy-sample-dropdown", "value"),
     State("base-outdir", "value")]
)
def refresh_taxonomy_samples(_, run_name, current_value, base_outdir):
    current_outdir = get_run_outdir(run_name, base_outdir)
    res_dir = os.path.join(current_outdir, "bracken")
    options = []
    
    if os.path.exists(res_dir):
        files = [f.replace(".bracken.tsv", "") for f in os.listdir(res_dir) if f.endswith(".bracken.tsv")]
        if files:
            options = [{"label": s, "value": s} for s in sorted(files)]
    
    # If options are empty but kraken2 directory exists, try checking there too
    if not options:
        k2_dir = os.path.join(current_outdir, "kraken2")
        if os.path.exists(k2_dir):
            files = [f.replace(".kraken2.report.txt", "").replace(".report", "") 
                     for f in os.listdir(k2_dir) if f.endswith(".kraken2.report.txt") or f.endswith(".report")]
            if files:
                options = [{"label": s, "value": s} for s in sorted(files)]

    new_value = current_value
    if options:
        val_list = [o["value"] for o in options]
        if not current_value or current_value not in val_list:
            new_value = val_list[0]
    else:
        new_value = None

    if new_value == current_value:
        return options, no_update
    return options, new_value


# Bracken bar plot callback
@app.callback(
    Output("bracken-bar-plot", "figure"),
    [Input("taxonomy-sample-dropdown", "value"),
     Input("output-run-name", "value"),
     Input("mt-update", "n_intervals")],
    [State("base-outdir", "value")]
)
def update_bracken_plot(sample, run_name, _, base_outdir):
    if not sample:
        return go.Figure().add_annotation(text="Select a sample to view data", showarrow=False)

    try:
        current_outdir = get_run_outdir(run_name, base_outdir)
        # Try both direct and with .bracken.tsv extension if sample was cleaned differently
        bracken_file = os.path.join(current_outdir, "bracken", f"{sample}.bracken.tsv")
        if not os.path.exists(bracken_file):
            # Try to see if sample already has the extension or if it's missing
            bracken_file = os.path.join(current_outdir, "bracken", f"{sample}")

        # Fallback to Kraken2 report if Bracken is missing
        if not os.path.exists(bracken_file):
            k2_dir = os.path.join(current_outdir, "kraken2")
            k2_report = None
            for ext in ['.kraken2.report.txt', '.report.txt', '.k2report', '.report']:
                temp_path = os.path.join(k2_dir, f"{sample}{ext}")
                if os.path.exists(temp_path):
                    k2_report = temp_path
                    break
                # Also try without adding extension if sample already has it
                temp_path = os.path.join(k2_dir, sample)
                if os.path.exists(temp_path) and (sample.endswith('.report') or sample.endswith('.txt')):
                    k2_report = temp_path
                    break
            
            if k2_report:
                k2_data = []
                with open(k2_report, 'r') as f:
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) < 6: continue
                        # Kraken2 report: 0:%, 1:clade_reads, 2:taxon_reads, 3:rank, 4:taxid, 5:name
                        if parts[3].strip() == 'S': # Species level
                            k2_data.append({
                                'name': parts[5].strip(), 
                                'fraction_total_reads': float(parts[0]) / 100.0
                            })
                if not k2_data:
                    return go.Figure().add_annotation(text=f"No species-level data for {sample}", showarrow=False)
                df = pd.DataFrame(k2_data)
                title_suffix = "(Kraken2 - Pre-Bracken)"
            else:
                return go.Figure().add_annotation(text=f"No Bracken or Kraken2 report for {sample}", showarrow=False)
        else:
            df = pd.read_csv(bracken_file, sep='\t')
            title_suffix = "(Bracken)"

        if df.empty or 'fraction_total_reads' not in df.columns:
            return go.Figure().add_annotation(text=f"No data available for {sample}", showarrow=False)

        top_10 = df.sort_values(by='fraction_total_reads', ascending=False).head(10)

        fig = px.bar(
            top_10,
            x='fraction_total_reads',
            y='name',
            orientation='h',
            color='fraction_total_reads',
            color_continuous_scale='Magma',
            title=f"Top 10 Species: {sample} {title_suffix}",
            labels={'fraction_total_reads': 'Fraction of Total Reads', 'name': 'Species'}
        )
        fig.update_layout(
            template="plotly_dark",
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            height=600,
            margin=dict(l=200) # Give space for species names
        )
        return fig
    except Exception as e:
        return go.Figure().add_annotation(text=f"Error loading taxonomy data: {str(e)[:50]}", showarrow=False)


# Krona plot callback
@app.callback(
    Output("krona-container", "children"),
    [Input("taxonomy-sample-dropdown", "value"),
     Input("output-run-name", "value"),
     Input("mt-update", "n_intervals")],
    [State("base-outdir", "value")]
)
def update_krona_plot(sample, run_name, _, base_outdir):
    if not sample:
        return dbc.Alert("Select a sample to view Krona plot", color="info")

    try:
        current_outdir = get_run_outdir(run_name, base_outdir)
        k2_dir = os.path.join(current_outdir, "kraken2")
        krona_dir = os.path.join(current_outdir, "krona")
        
        # Potential krona file names and paths
        search_paths = [
            os.path.join(k2_dir, f"{sample}.krona.html"),
            os.path.join(k2_dir, f"{sample}.kraken2.krona.html"),
            os.path.join(k2_dir, f"{sample}.html"),
            os.path.join(krona_dir, f"{sample}.krona.html"),
            os.path.join(krona_dir, f"{sample}.html"),
            # Also try without extensions if sample already has them
            os.path.join(k2_dir, sample if sample.endswith('.html') else '____none____')
        ]
        
        for kp in search_paths:
            if os.path.exists(kp) and os.path.isfile(kp):
                with open(kp, 'r', encoding='utf-8') as f:
                    html_content = f.read()
                return html.Iframe(srcDoc=html_content, style={"width": "100%", "height": "800px", "border": "none"})
        
        return dbc.Alert(f"Krona Plot HTML not found for {sample} (Searched in kraken2/ and krona/ folders)", color="info")
    except Exception as e:
        return dbc.Alert(f"Error loading Krona plot: {str(e)[:100]}", color="danger")

# Overall Abundance Plot callback (CTL vs CASES)
@app.callback(
    Output("overall-abundance-plot", "figure"),
    [Input("output-run-name", "value"),
     Input("mt-update", "n_intervals")],
    [State("base-outdir", "value")]
)
def update_overall_abundance_plot(run_name, _, base_outdir):
    if not run_name:
        return go.Figure()

    try:
        current_outdir = get_run_outdir(run_name, base_outdir)
        bracken_dir = os.path.join(current_outdir, "bracken")

        if not os.path.exists(bracken_dir):
            return go.Figure().add_annotation(text="Bracken results not available yet (check 'bracken' folder)", showarrow=False)

        bracken_files = [f for f in os.listdir(bracken_dir) if f.endswith(".bracken.tsv")]
        if not bracken_files:
            return go.Figure().add_annotation(text="No .bracken.tsv files found in bracken/ folder", showarrow=False)

        all_samples_data = {}
        for bf in bracken_files:
            sample_name = bf.replace(".bracken.tsv", "")
            try:
                df = pd.read_csv(os.path.join(bracken_dir, bf), sep='\t')
                if 'name' in df.columns and 'new_est_reads' in df.columns:
                    all_samples_data[sample_name] = df.groupby('name')['new_est_reads'].sum()
            except Exception:
                continue

        if not all_samples_data:
            return go.Figure().add_annotation(text="Could not extract data from Bracken files", showarrow=False)

        full_df = pd.DataFrame(all_samples_data).fillna(0)
        all_samples = sorted(full_df.columns.tolist())
        
        # Optionally sort by CTL, then CASES, then others to maintain some grouping
        ctl_samples = sorted([s for s in all_samples if "CTL" in s.upper()])
        cases_samples = sorted([s for s in all_samples if "CASES" in s.upper()])
        other_samples = sorted([s for s in all_samples if s not in ctl_samples and s not in cases_samples])
        all_ordered_samples = ctl_samples + cases_samples + other_samples

        if not all_ordered_samples:
            return go.Figure().add_annotation(text="No samples available.", showarrow=False)

        # Exclude Human
        df_no_human = full_df.drop(index="Homo sapiens", errors='ignore')
        total_non_human = df_no_human.sum()
        
        # Select Top 10 Taxa by total abundance across ALL samples
        top_taxa = df_no_human.sum(axis=1).sort_values(ascending=False).head(10).index.tolist()
        
        # Calculate Relative Abundance (%)
        total_non_human_safe = total_non_human.replace(0, 1)
        rel_abund = df_no_human.loc[top_taxa].div(total_non_human_safe) * 100

        fig = make_subplots(
            rows=2, cols=1,
            shared_xaxes=True,
            vertical_spacing=0.08,
            subplot_titles=("Filtered Reads", "Relative Abundance (%)"),
            row_heights=[0.3, 0.7]
        )

        colors = px.colors.qualitative.Alphabet + px.colors.qualitative.Dark24
        taxa_colors = {taxa: colors[i % len(colors)] for i, taxa in enumerate(top_taxa)}

        # Top: Reads
        fig.add_trace(
            go.Bar(x=all_ordered_samples, y=total_non_human[all_ordered_samples] / 1e6, 
                   marker_color='black', showlegend=False, hovertemplate="Sample: %{x}<br>Reads: %{y:.2f}M<extra></extra>"),
            row=1, col=1
        )
        
        # Bottom: Stacked Abundance
        for taxa in top_taxa:
            fig.add_trace(
                go.Bar(x=all_ordered_samples, y=rel_abund.loc[taxa, all_ordered_samples], name=taxa,
                       marker_color=taxa_colors[taxa], legendgroup=taxa, showlegend=True,
                       hovertemplate="Taxa: " + taxa + "<br>Sample: %{x}<br>Abundance: %{y:.2f}%<extra></extra>"),
                row=2, col=1
            )

        fig.update_layout(
            barmode='stack',
            template="plotly_dark",
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            height=700,
            margin=dict(t=50, b=50, l=80, r=20),
            legend=dict(orientation="v", y=0.5, x=1.02, xanchor="left", font=dict(size=10))
        )

        fig.update_yaxes(title_text="Filtered Reads (10⁶)", row=1, col=1)
        fig.update_yaxes(title_text="Relative Abundance (%)", range=[0, 100], row=2, col=1)
        
        max_reads = total_non_human.max() / 1e6 * 1.1 if not total_non_human.empty else 1
        fig.update_yaxes(range=[0, max_reads], row=1, col=1)

        return fig

    except Exception as e:
        return go.Figure().add_annotation(text=f"Error: {str(e)}", showarrow=False)


# Metatranscriptomics update indicator
@app.callback(
    Output("mt-update-indicator", "children"),
    [Input("mt-update", "n_intervals")]
)
def update_mt_indicator(n):
    import datetime
    now = datetime.datetime.now().strftime("%H:%M:%S")
    return f"Last background refresh: {now}"


# Unified Results Viewer callback - Fixed to prevent nextflow inception
@app.callback(
    [Output("results-files-table", "data"),
     Output("mqc-frame-v2", "src"),
     Output("mqc-frame-v2", "style"),
     Output("results-files-container", "style")],
    [Input("results-folder-dropdown", "value"),
     Input("btn-mqc-raw", "n_clicks"),
     Input("btn-mqc-trim", "n_clicks"),
     Input("btn-mqc-ribo", "n_clicks"),
     Input("btn-mqc-smr", "n_clicks"),
     Input("output-run-name", "value")],
    [State("mqc-frame-v2", "src"),
     State("base-outdir", "value")]
)
def unified_results_handler(folder, n1, n2, n3, n4, run_name, current_mqc_src, base_outdir):
    ctx = callback_context
    default_table_style = {"display": "block"}
    default_iframe_style = {"display": "none", "width": "100%", "height": "800px", "border": "none"}

    if not ctx.triggered:
        return [], None, default_iframe_style, default_table_style

    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
    current_outdir = get_run_outdir(run_name, base_outdir)

    # If run_name changed, decide what to refresh
    if trigger_id == "output-run-name":
        if folder:
            # We will fall through to Folder Selection logic
            trigger_id = "results-folder-dropdown"
        elif current_mqc_src and "/reports/serve" in current_mqc_src:
            # Try to identify which button was likely clicked by looking at the path
            if "multiqc_raw" in current_mqc_src: trigger_id = "btn-mqc-raw"
            elif "multiqc_trimmed" in current_mqc_src: trigger_id = "btn-mqc-trim"
            elif "multiqc_ribo" in current_mqc_src: trigger_id = "btn-mqc-ribo"
            elif "multiqc_sortmerna" in current_mqc_src: trigger_id = "btn-mqc-smr"
            else:
                return [], None, default_iframe_style, default_table_style

    # Handle MultiQC Buttons - Use src instead of srcDoc to prevent inception
    if "btn-mqc" in trigger_id:
        report_map = {
            "btn-mqc-raw": "multiqc_raw/multiqc_report_raw.html",
            "btn-mqc-trim": "multiqc_trimmed/multiqc_report_trimmed.html",
            "btn-mqc-ribo": "multiqc_ribo/multiqc_report_ribo.html",
            "btn-mqc-smr": "multiqc_sortmerna/multiqc_report_sortmerna.html"
        }
        rel_path = report_map.get(trigger_id)
        full_path = os.path.join(current_outdir, rel_path)

        if full_path and os.path.exists(full_path):
            params = urllib.parse.urlencode({"run_name": run_name, "path": rel_path, "base_outdir": base_outdir})
            return [], f"/reports/serve?{params}", {**default_iframe_style, "display": "block"}, {"display": "none"}

        return [], None, {**default_iframe_style, "display": "block"}, {"display": "none"}

    # Handle Folder Selection (File Table)
    if trigger_id == "results-folder-dropdown" and folder:
        path = os.path.join(current_outdir, folder)
        files = []
        if os.path.exists(path):
            for f in sorted(os.listdir(path)):
                full_p = os.path.join(path, f)
                if os.path.isfile(full_p):
                    size = os.path.getsize(full_p) / (1024 * 1024)
                    files.append({"name": f, "size": round(size, 2)})

        return files, None, default_iframe_style, default_table_style

    return [], None, default_iframe_style, default_table_style


# Metadata upload callback
@app.callback(
    [Output("metadata-datatable", "data"),
     Output("metadata-datatable", "columns")],
    [Input("upload-metadata", "contents")],
    [State("upload-metadata", "filename")]
)
def update_metadata_table(contents, filename):
    if contents is None:
        return [], []

    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)

    try:
        if filename.endswith('.csv'):
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
        elif filename.endswith('.tsv') or filename.endswith('.txt'):
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep='\t')
        else:
            return [], []

        columns = [{"name": i, "id": i} for i in df.columns]
        data = df.to_dict('records')
        return data, columns
    except Exception as e:
        print(f"Error reading file: {e}")
        return [], []


# Metadata statistics callback
@app.callback(
    Output("metadata-stats-container", "children"),
    Input("show-stats-check", "value"),
    State("metadata-datatable", "data")
)
def update_metadata_stats(show_stats, data):
    if not show_stats or 1 not in show_stats or not data:
        return ""

    df = pd.DataFrame(data)
    stats_df = df.describe(include='all').fillna("N/A").T

    table_header = [html.Thead(html.Tr([html.Th(col) for col in stats_df.columns]))]
    table_body = [html.Tbody([
        html.Tr([html.Td(stats_df.iloc[i][col]) for col in stats_df.columns])
        for i in range(len(stats_df))
    ])]

    return dbc.Table(table_header + table_body, bordered=True, hover=True, striped=True, responsive=True)


# Metadata download callback
@app.callback(
    Output("download-metadata-csv", "data"),
    Input("btn-csv-download", "n_clicks"),
    State("metadata-datatable", "data"),
    prevent_initial_call=True,
)
def download_metadata(n_clicks, data):
    if not data:
        return None

    df = pd.DataFrame(data)
    csv_string = df.to_csv(index=False, encoding='utf-8')
    return dict(content=csv_string, filename="metadata.csv")


# Host Transcriptomics Coverage callback
@app.callback(
    [Output("host-coverage-bar-plot", "figure"),
     Output("host-coverage-dist-plot", "figure"),
     Output("host-analysis-info", "children")],
    [Input("btn-refresh-host", "n_clicks"),
     Input("output-run-name", "value"),
     Input("host-min-reads", "value")],
    [State("rd-len", "value"),
     State("base-outdir", "value")]
)
def update_host_coverage(n_clicks, run_name, min_reads, read_len, base_outdir):
    if not run_name:
        return go.Figure(), go.Figure(), "Select a run name to begin."
    
    if read_len is None:
        read_len = 150 # Default
        
    try:
        current_outdir = get_run_outdir(run_name, base_outdir)
        # Search for results in various possible locations
        res_dirs = [
            os.path.join(current_outdir, "parabricks"),
            os.path.join(current_outdir, "featurecounts"),
            os.path.join(current_outdir, "salmon")
        ]
        
        fc_files = []
        for d in res_dirs:
            if not os.path.exists(d): continue
            
            if d.endswith("salmon"):
                # Salmon files are in subfolders
                for s in os.listdir(d):
                    s_path = os.path.join(d, s)
                    if os.path.isdir(s_path):
                        qf = os.path.join(s_path, "quant.sf")
                        if os.path.exists(qf): fc_files.append(qf)
            else:
                for f in os.listdir(d):
                    if f.endswith('.featureCounts.txt'):
                        fc_files.append(os.path.join(d, f))
        
        if not fc_files:
            return go.Figure(), go.Figure(), "No featureCounts or Salmon quantification files found."

        all_stats = []
        all_coverage_data = []

        for fc_file in sorted(fc_files):
            try:
                if fc_file.endswith('.featureCounts.txt'):
                    df = pd.read_csv(fc_file, sep='\t', comment='#')
                    if 'Length' in df.columns:
                        counts = df.iloc[:, -1]
                        lengths = df['Length']
                        sample_name = os.path.basename(fc_file).replace('.featureCounts.txt', '').replace('_featurecounts', '')
                    else: continue
                elif fc_file.endswith('quant.sf'):
                    df = pd.read_csv(fc_file, sep='\t')
                    counts = df['NumReads']
                    lengths = df['Length']
                    sample_name = os.path.basename(os.path.dirname(fc_file)).replace('_salmon', '')
                else: continue

                mask = counts >= min_reads
                valid_counts = counts[mask]
                valid_lengths = lengths[mask]
                
                if len(valid_counts) == 0: continue

                # Coverage estimate = (Count * ReadLen) / GeneLength
                coverage = (valid_counts * read_len) / valid_lengths
                
                avg_cov = coverage.mean()
                med_cov = coverage.median()
                
                all_stats.append({
                    "Sample": sample_name,
                    "Average Coverage": avg_cov,
                    "Median Coverage": med_cov,
                    "Genes Counted": len(valid_counts)
                })
                
                # Subset for distribution plot
                if len(coverage) > 3000:
                    subset = coverage.sample(3000, random_state=42)
                else:
                    subset = coverage
                
                for val in subset:
                    all_coverage_data.append({"Sample": sample_name, "Coverage": val})
                    
            except Exception as e:
                print(f"Error processing {fc_file}: {e}")
                continue

        if not all_stats:
            return go.Figure(), go.Figure(), "Could not calculate coverage from available files."

        df_stats = pd.DataFrame(all_stats)
        
        fig_bar = go.Figure()
        fig_bar.add_trace(go.Bar(x=df_stats["Sample"], y=df_stats["Average Coverage"], name="Average Coverage", marker_color='#00d9ff'))
        fig_bar.add_trace(go.Bar(x=df_stats["Sample"], y=df_stats["Median Coverage"], name="Median Coverage", marker_color='#00ff88'))
        
        fig_bar.update_layout(
            title="Average and Median Transcriptome Coverage",
            template="plotly_dark",
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            barmode='group',
            xaxis_title="Sample",
            yaxis_title="Coverage (Fold)"
        )

        df_dist = pd.DataFrame(all_coverage_data)
        fig_dist = px.box(df_dist, x="Sample", y="Coverage", color="Sample", 
                          title=f"Coverage Distribution (subset, min_reads={min_reads})",
                          template="plotly_dark",
                          points=False)
        
        fig_dist.update_layout(
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            yaxis_type="log" if df_dist["Coverage"].max() > 10 else "linear"
        )

        info_text = f"Analyzed {len(all_stats)} samples. Found {len(fc_files)} quantification files."
        return fig_bar, fig_dist, info_text

    except Exception as e:
        return go.Figure(), go.Figure(), f"Error during analysis: {e}"


# Metatranscriptomics pie charts callback
@app.callback(
    Output("mt-pie-charts-container", "children"),
    [Input("mt-update", "n_intervals"),
     Input("output-run-name", "value")],
    [State("base-outdir", "value")]
)
def update_mt_pie_charts(_, run_name, base_outdir):
    try:
        current_outdir = get_run_outdir(run_name, base_outdir)
        k2_results_dir = os.path.join(current_outdir, "kraken2")

        if not os.path.exists(k2_results_dir):
            return dbc.Alert("Kraken2 results not available yet", color="info")

        # Find all Kraken2 report files - prioritized extensions
        report_extensions = ['.kraken2.report.txt', '.report.txt', '.k2report', '.report']
        report_files = [f for f in os.listdir(k2_results_dir) if any(f.endswith(ext) for ext in report_extensions)]

        if not report_files:
            return dbc.Alert("Waiting for Kraken2 reports...", color="info")

        all_samples_data = []
        pie_charts = []

        # Limit individual pie charts if there are many samples to prevent browser lag
        MAX_PIES = 12
        sorted_files = sorted(report_files)
        
        for i, report_file in enumerate(sorted_files):
            # Consistent sample name extraction
            sample_name = report_file
            for ext in sorted(report_extensions, key=len, reverse=True):
                if sample_name.endswith(ext):
                    sample_name = sample_name[:-len(ext)]
                    break
            
            report_path = os.path.join(k2_results_dir, report_file)
            
            # Skip if it's not a file (could be krona dir if we are not careful)
            if not os.path.isfile(report_path): continue

            host_reads = 0
            bacterial_reads = 0
            viral_reads = 0
            unclassified_reads = 0
            total_classified = 0

            try:
                with open(report_path, 'r') as f:
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) < 6: continue
                        
                        num_reads = int(parts[1])
                        rank_code = parts[3].strip()
                        taxon_id = parts[4].strip()
                        
                        if rank_code == 'U':
                            unclassified_reads = num_reads
                        elif taxon_id == '1': # Root
                            total_classified = num_reads
                        elif taxon_id == '2': # Bacteria
                            bacterial_reads = num_reads
                        elif taxon_id == '10239': # Viruses
                            viral_reads = num_reads
                        elif taxon_id == '2759' or taxon_id == '9606': # Eukaryota or Human
                            host_reads = num_reads
                
                # Archaea and others
                other_classified = total_classified - (bacterial_reads + viral_reads + host_reads)
                if other_classified < 0: other_classified = 0
                
                total_reads = total_classified + unclassified_reads
                if total_reads == 0: continue

                all_samples_data.append({
                    "Sample": sample_name,
                    "Host": host_reads,
                    "Bacterial": bacterial_reads,
                    "Viral": viral_reads,
                    "Other Classified": other_classified,
                    "Unclassified": unclassified_reads,
                    "Total": total_reads
                })

                # Only create pie charts for the first MAX_PIES samples
                if i < MAX_PIES:
                    fig = go.Figure(data=[go.Pie(
                        labels=['Host', 'Bacterial', 'Viral', 'Other', 'Unclassified'],
                        values=[host_reads, bacterial_reads, viral_reads, other_classified, unclassified_reads],
                        marker=dict(colors=['#00d9ff', '#00ff88', '#ff6b6b', '#ff9800', '#888888']),
                        hole=0.4,
                        textinfo='percent',
                        showlegend=True
                    )])

                    fig.update_layout(
                        title=dict(text=f"<b>{sample_name}</b>", font=dict(size=12, color='#00d9ff'), x=0.5, xanchor='center'),
                        template="plotly_dark",
                        paper_bgcolor='rgba(0,0,0,0)',
                        plot_bgcolor='rgba(0,0,0,0)',
                        margin=dict(t=40, b=10, l=10, r=10),
                        height=300,
                        legend=dict(orientation="h", yanchor="bottom", y=-0.3, xanchor="center", x=0.5, font=dict(size=9))
                    )
                    pie_charts.append(dbc.Col(dcc.Graph(figure=fig, config={'displayModeBar': False}), width=4, className="mb-4"))
            except Exception as e:
                print(f"Error parsing Kraken2 report {report_file}: {e}")
                continue

        if not all_samples_data:
            return dbc.Alert("Could not parse any Kraken2 reports.", color="warning")

        # Create Summary Charts
        df = pd.DataFrame(all_samples_data)
        
        # 1. Percentage Composition Bar Chart
        df_pct = df.copy()
        for col in ['Host', 'Bacterial', 'Viral', 'Other Classified', 'Unclassified']:
            df_pct[col] = (df_pct[col] / df_pct['Total']) * 100
            
        fig_pct = px.bar(df_pct, x="Sample", y=['Host', 'Bacterial', 'Viral', 'Other Classified', 'Unclassified'],
                         title="Taxonomic Composition (%)",
                         color_discrete_map={'Host': '#00d9ff', 'Bacterial': '#00ff88', 'Viral': '#ff6b6b', 
                                            'Other Classified': '#ff9800', 'Unclassified': '#888888'},
                         template="plotly_dark")
        fig_pct.update_layout(paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)', barmode='stack', legend_title_text='')

        # 2. Absolute Reads Bar Chart
        fig_abs = px.bar(df, x="Sample", y=['Host', 'Bacterial', 'Viral', 'Other Classified', 'Unclassified'],
                         title="Absolute Read Distribution",
                         color_discrete_map={'Host': '#00d9ff', 'Bacterial': '#00ff88', 'Viral': '#ff6b6b', 
                                            'Other Classified': '#ff9800', 'Unclassified': '#888888'},
                         template="plotly_dark")
        fig_abs.update_layout(paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)', barmode='group', legend_title_text='')

        summary_row = dbc.Row([
            dbc.Col(dcc.Graph(figure=fig_pct), lg=6, md=12),
            dbc.Col(dcc.Graph(figure=fig_abs), lg=6, md=12)
        ], className="mb-5 border-bottom pb-4")

        # Sample Detail Header
        detail_msg = f"Individual Sample Detail (showing first {MAX_PIES} of {len(all_samples_data)})" if len(all_samples_data) > MAX_PIES else "Individual Sample Detail"
        detail_header = html.H5(detail_msg, className="mt-4 mb-3", style={"color": "#00d9ff"})

        pie_rows = []
        for i in range(0, len(pie_charts), 3):
            pie_rows.append(dbc.Row(pie_charts[i:i+3]))

        return html.Div([summary_row, detail_header, html.Div(pie_rows)])

    except Exception as e:
        return dbc.Alert(f"Error loading metatranscriptomics data: {str(e)[:100]}", color="danger")


# --- 7.5. DGE Analysis Callbacks ---
import time
import json

# DGE Results Update Callback
@app.callback(
    [Output("dge-comparison-dropdown", "options"),
     Output("dge-comparison-dropdown", "value"),
     Output("dge-analysis-info", "children", allow_duplicate=True),
     Output("dge-run-state", "data", allow_duplicate=True)],
    [Input("dge-refresh-interval", "n_intervals"),
     Input("btn-refresh-dge", "n_clicks"),
     Input("output-run-name", "value")],
    [State("dge-run-state", "data"),
     State("dge-comparison-dropdown", "value"),
     State("base-outdir", "value")],
    prevent_initial_call="initial_duplicate"
)
def update_dge_results(_n, _btn, run_name, state, current_comparison, base_outdir):
    """
    Auto-populate comparison dropdown and poll for run completion.
    """
    current_outdir = get_run_outdir(run_name, base_outdir)
    if not run_name or not os.path.exists(current_outdir):
        return [], None, no_update, {"running": False}

    dge_root = os.path.join(current_outdir, "dge_analysis")

    if not os.path.exists(dge_root):
        os.makedirs(dge_root, exist_ok=True)
        return [], None, no_update, {"running": False}

    # List comparisons that exist
    comparisons = [
        d for d in os.listdir(dge_root)
        if os.path.isdir(os.path.join(dge_root, d)) and not d.startswith("_")
    ]
    options = [{"label": comp, "value": comp} for comp in sorted(comparisons)]
    
    # Logic to handle comparison selection
    ctx = callback_context
    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0] if ctx.triggered else ""
    
    new_comparison = current_comparison
    if trigger_id == "output-run-name" or not current_comparison or current_comparison not in comparisons:
        new_comparison = comparisons[0] if comparisons else None

    info = no_update
    new_state = state or {"running": False}

    # If we are polling for a running process
    if new_state.get("running"):
        # Check if the expected comparison directory and summary file exist
        comp_name = new_state.get("comparison_name")
        if comp_name:
            summary_file = os.path.join(dge_root, comp_name, "dge_summary.json")
            if os.path.exists(summary_file):
                info = dbc.Alert(
                    f"Analysis for '{comp_name}' completed successfully!",
                    color="success",
                    dismissable=True
                )
                new_state["running"] = False
                # If this was the one running, select it
                new_comparison = comp_name
            else:
                info = dbc.Alert(
                    [dbc.Spinner(size="sm"), f" Running DGE analysis: {comp_name}..."],
                    color="info"
                )

    return options, new_comparison, info, new_state


@app.callback(
    Output("dge-summary-stats", "children"),
    [Input("dge-comparison-dropdown", "value"),
     Input("output-run-name", "value")],
    [State("base-outdir", "value")]
)
def display_dge_summary(comparison, run_name, base_outdir):
    if not comparison:
        return dbc.Alert("Select a comparison to view results", color="info")

    current_outdir = get_run_outdir(run_name, base_outdir)
    summary_file = os.path.join(current_outdir, "dge_analysis", comparison, "dge_summary.json")

    if not os.path.exists(summary_file):
        return dbc.Alert("Summary file (dge_summary.json) not found for this comparison.", color="warning")

    try:
        with open(summary_file, 'r') as f:
            summary = json.load(f)

        params = summary['analysis_parameters']
        stats = summary['summary_statistics']

        # Create summary cards
        cards = html.Div([
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H4(stats['total_genes'], className="text-info"),
                            html.P("Total Genes", className="mb-0")
                        ])
                    ], className="mb-2 text-center")
                ], width=3),
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H4(stats['significant_genes'], className="text-success"),
                            html.P("Significant (adj.P < thresh)", className="mb-0")
                        ])
                    ], className="mb-2 text-center")
                ], width=3),
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H4(stats['upregulated'], className="text-danger"),
                            html.P("Upregulated", className="mb-0")
                        ])
                    ], className="mb-2 text-center")
                ], width=3),
                dbc.Col([
                    dbc.Card([
                        dbc.CardBody([
                            html.H4(stats['downregulated'], className="text-primary"),
                            html.P("Downregulated", className="mb-0")
                        ])
                    ], className="mb-2 text-center")
                ], width=3),
            ]),
            html.Div([
                html.P([html.Strong("Control: "), params['control_group']], className="mb-1"),
                html.P([html.Strong("Treatment: "), params['treatment_group']], className="mb-1"),
                html.P([html.Strong("P-value Threshold: "), f"{params['p_threshold']}"], className="mb-1"),
                html.P([html.Strong("Log2FC Threshold: "), f"{params['fc_threshold']}"], className="mb-0")
            ], className="mt-3 p-3 bg-light border rounded", style={"fontSize": "0.9rem"})
        ])

        return cards

    except Exception as e:
        return dbc.Alert(f"Error loading summary: {str(e)}", color="danger")
# Load volcano plot
@app.callback(
    Output("dge-volcano-iframe", "src"),
    [Input("dge-comparison-dropdown", "value"),
     Input("output-run-name", "value")],
    [State("base-outdir", "value")]
)
def load_volcano_plot(comparison, run_name, base_outdir):
    if not comparison:
        return None

    volcano_file = f"dge_analysis/{comparison}/volcano_plot.html"
    current_outdir = get_run_outdir(run_name, base_outdir)
    full_path = os.path.join(current_outdir, volcano_file)

    if os.path.exists(full_path):
        params = urllib.parse.urlencode({"run_name": run_name, "path": volcano_file, "base_outdir": base_outdir})
        return f"/reports/serve?{params}"
    return None


# Load MA plot
@app.callback(
    Output("dge-ma-iframe", "src"),
    [Input("dge-comparison-dropdown", "value"),
     Input("output-run-name", "value")],
    [State("base-outdir", "value")]
)
def load_ma_plot(comparison, run_name, base_outdir):
    if not comparison:
        return None

    ma_file = f"dge_analysis/{comparison}/ma_plot.html"
    current_outdir = get_run_outdir(run_name, base_outdir)
    full_path = os.path.join(current_outdir, ma_file)

    if os.path.exists(full_path):
        params = urllib.parse.urlencode({"run_name": run_name, "path": ma_file, "base_outdir": base_outdir})
        return f"/reports/serve?{params}"
    return None


# Load PCA plot
@app.callback(
    Output("dge-pca-iframe", "src"),
    [Input("dge-comparison-dropdown", "value"),
     Input("output-run-name", "value")],
    [State("base-outdir", "value")]
)
def load_pca_plot(comparison, run_name, base_outdir):
    if not comparison:
        return None

    pca_file = f"dge_analysis/{comparison}/pca_plot.html"
    current_outdir = get_run_outdir(run_name, base_outdir)
    full_path = os.path.join(current_outdir, pca_file)

    if os.path.exists(full_path):
        params = urllib.parse.urlencode({"run_name": run_name, "path": pca_file, "base_outdir": base_outdir})
        return f"/reports/serve?{params}"
    return None


# Load scree plot
@app.callback(
    Output("dge-scree-iframe", "src"),
    [Input("dge-comparison-dropdown", "value"),
     Input("output-run-name", "value")],
    [State("base-outdir", "value")]
)
def load_scree_plot(comparison, run_name, base_outdir):
    if not comparison:
        return None

    scree_file = f"dge_analysis/{comparison}/pca_scree_plot.html"
    current_outdir = get_run_outdir(run_name, base_outdir)
    full_path = os.path.join(current_outdir, scree_file)

    if os.path.exists(full_path):
        params = urllib.parse.urlencode({"run_name": run_name, "path": scree_file, "base_outdir": base_outdir})
        return f"/reports/serve?{params}"
    return None


# Load heatmap
@app.callback(
    Output("dge-heatmap-iframe", "src"),
    [Input("dge-comparison-dropdown", "value"),
     Input("output-run-name", "value")],
    [State("base-outdir", "value")]
)
def load_heatmap(comparison, run_name, base_outdir):
    if not comparison:
        return None

    heatmap_file = f"dge_analysis/{comparison}/heatmap.html"
    current_outdir = get_run_outdir(run_name, base_outdir)
    full_path = os.path.join(current_outdir, heatmap_file)

    if os.path.exists(full_path):
        params = urllib.parse.urlencode({"run_name": run_name, "path": heatmap_file, "base_outdir": base_outdir})
        return f"/reports/serve?{params}"
    return None


# Load top 50 DEGs table
@app.callback(
    [Output("dge-top50-table", "data"),
     Output("dge-top50-table", "columns")],
    [Input("dge-comparison-dropdown", "value"),
     Input("output-run-name", "value")],
    [State("base-outdir", "value")]
)
def load_top50_table(comparison, run_name, base_outdir):
    if not comparison:
        return [], []

    current_outdir = get_run_outdir(run_name, base_outdir)
    top50_file = os.path.join(current_outdir, "dge_analysis", comparison, "top50_degs.csv")

    if not os.path.exists(top50_file):
        return [], []

    try:
        df = pd.read_csv(top50_file)
        # Round numeric columns
        numeric_cols = ['log2FoldChange', 'pvalue', 'padj', 'baseMean_control', 'baseMean_treatment']
        for col in numeric_cols:
            if col in df.columns:
                df[col] = df[col].round(4)

        columns = [{"name": c, "id": c} for c in df.columns]
        data = df.to_dict('records')
        return data, columns

    except Exception as e:
        print(f"Error loading top50 table: {e}")
        return [], []


# Load markdown report
@app.callback(
    [Output("dge-markdown-iframe", "src"),
     Output("dge-markdown-missing", "style"),
     Output("dge-markdown-iframe", "style")],
    [Input("dge-comparison-dropdown", "value"),
     Input("output-run-name", "value")],
    [State("base-outdir", "value")]
)
def load_markdown_report(comparison, run_name, base_outdir):
    base_style_iframe = {"width": "100%", "height": "1000px", "border": "none"}
    if not comparison:
        return None, {"display": "block"}, {**base_style_iframe, "display": "none"}

    rel_path = f"dge_analysis/{comparison}/dge_report.html"
    current_outdir = get_run_outdir(run_name, base_outdir)
    full_path = os.path.join(current_outdir, rel_path)

    if os.path.exists(full_path):
        params = urllib.parse.urlencode({"run_name": run_name, "path": rel_path, "base_outdir": base_outdir})
        return f"/reports/serve?{params}", {"display": "none"}, {**base_style_iframe, "display": "block"}

    return None, {"display": "block"}, {**base_style_iframe, "display": "none"}


# Trigger DGE Analysis manually
@app.callback(
    [Output("dge-analysis-info", "children", allow_duplicate=True),
     Output("dge-run-state", "data", allow_duplicate=True)],
    [Input("btn-run-dge", "n_clicks")],
    [State("output-run-name", "value"),
     State("dge-tool", "value"),
     State("dge-group-col", "value"),
     State("dge-control", "value"),
     State("dge-treatment", "value"),
     State("dge-covariates", "value"),
     State("dge-comparison", "value"),
     State("batch-method", "value"),
     State("batch-col", "value"),
     State("dge-p-thresh", "value"),
     State("dge-fc-thresh", "value"),
     State("dge-enrichment-checks", "value"),
     State("dge-ontology", "value"),
     State("dge-organism-db", "value"),
     State("dge-keytype", "value"),
     State("metadata-datatable", "data"),
     State("base-outdir", "value"),
     State("dge-use-biomart", "value"),
     State("dge-auto-install", "value"),
     State("dge-proxy", "value")],
    prevent_initial_call=True
)
def trigger_dge_analysis(n_clicks, run_name, tool, group_col, control, treatment,
                         covariates, comparison, batch_method, batch_col,
                         p_thresh, fc_thresh, enrichment_checks, ontology, organism_db, keytype,
                         metadata_data, base_outdir, use_biomart, auto_install, proxy):
    if not n_clicks:
        return no_update, no_update

    if not run_name or not comparison:
        return dbc.Alert("Run Name and Comparison Name are required.", color="danger"), no_update

    if not group_col:
        return dbc.Alert("Please select a Group Column first.", color="danger"), no_update

    if control is None or treatment is None:
        return dbc.Alert("Please select BOTH Control and Treatment group values.", color="danger"), no_update

    if str(control).strip() == str(treatment).strip():
        return dbc.Alert("Control and Treatment must be different values.", color="danger"), no_update

    current_outdir = get_run_outdir(run_name, base_outdir)
    outdir = os.path.join(current_outdir, "dge_analysis", comparison)
    os.makedirs(outdir, exist_ok=True)

    # Locate counts file using robust finder
    counts_file = find_counts_file(current_outdir)
    if not counts_file:
         # Fallback to display the attempted path in error
         attempted = os.path.join(current_outdir, "counts", "gene_counts_matrix.tsv")
         return dbc.Alert(f"Counts file not found at {attempted}. Please run the quantification pipeline first.", color="danger"), no_update

    # Save metadata to temp file
    if not metadata_data:
        return dbc.Alert("Metadata is missing. Please upload metadata first.", color="danger"), no_update

    df_meta = pd.DataFrame(metadata_data)
    meta_path = os.path.join(outdir, "metadata_dge.csv")
    df_meta.to_csv(meta_path, index=False)

    # Build command
    r_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "bin", "run_dge.R")
    cmd = [
        "Rscript", r_script,
        "--counts", counts_file,
        "--metadata", meta_path,
        "--control", str(control),
        "--treatment", str(treatment),
        "--outdir", outdir,
        "--tool", tool or "deseq2",
        "--group_col", group_col or "group",
        "--p_threshold", str(p_thresh or 0.05),
        "--fc_threshold", str(fc_thresh or 1.0)
    ]
    if batch_method and batch_method != "none":
        cmd.extend(["--batch_method", batch_method])
        if batch_col:
            cmd.extend(["--batch_col", batch_col])
    if covariates:
        cmd.extend(["--covariates", ",".join(covariates)])
    # Enrichment flags
    if enrichment_checks and ("go" in enrichment_checks or "gsea" in enrichment_checks):
        if "go" in enrichment_checks:
            cmd.append("--run_go")
        if "gsea" in enrichment_checks:
            cmd.append("--run_gsea")
        if ontology:
            cmd.extend(["--ontology", ontology])
        if organism_db:
            cmd.extend(["--organism_db", organism_db])
        if keytype:
            cmd.extend(["--keytype", keytype])
    # biomaRt for plot labels
    if use_biomart and 1 in use_biomart:
        cmd.append("--use_biomart")
    
    # Environment variables for R auto-install and Proxy
    env = os.environ.copy()
    if not auto_install or 1 not in auto_install:
        env["METALL_NO_R_AUTOINSTALL"] = "1"
    
    if proxy:
        env["http_proxy"] = proxy
        env["https_proxy"] = proxy

    try:
        log_path = os.path.join(outdir, "dge_run.log")
        with open(log_path, "w") as f:
            subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT, env=env)

        return dbc.Alert([dbc.Spinner(size="sm"), f" Started DGE analysis: {comparison}"], color="info"), \
               {"running": True, "comparison_name": comparison, "started_at": time.time(),
                "run_name": run_name, "base_outdir": base_outdir, "outdir": outdir, "log_path": log_path}
    except Exception as e:
        return dbc.Alert(f"Error starting DGE: {e}", color="danger"), {"running": False}

# Enable/disable DGE log polling based on run state
@app.callback(
    Output("dge-log-interval", "disabled"),
    Input("dge-run-state", "data")
)
def toggle_dge_log_interval(state):
    try:
        if not state:
            return True
        return not bool(state.get("running", False))
    except Exception:
        return True

# Stream DGE log content
@app.callback(
    [Output("dge-live-log", "children"),
     Output("dge-run-state", "data", allow_duplicate=True)],
    Input("dge-log-interval", "n_intervals"),
    State("dge-run-state", "data"),
    prevent_initial_call=True
)
def stream_dge_log(_n, state):
    if not state:
        return "", {"running": False}
    log_path = state.get("log_path")
    outdir = state.get("outdir")
    # Default: keep current running flag
    running = bool(state.get("running", False))
    text = ""
    try:
        if log_path and os.path.exists(log_path):
            # Read last ~4000 chars (approx 100 lines)
            with open(log_path, "rb") as f:
                try:
                    f.seek(-4000, os.SEEK_END)
                except Exception:
                    f.seek(0)
                text = f.read().decode(errors="replace")
        # Stop when summary exists or completion message appears
        done = False
        if outdir and os.path.exists(os.path.join(outdir, "dge_summary.json")):
            done = True
        if "DGE analysis complete" in text:
            done = True
        if done:
            running = False
    except Exception as _e:
        text = f"[log unavailable] {_e}"
        running = False
    # Update state while preserving keys
    new_state = dict(state)
    new_state["running"] = running
    return text, new_state

# Populate Enrichment figures
@app.callback(
    [Output("dge-go-figs", "children"),
     Output("dge-gsea-figs", "children")],
    [Input("dge-comparison-dropdown", "value"),
     Input("output-run-name", "value")],
    State("base-outdir", "value")
)
def load_enrichment_figs(comparison, run_name, base_outdir):
    if not comparison or not run_name:
        return "", ""
    current_outdir = get_run_outdir(run_name, base_outdir)
    enrich_dir = os.path.join(current_outdir, "dge_analysis", comparison, "enrichment")
    if not os.path.isdir(enrich_dir):
        return html.Div("No GO enrichment results found yet."), html.Div("No GSEA results found yet.")

    def build_img(path_rel):
        params = urllib.parse.urlencode({"run_name": run_name, "path": path_rel, "base_outdir": base_outdir})
        return html.Img(src=f"/reports/serve?{params}", style={"maxWidth": "100%", "maxHeight": "650px", "objectFit": "contain", "marginBottom": "10px", "border": "1px solid #333"})

    # GO images
    go_imgs = []
    for patt in ("GO_*_UP_dotplot.png", "GO_*_DOWN_dotplot.png"):
        for fp in sorted(glob.glob(os.path.join(enrich_dir, patt))):
            rel = os.path.relpath(fp, current_outdir)
            go_imgs.append(build_img(rel))
    if not go_imgs:
        go_div = html.Div("No GO enrichment plots found.")
    else:
        go_div = html.Div(go_imgs)

    # GSEA images (dotplots + top curves if present)
    gsea_imgs = []
    for patt in ("GSEA_GO_*_dotplot.png", "GSEA_GO_*_curve_*.png"):
        for fp in sorted(glob.glob(os.path.join(enrich_dir, patt))):
            rel = os.path.relpath(fp, current_outdir)
            gsea_imgs.append(build_img(rel))
    if not gsea_imgs:
        gsea_div = html.Div("No GSEA plots found.")
    else:
        gsea_div = html.Div(gsea_imgs)

    return go_div, gsea_div

# Load Rmd into editor
@app.callback(
    Output("dge-markdown-editor", "value"),
    Input("dge-comparison-dropdown", "value")
)
def load_markdown_editor(comparison):
    rmd_path = os.path.join(os.path.dirname(__file__), "bin", "dge_report.Rmd")
    if os.path.exists(rmd_path):
        with open(rmd_path, 'r') as f:
            return f.read()
    return "# R Markdown Analysis\n\nSelect a comparison or edit this code."

# Run Markdown Analysis
@app.callback(
    [Output("markdown-run-status", "children"),
     Output("dge-markdown-iframe", "src", allow_duplicate=True)],
    Input("btn-run-markdown", "n_clicks"),
    [State("dge-markdown-editor", "value"),
     State("dge-comparison-dropdown", "value"),
     State("output-run-name", "value"),
     State("dge-tool", "value"),
     State("dge-group-col", "value"),
     State("dge-control", "value"),
     State("dge-treatment", "value"),
     State("dge-covariates", "value"),
     State("batch-method", "value"),
     State("batch-col", "value"),
     State("dge-p-thresh", "value"),
     State("dge-fc-thresh", "value"),
    State("base-outdir", "value")],
    prevent_initial_call=True
)
def run_markdown_analysis(n_clicks, code, comparison, run_name, tool, group_col, control,
                          treatment, covariates, batch_method, batch_col, p_thresh, fc_thresh, base_outdir):
    if not n_clicks or not comparison or not run_name:
        return no_update, no_update

    current_outdir = get_run_outdir(run_name, base_outdir)
    outdir = os.path.join(current_outdir, "dge_analysis", comparison)
    os.makedirs(outdir, exist_ok=True)

    rmd_path = os.path.join(outdir, "custom_report.Rmd")
    html_path = os.path.join(outdir, "dge_report.html")

    with open(rmd_path, "w") as f:
        f.write(code)

    counts_file = find_counts_file(current_outdir)
    if not counts_file:
        counts_file = os.path.join(current_outdir, "counts", "gene_counts_matrix.tsv")

    meta_file = os.path.join(outdir, "metadata_dge.csv")

    # Build R command for rendering
    params = f"list(counts='{counts_file}', metadata='{meta_file}', control='{control}', treatment='{treatment}', tool='{tool}', group_col='{group_col}', p_threshold={p_thresh}, fc_threshold={fc_thresh}"
    if batch_method: params += f", batch_method='{batch_method}'"
    if batch_col: params += f", batch_col='{batch_col}'"
    if covariates: params += f", covariates='{','.join(covariates)}'"
    params += ")"

    render_cmd = f"rmarkdown::render('{rmd_path}', output_file='{html_path}', params={params}, quiet=TRUE)"

    try:
        subprocess.check_call(["Rscript", "-e", render_cmd])
        rel_path = f"dge_analysis/{comparison}/dge_report.html"
        params = urllib.parse.urlencode({"run_name": run_name, "path": rel_path, "base_outdir": base_outdir})
        return dbc.Alert("Report rendered successfully!", color="success"), f"/reports/serve?{params}"
    except Exception as e:
        return dbc.Alert(f"Error rendering report: {e}", color="danger"), no_update


# --- 8. Run Server ---
if __name__ == "__main__":
    app.run(debug=False, host="0.0.0.0", port=8502)
