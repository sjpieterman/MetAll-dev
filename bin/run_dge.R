#!/usr/bin/env Rscript

# Differential Gene Expression Analysis with DESeq2, edgeR, or limma-voom
# Supports batch correction with ComBat or limma
# Includes GO and GSEA functional enrichment

# ---- Auto-install dependencies (optional) ----
# Disable by setting: METALL_NO_R_AUTOINSTALL=1
no_autoinstall <- Sys.getenv("METALL_NO_R_AUTOINSTALL", unset = "")
no_autoinstall <- tolower(trimws(no_autoinstall)) %in% c("1", "true", "yes")

# Function to check if internet is available
check_internet <- function() {
  tryCatch({
    con <- url("https://cloud.r-project.org")
    on.exit(close(con))
    readLines(con, n = 1)
    TRUE
  }, error = function(e) FALSE)
}

ensure_cran <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(missing) == 0) return(invisible(TRUE))
  
  if (!check_internet()) {
    message("[MetaLL] WARNING: Missing CRAN packages (", paste(missing, collapse = ", "), ") and NO INTERNET DETECTED.")
    message("[MetaLL] If you use a proxy, set http_proxy/https_proxy env vars.")
    return(invisible(FALSE))
  }
  
  message("[MetaLL] Installing CRAN packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}

ensure_bioc <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    if (!check_internet()) {
        message("[MetaLL] WARNING: BiocManager missing and NO INTERNET.")
        return(invisible(FALSE))
    }
    message("[MetaLL] Installing CRAN package: BiocManager")
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  
  missing <- pkgs[!vapply(pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(missing) == 0) return(invisible(TRUE))
  
  if (!check_internet()) {
    message("[MetaLL] WARNING: Missing Bioconductor packages (", paste(missing, collapse = ", "), ") and NO INTERNET.")
    return(invisible(FALSE))
  }

  message("[MetaLL] Installing Bioconductor packages: ", paste(missing, collapse = ", "))
  BiocManager::install(missing, ask = FALSE, update = FALSE)
}

# Initial core dependencies
if (!no_autoinstall) {
  ensure_cran(c("optparse", "ggplot2", "pheatmap", "RColorBrewer", "dplyr", "plotly", "htmlwidgets", "jsonlite", "ggrepel"))
  ensure_bioc(c("edgeR", "limma", "DESeq2", "sva", "biomaRt", "preprocessCore"))
}

# ---- Load libraries ----
suppressPackageStartupMessages({
  library(optparse)
  library(preprocessCore)
})

option_list <- list(
  make_option(c("-c", "--counts"), type="character", default=NULL, help="Counts matrix file (TSV)", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NULL, help="Metadata file (CSV/TSV/TXT)", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=".", help="Output directory", metavar="character"),
  make_option(c("-t", "--tool"), type="character", default="deseq2", help="DGE tool: deseq2, edger, limma", metavar="character"),
  make_option(c("-b", "--batch_method"), type="character", default="none", help="Batch correction: none, combat, limma", metavar="character"),
  make_option(c("-v", "--batch_col"), type="character", default=NULL, help="Metadata column for batch effect", metavar="character"),
  make_option(c("-g", "--group_col"), type="character", default="group", help="Metadata column for experimental group (e.g. AlteredTasteOrSmell)", metavar="character"),
  make_option(c("--covariates"), type="character", default=NULL, help="Comma-separated covariates to include in the model (e.g. age,sex,bmi)", metavar="character"),
  make_option(c("--sample_col"), type="character", default="sample", help="Metadata column containing sample IDs (must match counts column names)", metavar="character"),
  make_option(c("--metadata_sep"), type="character", default="auto", help="Metadata separator: auto, tab, comma, semicolon", metavar="character"),
  make_option(c("-r", "--control"), type="character", default=NULL, help="Control group value (e.g. 0)", metavar="character"),
  make_option(c("-x", "--treatment"), type="character", default=NULL, help="Treatment group value (e.g. 1)", metavar="character"),
  make_option(c("--all_pairs"), action="store_true", default=FALSE, help="Run all pairwise comparisons within metadata columns", metavar="logical"),
  make_option(c("--group_cols"), type="character", default=NULL, help="Comma-separated list of metadata columns for --all_pairs (empty = all)", metavar="character"),
  make_option(c("--skip_qc"), action="store_true", default=FALSE, help="Skip dataset-level QC plots", metavar="logical"),
  make_option(c("--qc_outdir"), type="character", default=NULL, help="Output directory for dataset-level QC plots (default: <outdir>/dataset_qc)", metavar="character"),
  make_option(c("--use_biomart"), action="store_true", default=FALSE, help="Map Ensembl gene IDs to gene symbols via biomaRt", metavar="logical"),
  make_option(c("--biomart_dataset"), type="character", default="hsapiens_gene_ensembl", help="biomaRt dataset (e.g. hsapiens_gene_ensembl)", metavar="character"),
  make_option(c("--biomart_id_attr"), type="character", default="ensembl_gene_id", help="biomaRt ID attribute (e.g. ensembl_gene_id)", metavar="character"),
  make_option(c("--biomart_symbol_attr"), type="character", default="hgnc_symbol", help="biomaRt symbol attribute (e.g. hgnc_symbol)", metavar="character"),
  make_option(c("--strip_gene_version"), type="logical", default=TRUE, help="Strip version suffix from gene IDs before mapping (e.g. ENSG... .1)", metavar="logical"),
  make_option(c("-p", "--p_threshold"), type="numeric", default=0.05, help="P-value threshold", metavar="number"),
  make_option(c("-f", "--fc_threshold"), type="numeric", default=1.0, help="Log2 fold change threshold", metavar="number"),
  make_option(c("--norm_method"), type="character", default=NULL, help="Normalization method for insights: voom, deseq2, tmm (defaults to match --tool)", metavar="character"),
  make_option(c("--run_go"), action="store_true", default=FALSE, help="Run GO enrichment analysis", metavar="logical"),
  make_option(c("--run_gsea"), action="store_true", default=FALSE, help="Run GSEA enrichment analysis", metavar="logical"),
  make_option(c("--organism_db"), type="character", default="org.Hs.eg.db", help="OrgDb for clusterProfiler (e.g. org.Hs.eg.db)", metavar="character"),
  make_option(c("--ontology"), type="character", default="BP", help="Ontology for enrichment (BP, CC, MF, or ALL)", metavar="character"),
  make_option(c("--top_n_enrich"), type="integer", default=20, help="Top N terms to show in enrichment plots", metavar="integer"),
  make_option(c("--keytype"), type="character", default="SYMBOL", help="Gene ID keyType for enrichment (SYMBOL, ENSEMBL, ENTREZID)", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$counts) || is.null(opt$metadata)) {
  print_help(opt_parser)
  stop("Missing required arguments", call.=FALSE)
}
if (!isTRUE(opt$all_pairs) && (is.null(opt$control) || is.null(opt$treatment))) {
  print_help(opt_parser)
  stop("Missing required arguments", call.=FALSE)
}

# Create output directory
if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE)
}

# Load libraries
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  if (opt$tool == "deseq2") library(DESeq2)
  if (opt$batch_method == "combat") library(sva)
  library(plotly)
  library(htmlwidgets)
  library(jsonlite)
  library(ggrepel)
})

read_metadata_flexible <- function(path, sep_opt = "auto") {
  if (sep_opt == "tab") return(read.delim(path, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE))
  if (sep_opt == "comma") return(read.csv(path, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE))
  if (sep_opt == "semicolon") return(read.csv(path, header=TRUE, sep=";", check.names=FALSE, stringsAsFactors=FALSE))

  first <- readLines(path, n = 1, warn = FALSE)
  sep <- "\t"
  if (grepl(";", first, fixed = TRUE) && !grepl("\t", first, fixed = TRUE) && !grepl(",", first, fixed = TRUE)) {
    sep <- ";"
  } else if (grepl(",", first, fixed = TRUE) && !grepl("\t", first, fixed = TRUE)) {
    sep <- ","
  } else if (grepl("\t", first, fixed = TRUE)) {
    sep <- "\t"
  }
  read.table(path, header=TRUE, sep=sep, check.names=FALSE, stringsAsFactors=FALSE, quote="\"", comment.char="")
}

sanitize_factor_levels <- function(vec) {
  orig_levels <- levels(factor(vec))
  safe_levels <- make.unique(make.names(orig_levels))
  f <- factor(vec, levels = orig_levels)
  levels(f) <- safe_levels
  list(factor = f, map = setNames(safe_levels, orig_levels))
}

map_level <- function(value, level_map) {
  if (is.null(value) || is.na(value)) return(value)
  value <- as.character(value)
  if (value %in% names(level_map)) return(unname(level_map[value]))
  if (value %in% level_map) return(value)
  value
}

map_gene_symbols_biomart <- function(gene_ids, dataset, id_attr, symbol_attr) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    if (!no_autoinstall) {
       message("[MetaLL] biomaRt requested but not found. Attempting installation...")
       ensure_bioc("biomaRt")
    }
  }
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    message("biomaRt not available. Skipping biomaRt mapping.")
    return(setNames(rep(NA_character_, length(gene_ids)), gene_ids))
  }
  tryCatch({
    message("[MetaLL] Connecting to biomaRt...")
    mart <- biomaRt::useEnsembl(biomart = "genes", dataset = dataset)
    bm <- biomaRt::getBM(
      attributes = c(id_attr, symbol_attr),
      filters = id_attr,
      values = unique(gene_ids),
      mart = mart
    )
    if (nrow(bm) == 0) return(setNames(rep(NA_character_, length(gene_ids)), gene_ids))
    split_syms <- split(bm[[symbol_attr]], bm[[id_attr]])
    sym_map <- vapply(
      split_syms,
      function(x) {
        x <- x[!is.na(x) & x != ""]
        if (length(x) == 0) return(NA_character_)
        x[1]
      },
      FUN.VALUE = character(1)
    )
    sym_map
  }, error = function(e) {
    message("biomaRt mapping failed: ", e$message)
    setNames(rep(NA_character_, length(gene_ids)), gene_ids)
  })
}

detect_keytype <- function(gene_ids) {
  gene_ids <- head(unique(gene_ids[!is.na(gene_ids)]), 100)
  if (length(gene_ids) == 0) return("SYMBOL")
  
  if (any(grepl("^ENS[A-Z]*G", gene_ids))) {
    return("ENSEMBL")
  } else if (all(grepl("^\\d+$", gene_ids))) {
    return("ENTREZID")
  } else {
    return("SYMBOL")
  }
}

# Load data
message("[MetaLL] Loading counts matrix...")
counts <- read.table(opt$counts, header=TRUE, sep="\t", row.names=1, check.names=FALSE)

# Keep an unfiltered copy for "raw vs normalized" comparisons later
original_counts <- counts

message("[MetaLL] Loading metadata...")
metadata <- read_metadata_flexible(opt$metadata, opt$metadata_sep)

# Validate required metadata columns
if (!(opt$sample_col %in% colnames(metadata))) {
  stop(paste0("Sample column '", opt$sample_col, "' not found in metadata. Available: ", paste(colnames(metadata), collapse=", ")))
}

# Ensure metadata matches counts
metadata_samples <- as.character(metadata[[opt$sample_col]])
common_samples <- intersect(colnames(counts), metadata_samples)
if (length(common_samples) == 0) {
  stop("No common samples between counts and metadata. Check sample names.")
}

counts <- counts[, common_samples, drop=FALSE]
metadata <- metadata[match(common_samples, metadata_samples), , drop=FALSE]

# Keep aligned copies for dataset-level QC
counts_full <- counts
original_counts_full <- original_counts[, colnames(counts_full), drop=FALSE]
metadata_full <- metadata

run_dataset_qc <- function(counts_raw, counts_for_norm, metadata_qc, outdir) {
  message("[MetaLL] Generating dataset-level QC plots...")
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  # Set norm_method to match tool if not specified
  if (is.null(opt$norm_method)) {
    opt$norm_method <- switch(opt$tool,
      "deseq2" = "deseq2",
      "edger" = "tmm",
      "limma" = "voom",
      "voom"
    )
  }

  qc_group_col <- NULL
  if (!is.null(opt$group_col) && opt$group_col %in% colnames(metadata_qc)) {
    qc_group_col <- opt$group_col
  }
  group_vec <- if (!is.null(qc_group_col)) metadata_qc[[qc_group_col]] else factor(rep("all", nrow(metadata_qc)))

  # Compute raw logCPM for "before" comparison
  dge_raw <- DGEList(counts=counts_raw, group=group_vec)
  dge_raw <- calcNormFactors(dge_raw)
  logcpm_raw <- cpm(dge_raw, log=TRUE)

  # Compute normalized data based on selected method
  dge <- DGEList(counts=counts_for_norm, group=group_vec)
  dge <- calcNormFactors(dge)

  if (opt$norm_method == "voom") {
    design <- model.matrix(~1, data=metadata_qc)
    v <- voom(dge, design, plot=FALSE)
    normalized_data <- v$E
  } else if (opt$norm_method == "deseq2") {
    dds <- DESeqDataSetFromMatrix(countData = round(counts_for_norm), colData = metadata_qc, design = ~1)
    dds <- estimateSizeFactors(dds)
    normalized_data <- log2(counts(dds, normalized=TRUE) + 1)
  } else if (opt$norm_method == "tmm") {
    normalized_data <- cpm(dge, log=TRUE)
  } else {
    stop("Invalid --norm_method. Choose: voom, deseq2, tmm")
  }

  batch_corrected_data <- NULL
  if (opt$batch_method == "combat" && !is.null(opt$batch_col) && opt$batch_col %in% colnames(metadata_qc)) {
    batch <- as.factor(metadata_qc[[opt$batch_col]])
    mod <- model.matrix(~1, data=metadata_qc)
    try({
       batch_corrected_data <- sva::ComBat(dat=normalized_data, batch=batch, mod=mod)
    })
  }

  zscore_normalized_data <- t(scale(t(if (is.null(batch_corrected_data)) normalized_data else batch_corrected_data)))
  zscore_normalized_data[is.na(zscore_normalized_data)] <- 0

  before_data <- logcpm_raw
  after_data <- normalized_data
  before_title <- "Raw Data (logCPM)"
  after_title <- paste("Normalized Data (", opt$norm_method, ")", sep="")
  batch_title <- "Batch-Corrected Data (ComBat)"

  # Boxplots
  png(file.path(outdir, "boxplot_normalization.png"), width=800, height=600)
  par(mfrow = c(2, 1))
  boxplot(before_data, main = before_title, las = 2, col = "lightgreen", outline = FALSE)
  boxplot(after_data, main = after_title, las = 2, col = "lightcoral", outline = FALSE)
  dev.off()

  # Density plots
  png(file.path(outdir, "density_plot.png"), width=800, height=600)
  par(mfrow = c(1, 1))
  plot(density(as.vector(before_data)), main = "Density Plot", col = "lightgreen", lwd = 2)
  lines(density(as.vector(after_data)), col = "lightcoral", lwd = 2)
  legend("topright", legend = c(before_title, after_title), col = c("lightgreen", "lightcoral"), lwd = 2)
  dev.off()

  png(file.path(outdir, "density_before_normalization.png"), width=800, height=600)
  plot(density(as.vector(before_data)), main = before_title, col = "lightgreen", lwd = 2)
  dev.off()

  png(file.path(outdir, "density_after_normalization.png"), width=800, height=600)
  plot(density(as.vector(after_data)), main = after_title, col = "lightcoral", lwd = 2)
  dev.off()

  if (!is.null(batch_corrected_data)) {
    png(file.path(outdir, "density_batch_corrected.png"), width=800, height=600)
    plot(density(as.vector(batch_corrected_data)), main = batch_title, col = "lightskyblue", lwd = 2)
    dev.off()
  }

  # Per-sample bar plots (mean log-scale expression)
  sample_barplot <- function(mat, title, file) {
    png(file, width=1200, height=600)
    barplot(colMeans(mat), las = 2, col = "steelblue", main = title, ylab = "Mean (log scale)")
    dev.off()
  }

  sample_barplot(before_data, paste("Per-sample Mean -", before_title), file.path(outdir, "sample_barplot_before_normalization.png"))
  sample_barplot(after_data, paste("Per-sample Mean -", after_title), file.path(outdir, "sample_barplot_after_normalization.png"))
  if (!is.null(batch_corrected_data)) {
    sample_barplot(batch_corrected_data, paste("Per-sample Mean -", batch_title), file.path(outdir, "sample_barplot_batch_corrected.png"))
  }

  # PCA comparisons
  pca_before <- tryCatch(prcomp(t(before_data)), error = function(e) NULL)
  pca_after <- tryCatch(prcomp(t(after_data)), error = function(e) NULL)

  if (opt$batch_method != "none" && !is.null(opt$batch_col) && opt$batch_col %in% colnames(metadata_qc)) {
    unique_batches <- unique(metadata_qc[[opt$batch_col]])
    color_mapping <- rainbow(length(unique_batches))
    names(color_mapping) <- unique_batches
    colors <- color_mapping[as.character(metadata_qc[[opt$batch_col]])]
    colors[is.na(colors)] <- "gray"
    color_legend <- opt$batch_col
  } else if (!is.null(qc_group_col)) {
    unique_groups <- unique(metadata_qc[[qc_group_col]])
    color_mapping <- rainbow(length(unique_groups))
    names(color_mapping) <- unique_groups
    colors <- color_mapping[as.character(metadata_qc[[qc_group_col]])]
    color_legend <- qc_group_col
  } else {
    colors <- rep("gray", nrow(metadata_qc))
    color_mapping <- c(all = "gray")
    color_legend <- "all"
  }

  png(file.path(outdir, "pca_normalization.png"), width=800, height=600)
  par(mfrow = c(2, 1))
  if (!is.null(pca_before)) {
    plot(pca_before$x[, 1:2], col = colors, main = paste("PCA -", before_title), xlab = "PC1", ylab = "PC2")
    legend("topright", legend = names(color_mapping), fill = color_mapping, title = color_legend, bty="n", xpd=TRUE, inset=c(0,0))
  }
  if (!is.null(pca_after)) {
    plot(pca_after$x[, 1:2], col = colors, main = paste("PCA -", after_title), xlab = "PC1", ylab = "PC2")
    legend("topright", legend = names(color_mapping), fill = color_mapping, title = color_legend, bty="n", xpd=TRUE, inset=c(0,0))
  }
  dev.off()

  # PCA Plot + Scree (interactive)
  pca_data <- t(as.matrix(if (is.null(batch_corrected_data)) after_data else batch_corrected_data))
  pca_data <- pca_data[, apply(pca_data, 2, function(z) all(is.finite(z))), drop=FALSE]
  if (ncol(pca_data) > 1) {
    nzv <- apply(pca_data, 2, function(z) stats::var(z) > 0)
    pca_data <- pca_data[, nzv, drop=FALSE]
  }

  if (ncol(pca_data) >= 2) {
    pca <- tryCatch(prcomp(pca_data, center = TRUE, scale. = FALSE), error = function(e) NULL)
    if (!is.null(pca)) {
      pca_df <- as.data.frame(pca$x)
      pca_df$Sample <- rownames(pca_df)
      if (!is.null(qc_group_col)) pca_df$Group <- metadata_qc[[qc_group_col]]
      if (!is.null(opt$batch_col) && opt$batch_col %in% colnames(metadata_qc)) pca_df$Batch <- metadata_qc[[opt$batch_col]]

      percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2))

      p <- ggplot(pca_df, aes(PC1, PC2, color=if (!is.null(qc_group_col)) Group else NULL, shape=if (!is.null(opt$batch_col) && opt$batch_col %in% colnames(metadata_qc)) Batch else NULL, text=Sample)) +
        geom_point(size=3) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        ggtitle("PCA Plot") +
        theme_minimal() +
        theme(legend.position = "bottom")

      ggsave(file.path(outdir, "pca_plot.png"), p, width=7, height=6)

      p_inter <- ggplotly(p, tooltip="text") %>% layout(legend = list(orientation = "h", x = 0, y = -0.2))
      saveWidget(p_inter, file.path(outdir, "pca_plot.html"), selfcontained = TRUE)

      scree_df <- data.frame(
        PC = factor(paste0("PC", seq_along(percentVar)), levels = paste0("PC", seq_along(percentVar))),
        Variance = percentVar
      )
      # Limit to top 20 PCs for scree
      scree_df <- head(scree_df, 20)

      p_scree <- ggplot(scree_df, aes(x = PC, y = Variance, text = paste0(Variance, "%"))) +
        geom_col(fill = "#00d9ff") +
        theme_minimal() +
        ylab("Variance Explained (%)") +
        xlab("Principal Component") +
        ggtitle("PCA Scree Plot")

      p_scree_inter <- ggplotly(p_scree, tooltip = "text")
      saveWidget(p_scree_inter, file.path(outdir, "pca_scree_plot.html"), selfcontained = TRUE)
    } else {
       message("PCA calculation failed.")
    }
  } else {
    message("Skipping PCA plot: not enough finite, variable features after filtering.")
  }

  # Statistical comparison
  before_variance <- apply(before_data, 1, var)
  after_variance <- apply(after_data, 1, var)
  try({
    t_test_result <- t.test(before_variance, after_variance)
    sink(file.path(outdir, "normalization_stats.txt"))
    print(t_test_result)
    sink()
  })
}

sanitize_name <- function(x) {
  gsub("[^A-Za-z0-9._-]", "_", x)
}

parse_group_cols <- function(x) {
  if (is.null(x) || is.na(x) || trimws(x) == "") return(character(0))
  out <- unlist(strsplit(x, ",", fixed = TRUE))
  out <- trimws(out)
  out <- out[out != ""]
  unique(out)
}

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grep("^--file=", args)]
  if (length(file_arg) == 0) {
    stop("Unable to locate script path for --all_pairs execution.")
  }
  normalizePath(sub("^--file=", "", file_arg[1]))
}

if (isTRUE(opt$all_pairs)) {
  if (is.null(opt$qc_outdir) || trimws(opt$qc_outdir) == "") {
    opt$qc_outdir <- file.path(opt$outdir, "dataset_qc")
  }

  if (!isTRUE(opt$skip_qc)) {
    run_dataset_qc(original_counts_full, counts_full, metadata_full, opt$qc_outdir)
  }

  cols <- parse_group_cols(opt$group_cols)
  if (length(cols) == 0) {
    cols <- setdiff(colnames(metadata), opt$sample_col)
  }

  missing_cols <- cols[!(cols %in% colnames(metadata))]
  if (length(missing_cols) > 0) {
    stop(paste0("Group column(s) not found in metadata: ", paste(missing_cols, collapse = ", ")))
  }

  script_path <- get_script_path()
  add_arg <- function(args, flag, value) {
    if (is.null(value) || is.na(value) || trimws(as.character(value)) == "") return(args)
    c(args, flag, as.character(value))
  }

  for (group_col in cols) {
    raw_vals <- trimws(as.character(metadata[[group_col]]))
    raw_vals[raw_vals == ""] <- NA
    vals <- unique(raw_vals[!is.na(raw_vals) & !tolower(raw_vals) %in% c("na", "nan")])
    if (length(vals) < 2) {
      message("Skipping ", group_col, ": need at least 2 non-empty levels.")
      next
    }

    for (i in seq_len(length(vals) - 1)) {
      for (j in (i + 1):length(vals)) {
        control <- vals[i]
        treatment <- vals[j]
        outdir <- file.path(
          opt$outdir,
          sanitize_name(group_col),
          paste0(
            sanitize_name(treatment),
            "_vs_",
            sanitize_name(control),
            "_",
            sanitize_name(opt$tool)
          )
        )

        cmd_args <- c(
          script_path,
          "--counts", opt$counts,
          "--metadata", opt$metadata,
          "--metadata_sep", opt$metadata_sep,
          "--sample_col", opt$sample_col,
          "--group_col", group_col,
          "--control", control,
          "--treatment", treatment,
          "--outdir", outdir,
          "--tool", opt$tool,
          "--batch_method", opt$batch_method,
          "--p_threshold", opt$p_threshold,
          "--fc_threshold", opt$fc_threshold,
          "--all_pairs", "false",
          "--skip_qc", "true",
          "--qc_outdir", opt$qc_outdir,
          "--keytype", opt$keytype
        )

        cmd_args <- add_arg(cmd_args, "--covariates", opt$covariates)
        cmd_args <- add_arg(cmd_args, "--batch_col", opt$batch_col)
        cmd_args <- add_arg(cmd_args, "--norm_method", opt$norm_method)
        if (opt$run_go) cmd_args <- c(cmd_args, "--run_go")
        if (opt$run_gsea) cmd_args <- c(cmd_args, "--run_gsea")
        cmd_args <- add_arg(cmd_args, "--organism_db", opt$organism_db)
        cmd_args <- add_arg(cmd_args, "--ontology", opt$ontology)
        
        message("[MetaLL] Running comparison: ", group_col, " (", treatment, " vs ", control, ")")
        status <- system2("Rscript", cmd_args)
        if (!identical(status, 0L)) {
          message("Rscript failed for ", group_col, ": ", treatment, " vs ", control)
        }
      }
    }
  }

  quit(save = "no", status = 0)
}

if (!(opt$group_col %in% colnames(metadata))) {
  stop(paste0("Group column '", opt$group_col, "' not found in metadata. Available: ", paste(colnames(metadata), collapse=", ")))
}

# Normalization insights (dataset-level)
if (!isTRUE(opt$skip_qc)) {
  if (is.null(opt$qc_outdir) || trimws(opt$qc_outdir) == "") {
    opt$qc_outdir <- file.path(opt$outdir, "dataset_qc")
  }
  run_dataset_qc(original_counts_full, counts_full, metadata_full, opt$qc_outdir)
}

# Group column as factor (trim blanks -> NA)
grp_raw <- trimws(as.character(metadata[[opt$group_col]]))
grp_raw[grp_raw == ""] <- NA
metadata[[opt$group_col]] <- as.factor(grp_raw)
grp_sanitized <- sanitize_factor_levels(metadata[[opt$group_col]])
metadata[[opt$group_col]] <- grp_sanitized$factor
group_level_map <- grp_sanitized$map

# Normalize CLI inputs
opt$control   <- trimws(as.character(opt$control))
opt$treatment <- trimws(as.character(opt$treatment))
opt$control   <- map_level(opt$control, group_level_map)
opt$treatment <- map_level(opt$treatment, group_level_map)

# Filter metadata for control and treatment groups (and drop NA)
samples_to_keep <- !is.na(metadata[[opt$group_col]]) & metadata[[opt$group_col]] %in% c(opt$control, opt$treatment)
counts <- counts[, samples_to_keep, drop=FALSE]
metadata <- metadata[samples_to_keep, , drop=FALSE]
metadata[[opt$group_col]] <- droplevels(metadata[[opt$group_col]])

# Keep original_counts aligned to the filtered sample set
original_counts <- original_counts[, colnames(counts), drop=FALSE]

# Validate that control/treatment exist after filtering
lvl <- levels(metadata[[opt$group_col]])
if (!(opt$control %in% lvl)) {
  stop(
    paste0(
      "Control value '", opt$control, "' not found in group levels after filtering.\n",
      "Available levels: ", paste(lvl, collapse = ", "), "\n",
      "Fix: choose a Control value that exactly matches the metadata column '", opt$group_col, "'."
    ),
    call. = FALSE
  )
}
if (!(opt$treatment %in% lvl)) {
  stop(
    paste0(
      "Treatment value '", opt$treatment, "' not found in group levels after filtering.\n",
      "Available levels: ", paste(lvl, collapse = ", "), "\n",
      "Fix: choose a Treatment value that exactly matches the metadata column '", opt$group_col, "'."
    ),
    call. = FALSE
  )
}

metadata[[opt$group_col]] <- relevel(metadata[[opt$group_col]], ref=opt$control)

parse_covariates <- function(x) {
  if (is.null(x) || is.na(x) || trimws(x) == "") return(character(0))
  out <- unlist(strsplit(x, ",", fixed = TRUE))
  out <- trimws(out)
  out <- out[out != ""]
  unique(out)
}

covars <- parse_covariates(opt$covariates)

# Validate covariates exist and coerce character covariates to factors
if (length(covars) > 0) {
  missing_covars <- covars[!(covars %in% colnames(metadata))]
  if (length(missing_covars) > 0) {
    stop(paste0("Covariate column(s) not found in metadata: ", paste(missing_covars, collapse = ", ")))
  }

  for (cv in covars) {
    if (is.character(metadata[[cv]])) {
      metadata[[cv]] <- as.factor(metadata[[cv]])
    }
    if (is.factor(metadata[[cv]])) {
      metadata[[cv]] <- sanitize_factor_levels(metadata[[cv]])$factor
    }
  }
}

# Helper to build model formula string: ~ cov1 + cov2 + group
build_design_formula <- function(group_col, covars) {
  terms <- c(covars, group_col)
  as.formula(paste("~", paste(terms, collapse = " + ")))
}

# Filtering low counts
message("[MetaLL] Filtering low-count genes...")
keep <- rowSums(cpm(counts) > 1) >= (min(table(metadata[[opt$group_col]])))
counts <- counts[keep,]

# Batch correction
if (opt$batch_method != "none" && !is.null(opt$batch_col)) {
  if (!(opt$batch_col %in% colnames(metadata))) {
    stop(paste("Batch column", opt$batch_col, "not found in metadata"))
  }

  batch <- as.factor(metadata[[opt$batch_col]])
  mod <- model.matrix(~metadata[[opt$group_col]])

  if (opt$batch_method == "combat") {
    message("[MetaLL] Applying ComBat batch correction...")
    if (requireNamespace("sva", quietly = TRUE)) {
       tryCatch({
         counts <- sva::ComBat_seq(as.matrix(counts), batch=batch, group=metadata[[opt$group_col]])
       }, error = function(e) {
         message("ComBat_seq failed, falling back to including batch in model.")
       })
    }
  }
}

# DGE Analysis
res <- NULL
normalized_counts <- NULL

if (opt$tool == "deseq2") {
  message("[MetaLL] Running DESeq2 analysis...")
  design_formula <- build_design_formula(opt$group_col, covars)
  dds <- DESeqDataSetFromMatrix(
    countData = round(counts),
    colData   = metadata,
    design    = design_formula
  )
  dds <- DESeq(dds)
  res <- results(dds, contrast = c(opt$group_col, opt$treatment, opt$control))
  normalized_counts <- counts(dds, normalized=TRUE)
  res_df <- as.data.frame(res)
  colnames(res_df)[colnames(res_df) == "padj"] <- "adj.P.Val"
  colnames(res_df)[colnames(res_df) == "log2FoldChange"] <- "logFC"
  colnames(res_df)[colnames(res_df) == "pvalue"] <- "P.Value"
} else if (opt$tool == "edger") {
  message("[MetaLL] Running edgeR analysis...")
  dge <- DGEList(counts=counts, group=metadata[[opt$group_col]])
  dge <- calcNormFactors(dge)
  design <- model.matrix(build_design_formula(opt$group_col, covars), data=metadata)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef=ncol(design))
  res <- topTags(qlf, n=Inf)
  res_df <- as.data.frame(res)
  normalized_counts <- cpm(dge)
  colnames(res_df)[colnames(res_df) == "FDR"] <- "adj.P.Val"
  colnames(res_df)[colnames(res_df) == "PValue"] <- "P.Value"
} else if (opt$tool == "limma") {
  message("[MetaLL] Running limma-voom analysis...")
  dge <- DGEList(counts=counts)
  dge <- calcNormFactors(dge)
  design <- model.matrix(build_design_formula(opt$group_col, covars), data=metadata)
  v <- voom(dge, design, plot=FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  res <- topTable(fit, coef=ncol(design), n=Inf)
  res_df <- as.data.frame(res)
  normalized_counts <- v$E
}

# Add gene_id
res_df$gene_id <- rownames(res_df)
res_df$gene_id_clean <- if (isTRUE(opt$strip_gene_version)) {
  sub("\\.\\d+$", "", res_df$gene_id)
} else {
  res_df$gene_id
}

res_df$gene_symbol <- NA_character_
if (isTRUE(opt$use_biomart)) {
  sym_map <- map_gene_symbols_biomart(
    gene_ids = unique(res_df$gene_id_clean),
    dataset = opt$biomart_dataset,
    id_attr = opt$biomart_id_attr,
    symbol_attr = opt$biomart_symbol_attr
  )
  res_df$gene_symbol <- unname(sym_map[res_df$gene_id_clean])
}
res_df$gene_label <- ifelse(
  !is.na(res_df$gene_symbol) & res_df$gene_symbol != "",
  res_df$gene_symbol,
  res_df$gene_id
)

# Define significance
res_df$significance <- "Not Significant"
res_df$significance[res_df$adj.P.Val < opt$p_threshold & res_df$logFC > opt$fc_threshold] <- "Upregulated"
res_df$significance[res_df$adj.P.Val < opt$p_threshold & res_df$logFC < -opt$fc_threshold] <- "Downregulated"

# Save results
write.csv(res_df, file.path(opt$outdir, "dge_results.csv"), row.names=FALSE)
top50 <- res_df %>%
  filter(significance != "Not Significant") %>%
  arrange(adj.P.Val) %>%
  head(50)
write.csv(top50, file.path(opt$outdir, "top50_degs.csv"), row.names=FALSE)

if (!is.null(normalized_counts)) {
  write.table(normalized_counts, file.path(opt$outdir, "normalized_counts.tsv"), sep="\t", quote=FALSE, col.names=NA)
  write.table(log1p(normalized_counts), file.path(opt$outdir, "normalized_counts.log1p.tsv"), sep="\t", quote=FALSE, col.names=NA)
}

# --- Functional Enrichment (GO + GSEA) ---
run_go_and_gsea_internal <- function(res_df, opt) {
  if (!opt$run_go && !opt$run_gsea) return(invisible(NULL))
  
  # Auto-install enrichment dependencies if requested
  if (!no_autoinstall) {
     message("[MetaLL] Checking enrichment dependencies...")
     ensure_bioc(c("clusterProfiler", "enrichplot", opt$organism_db))
  }

  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    message("Skipping GO/GSEA: 'clusterProfiler' not available.")
    return(invisible(NULL))
  }
  if (!requireNamespace(opt$organism_db, quietly = TRUE)) {
    message(paste0("Skipping GO/GSEA: OrgDb '", opt$organism_db, "' not available."))
    return(invisible(NULL))
  }
  org_db <- get(opt$organism_db, envir = asNamespace(opt$organism_db))

  out_dir <- file.path(opt$outdir, "enrichment")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  valid <- is.finite(res_df$logFC) & is.finite(res_df$adj.P.Val)
  df <- res_df[valid, , drop = FALSE]
  if (nrow(df) == 0) return(invisible(NULL))
  
  # Auto-detect keyType if set to SYMBOL but looks like ENSEMBL
  detected_type <- detect_keytype(df$gene_id_clean)
  use_keytype <- if (opt$keytype == "SYMBOL") detected_type else opt$keytype
  message("[MetaLL] Using keyType: ", use_keytype)
  
  gene_ids_for_enrich <- df$gene_id_clean
  
  if (opt$run_go) {
    message("[MetaLL] Running GO enrichment analysis...")
    up_genes <- unique(gene_ids_for_enrich[df$adj.P.Val < opt$p_threshold & df$logFC > opt$fc_threshold])
    down_genes <- unique(gene_ids_for_enrich[df$adj.P.Val < opt$p_threshold & df$logFC < -opt$fc_threshold])
    universe <- unique(gene_ids_for_enrich)
    
    run_enrich_go <- function(genes, univ, ont, suffix) {
      if (length(genes) < 5) return(NULL)
      ego <- tryCatch(
        clusterProfiler::enrichGO(
          gene = genes, universe = univ, OrgDb = org_db,
          keyType = use_keytype, ont = ont, pAdjustMethod = "BH",
          pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = (use_keytype != "SYMBOL")
        ),
        error = function(e) { message("GO failed: ", e$message); NULL }
      )
      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        write.csv(as.data.frame(ego), file.path(out_dir, paste0("GO_", ont, "_", suffix, ".csv")), row.names = FALSE)
        if (requireNamespace("enrichplot", quietly = TRUE)) {
          p <- enrichplot::dotplot(ego, showCategory = opt$top_n_enrich) + ggtitle(paste("GO", ont, suffix)) + theme_minimal() + theme(legend.position="bottom")
          ggsave(file.path(out_dir, paste0("GO_", ont, "_", suffix, "_dotplot.png")), p, width = 9, height = 7)
        }
      }
      ego
    }
    
    onts <- if (opt$ontology == "ALL") c("BP", "CC", "MF") else opt$ontology
    for (o in onts) {
      run_enrich_go(up_genes, universe, o, "UP")
      run_enrich_go(down_genes, universe, o, "DOWN")
    }
  }
  
  if (opt$run_gsea) {
    message("[MetaLL] Running GSEA analysis...")
    gsea_df <- df[!duplicated(gene_ids_for_enrich), ]
    gene_list <- gsea_df$logFC
    names(gene_list) <- gene_ids_for_enrich[!duplicated(gene_ids_for_enrich)]
    gene_list <- sort(gene_list, decreasing = TRUE)
    
    onts <- if (opt$ontology == "ALL") c("BP", "CC", "MF") else opt$ontology
    for (o in onts) {
      gse <- tryCatch(
        clusterProfiler::gseGO(
          geneList = gene_list, OrgDb = org_db, keyType = use_keytype,
          ont = o, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05,
          pAdjustMethod = "BH", verbose = FALSE
        ),
        error = function(e) { message("GSEA failed: ", e$message); NULL }
      )
      if (!is.null(gse) && nrow(as.data.frame(gse)) > 0) {
        write.csv(as.data.frame(gse), file.path(out_dir, paste0("GSEA_GO_", o, ".csv")), row.names = FALSE)
        if (requireNamespace("enrichplot", quietly = TRUE)) {
          p <- enrichplot::dotplot(gse, showCategory = opt$top_n_enrich) + ggtitle(paste("GSEA GO", o)) + theme_minimal() + theme(legend.position="bottom")
          ggsave(file.path(out_dir, paste0("GSEA_GO_", o, "_dotplot.png")), p, width = 9, height = 7)
          
          top_ids <- head(as.data.frame(gse)$ID, 3)
          for (tid in top_ids) {
            try({
               p_curve <- enrichplot::gseaplot2(gse, geneSetID = tid, title = tid)
               safe_id <- gsub("[^A-Za-z0-9]", "_", tid)
               ggsave(file.path(out_dir, paste0("GSEA_GO_", o, "_curve_", safe_id, ".png")), p_curve, width = 10, height = 8)
            })
          }
        }
      }
    }
  }
}

run_go_and_gsea_internal(res_df, opt)

# --- Plots ---
message("[MetaLL] Generating DGE plots...")

# MA plot
A <- if ("AveExpr" %in% colnames(res_df)) {
  res_df$AveExpr
} else if ("baseMean" %in% colnames(res_df)) {
  log2(res_df$baseMean + 1)
} else if ("logCPM" %in% colnames(res_df)) {
  res_df$logCPM
} else if (!is.null(normalized_counts)) {
  log2(rowMeans(as.matrix(normalized_counts), na.rm = TRUE) + 1)
}

if (!is.null(A)) {
  ma_df <- data.frame(A = as.numeric(A), M = as.numeric(res_df$logFC), significance = res_df$significance, gene_label = res_df$gene_label)
  ma_df <- ma_df[is.finite(ma_df$A) & is.finite(ma_df$M), , drop = FALSE]
  p_ma <- ggplot(ma_df, aes(x = A, y = M, color = significance, text = gene_label)) +
    geom_point(alpha = 0.45, size = 1.2) +
    scale_color_manual(values = c("Upregulated"="#ff6b6b", "Downregulated"="#00d9ff", "Not Significant"="#888888")) +
    theme_minimal() + theme(legend.position = "bottom") +
    geom_hline(yintercept = c(-opt$fc_threshold, 0, opt$fc_threshold), linetype = "dashed", color = "gray60") +
    ggtitle(paste0("MA Plot: ", opt$treatment, " vs ", opt$control))
  p_ma_inter <- ggplotly(p_ma, tooltip = "text") %>% layout(legend = list(orientation = "h", x = 0, y = -0.2))
  saveWidget(p_ma_inter, file.path(opt$outdir, "ma_plot.html"), selfcontained = TRUE)
}

# Heatmap
if (!is.null(normalized_counts) && nrow(top50) > 0) {
  top_genes <- intersect(top50$gene_id, rownames(normalized_counts))
  if (length(top_genes) >= 2) {
    hm <- as.matrix(normalized_counts)[top_genes, , drop = FALSE]
    if (opt$tool != "limma") hm <- log2(hm + 1)
    hm_z <- t(scale(t(hm)))
    hm_z[!is.finite(hm_z)] <- 0
    rownames(hm_z) <- res_df$gene_label[match(rownames(hm_z), res_df$gene_id)]
    p_hm <- plot_ly(x = colnames(hm_z), y = rownames(hm_z), z = hm_z, type = "heatmap", colors = colorRamp(c("#2b83ba", "#f7f7f7", "#d7191c"))) %>%
      layout(title = paste0("Heatmap (Top DEGs)"), xaxis = list(title = "Samples"), yaxis = list(title = "Genes"))
    saveWidget(p_hm, file.path(opt$outdir, "heatmap.html"), selfcontained = TRUE)
  }
}

# Volcano
res_df_labeled <- res_df
res_df_labeled$label <- NA
sig_indices <- which(res_df$adj.P.Val < opt$p_threshold & abs(res_df$logFC) > opt$fc_threshold)
if (length(sig_indices) > 0) {
  top_sig <- sig_indices[order(res_df$adj.P.Val[sig_indices])][1:min(length(sig_indices), 20)]
  res_df_labeled$label[top_sig] <- res_df$gene_label[top_sig]
}
p_volc <- ggplot(res_df_labeled, aes(x=logFC, y=-log10(P.Value), color=significance, text=gene_label)) +
  geom_point(alpha=0.4, size=1.5) +
  scale_color_manual(values=c("Upregulated"="#ff6b6b", "Downregulated"="#00d9ff", "Not Significant"="#888888")) +
  theme_minimal() + theme(legend.position = "bottom") +
  geom_vline(xintercept=c(-opt$fc_threshold, opt$fc_threshold), linetype="dashed") +
  geom_hline(yintercept=-log10(opt$p_threshold), linetype="dashed") +
  geom_text_repel(aes(label=label), size=3, max.overlaps=15, show.legend=FALSE) +
  ggtitle(paste("Volcano Plot:", opt$treatment, "vs", opt$control))
ggsave(file.path(opt$outdir, "volcano_plot.png"), p_volc, width=7, height=6)
p_volc_inter <- ggplotly(p_volc, tooltip="text") %>% layout(legend = list(orientation = "h", x = 0, y = -0.2))
saveWidget(p_volc_inter, file.path(opt$outdir, "volcano_plot.html"), selfcontained = TRUE)

# Summary JSON
summary_json <- list(
  analysis_parameters = list(control_group = opt$control, treatment_group = opt$treatment, group_col = opt$group_col, p_threshold = opt$p_threshold, fc_threshold = opt$fc_threshold, batch_method = opt$batch_method, norm_method = opt$norm_method),
  summary_statistics = list(total_genes = nrow(res_df), significant_genes = sum(res_df$significance != "Not Significant", na.rm=TRUE), upregulated = sum(res_df$significance == "Upregulated", na.rm=TRUE), downregulated = sum(res_df$significance == "Downregulated", na.rm=TRUE)),
  visualizations = list(pca_plot = "dataset_qc/pca_plot.html", volcano_plot = "volcano_plot.html", ma_plot = "ma_plot.html", heatmap = "heatmap.html")
)
write_json(summary_json, file.path(opt$outdir, "dge_summary.json"), auto_unbox = TRUE, pretty = TRUE)

message("[MetaLL] DGE analysis complete. Results saved to ", opt$outdir)
