# ---------------------------------------------------------
# RNA-seq Workflow: P29534 Analysis (DESeq2 version)
# Mirrors test_2.R workflow but uses DESeq2 instead of voom/limma
# ---------------------------------------------------------

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(biomaRt)
library(ggrepel)

# -----------------------------
# 1. Configuration & Paths
# -----------------------------
config <- list(
  annotation_path = "/Users/stephen.pieterman/Documents/PhD/P29453/DownStreamAnalysis/Clean_20251116/annotationFile_P29453_B22GHKVLT3_basic.csv",
  count_data_path = "/Users/stephen.pieterman/Documents/PhD/P29453/DownStreamAnalysis/Clean_20251116/salmon.merged.gene_counts_270924.tsv",
  min_count = 10,
  min_prop = 0.5,
  ensembl_host = "https://useast.ensembl.org",
  go_gsea_output_dir = "/Users/stephen.pieterman/handyTools/RNAlysis/P29453/enrichment_results",
  enrichment_fdr = 0.05,
  gsea_min_gs = 10,
  gsea_max_gs = 500,
  enrichment_plot_top_n = 20,
  gsea_plot_top_n = 3
)

# -----------------------------
# 2. Data Loading & Preprocessing
# -----------------------------
load_and_clean_data <- function(cfg) {
  message("Loading data...")

  ann <- read.csv(cfg$annotation_path, header = TRUE, sep = ";")
  counts <- read.csv(cfg$count_data_path, header = TRUE, row.names = 1, sep = "\t")
  print(ann)

  ann <- ann[ann$lane %in% c("L001", "L002") & ann$sampleType != "NEG", ]
  ann <- ann[!is.na(ann$AlteredTasteOrSmell), ]

  ann$SensoryStatus <- ifelse(ann$Hyposmia == "Hyposmia/anosmia", "SensoryLoss",
                              ifelse(ann$Hyposmia == "Normosmia", "Normal",
                                     ifelse(ann$AlteredTasteOrSmell == 1, "SensoryLoss", "Normal")))

  ann$SensoryStatus <- factor(ann$SensoryStatus, levels = c("Normal", "SensoryLoss"))
  ann$timepoint <- factor(ann$timepoint, levels = c("Inclusion", "3_weeks", "2_months"))
  ann$patientID <- factor(ann$patientID)
  ann$lane <- factor(ann$lane)

  # Adjust for sex and bmi (ensure these exist in your CSV)
  if ("sex" %in% colnames(ann)) ann$sex <- factor(ann$sex)
  if ("bmi" %in% colnames(ann)) ann$bmi <- as.numeric(gsub(",", ".", as.character(ann$bmi)))

  common_samples <- intersect(ann$sample, colnames(counts))
  counts <- counts[, common_samples, drop = FALSE]
  ann <- ann[match(common_samples, ann$sample), ]
  rownames(ann) <- ann$sample

  message(sprintf("Loaded %d samples and %d genes.", ncol(counts), nrow(counts)))
  return(list(counts = counts, ann = ann))
}

# -----------------------------
# 3. Design Selection
# -----------------------------
choose_full_rank_formula <- function(ann) {
  # Standardize covariates
  covs <- intersect(c("sex", "bmi"), colnames(ann))
  cov_str <- if (length(covs) > 0) paste(covs, collapse = " + ") else NULL

  candidates <- list()
  if (!is.null(cov_str)) {
    candidates <- list(
      as.formula(paste("~", cov_str, "+ patientID + lane + SensoryStatus * timepoint")),
      as.formula(paste("~", cov_str, "+ patientID + SensoryStatus * timepoint")),
      as.formula(paste("~", cov_str, "+ lane + SensoryStatus * timepoint")),
      as.formula(paste("~", cov_str, "+ SensoryStatus * timepoint"))
    )
  }

  candidates <- c(candidates, list(
    ~ patientID + lane + SensoryStatus * timepoint,
    ~ patientID + SensoryStatus * timepoint,
    ~ lane + SensoryStatus * timepoint,
    ~ SensoryStatus * timepoint
  ))

  for (f in candidates) {
    mm <- tryCatch(model.matrix(f, data = ann), error = function(e) NULL)
    if (!is.null(mm) && qr(mm)$rank == ncol(mm)) {
      return(f)
    }
  }

  stop("No candidate design formula produced a full-rank model matrix.")
}

# -----------------------------
# 4. DESeq2 Modeling
# -----------------------------
run_de_analysis <- function(counts, ann, cfg) {
  message("Running DESeq2 analysis...")

  design_formula <- choose_full_rank_formula(ann)
  message("Using design: ", paste(deparse(design_formula), collapse = " "))

  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(counts)),
    colData = ann,
    design = design_formula
  )

  keep <- rowSums(counts(dds) >= cfg$min_count) >= ceiling(ncol(dds) * cfg$min_prop)
  dds <- dds[keep, ]

  dds <- DESeq(dds)

  rn <- resultsNames(dds)
  message("Model coefficients: ", paste(rn, collapse = ", "))

  main_coef <- rn[grepl("^SensoryStatus_.*_vs_.*$", rn)][1]
  if (is.na(main_coef)) {
    main_coef <- rn[grepl("^SensoryStatus.*$", rn) & !grepl("timepoint", rn)][1]
  }
  int_t1 <- rn[grepl("SensoryStatus.*timepoint.*3_weeks", rn)][1]
  int_t2 <- rn[grepl("SensoryStatus.*timepoint.*2_months", rn)][1]

  if (is.na(main_coef)) {
    stop("Could not locate the SensoryStatus main-effect coefficient in the DESeq2 model.")
  }

  res_t0 <- results(dds, name = main_coef)
  res_t1 <- if (!is.na(int_t1)) {
    results(dds, contrast = list(c(main_coef, int_t1)))
  } else {
    warning("No SensoryStatus:timepoint(3_weeks) interaction coefficient found; using SensoryStatus main effect for T1.")
    results(dds, name = main_coef)
  }
  res_t2 <- if (!is.na(int_t2)) {
    results(dds, contrast = list(c(main_coef, int_t2)))
  } else {
    warning("No SensoryStatus:timepoint(2_months) interaction coefficient found; using SensoryStatus main effect for T2.")
    results(dds, name = main_coef)
  }

  to_df <- function(res) {
    df <- as.data.frame(res)
    df$P.Value <- df$pvalue
    df$adj.P.Val <- df$padj
    df$logFC <- df$log2FoldChange
    df$logP <- -log10(df$P.Value)
    df
  }

  results_list <- list(
    SensoryLoss_vs_Normal_T0 = to_df(res_t0),
    SensoryLoss_vs_Normal_T1 = to_df(res_t1),
    SensoryLoss_vs_Normal_T2 = to_df(res_t2)
  )

  # Use log2-normalized counts for visualization to avoid VST errors with small gene sets
  ntd <- normTransform(dds)
  vst_mat <- assay(ntd)

  return(list(results = results_list, vst_mat = vst_mat))
}

# -----------------------------
# 5. Gene Annotation Mapping
# -----------------------------
map_gene_symbols <- function(gene_ids, cfg) {
  message("Mapping Ensembl IDs to HGNC symbols...")

  tryCatch({
    mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = cfg$ensembl_host)
    mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                     filters = "ensembl_gene_id", values = gene_ids, mart = mart)
    mapping <- mapping[mapping$hgnc_symbol != "", ]

    symbol_map <- setNames(mapping$hgnc_symbol, mapping$ensembl_gene_id)
    mapped_names <- ifelse(gene_ids %in% names(symbol_map), symbol_map[gene_ids], gene_ids)
    make.unique(mapped_names)

  }, error = function(e) {
    warning("biomaRt mapping failed. Returning original IDs. Error: ", e$message)
    return(gene_ids)
  })
}

# -----------------------------
# 6. Functional Enrichment (GO + GSEA)
# -----------------------------
run_go_and_gsea <- function(results, method_name, cfg) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    warning("Skipping GO/GSEA: install 'clusterProfiler' and 'org.Hs.eg.db'.")
    return(invisible(NULL))
  }
  has_enrichplot <- requireNamespace("enrichplot", quietly = TRUE)
  if (!has_enrichplot) {
    warning("Package 'enrichplot' not installed; writing GO/GSEA tables only (no enrichment plots).")
  }

  out_dir <- file.path(cfg$go_gsea_output_dir, method_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  plot_dir <- file.path(out_dir, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  run_go <- function(genes, universe, ontology = "CC") {
    if (length(genes) < 5) return(NULL)
    tryCatch(
      clusterProfiler::enrichGO(
        gene = genes,
        universe = universe,
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = ontology,
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        readable = TRUE
      ),
      error = function(e) {
        warning("GO failed: ", e$message)
        NULL
      }
    )
  }

  run_gsea <- function(df, ontology = "BP") {
    keep <- is.finite(df$logFC) & !is.na(rownames(df)) & nzchar(rownames(df))
    df <- df[keep, , drop = FALSE]
    if (nrow(df) < 20) return(NULL)

    ord <- order(-abs(df$logFC))
    df <- df[ord, , drop = FALSE]
    df <- df[!duplicated(rownames(df)), , drop = FALSE]

    gene_list <- df$logFC
    names(gene_list) <- rownames(df)
    gene_list <- sort(gene_list, decreasing = TRUE)

    tryCatch(
      clusterProfiler::gseGO(
        geneList = gene_list,
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = ontology,
        minGSSize = cfg$gsea_min_gs,
        maxGSSize = cfg$gsea_max_gs,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        verbose = FALSE
      ),
      error = function(e) {
        warning("GSEA failed: ", e$message)
        NULL
      }
    )
  }

  save_enrichment_plot <- function(plot_obj, out_file, width = 9, height = 6) {
    if (is.null(plot_obj)) return(invisible(NULL))
    tryCatch(
      ggplot2::ggsave(filename = out_file, plot = plot_obj, width = width, height = height, dpi = 300),
      error = function(e) warning("Plot save failed for ", basename(out_file), ": ", e$message)
    )
    invisible(NULL)
  }

  enrichment <- list()
  for (contrast in names(results)) {
    df <- results[[contrast]]
    valid <- is.finite(df$logFC) & is.finite(df$adj.P.Val) & !is.na(rownames(df)) & nzchar(rownames(df))
    df <- df[valid, , drop = FALSE]
    if (nrow(df) == 0) next

    universe <- unique(rownames(df))
    up_genes <- unique(rownames(df)[df$adj.P.Val < cfg$enrichment_fdr & df$logFC > 0])
    down_genes <- unique(rownames(df)[df$adj.P.Val < cfg$enrichment_fdr & df$logFC < 0])

    go_up <- run_go(up_genes, universe, ontology = "BP")
    go_down <- run_go(down_genes, universe, ontology = "BP")
    gsea <- run_gsea(df, ontology = "BP")

    if (!is.null(go_up)) {
      go_up_df <- as.data.frame(go_up)
      utils::write.csv(go_up_df, file.path(out_dir, paste0(contrast, "_GO_BP_UP.csv")), row.names = FALSE)
      if (has_enrichplot && nrow(go_up_df) > 0) {
        p_go_up <- enrichplot::dotplot(go_up, showCategory = min(cfg$enrichment_plot_top_n, nrow(go_up_df))) +
          ggplot2::ggtitle(paste(contrast, "GO BP UP"))
        save_enrichment_plot(p_go_up, file.path(plot_dir, paste0(contrast, "_GO_BP_UP_dotplot.png")))
      }
    }
    if (!is.null(go_down)) {
      go_down_df <- as.data.frame(go_down)
      utils::write.csv(go_down_df, file.path(out_dir, paste0(contrast, "_GO_BP_DOWN.csv")), row.names = FALSE)
      if (has_enrichplot && nrow(go_down_df) > 0) {
        p_go_down <- enrichplot::dotplot(go_down, showCategory = min(cfg$enrichment_plot_top_n, nrow(go_down_df))) +
          ggplot2::ggtitle(paste(contrast, "GO BP DOWN"))
        save_enrichment_plot(p_go_down, file.path(plot_dir, paste0(contrast, "_GO_BP_DOWN_dotplot.png")))
      }
    }
    if (!is.null(gsea)) {
      gsea_df <- as.data.frame(gsea)
      utils::write.csv(gsea_df, file.path(out_dir, paste0(contrast, "_GSEA_GO_BP.csv")), row.names = FALSE)
      if (has_enrichplot && nrow(gsea_df) > 0) {
        p_gsea <- enrichplot::dotplot(gsea, showCategory = min(cfg$enrichment_plot_top_n, nrow(gsea_df))) +
          ggplot2::ggtitle(paste(contrast, "GSEA GO BP"))
        save_enrichment_plot(p_gsea, file.path(plot_dir, paste0(contrast, "_GSEA_GO_BP_dotplot.png")))

        gsea_ranked <- gsea_df[order(gsea_df$p.adjust, na.last = NA), , drop = FALSE]
        top_ids <- utils::head(gsea_ranked$ID, cfg$gsea_plot_top_n)
        for (term_id in top_ids) {
          p_curve <- tryCatch(
            enrichplot::gseaplot2(gsea, geneSetID = term_id, title = paste(contrast, term_id)),
            error = function(e) {
              warning("GSEA curve plot failed for ", term_id, ": ", e$message)
              NULL
            }
          )
          safe_id <- gsub("[^A-Za-z0-9_\\-]", "_", term_id)
          save_enrichment_plot(p_curve, file.path(plot_dir, paste0(contrast, "_GSEA_curve_", safe_id, ".png")), width = 10, height = 7)
        }
      }
    }

    enrichment[[contrast]] <- list(go_up = go_up, go_down = go_down, gsea = gsea)
  }

  message("GO/GSEA tables written to: ", out_dir)
  if (has_enrichplot) message("GO/GSEA plots written to: ", plot_dir)
  invisible(enrichment)
}

# -----------------------------
# 7. Visualization Functions
# -----------------------------
plot_volcano_enhanced <- function(df, title) {
  df <- df[is.finite(df$logFC) & is.finite(df$logP), , drop = FALSE]

  if (nrow(df) == 0) {
    return(
      ggplot() +
        theme_minimal(base_size = 12) +
        labs(title = title, x = "log2 Fold Change", y = "-log10(p-value)") +
        annotate("text", x = 0, y = 0, label = "No data to plot", size = 4)
    )
  }

  top_n <- min(30, nrow(df))
  top_genes <- rownames(df)[order(df$adj.P.Val)][seq_len(top_n)]

  df$label <- ifelse(
    rownames(df) %in% top_genes | grepl("^MT-", rownames(df)) | rownames(df) == "TPPP3",
    rownames(df),
    NA
  )

  ggplot(df, aes(x = logFC, y = logP)) +
    geom_point(aes(color = adj.P.Val < 0.05), alpha = 0.6) +
    scale_color_manual(values = c("grey70", "firebrick3"), name = "FDR < 0.05") +
    theme_minimal(base_size = 12) +
    labs(title = title, x = "log2 Fold Change", y = "-log10(p-value)") +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
    geom_text_repel(aes(label = label), size = 3, max.overlaps = 20, na.rm = TRUE)
}

plot_heatmap_enhanced <- function(df, expr_data, ann_subset, title) {
  if (is.null(ann_subset) || nrow(ann_subset) == 0) {
    message("Skipping heatmap (no samples): ", title)
    return(invisible(NULL))
  }

  zscore_mat <- t(scale(t(expr_data)))
  zscore_mat[is.na(zscore_mat)] <- 0

  ann_subset$sample_label <- paste(ann_subset$sample, ann_subset$patientID, sep = "_")

  top_n <- min(30, nrow(df))
  top_genes <- rownames(df)[order(df$adj.P.Val)][seq_len(top_n)]
  selected_genes <- unique(c(top_genes, rownames(df)[grepl("^TPP-", rownames(df))]))

  selected_genes <- intersect(selected_genes, rownames(zscore_mat))
  if (length(selected_genes) == 0) {
    message("Skipping heatmap (no matching genes): ", title)
    return(invisible(NULL))
  }

  keep_samples <- intersect(ann_subset$sample, colnames(zscore_mat))
  if (length(keep_samples) == 0) {
    message("Skipping heatmap (no matching samples): ", title)
    return(invisible(NULL))
  }

  mat <- zscore_mat[selected_genes, keep_samples, drop = FALSE]
  colnames(mat) <- ann_subset$sample_label[match(keep_samples, ann_subset$sample)]

  ann_col <- ann_subset[match(keep_samples, ann_subset$sample), c("SensoryStatus", "timepoint"), drop = FALSE]
  rownames(ann_col) <- colnames(mat)

  if (is.null(ann_col) || nrow(ann_col) == 0 || nrow(mat) == 0 || ncol(mat) == 0) {
    message("Skipping heatmap (insufficient data after filtering): ", title)
    return(invisible(NULL))
  }

  ann_col <- ann_col[, colSums(!is.na(ann_col)) > 0, drop = FALSE]
  if (ncol(ann_col) == 0) {
    message("Skipping heatmap (no valid annotation columns): ", title)
    return(invisible(NULL))
  }

  tryCatch({
    pheatmap(mat,
             annotation_col = ann_col,
             show_rownames = TRUE,
             show_colnames = TRUE,
             cluster_cols = TRUE,
             cluster_rows = TRUE,
             main = title,
             fontsize_row = 8,
             fontsize_col = 6,
             color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
  }, error = function(e) {
    message("Skipping heatmap due to plotting error: ", e$message)
  })
}

# -----------------------------
# 8. Execution Workflow
# -----------------------------
main <- function() {
  data <- load_and_clean_data(config)

  de_out <- run_de_analysis(data$counts, data$ann, config)
  results <- de_out$results
  expr_mat <- de_out$vst_mat

  mapped_names <- map_gene_symbols(rownames(expr_mat), config)
  rownames(expr_mat) <- mapped_names

  for (nm in names(results)) {
    rownames(results[[nm]]) <- mapped_names[match(rownames(results[[nm]]), rownames(de_out$vst_mat))]
  }
  enrichment <- run_go_and_gsea(results, method_name = "DESeq2", cfg = config)
  attr(results, "enrichment") <- enrichment

  for (tp in c("T0", "T1", "T2")) {
    res_name <- paste0("SensoryLoss_vs_Normal_", tp)
    tp_label <- switch(tp, T0 = "Inclusion", T1 = "3_weeks", T2 = "2_months")

    print(plot_volcano_enhanced(results[[res_name]], paste("Volcano:", res_name)))

    ann_sub <- data$ann[data$ann$timepoint == tp_label, , drop = FALSE]
    plot_heatmap_enhanced(results[[res_name]], expr_mat, ann_sub, paste("Heatmap:", res_name))
  }

  message("Workflow complete.")
  return(results)
}

# Uncomment to run:
results <- main()
