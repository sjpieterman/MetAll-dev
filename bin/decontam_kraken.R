#!/usr/bin/env Rscript
library(decontam)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: decontam_kraken.R <metadata.csv> <neg_col> <method> <threshold> <report1> [report2 ...]")
}

metadata_file <- args[1]
neg_col <- args[2]
method <- args[3]
threshold <- as.numeric(args[4])
report_files <- args[5:length(args)]

# Load metadata
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
# Assume sample IDs match the filenames (removing .kraken2.report)
# We need to match metadata to the report files.
# Usually, the process tag $sample_id is used in the filename.

# Read reports and build a count matrix of directly assigned reads (col 3)
taxa_counts <- list()
all_taxids <- c()

for (f in report_files) {
  sample_id <- basename(f)
  sample_id <- gsub("\\.kraken2\\.report\\.txt$", "", sample_id)
  sample_id <- gsub("\\.kraken2\\.report$", "", sample_id)
  report <- read.delim(f, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
  # col 3: directly assigned, col 5: taxid
  counts <- report[[3]]
  names(counts) <- report[[5]]
  taxa_counts[[sample_id]] <- counts
  all_taxids <- unique(c(all_taxids, report[[5]]))
}

# Build matrix
mat <- matrix(0, nrow = length(all_taxids), ncol = length(taxa_counts),
              dimnames = list(all_taxids, names(taxa_counts)))

for (sample_id in names(taxa_counts)) {
  counts <- taxa_counts[[sample_id]]
  mat[names(counts), sample_id] <- counts
}

# Filter metadata to only include samples we have reports for
samples_in_reports <- colnames(mat)
metadata_filtered <- metadata[metadata$sample %in% samples_in_reports, ]
# Sort metadata to match matrix columns
metadata_filtered <- metadata_filtered[match(samples_in_reports, metadata_filtered$sample), ]

# Logical vector for negative controls
# Handle various ways it might be specified (TRUE, "TRUE", "1", "Negative Control")
is_neg <- as.logical(metadata_filtered[[neg_col]])
if (any(is.na(is_neg))) {
    # Try string matching if logical conversion failed
    val <- as.character(metadata_filtered[[neg_col]])
    is_neg <- grepl("neg|control|true|yes|1", val, ignore.case = TRUE)
}

if (sum(is_neg) == 0) {
  message("No negative controls found in column ", neg_col)
  # Write empty contaminants file
  writeLines(character(0), "contaminants.txt")
  quit(save="no")
}

# Run decontam
# Transpose matrix for decontam (samples as rows)
mat_t <- t(mat)
contam_df <- isContaminant(mat_t, neg=is_neg, method=method, threshold=threshold)

# Extract contaminant taxids
contaminants <- rownames(contam_df)[contam_df$contaminant]

# Write out contaminant taxids
writeLines(contaminants, "contaminants.txt")

message("Identified ", length(contaminants), " contaminant taxids.")
