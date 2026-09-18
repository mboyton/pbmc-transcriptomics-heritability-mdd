# ============================================================
# Add combined cell-type and stimulation-condition labels
# ------------------------------------------------------------
# Purpose:
#   Create a combined cell-state label for CELLEX by joining
#   the Lawlor cell-type annotation and stimulation condition.
#
# Input:
#   lawlor_cells_for_cellex.csv.gz
#
# Output:
#   The same metadata file, updated with:
#     cell_type_condition
#
# Example:
#   NK + LPS -> NK_LPS
# ============================================================

# ============================================================
# Paths
# ============================================================

base_dir <- "/path/to/Lawlor_PBMC_CITEseq"

meta_file <- file.path(
  base_dir,
  "02_processed",
  "cellex_inputs",
  "lawlor_cells_for_cellex.csv.gz"
)

# ============================================================
# Load CELLEX cell metadata
# ============================================================

meta <- read.csv(
  gzfile(meta_file),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ============================================================
# Check required columns
# ============================================================

required_cols <- c(
  "Celltype_Annotation",
  "HTO_Barcodes"
)

missing_cols <- setdiff(
  required_cols,
  colnames(meta)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing expected metadata columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# ============================================================
# Create combined cell-state label
# ============================================================

meta$cell_type_condition <- paste(
  meta$Celltype_Annotation,
  meta$HTO_Barcodes,
  sep = "_"
)

# ============================================================
# Inspect resulting states
# ============================================================

cat("Cell-type / condition labels:\n")
print(unique(meta$cell_type_condition))

cat("\nCell counts per state:\n")
print(table(meta$cell_type_condition))

# ============================================================
# Overwrite CELLEX metadata file
# ============================================================

write.csv(
  meta,
  gzfile(meta_file),
  row.names = FALSE,
  quote = FALSE
)

cat(
  "\nUpdated CELLEX metadata written to:\n",
  meta_file,
  "\n"
)