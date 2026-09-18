# ============================================================
# Export Lawlor PBMC data in CELLEX-compatible format
# ------------------------------------------------------------
# Purpose:
#   Convert the filtered Lawlor RNA count matrix and aligned
#   metadata into MatrixMarket and compressed CSV files for
#   downstream CELLEX analysis.
#
# Input:
#   lawlor_rna_counts_and_metadata.filtered.rds
#
# Output:
#   lawlor_rna.filtered.counts.mtx.gz
#   lawlor_genes.csv.gz
#   lawlor_cells_for_cellex.csv.gz
# ============================================================

suppressPackageStartupMessages({
  library(Matrix)
})

# ============================================================
# Paths
# ============================================================

base_dir <- "/path/to/Lawlor_PBMC_CITEseq"

in_file <- file.path(
  base_dir,
  "02_processed",
  "lawlor_preprocessed",
  "lawlor_rna_counts_and_metadata.filtered.rds"
)

out_dir <- file.path(
  base_dir,
  "02_processed",
  "cellex_inputs"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load filtered Lawlor object
# ============================================================

obj <- readRDS(in_file)

mat <- obj$counts
meta <- obj$cell_metadata

cat(
  "Expression matrix dimensions (genes x cells): ",
  nrow(mat), " x ", ncol(mat), "\n",
  sep = ""
)

cat(
  "Metadata dimensions (cells x variables): ",
  nrow(meta), " x ", ncol(meta), "\n",
  sep = ""
)

# ============================================================
# Confirm alignment
# ============================================================

stopifnot(
  identical(colnames(mat), rownames(meta))
)

# ============================================================
# Prepare gene metadata
# ============================================================

gene_meta <- data.frame(
  gene_id = rownames(mat),
  stringsAsFactors = FALSE
)

# ============================================================
# Prepare cell metadata for CELLEX
# ============================================================

required_cols <- c(
  "HTO_Barcodes",
  "HTO_Classification",
  "Run_Identifier",
  "Demuxlet_Classification",
  "Donor_of_Origin",
  "Celltype_Annotation"
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

cell_meta <- data.frame(
  cell_id = rownames(meta),
  HTO_Barcodes = meta$HTO_Barcodes,
  HTO_Classification = meta$HTO_Classification,
  Run_Identifier = meta$Run_Identifier,
  Demuxlet_Classification = meta$Demuxlet_Classification,
  Donor_of_Origin = meta$Donor_of_Origin,
  Celltype_Annotation = meta$Celltype_Annotation,
  stringsAsFactors = FALSE
)

# ============================================================
# Output paths
# ============================================================

mtx_file <- file.path(
  out_dir,
  "lawlor_rna.filtered.counts.mtx"
)

genes_file <- file.path(
  out_dir,
  "lawlor_genes.csv.gz"
)

cells_file <- file.path(
  out_dir,
  "lawlor_cells_for_cellex.csv.gz"
)

# ============================================================
# Export count matrix
# ============================================================

writeMM(
  mat,
  file = mtx_file
)

# Compress MatrixMarket file
system2(
  "gzip",
  args = c("-f", shQuote(mtx_file))
)

# ============================================================
# Export gene and cell metadata
# ============================================================

write.csv(
  gene_meta,
  gzfile(genes_file),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  cell_meta,
  gzfile(cells_file),
  row.names = FALSE,
  quote = FALSE
)

# ============================================================
# Summary
# ============================================================

cat("CELLEX input export complete.\n")
cat("Files written:\n")
cat(" - ", paste0(mtx_file, ".gz"), "\n", sep = "")
cat(" - ", genes_file, "\n", sep = "")
cat(" - ", cells_file, "\n", sep = "")