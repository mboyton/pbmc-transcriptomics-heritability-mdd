# ============================================================
# Prepare Lawlor PBMC CITE-seq data for downstream analyses
# ------------------------------------------------------------
# Purpose:
#   Load the Lawlor PBMC RNA count matrix and cell annotations,
#   align metadata to the expression matrix, remove cells without
#   a cell-type annotation, and save a filtered analysis object.
#
# Input:
#   - CZI.PBMC.RNA.matrix.Rds
#   - CZI.PBMC.cell.annotations.csv
#
# Output:
#   lawlor_rna_counts_and_metadata.filtered.rds
#
# The saved object contains:
#   counts        : genes x cells sparse count matrix
#   cell_metadata : cell-level metadata aligned to counts
# ============================================================

suppressPackageStartupMessages({
  library(Matrix)
})

# ============================================================
# Paths
# ============================================================

base_dir <- "/path/to/Lawlor_PBMC_CITEseq"

raw_dir <- file.path(
  base_dir,
  "01_raw"
)

out_dir <- file.path(
  base_dir,
  "02_processed",
  "lawlor_preprocessed"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rna_file <- file.path(
  raw_dir,
  "CZI.PBMC.RNA.matrix.Rds"
)

metadata_file <- file.path(
  raw_dir,
  "CZI.PBMC.cell.annotations.csv"
)

out_file <- file.path(
  out_dir,
  "lawlor_rna_counts_and_metadata.filtered.rds"
)

# ============================================================
# Load expression matrix and metadata
# ============================================================

counts <- readRDS(rna_file)

cell_metadata <- read.csv(
  metadata_file,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat(
  "RNA matrix dimensions (genes x cells): ",
  nrow(counts), " x ", ncol(counts), "\n",
  sep = ""
)

cat(
  "Metadata dimensions (cells x variables): ",
  nrow(cell_metadata), " x ", ncol(cell_metadata), "\n",
  sep = ""
)

# ============================================================
# Align RNA matrix and metadata
# ============================================================

missing_metadata <- setdiff(
  colnames(counts),
  rownames(cell_metadata)
)

if (length(missing_metadata) > 0) {
  stop(
    length(missing_metadata),
    " RNA cells are missing from the metadata."
  )
}

# Reorder metadata to exactly match expression-matrix columns
cell_metadata <- cell_metadata[
  colnames(counts),
  ,
  drop = FALSE
]

stopifnot(
  identical(colnames(counts), rownames(cell_metadata))
)

# ============================================================
# Remove cells without a cell-type annotation
# ============================================================

if (!"Celltype_Annotation" %in% colnames(cell_metadata)) {
  stop("Metadata does not contain 'Celltype_Annotation'.")
}

keep_cells <- !is.na(cell_metadata$Celltype_Annotation)

counts <- counts[, keep_cells, drop = FALSE]
cell_metadata <- cell_metadata[keep_cells, , drop = FALSE]

stopifnot(
  identical(colnames(counts), rownames(cell_metadata))
)

cat(
  "Cells retained after filtering: ",
  ncol(counts), "\n",
  sep = ""
)

# ============================================================
# Save processed object
# ============================================================

lawlor_data <- list(
  counts = counts,
  cell_metadata = cell_metadata
)

saveRDS(
  lawlor_data,
  out_file
)

cat(
  "Saved processed Lawlor object to:\n",
  out_file,
  "\n"
)