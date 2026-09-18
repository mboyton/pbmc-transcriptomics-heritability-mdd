# ============================================================
# Prepare lineage-specific activation pseudobulk data for WGCNA
# ------------------------------------------------------------
# This script generalises the pseudobulk preparation used for
# the Lawlor activation analyses.
#
# Raw counts are summed within donor x condition, then converted
# to logCPM for input to the downstream WGCNA analysis.
#
# Output:
#   pb$expr : genes x pseudobulk samples (logCPM)
#   pb$meta : pseudobulk samples x metadata
# ============================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
})

# -----------------------------
# User settings
# -----------------------------

focal_celltype <- "<FOCAL_CELLTYPE>"

baseline_condition <- "Baseline"
stim_condition <- "<STIM_CONDITION>"

# Examples used in the study:
#
# CD4T_Mem   : Baseline vs CD3_CD28
# CD4T_Naive : Baseline vs CD3_CD28
# CD8T_Mem   : Baseline vs CD3_CD28
# CD8T_Naive : Baseline vs CD3_CD28
# NK         : Baseline vs CD3_CD28
# NK         : Baseline vs LPS

in_rds <- "/path/to/lawlor_rna_counts_and_metadata.filtered.rds"

out_dir <- "/path/to/wgcna_ready_data"

# -----------------------------
# Load Lawlor data
# -----------------------------

x <- readRDS(in_rds)

counts <- x$counts
meta <- x$cell_metadata

stopifnot(
  identical(colnames(counts), rownames(meta))
)

# -----------------------------
# Define analysis variables
# -----------------------------

meta$celltype  <- meta$Celltype_Annotation
meta$condition <- meta$HTO_Classification
meta$donor     <- meta$Donor_of_Origin

# -----------------------------
# Restrict to focal lineage and
# selected conditions
# -----------------------------

keep_cells <- rownames(meta)[
  meta$celltype == focal_celltype &
    meta$condition %in% c(
      baseline_condition,
      stim_condition
    )
]

counts_sub <- counts[
  ,
  keep_cells,
  drop = FALSE
]

meta_sub <- meta[
  keep_cells,
  ,
  drop = FALSE
]

stopifnot(
  identical(colnames(counts_sub), rownames(meta_sub))
)

cat(
  "Cell type: ",
  focal_celltype,
  "\n",
  sep = ""
)

cat("Cells per condition:\n")
print(table(meta_sub$condition))

cat("\nCells per donor x condition:\n")
print(
  table(
    meta_sub$donor,
    meta_sub$condition
  )
)

# -----------------------------
# Create donor x condition groups
# -----------------------------

meta_sub$sample_id <- paste(
  meta_sub$donor,
  meta_sub$condition,
  sep = "__"
)

sample_ids <- unique(
  meta_sub$sample_id
)

# -----------------------------
# Sum raw counts within each
# donor x condition
# -----------------------------

pb_counts <- sapply(
  sample_ids,
  function(s) {

    cells <- rownames(meta_sub)[
      meta_sub$sample_id == s
    ]

    Matrix::rowSums(
      counts_sub[
        ,
        cells,
        drop = FALSE
      ]
    )
  }
)

pb_counts <- as.matrix(pb_counts)

# -----------------------------
# Pseudobulk metadata
# -----------------------------

pb_meta <- unique(
  meta_sub[
    ,
    c(
      "sample_id",
      "donor",
      "condition"
    )
  ]
)

pb_meta <- pb_meta[
  match(
    colnames(pb_counts),
    pb_meta$sample_id
  ),
  ,
  drop = FALSE
]

rownames(pb_meta) <- pb_meta$sample_id

# -----------------------------
# Convert pseudobulk counts to logCPM
# -----------------------------

dge <- DGEList(
  counts = pb_counts
)

logcpm <- cpm(
  dge,
  log = TRUE,
  prior.count = 1
)

# -----------------------------
# Create WGCNA input object
# -----------------------------

pb <- list(
  expr = logcpm,
  meta = pb_meta
)

# -----------------------------
# Sanity checks
# -----------------------------

cat(
  "\nExpression dimensions (genes x samples): ",
  nrow(pb$expr),
  " x ",
  ncol(pb$expr),
  "\n",
  sep = ""
)

cat(
  "Metadata dimensions (samples x variables): ",
  nrow(pb$meta),
  " x ",
  ncol(pb$meta),
  "\n",
  sep = ""
)

stopifnot(
  identical(
    colnames(pb$expr),
    rownames(pb$meta)
  )
)

cat("\nPseudobulk samples per condition:\n")
print(table(pb$meta$condition))

cat("\nPseudobulk samples per donor:\n")
print(table(pb$meta$donor))

# -----------------------------
# Save
# -----------------------------

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

celltype_slug <- tolower(
  gsub("_", "-", focal_celltype)
)

stim_slug <- tolower(
  gsub("_", "-", stim_condition)
)

timestamp <- format(
  Sys.time(),
  "%Y%m%d_%H%M%S"
)

out_file <- file.path(
  out_dir,
  paste0(
    celltype_slug,
    "_rest_vs_",
    stim_slug,
    "_for_wgcna_",
    timestamp,
    ".rda"
  )
)

params <- list(
  focal_celltype = focal_celltype,
  baseline_condition = baseline_condition,
  stim_condition = stim_condition,
  pseudobulk = "sum raw counts per donor x condition",
  expression = "logCPM with prior.count = 1"
)

save(
  pb,
  params,
  file = out_file
)

cat(
  "\nSaved WGCNA input object to:\n",
  out_file,
  "\n"
)