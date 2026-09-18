# ============================================================
# Identify resting cell-type markers
# ------------------------------------------------------------
# Focal resting cell type vs all other Baseline cells
# Lawlor PBMC CITE-seq pseudobulk + edgeR QL-GLM
#
# Set focal_celltype to one of:
#   "B"
#   "CD14_Mono"
#   "CD4T_Mem"
#   "CD4T_Naive"
#   "CD8T_Mem"
#   "CD8T_Naive"
#   "NK"
#
# Output:
#   timestamped .rda containing:
#     rest_vs_other_markers
#     params
#     session_info
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

in_rds <- "/path/to/lawlor_rna_counts_and_metadata.filtered.rds"

base_out_dir <- "/path/to/resting_transcriptional_signatures"

# -----------------------------
# Output directory
# -----------------------------

analysis_date <- format(Sys.Date(), "%Y%m%d")
celltype_slug <- tolower(gsub("_", "-", focal_celltype))

out_dir <- file.path(
  base_out_dir,
  paste0(
    analysis_date,
    "_identify-cluster-markers-",
    celltype_slug,
    "-rest-vs-other-rest"
  )
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

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
# Define variables
# -----------------------------

meta$condition <- meta$HTO_Classification
meta$celltype  <- meta$Celltype_Annotation
meta$donor     <- meta$Donor_of_Origin

available_celltypes <- sort(unique(meta$celltype))

cat("Available cell types:\n")
print(available_celltypes)

if (!focal_celltype %in% available_celltypes) {
  stop(
    "focal_celltype not found in Celltype_Annotation.\n",
    "Requested: ", focal_celltype, "\n",
    "Available: ",
    paste(available_celltypes, collapse = ", ")
  )
}

# -----------------------------
# Define resting comparison
# -----------------------------

cond_keep <- baseline_condition

meta$group <- ifelse(
  meta$celltype == focal_celltype,
  focal_celltype,
  "Other"
)

meta$group <- factor(
  meta$group,
  levels = c("Other", focal_celltype)
)

# -----------------------------
# Restrict to Baseline cells
# -----------------------------

meta_sub <- meta[
  meta$condition %in% cond_keep,
  ,
  drop = FALSE
]

if (nrow(meta_sub) == 0) {
  stop("No cells found for condition: ", cond_keep)
}

counts_sub <- counts[
  ,
  rownames(meta_sub),
  drop = FALSE
]

stopifnot(
  identical(colnames(counts_sub), rownames(meta_sub))
)

if (!all(c(focal_celltype, "Other") %in% unique(meta_sub$group))) {
  stop(
    "Both focal cell type and Other groups must be present within ",
    cond_keep,
    " cells.\n",
    "Focal cell type: ",
    focal_celltype
  )
}

cat(
  "\nCells per group within ",
  cond_keep,
  ":\n",
  sep = ""
)

print(table(meta_sub$group))

cat("\nCells per donor x group:\n")
print(table(meta_sub$donor, meta_sub$group))

# -----------------------------
# Pseudobulk:
# sum counts per donor x group
# -----------------------------

sample_id <- paste(
  meta_sub$donor,
  meta_sub$group,
  sep = "___"
)

grp <- factor(sample_id)

# cells x genes
counts_cxg <- t(counts_sub)

# Map each cell to its pseudobulk sample
M <- sparse.model.matrix(~ 0 + grp)

# genes x pseudobulk samples
pb_counts <- t(counts_cxg) %*% M

colnames(pb_counts) <- sub(
  "^grp",
  "",
  colnames(M)
)

# Pseudobulk metadata
pb_meta <- data.frame(
  sample_id = colnames(pb_counts),
  donor = sub("___.*$", "", colnames(pb_counts)),
  group = sub("^.*___", "", colnames(pb_counts)),
  stringsAsFactors = FALSE,
  row.names = colnames(pb_counts)
)

pb_meta$donor <- factor(pb_meta$donor)

pb_meta$group <- factor(
  pb_meta$group,
  levels = c("Other", focal_celltype)
)

if (any(is.na(pb_meta$group))) {
  stop(
    "Unexpected group labels found after parsing sample_id."
  )
}

if (length(unique(pb_meta$group)) < 2) {
  stop(
    "Only one group present after subsetting; ",
    "cannot run focal vs Other contrast."
  )
}

cat("\nPseudobulk samples per group:\n")
print(table(pb_meta$group))

# -----------------------------
# edgeR QL-GLM
# donor-blocked
# -----------------------------

y <- DGEList(
  counts = pb_counts
)

keep <- filterByExpr(
  y,
  group = pb_meta$group
)

y <- y[
  keep,
  ,
  keep.lib.sizes = FALSE
]

# TMM normalisation
y <- calcNormFactors(y)

# Group effect with donor blocking
design <- model.matrix(
  ~ 0 + group + donor,
  data = pb_meta
)

colnames(design) <- make.names(
  colnames(design)
)

y <- estimateDisp(
  y,
  design
)

fit <- glmQLFit(
  y,
  design
)

# -----------------------------
# Contrast:
# focal cell type vs Other
# -----------------------------

coef_other <- make.names(
  "groupOther"
)

coef_focal <- make.names(
  paste0("group", focal_celltype)
)

if (!all(
  c(coef_other, coef_focal) %in% colnames(design)
)) {
  stop(
    "Could not find expected group coefficients in design.\n",
    "Expected: ",
    coef_other,
    " and ",
    coef_focal,
    "\n",
    "Available: ",
    paste(colnames(design), collapse = ", ")
  )
}

contrast <- rep(
  0,
  ncol(design)
)

names(contrast) <- colnames(design)

contrast[coef_focal] <- 1
contrast[coef_other] <- -1

res <- glmQLFTest(
  fit,
  contrast = contrast
)

tt <- topTags(
  res,
  n = Inf
)$table

rest_vs_other_markers <- as.data.frame(tt)

rest_vs_other_markers$gene <- rownames(
  rest_vs_other_markers
)

rest_vs_other_markers$contrast <- paste0(
  cond_keep,
  ": ",
  focal_celltype,
  " vs Other"
)

# -----------------------------
# Save results
# -----------------------------

timestamp <- format(
  Sys.time(),
  "%Y%m%d_%H%M%S"
)

out_prefix <- paste0(
  celltype_slug,
  "_rest_vs_other_rest_markers_"
)

out_file <- file.path(
  out_dir,
  paste0(
    out_prefix,
    timestamp,
    ".rda"
  )
)

params <- list(
  in_rds = in_rds,
  out_dir = out_dir,
  condition = cond_keep,
  focal_celltype = focal_celltype,
  comparison = paste0(
    focal_celltype,
    " vs all other ",
    cond_keep,
    " cells"
  ),
  grouping = paste0(
    "group = ifelse(Celltype_Annotation == '",
    focal_celltype,
    "', '",
    focal_celltype,
    "', 'Other')"
  ),
  pseudobulk = paste0(
    "sum counts per donor x group within ",
    cond_keep
  ),
  model = "~ 0 + group + donor",
  contrast = paste0(
    coef_focal,
    " - ",
    coef_other
  )
)

session_info <- sessionInfo()

save(
  rest_vs_other_markers,
  params,
  session_info,
  file = out_file
)

cat(
  "\nSaved: ",
  out_file,
  "\n",
  sep = ""
)

cat(
  "Rows (genes tested): ",
  nrow(rest_vs_other_markers),
  "\n",
  sep = ""
)

cat("\nTop 10 genes by FDR:\n")

print(
  head(
    rest_vs_other_markers[
      order(rest_vs_other_markers$FDR),
    ],
    10
  )
)