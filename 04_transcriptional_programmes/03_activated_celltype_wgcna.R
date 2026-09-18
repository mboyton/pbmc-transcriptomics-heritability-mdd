# ============================================================
# WGCNA of lineage-specific activation pseudobulk expression
# ------------------------------------------------------------
# Input:
#   .rda file containing:
#     pb$expr : genes x samples
#     pb$meta : samples x metadata
#
# Analysis:
#   - variance filtering
#   - automatic soft-threshold selection
#   - signed WGCNA network
#   - module association with activation condition
#   - gene-level kME and expression summaries
# ============================================================

suppressPackageStartupMessages({
  library(WGCNA)
  library(dplyr)
  library(data.table)
})

enableWGCNAThreads()

# ============================================================
# User inputs
# ============================================================

in_rda <- "/path/to/pseudobulk_for_wgcna.rda"
out_dir <- "/path/to/wgcna_results"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# -----------------------------
# Analysis parameters
# -----------------------------

variance_quantile_keep  <- 0.50
soft_power_candidates   <- c(1:10, 12, 14, 16, 18, 20)
target_sft_r2           <- 0.80
network_type            <- "signed"
tom_type                <- "signed"
min_module_size         <- 50
merge_cut_height        <- 0.25
module_trait_fdr_cutoff <- 0.05

# ============================================================
# Output naming
# ============================================================

timestamp <- format(
  Sys.time(),
  "%Y%m%d_%H%M%S"
)

run_name <- tools::file_path_sans_ext(
  basename(in_rda)
)

run_dir <- file.path(
  out_dir,
  paste0(
    "wgcna_",
    run_name,
    "_",
    timestamp
  )
)

dir.create(
  run_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Load input
# ============================================================

load(in_rda)

if (!exists("pb")) {
  stop("Input .rda does not contain object 'pb'.")
}

if (!all(c("expr", "meta") %in% names(pb))) {
  stop("Object 'pb' must contain 'expr' and 'meta'.")
}

expr <- pb$expr
meta <- pb$meta

cat("Loaded input.\n")

cat(
  "Expression dim (genes x samples): ",
  nrow(expr),
  " x ",
  ncol(expr),
  "\n",
  sep = ""
)

cat(
  "Metadata dim (samples x vars): ",
  nrow(meta),
  " x ",
  ncol(meta),
  "\n",
  sep = ""
)

if (!identical(colnames(expr), rownames(meta))) {
  stop("colnames(pb$expr) must match rownames(pb$meta).")
}

if (!"condition" %in% colnames(meta)) {
  stop("pb$meta must contain a 'condition' column.")
}

# ============================================================
# Condition handling
# ============================================================

conditions_present <- unique(
  as.character(meta$condition)
)

if (length(conditions_present) != 2) {
  stop(
    "This script expects exactly 2 conditions in pb$meta$condition.\n",
    "Found: ",
    paste(conditions_present, collapse = ", ")
  )
}

# Preserve condition order from the input metadata
condition_levels <- unique(
  as.character(meta$condition)
)

condition <- factor(
  meta$condition,
  levels = condition_levels
)

baseline_condition <- condition_levels[1]
stim_condition <- condition_levels[2]

cat("Conditions detected:\n")
print(condition_levels)

cat(
  "Reference/baseline condition: ",
  baseline_condition,
  "\n",
  sep = ""
)

cat(
  "Stim/contrast condition: ",
  stim_condition,
  "\n",
  sep = ""
)

# ============================================================
# Prepare WGCNA input
# WGCNA expects samples x genes
# ============================================================

datExpr <- t(expr)

cat(
  "WGCNA input dim (samples x genes): ",
  nrow(datExpr),
  " x ",
  ncol(datExpr),
  "\n",
  sep = ""
)

if (!identical(rownames(datExpr), rownames(meta))) {
  stop("rownames(datExpr) must match rownames(meta).")
}

# ============================================================
# Remove problematic samples / genes
# ============================================================

gsg <- goodSamplesGenes(
  datExpr,
  verbose = 3
)

if (!gsg$allOK) {

  cat(
    "Removing bad samples / genes identified by ",
    "goodSamplesGenes().\n"
  )

  datExpr <- datExpr[
    gsg$goodSamples,
    gsg$goodGenes,
    drop = FALSE
  ]

  meta <- meta[
    rownames(datExpr),
    ,
    drop = FALSE
  ]

  condition <- factor(
    meta$condition,
    levels = condition_levels
  )
}

cat(
  "After goodSamplesGenes: ",
  nrow(datExpr),
  " samples x ",
  ncol(datExpr),
  " genes\n",
  sep = ""
)

# ============================================================
# Variance filtering
# ============================================================

gene_sd <- apply(
  datExpr,
  2,
  sd,
  na.rm = TRUE
)

var_cut <- quantile(
  gene_sd,
  probs = 1 - variance_quantile_keep,
  na.rm = TRUE
)

datExpr_filt <- datExpr[
  ,
  gene_sd > var_cut,
  drop = FALSE
]

cat(
  "After variance filtering: ",
  nrow(datExpr_filt),
  " samples x ",
  ncol(datExpr_filt),
  " genes\n",
  sep = ""
)

cat(
  "Variance filter kept top ",
  variance_quantile_keep * 100,
  "% most variable genes\n",
  sep = ""
)

if (ncol(datExpr_filt) < min_module_size) {
  stop(
    "Too few genes remain after variance filtering for WGCNA."
  )
}

# ============================================================
# Soft-threshold selection
# ============================================================

sft <- pickSoftThreshold(
  datExpr_filt,
  powerVector = soft_power_candidates,
  networkType = network_type,
  verbose = 5
)

fit_tbl <- as.data.frame(
  sft$fitIndices
)

fit_tbl$signed_R2 <-
  -sign(fit_tbl$slope) * fit_tbl$SFT.R.sq

candidate_ok <- fit_tbl %>%
  filter(
    !is.na(signed_R2),
    signed_R2 >= target_sft_r2
  ) %>%
  arrange(Power)

if (nrow(candidate_ok) > 0) {

  softPower <- candidate_ok$Power[1]

  softpower_rule <- paste0(
    "lowest power with signed_R2 >= ",
    target_sft_r2
  )

} else {

  best_row <- fit_tbl %>%
    filter(!is.na(signed_R2)) %>%
    arrange(
      desc(signed_R2),
      Power
    ) %>%
    slice(1)

  softPower <- best_row$Power

  softpower_rule <- paste0(
    "max signed_R2 because no power reached ",
    target_sft_r2
  )
}

cat(
  "Selected softPower: ",
  softPower,
  "\n",
  sep = ""
)

cat(
  "Selection rule: ",
  softpower_rule,
  "\n",
  sep = ""
)

# ============================================================
# Soft-threshold diagnostic plots
# ============================================================

png(
  file.path(
    run_dir,
    "soft_threshold_diagnostics.png"
  ),
  width = 1600,
  height = 700,
  res = 120
)

par(mfrow = c(1, 2))

plot(
  fit_tbl$Power,
  fit_tbl$signed_R2,
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit (signed R^2)",
  type = "n",
  main = "Scale independence"
)

text(
  fit_tbl$Power,
  fit_tbl$signed_R2,
  labels = fit_tbl$Power,
  col = "red"
)

abline(
  h = target_sft_r2,
  lty = 2
)

abline(
  v = softPower,
  lty = 3
)

plot(
  fit_tbl$Power,
  fit_tbl$mean.k.,
  xlab = "Soft Threshold (power)",
  ylab = "Mean connectivity",
  type = "n",
  main = "Mean connectivity"
)

text(
  fit_tbl$Power,
  fit_tbl$mean.k.,
  labels = fit_tbl$Power,
  col = "red"
)

abline(
  v = softPower,
  lty = 3
)

par(mfrow = c(1, 1))
dev.off()

# ============================================================
# Network construction and module detection
# ============================================================

net <- blockwiseModules(
  datExpr_filt,
  power = softPower,
  networkType = network_type,
  TOMType = tom_type,
  minModuleSize = min_module_size,
  reassignThreshold = 0,
  mergeCutHeight = merge_cut_height,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

moduleColors <- labels2colors(
  net$colors
)

cat("Detected modules:\n")
print(table(moduleColors))

# ============================================================
# Module dendrogram
# ============================================================

png(
  file.path(
    run_dir,
    "module_dendrogram.png"
  ),
  width = 1800,
  height = 900,
  res = 120
)

plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)

dev.off()

# ============================================================
# Module eigengenes
# ============================================================

MEs0 <- moduleEigengenes(
  datExpr_filt,
  colors = moduleColors
)$eigengenes

MEs <- orderMEs(MEs0)

cat(
  "Module eigengenes dim: ",
  nrow(MEs),
  " x ",
  ncol(MEs),
  "\n",
  sep = ""
)

# ============================================================
# Association with activation
# Baseline = 0; stimulated = 1
# ============================================================

trait_activation <- as.numeric(
  condition == stim_condition
)

names(trait_activation) <- rownames(meta)

if (!identical(
  rownames(MEs),
  names(trait_activation)
)) {
  stop(
    "Module eigengenes and trait vector are misaligned."
  )
}

moduleTraitCor <- cor(
  MEs,
  trait_activation,
  use = "p"
)

moduleTraitP <- corPvalueStudent(
  moduleTraitCor,
  nSamples = nrow(datExpr_filt)
)

module_trait_table <- data.frame(
  module = rownames(moduleTraitCor),
  cor = as.numeric(moduleTraitCor),
  pval = as.numeric(moduleTraitP),
  stringsAsFactors = FALSE
) %>%
  mutate(
    fdr = p.adjust(
      pval,
      method = "fdr"
    ),
    direction = case_when(
      cor > 0 ~ paste0(
        "higher_in_",
        stim_condition
      ),
      cor < 0 ~ paste0(
        "lower_in_",
        stim_condition
      ),
      TRUE ~ "no_direction"
    )
  ) %>%
  arrange(
    fdr,
    desc(abs(cor))
  )

cat("Module-trait table:\n")
print(module_trait_table)

# ============================================================
# Module-trait heatmap
# ============================================================

textMatrix <- paste0(
  signif(moduleTraitCor, 2),
  "\n(",
  signif(moduleTraitP, 1),
  ")"
)

dim(textMatrix) <- dim(
  moduleTraitCor
)

png(
  file.path(
    run_dir,
    "module_trait_heatmap.png"
  ),
  width = 900,
  height = 900,
  res = 120
)

par(
  mar = c(
    6,
    8.5,
    3,
    3
  )
)

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = "Activation",
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.8,
  zlim = c(-1, 1),
  main = "Module-trait relationships"
)

dev.off()

# ============================================================
# Significant activation-associated modules
# ============================================================

significant_modules <- module_trait_table %>%
  filter(
    fdr < module_trait_fdr_cutoff
  ) %>%
  mutate(
    module_color = sub(
      "^ME",
      "",
      module
    )
  )

cat(
  "Significant modules (FDR < ",
  module_trait_fdr_cutoff,
  "):\n",
  sep = ""
)

print(significant_modules)

# ============================================================
# Gene-level table
# ============================================================

gene_ids <- colnames(
  datExpr_filt
)

if (length(moduleColors) != length(gene_ids)) {
  stop(
    "Length of moduleColors does not match number of genes."
  )
}

stim_idx <- condition == stim_condition
base_idx <- condition == baseline_condition

gene_mean_base <- colMeans(
  datExpr_filt[
    base_idx,
    ,
    drop = FALSE
  ],
  na.rm = TRUE
)

gene_mean_stim <- colMeans(
  datExpr_filt[
    stim_idx,
    ,
    drop = FALSE
  ],
  na.rm = TRUE
)

gene_mean_diff <-
  gene_mean_stim - gene_mean_base

gene_abs_mean_diff <- abs(
  gene_mean_diff
)

gene_module_table <- data.frame(
  gene = gene_ids,
  module = moduleColors,
  mean_logcpm_baseline = gene_mean_base[gene_ids],
  mean_logcpm_stim = gene_mean_stim[gene_ids],
  mean_logcpm_diff_stim_minus_base =
    gene_mean_diff[gene_ids],
  abs_mean_logcpm_diff_stim_minus_base =
    gene_abs_mean_diff[gene_ids],
  stringsAsFactors = FALSE
)

# ============================================================
# kME
# ============================================================

kME_all <- as.data.frame(
  cor(
    datExpr_filt,
    MEs,
    use = "p"
  )
)

kME_all$gene <- rownames(
  kME_all
)

kME_long <- as.data.table(
  kME_all
)

kME_long <- melt(
  kME_long,
  id.vars = "gene",
  variable.name = "ME",
  value.name = "kME"
)

kME_long <- as.data.frame(
  kME_long
) %>%
  mutate(
    module = sub(
      "^ME",
      "",
      ME
    )
  )

# Keep kME for each gene's assigned module
gene_module_table <- gene_module_table %>%
  left_join(
    kME_long %>%
      select(
        gene,
        module,
        kME
      ),
    by = c(
      "gene",
      "module"
    )
  )

# Add activation association statistics
gene_module_table <- gene_module_table %>%
  left_join(
    significant_modules %>%
      select(
        module_color,
        cor,
        pval,
        fdr,
        direction
      ) %>%
      rename(
        module = module_color,
        module_trait_cor = cor,
        module_trait_pval = pval,
        module_trait_fdr = fdr,
        module_trait_direction = direction
      ),
    by = "module"
  )

# ============================================================
# Significant-module gene tables
# ============================================================

significant_module_gene_tables <- list()

if (nrow(significant_modules) > 0) {

  for (mod in significant_modules$module_color) {

    tbl_mod <- gene_module_table %>%
      filter(
        module == mod
      ) %>%
      arrange(
        desc(abs(kME)),
        desc(
          abs_mean_logcpm_diff_stim_minus_base
        )
      )

    significant_module_gene_tables[[mod]] <-
      tbl_mod

    fwrite(
      tbl_mod,
      file = file.path(
        run_dir,
        paste0(
          "genes_in_module_",
          mod,
          ".tsv"
        )
      ),
      sep = "\t"
    )
  }
}

# ============================================================
# Save summary tables
# ============================================================

fwrite(
  fit_tbl,
  file = file.path(
    run_dir,
    "soft_threshold_fit_table.tsv"
  ),
  sep = "\t"
)

fwrite(
  module_trait_table,
  file = file.path(
    run_dir,
    "module_trait_table.tsv"
  ),
  sep = "\t"
)

fwrite(
  gene_module_table,
  file = file.path(
    run_dir,
    "gene_module_table_all.tsv"
  ),
  sep = "\t"
)

fwrite(
  meta %>%
    as.data.frame() %>%
    tibble::rownames_to_column(
      "sample_id"
    ),
  file = file.path(
    run_dir,
    "sample_metadata.tsv"
  ),
  sep = "\t"
)

fwrite(
  as.data.frame(MEs) %>%
    tibble::rownames_to_column(
      "sample_id"
    ),
  file = file.path(
    run_dir,
    "module_eigengenes.tsv"
  ),
  sep = "\t"
)

# ============================================================
# Save full R object bundle
# ============================================================

wgcna_params <- list(
  in_rda = in_rda,
  run_dir = run_dir,
  variance_quantile_keep = variance_quantile_keep,
  soft_power_candidates = soft_power_candidates,
  target_sft_r2 = target_sft_r2,
  softPower = softPower,
  softpower_rule = softpower_rule,
  network_type = network_type,
  tom_type = tom_type,
  min_module_size = min_module_size,
  merge_cut_height = merge_cut_height,
  module_trait_fdr_cutoff =
    module_trait_fdr_cutoff,
  baseline_condition =
    baseline_condition,
  stim_condition =
    stim_condition
)

session_info <- sessionInfo()

save(
  pb,
  meta,
  datExpr,
  datExpr_filt,
  condition,
  sft,
  fit_tbl,
  net,
  moduleColors,
  MEs,
  module_trait_table,
  significant_modules,
  gene_module_table,
  significant_module_gene_tables,
  wgcna_params,
  session_info,
  file = file.path(
    run_dir,
    "wgcna_results_bundle.rda"
  )
)

cat(
  "\n============================================================\n"
)

cat("WGCNA complete.\n")

cat(
  "Run directory:\n",
  run_dir,
  "\n",
  sep = ""
)

cat(
  "Selected softPower: ",
  softPower,
  "\n",
  sep = ""
)

cat(
  "Significant modules: ",
  nrow(significant_modules),
  "\n",
  sep = ""
)

if (nrow(significant_modules) > 0) {
  cat(
    "Module names: ",
    paste(
      significant_modules$module_color,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
}

cat(
  "============================================================\n"
)