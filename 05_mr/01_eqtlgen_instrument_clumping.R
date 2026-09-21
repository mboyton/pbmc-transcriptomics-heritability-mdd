library(data.table)
library(ieugwasr)
library(genetics.binaRies)

# =========================
# Paths
# =========================

eqtl_file <- "path/to/eQTLGen_cis_eQTLs.txt"
ref_dir <- "path/to/1000G_EUR_Phase3_plink"
out_dir <- "path/to/output"
log_dir <- "path/to/logs"

# =========================
# Timestamp
# =========================

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# =========================
# Output files
# =========================

out_file <- file.path(
 out_dir,
 paste0("eQTLGen_cis_eqtl_clumped_all_genes_", timestamp, ".txt")
)

log_file <- file.path(
 log_dir,
 paste0("eQTLGen_cis_eqtl_clumping_log_", timestamp, ".txt")
)

# =========================
# Settings
# =========================

clump_kb <- 10000
clump_r2 <- 0.001
clump_p  <- 5e-8

plink_bin <- genetics.binaRies::get_plink_binary()

# =========================
# Setup
# =========================

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

if (file.exists(out_file)) file.remove(out_file)
if (file.exists(log_file)) file.remove(log_file)

log_message <- function(msg) {
 cat(
  sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg),
  file = log_file,
  append = TRUE
 )
}

log_message("Starting eQTLGen per-gene LD clumping")

# =========================
# Read full dataset
# =========================

eqtl <- fread(eqtl_file)

# Basic cleanup and genome-wide significance filter
eqtl <- eqtl[
 !is.na(Gene) &
  !is.na(GeneSymbol) &
  !is.na(SNP) &
  !is.na(Pvalue) &
  !is.na(SNPChr) &
  Pvalue < 5e-8
]

# Unique genes
genes <- unique(eqtl$Gene)

log_message(
 sprintf(
  "Loaded %s rows across %s unique genes",
  nrow(eqtl),
  length(genes)
 )
)

# =========================
# Loop across genes
# =========================

first_write <- TRUE

for (i in seq_along(genes)) {
 g <- genes[i]
 
 gene_dat <- eqtl[Gene == g]
 
 if (nrow(gene_dat) == 0) {
  log_message(sprintf("SKIP: %s | no rows", g))
  next
 }
 
 gene_symbol <- unique(gene_dat$GeneSymbol)
 gene_symbol <- gene_symbol[!is.na(gene_symbol)]
 gene_symbol <- if (length(gene_symbol) > 0) gene_symbol[1] else g
 
 chr_vals <- unique(gene_dat$SNPChr)
 chr_vals <- chr_vals[!is.na(chr_vals)]
 
 if (length(chr_vals) == 0) {
  log_message(
   sprintf(
    "SKIP: %s (%s) | missing chromosome",
    g,
    gene_symbol
   )
  )
  next
 }
 
 if (length(chr_vals) > 1) {
  log_message(
   sprintf(
    "WARN: %s (%s) | multiple SNP chromosomes found: %s | using first",
    g,
    gene_symbol,
    paste(chr_vals, collapse = ",")
   )
  )
 }
 
 chr_use <- chr_vals[1]
 
 bfile_prefix <- file.path(
  ref_dir,
  paste0("1000G.EUR.QC.", chr_use)
 )
 
 if (!all(file.exists(
  paste0(bfile_prefix, c(".bed", ".bim", ".fam"))
 ))) {
  log_message(
   sprintf(
    "SKIP: %s (%s) | missing reference files for chr%s",
    g,
    gene_symbol,
    chr_use
   )
  )
  next
 }
 
 clump_input <- unique(
  gene_dat[, .(
   rsid = SNP,
   pval = Pvalue,
   id = gene_symbol
  )]
 )
 
 if (nrow(clump_input) == 0) {
  log_message(
   sprintf(
    "SKIP: %s (%s) | no clumpable SNPs",
    g,
    gene_symbol
   )
  )
  next
 }
 
 clumped <- tryCatch(
  {
   ld_clump(
    d = clump_input,
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    clump_p = clump_p,
    plink_bin = plink_bin,
    bfile = bfile_prefix
   )
  },
  error = function(e) {
   log_message(
    sprintf(
     "ERROR: %s (%s) | ld_clump failed | %s",
     g,
     gene_symbol,
     e$message
    )
   )
   return(NULL)
  }
 )
 
 if (is.null(clumped) || nrow(clumped) == 0) {
  log_message(
   sprintf(
    "SKIP: %s (%s) | no SNPs retained after clumping",
    g,
    gene_symbol
   )
  )
  next
 }
 
 gene_dat_clumped <- gene_dat[
  SNP %in% clumped$rsid
 ]
 
 if (nrow(gene_dat_clumped) == 0) {
  log_message(
   sprintf(
    "SKIP: %s (%s) | clumped rsids not found back in source table",
    g,
    gene_symbol
   )
  )
  next
 }
 
 gene_dat_clumped[, clump_gene := g]
 gene_dat_clumped[, clump_gene_symbol := gene_symbol]
 
 fwrite(
  gene_dat_clumped,
  file = out_file,
  sep = "\t",
  quote = FALSE,
  append = !first_write,
  col.names = first_write
 )
 
 first_write <- FALSE
 
 if (i %% 100 == 0) {
  log_message(
   sprintf(
    "Progress: processed %s / %s genes",
    i,
    length(genes)
   )
  )
 }
}

log_message("Finished eQTLGen per-gene LD clumping")
