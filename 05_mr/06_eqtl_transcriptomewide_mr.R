library(data.table)
library(TwoSampleMR)

# =====================================
# File paths
# =====================================
exposure_filepath <- "path/to/eQTLGen_MR_ready.txt"

outcome_filepath <- "path/to/outcome_GWAS.txt"

output_dir <- "path/to/output_directory"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

outcome_phenotype <- "OUTCOME_NAME"

# =====================================
# Load exposure
# =====================================
message("Loading exposure data...")
exp <- fread(exposure_filepath)

required_exp_cols <- c(
 "SNP",
 "beta.exposure",
 "se.exposure",
 "effect_allele.exposure",
 "other_allele.exposure",
 "eaf.exposure",
 "pval.exposure",
 "samplesize.exposure",
 "Gene",
 "GeneSymbol",
 "id.exposure"
)

missing_exp <- setdiff(required_exp_cols, names(exp))
if (length(missing_exp) > 0) {
 stop("Missing exposure columns: ", paste(missing_exp, collapse = ", "))
}

exp <- exp[
 !is.na(SNP) &
  !is.na(beta.exposure) &
  !is.na(se.exposure) &
  !is.na(effect_allele.exposure) &
  !is.na(other_allele.exposure) &
  !is.na(pval.exposure) &
  !is.na(id.exposure) &
  !is.na(Gene)
]

# Make exposure label unique and readable
exp[, id.exposure := paste0("eQTLGen_", Gene)]
exp[, exposure := paste0("eQTLGen_", GeneSymbol, "_", Gene)]

# Optional rare-variant filter
# exp <- exp[eaf.exposure > 0.01 & eaf.exposure < 0.99]

exp <- unique(exp)

message("Exposure rows: ", nrow(exp))
message("Unique ENSID exposures: ", uniqueN(exp$id.exposure))

# =====================================
# Load outcome
# =====================================
message("Loading outcome data...")
out <- fread(outcome_filepath)

required_out_cols <- c("SNP", "A1", "A2", "BETA", "SE", "P", "N")
missing_out <- setdiff(required_out_cols, names(out))
if (length(missing_out) > 0) {
 stop("Missing outcome columns: ", paste(missing_out, collapse = ", "))
}

out <- out[, .(
 SNP = SNP,
 beta.outcome = BETA,
 se.outcome = SE,
 effect_allele.outcome = A1,
 other_allele.outcome = A2,
 eaf.outcome = NA_real_,
 pval.outcome = P,
 samplesize.outcome = N,
 outcome = outcome_phenotype,
 id.outcome = outcome_phenotype
)]

out <- out[
 !is.na(SNP) &
  !is.na(beta.outcome) &
  !is.na(se.outcome) &
  !is.na(effect_allele.outcome) &
  !is.na(other_allele.outcome) &
  !is.na(pval.outcome)
]

out <- unique(out, by = "SNP")
setkey(out, SNP)

message("Outcome rows: ", nrow(out))

# =====================================
# Prepare unique ENSID exposure list
# =====================================
exposure_ids <- unique(exp[, .(id.exposure, exposure, Gene, GeneSymbol)])
setorder(exposure_ids, id.exposure)

total_exposures <- nrow(exposure_ids)
message("Total ENSID exposures to analyse: ", total_exposures)

# =====================================
# Batch settings
# =====================================
save_interval <- 1000
batch_number <- 1
all_files <- character()
ind <- 1

# =====================================
# Main loop
# =====================================
while (ind <= total_exposures) {
 batch_start <- ind
 batch_end <- min(total_exposures, batch_start + save_interval - 1)
 
 message("Starting batch ", batch_number, ": exposures ", batch_start, " to ", batch_end)
 
 res_batch <- vector("list", length = 0)
 
 for (i in batch_start:batch_end) {
  current_id <- exposure_ids$id.exposure[i]
  current_exposure <- exposure_ids$exposure[i]
  current_gene <- exposure_ids$Gene[i]
  current_symbol <- exposure_ids$GeneSymbol[i]
  
  message("Analysing exposure ", i, " of ", total_exposures, ": ", current_id)
  
  tryCatch({
   exp_data <- exp[id.exposure == current_id]
   
   if (nrow(exp_data) == 0) next
   
   out_data <- out[exp_data, on = "SNP", nomatch = 0]
   if (nrow(out_data) == 0) next
   
   exp_data <- exp_data[SNP %in% out_data$SNP]
   
   exp_df <- as.data.frame(exp_data[, .(
    SNP,
    beta.exposure,
    se.exposure,
    effect_allele.exposure,
    other_allele.exposure,
    eaf.exposure,
    pval.exposure,
    samplesize.exposure,
    exposure,
    id.exposure
   )])
   
   out_df <- as.data.frame(out_data[, .(
    SNP,
    beta.outcome,
    se.outcome,
    effect_allele.outcome,
    other_allele.outcome,
    eaf.outcome,
    pval.outcome,
    samplesize.outcome,
    outcome,
    id.outcome
   )])
   
   mr_data <- harmonise_data(exp_df, out_df, action = 3)
   if (nrow(mr_data) == 0) next
   
   n_harmonised <- nrow(mr_data)
   mr_data <- mr_data[mr_data$mr_keep, ]
   n_mr_keep <- nrow(mr_data)
   
   if (n_mr_keep == 0) next
   
   mr_res <- mr(mr_data)
   if (nrow(mr_res) == 0) next
   
   mr_res_dt <- as.data.table(mr_res)
   
   # add useful metadata per exposure
   mr_res_dt[, gene := current_gene]
   mr_res_dt[, gene_symbol := current_symbol]
   mr_res_dt[, n_harmonised := n_harmonised]
   mr_res_dt[, n_mr_keep := n_mr_keep]
   
   # keep all MR methods in same output
   setcolorder(
    mr_res_dt,
    c(
     "exposure", "id.exposure", "gene", "gene_symbol",
     "outcome", "id.outcome", "method", "nsnp", "b", "se", "pval",
     "n_harmonised", "n_mr_keep"
    )
   )
   
   res_batch[[length(res_batch) + 1]] <- mr_res_dt
  }, error = function(e) {
   message("Error for exposure ", current_id, ": ", e$message)
  })
 }
 
 if (length(res_batch) > 0) {
  res_batch_dt <- rbindlist(res_batch, use.names = TRUE, fill = TRUE)
 } else {
  res_batch_dt <- data.table(
   exposure = character(),
   id.exposure = character(),
   gene = character(),
   gene_symbol = character(),
   outcome = character(),
   id.outcome = character(),
   method = character(),
   nsnp = integer(),
   b = numeric(),
   se = numeric(),
   pval = numeric(),
   n_harmonised = integer(),
   n_mr_keep = integer()
  )
 }
 
 timestamp <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
 batch_filename <- paste0(
  "2SMR_eqtlgen_vs_", outcome_phenotype,
  "_all_methods_batch", batch_number, "_", timestamp, ".txt.gz"
 )
 
 fwrite(
  res_batch_dt,
  file = file.path(output_dir, batch_filename),
  sep = "\t"
 )
 
 all_files <- c(all_files, file.path(output_dir, batch_filename))
 batch_number <- batch_number + 1
 ind <- batch_end + 1
}

# =====================================
# Merge all batch files
# =====================================
final_output_file <- paste0(
 "2SMR_eqtlgen_vs_", outcome_phenotype,
 "_all_methods_final_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".txt.gz"
)

system(
 paste(
  "zcat", paste(all_files, collapse = " "),
  "| awk 'NR==1 || FNR>1'",
  "| gzip >",
  file.path(output_dir, final_output_file)
 )
)

file.remove(all_files)

message("Final merged file saved as: ", file.path(output_dir, final_output_file))