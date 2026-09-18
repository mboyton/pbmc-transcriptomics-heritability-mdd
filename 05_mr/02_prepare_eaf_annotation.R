library(data.table)

# =========================
# Paths
# =========================

exp_file <- "path/to/eQTLGen_cis_eqtl_clumped_all_genes.txt"

work_dir <- "path/to/eaf_annotation_work_directory"

# =========================
# Setup
# =========================

dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(work_dir, "per_chr"), showWarnings = FALSE)
dir.create(file.path(work_dir, "logs"), showWarnings = FALSE)

# =========================
# Load exposure data
# =========================

exp <- fread(exp_file)

cat("Rows in input:", nrow(exp), "\n")

required_cols <- c(
 "Pvalue", "SNP", "SNPChr", "SNPPos",
 "AssessedAllele", "OtherAllele",
 "Zscore", "Gene", "GeneSymbol",
 "NrSamples"
)

missing_cols <- setdiff(required_cols, names(exp))

if (length(missing_cols) > 0) {
 stop(
  "Missing required columns: ",
  paste(missing_cols, collapse = ", ")
 )
}

# =========================
# Basic variant cleaning
# =========================

# Keep clean autosomal biallelic SNPs
exp_clean <- exp[
 SNPChr %in% 1:22 &
  AssessedAllele %in% c("A", "C", "G", "T") &
  OtherAllele %in% c("A", "C", "G", "T") &
  AssessedAllele != OtherAllele
]

cat("Rows after basic cleaning:", nrow(exp_clean), "\n")

# Save cleaned copy for downstream merge
fwrite(
 exp_clean,
 file.path(work_dir, "eqtlgen_input_cleaned.txt"),
 sep = "\t"
)

# =========================
# Prepare EAF annotation files
# =========================

# Unique SNP-level request table for EAF annotation
req <- unique(
 exp_clean[, .(
  SNP,
  SNPChr,
  AssessedAllele
 )]
)

cat("Unique SNPs for frequency annotation:", nrow(req), "\n")

for (chr in 1:22) {
 
 chr_req <- req[SNPChr == chr]
 
 extract_file <- file.path(
  work_dir,
  "per_chr",
  paste0("chr", chr, "_extract.txt")
 )
 
 allele_file <- file.path(
  work_dir,
  "per_chr",
  paste0("chr", chr, "_count_allele.txt")
 )
 
 if (nrow(chr_req) > 0) {
  
  fwrite(
   chr_req[, .(SNP)],
   extract_file,
   sep = "\t",
   col.names = FALSE
  )
  
  fwrite(
   chr_req[, .(SNP, AssessedAllele)],
   allele_file,
   sep = "\t",
   col.names = FALSE
  )
  
 } else {
  
  file.create(extract_file)
  file.create(allele_file)
  
 }
}

# =========================
# QC summary
# =========================

chr_counts <- req[
 , .N,
 by = SNPChr
][order(SNPChr)]

fwrite(
 chr_counts,
 file.path(work_dir, "per_chr_snp_counts.txt"),
 sep = "\t"
)

cat("Preparation complete.\n")