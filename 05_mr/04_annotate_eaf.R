library(data.table)

# =========================
# Paths
# =========================

work_dir <- "path/to/eaf_annotation_work_directory"

raw_dir <- file.path(work_dir, "plink_exportA")

# =========================
# Load PLINK dosage files
# =========================

raw_files <- list.files(
 raw_dir,
 pattern = "\\.raw$",
 full.names = TRUE
)

cat("Number of .raw files found:", length(raw_files), "\n")

# =========================
# Calculate EAF
# =========================

get_eaf_from_raw <- function(raw_file) {
 
 dt <- fread(raw_file)
 
 # First 6 columns are pedigree/sample metadata
 if (ncol(dt) <= 6) return(NULL)
 
 geno_cols <- names(dt)[7:ncol(dt)]
 
 eaf_vec <- sapply(geno_cols, function(x) {
  mean(dt[[x]], na.rm = TRUE) / 2
 })
 
 out <- data.table(
  raw_col = geno_cols,
  eaf.exposure = as.numeric(eaf_vec)
 )
 
 # Convert e.g. rs123_A -> rs123
 out[, SNP := sub("_[^_]+$", "", raw_col)]
 
 out[, .(SNP, eaf.exposure)]
}

# =========================
# Combine chromosomes
# =========================

eaf_list <- lapply(
 raw_files,
 get_eaf_from_raw
)

eaf_dt <- rbindlist(
 eaf_list,
 use.names = TRUE,
 fill = TRUE
)

eaf_dt <- unique(
 eaf_dt,
 by = "SNP"
)

cat("Annotated SNPs with EAF:", nrow(eaf_dt), "\n")

# =========================
# Save
# =========================

fwrite(
 eaf_dt,
 file.path(
  work_dir,
  "annotated_eaf_from_1000G_EUR.txt"
 ),
 sep = "\t"
)