module load languages/R/4.5.1

# paths
RAW_GWAS="path/to/psoriasis_gwas.tsv.gz"
DBSNP_REF="path/to/dbsnp_GRCh37_master_table.tsv.gz"
OUT_TSV="path/to/Dand_Psoriasis_2025_GRCh37_ready.tsv"

FIXED_N=56000

# prepare GWAS
Rscript - <<EOF
library(data.table)

raw_gwas  <- "${RAW_GWAS}"
dbsnp_ref <- "${DBSNP_REF}"
out_tsv   <- "${OUT_TSV}"
fixed_N   <- as.numeric("${FIXED_N}")

gwas <- fread(
 cmd = paste(
  "zcat",
  shQuote(raw_gwas)
 )
)

# standardise column names
orig_names <- names(gwas)
clean_names <- tolower(orig_names)

setnames(
 gwas,
 old = orig_names,
 new = clean_names
)

required_cols <- c(
 "chromosome",
 "base_pair_location",
 "effect_allele",
 "other_allele",
 "beta",
 "standard_error",
 "p_value",
 "cum_eff_sample_size"
)

missing <- setdiff(
 required_cols,
 names(gwas)
)

if (length(missing) > 0) {
 stop(
  "Missing expected columns in GWAS: ",
  paste(missing, collapse = ", ")
 )
}

# rename columns
setnames(
 gwas,
 c(
  "chromosome",
  "base_pair_location",
  "effect_allele",
  "other_allele",
  "beta",
  "standard_error",
  "p_value"
 ),
 c(
  "CHR_raw",
  "BP_raw",
  "A1_raw",
  "A2_raw",
  "BETA_raw",
  "SE_raw",
  "P_raw"
 )
)

# format GWAS
gwas[, CHR := gsub(
 "^chr",
 "",
 as.character(CHR_raw)
)]

gwas[, BP := as.integer(BP_raw)]

gwas[, A1 := toupper(A1_raw)]
gwas[, A2 := toupper(A2_raw)]

gwas[, BETA := as.numeric(BETA_raw)]
gwas[, SE := as.numeric(SE_raw)]
gwas[, P := as.numeric(P_raw)]


# read dbSNP reference
db <- fread(
 cmd = paste(
  "zcat",
  shQuote(dbsnp_ref)
 ),
 header = FALSE
)

setnames(
 db,
 c(
  "CHR",
  "BP",
  "REF",
  "ALT",
  "rsid"
 )
)

db[, CHR := as.character(CHR)]
db[, BP := as.integer(BP)]

# retain simple SNPs
db <- db[
 nchar(REF) == 1 &
 nchar(ALT) == 1
]


# annotate rsIDs
gwas[, CHR := as.character(CHR)]
gwas[, BP := as.integer(BP)]

setkey(
 db,
 CHR,
 BP
)

setkey(
 gwas,
 CHR,
 BP
)

g_annot <- db[gwas]


# construct SNP identifier
g_annot[, SNP := ifelse(
 !is.na(rsid),
 rsid,
 paste(
  CHR,
  BP,
  A1,
  A2,
  sep = ":"
 )
)]


# sample size and Z-score
g_annot[, N := fixed_N]
g_annot[, Z := BETA / SE]


# final GWAS
out <- g_annot[, .(
 SNP,
 A1,
 A2,
 BETA,
 SE,
 P,
 N,
 Z,
 CHR,
 BP,
 rsid,
 cum_eff_sample_size
)]

fwrite(
 out,
 out_tsv,
 sep = "\t",
 quote = FALSE,
 na = "NA"
)

EOF

# compress output
gzip -f "${OUT_TSV}"