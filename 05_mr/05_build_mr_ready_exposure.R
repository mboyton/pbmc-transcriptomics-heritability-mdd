library(data.table)

work_dir <- "path/to/eaf_annotation_work_directory"

exp <- fread(file.path(work_dir, "eqtlgen_input_cleaned.txt"))
eaf <- fread(file.path(work_dir, "annotated_eaf_from_1000G_EUR.txt"))

cat("Exposure rows:", nrow(exp), "\n")
cat("EAF rows:", nrow(eaf), "\n")

dat <- merge(exp, eaf, by = "SNP", all.x = TRUE)

cat("Rows after merge:", nrow(dat), "\n")
cat("Rows with missing EAF:", sum(is.na(dat$eaf.exposure)), "\n")

# keep only rows with valid EAF
dat <- dat[!is.na(eaf.exposure) & eaf.exposure > 0 & eaf.exposure < 1]

cat("Rows with valid EAF:", nrow(dat), "\n")

# derive beta and se
denom <- sqrt(2 * dat$eaf.exposure * (1 - dat$eaf.exposure) * (dat$NrSamples + dat$Zscore^2))
dat[, beta.exposure := Zscore / denom]
dat[, se.exposure   := 1 / denom]

# standardised MR columns
dat[, effect_allele.exposure := AssessedAllele]
dat[, other_allele.exposure  := OtherAllele]
dat[, pval.exposure := Pvalue]
dat[, samplesize.exposure := NrSamples]
dat[, chr.exposure := SNPChr]
dat[, pos.exposure := SNPPos]
dat[, exposure := paste0("eQTLGen_", GeneSymbol)]
dat[, id.exposure := paste0("eQTLGen_", Gene)]

# useful QC flag
dat[, palindromic := (
 (effect_allele.exposure == "A" & other_allele.exposure == "T") |
  (effect_allele.exposure == "T" & other_allele.exposure == "A") |
  (effect_allele.exposure == "C" & other_allele.exposure == "G") |
  (effect_allele.exposure == "G" & other_allele.exposure == "C")
)]

mr_ready <- dat[, .(
 SNP,
 chr.exposure,
 pos.exposure,
 effect_allele.exposure,
 other_allele.exposure,
 eaf.exposure,
 beta.exposure,
 se.exposure,
 pval.exposure,
 samplesize.exposure,
 exposure,
 id.exposure,
 Gene,
 GeneSymbol,
 Zscore,
 NrCohorts,
 FDR,
 BonferroniP,
 clump_gene,
 clump_gene_symbol,
 palindromic
)]

fwrite(
 mr_ready,
 file.path(
  work_dir,
  "eQTLGen_cis_eqtl_clumped_MR_ready.txt"
 ),
 sep = "\t"
)

cat("Final MR-ready rows:", nrow(mr_ready), "\n")