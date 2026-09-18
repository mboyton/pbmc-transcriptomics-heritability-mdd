#!/bin/bash

set -euo pipefail

########################################
# USER SETTINGS
########################################

# CELLECT repository root
CELLECT_DIR="/path/to/CELLECT"

# CELLEX expression-specificity file
SPEC_RAW="/path/to/cellex_output.esmu.csv.gz"

# Pre-standardised GWAS input
# Expected first eight columns:
# SNP A1 A2 BETA SE P N Z
GWAS_IN="/path/to/gwas_cellect_ready.tsv"

# Short label for this CELLECT run
RUN_TAG="example_run"

# Fixed global sample size used during LDSC munging
FIXED_N=""

# Number of Snakemake cores
SMK_CORES=4

########################################
# TOOL PATHS
########################################

SNAKEMAKE_BIN="/path/to/snakemake"
CONDA_SH="/path/to/conda.sh"
MAMBA_BIN="/path/to/mamba"
CONDA_PREFIX_DIR="/path/to/snakemake_envs"
MUNGE_ENV="/path/to/munge_ldsc"

########################################
# DERIVED PATHS
########################################

BASE_OUT="$CELLECT_DIR/CELLECT-${RUN_TAG}"
LDSC_OUT="$BASE_OUT/CELLECT-LDSC"

CUSTOM_IN="$CELLECT_DIR/custom_inputs/${RUN_TAG}"

mkdir -p "$CUSTOM_IN" "$BASE_OUT"

MUNGED_PREFIX="$CUSTOM_IN/${RUN_TAG}"
CONFIG_YAML="$CELLECT_DIR/config_${RUN_TAG}.yml"

TS=$(date +%Y%m%d_%H%M%S)
LOG_FILE="$CELLECT_DIR/cellect_ldsc_run_${RUN_TAG}_${TS}.log"

########################################
# SANITY CHECKS
########################################

[ -d "$CELLECT_DIR" ] || {
  echo "ERR: CELLECT_DIR not found: $CELLECT_DIR"
  exit 1
}

[ -f "$SPEC_RAW" ] || {
  echo "ERR: CELLEX specificity file not found: $SPEC_RAW"
  exit 1
}

[ -x "$SNAKEMAKE_BIN" ] || {
  echo "ERR: Snakemake not found at: $SNAKEMAKE_BIN"
  exit 1
}

[ -f "$CONDA_SH" ] || {
  echo "ERR: conda.sh not found at: $CONDA_SH"
  exit 1
}

[ -x "$MAMBA_BIN" ] || {
  echo "ERR: mamba not found at: $MAMBA_BIN"
  exit 1
}

[ -f "$GWAS_IN" ] || {
  echo "ERR: GWAS input not found: $GWAS_IN"
  exit 1
}

########################################
# LOAD CONDA
########################################

source "$CONDA_SH"

########################################
# ENSURE LDSC MUNGE ENVIRONMENT EXISTS
########################################

if [ ! -d "$MUNGE_ENV" ]; then
  echo "[INFO] Creating LDSC munge environment: $MUNGE_ENV"

  conda env create \
    -f "$CELLECT_DIR/ldsc/environment_munge_ldsc.yml" \
    -p "$MUNGE_ENV"
fi

echo "[OK] Munge environment ready: $MUNGE_ENV"

########################################
# CHECK GWAS HEADER
########################################

echo "[INFO] Checking GWAS header: $GWAS_IN"

awk 'BEGIN{FS=OFS="\t"}
NR==1{
  expected = "SNP\tA1\tA2\tBETA\tSE\tP\tN\tZ"
  got = $1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7 FS $8

  if (got != expected) {
    print "ERR: Header mismatch." > "/dev/stderr"
    print "ERR: Expected: " expected > "/dev/stderr"
    print "ERR: Got:      " got > "/dev/stderr"
    exit 1
  }

  exit 0
}' "$GWAS_IN"

echo "[OK] GWAS header looks good."

GWAS_FOR_MUNGE="$GWAS_IN"

########################################
# FIXED SAMPLE SIZE FOR LDSC
########################################

if [ -z "$FIXED_N" ]; then
  echo "ERR: FIXED_N is empty but a global N is required."
  exit 1
fi

echo "[OK] Using fixed global N = $FIXED_N for LDSC"

########################################
# LDSC MUNGING
########################################

echo "[INFO] Munging GWAS -> ${MUNGED_PREFIX}.sumstats.gz"

conda run -p "$MUNGE_ENV" \
  python "$CELLECT_DIR/ldsc/mtag_munge.py" \
  --sumstats "$GWAS_FOR_MUNGE" \
  --snp SNP \
  --a1 A1 \
  --a2 A2 \
  --p P \
  --merge-alleles "$CELLECT_DIR/data/ldsc/w_hm3.snplist" \
  --keep-pval \
  --out "$MUNGED_PREFIX" \
  --n-value "$FIXED_N" \
  --signed-sumstats Z,0

[ -s "${MUNGED_PREFIX}.sumstats.gz" ] || {
  echo "ERR: munged summary statistics were not produced."
  exit 1
}

echo "[OK] Munged summary statistics: ${MUNGED_PREFIX}.sumstats.gz"

########################################
# WRITE CELLECT CONFIG
########################################

echo "[INFO] Writing CELLECT configuration: $CONFIG_YAML"

cat > "$CONFIG_YAML" <<YML
---
BASE_OUTPUT_DIR: ${BASE_OUT}

SPECIFICITY_INPUT:
  - id: cellect_output
    path: ${SPEC_RAW}

GWAS_SUMSTATS:
  - id: ${RUN_TAG}
    path: ${MUNGED_PREFIX}.sumstats.gz

ANALYSIS_TYPE:
  prioritization: True
  conditional: False
  heritability: False
  heritability_intervals: False

WINDOW_DEFINITION:
  WINDOW_SIZE_KB:
    100

GENE_COORD_FILE:
  data/shared/gene_coordinates.GRCh37.ensembl_v91.txt

KEEP_ANNOTS:
  False

LDSC_CONST:
  DATA_DIR:
    data/ldsc
  LDSC_DIR:
    ldsc
  NUMPY_CORES:
    1
YML

########################################
# RUN CELLECT-LDSC
########################################

cd "$CELLECT_DIR"

"$SNAKEMAKE_BIN" \
  --unlock \
  -s cellect-ldsc.snakefile \
  --configfile "$CONFIG_YAML" || true

rm -rf ".snakemake/locks"

echo "[INFO] Starting CELLECT-LDSC."

"$SNAKEMAKE_BIN" \
  --use-conda \
  --conda-prefix "$CONDA_PREFIX_DIR" \
  --conda-frontend mamba \
  -j "$SMK_CORES" \
  -s cellect-ldsc.snakefile \
  --configfile "$CONFIG_YAML" \
  --rerun-incomplete \
  > "$LOG_FILE" 2>&1

echo
echo "=== CELLECT-LDSC completed ==="
echo "Log file: $LOG_FILE"
echo
echo "Key outputs:"
echo "  Prioritization CSV:"
echo "    $LDSC_OUT/results/prioritization.csv"
echo
echo "  Per-GWAS cell-type results:"
echo "    $LDSC_OUT/out/prioritization/cellect_output__${RUN_TAG}.cell_type_results.txt"