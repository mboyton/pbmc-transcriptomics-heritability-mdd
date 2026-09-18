PLINK2=plink2

REF_DIR="path/to/1000G_EUR_Phase3_plink"

WORK_DIR="path/to/eaf_annotation_work_directory"

# =========================
# Setup
# =========================

mkdir -p ${WORK_DIR}/plink_exportA
mkdir -p ${WORK_DIR}/logs

# =========================
# Extract allele dosages
# =========================

for chr in $(seq 1 22); do

  BFILE=${REF_DIR}/1000G.EUR.QC.${chr}
  EXTRACT=${WORK_DIR}/per_chr/chr${chr}_extract.txt
  ALLELES=${WORK_DIR}/per_chr/chr${chr}_count_allele.txt
  OUT=${WORK_DIR}/plink_exportA/chr${chr}

  if [ -s "${EXTRACT}" ]; then

    echo "Running chr${chr}"

    ${PLINK2} \
      --bfile ${BFILE} \
      --extract ${EXTRACT} \
      --export A \
      --export-allele ${ALLELES} \
      --out ${OUT}

  else

    echo "Skipping chr${chr}"

  fi

done