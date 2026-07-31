#!/bin/bash
# =============================================================================
# METAL Meta-Analysis — Neurodegenerative STRs
# =============================================================================
# Gabrielle N. Altman
#
# Fixed-effects meta-analysis (METAL) of TR expansion association results
# (REGENIE step 2 output) across two cohorts, e.g. AoU + UKB, for one
# phenotype x expansion cutoff combination.
#
# Usage: bash runMETAL.neuroDegen.sh <PHENO> <THRESHOLD>
#   PHENO      e.g. neuroDegen | neuroDegen_noPDorAD
#   THRESHOLD  e.g. Cutoff95 | Cutoff99 | Cutoff9995
#
# HPC: bsub -P acc_PROJECTID -L /bin/bash -q express -n 1 \
#           -R rusage[mem=8000] -W 4:00 \
#           bash runMETAL.neuroDegen.sh neuroDegen Cutoff99
#
# Dependencies: metal (loaded via environment modules)
# =============================================================================

set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "ERROR: Missing required arguments." >&2
  echo "Usage: bash $(basename "$0") <PHENO> <THRESHOLD>" >&2
  echo "  PHENO      e.g. neuroDegen | neuroDegen_noPDorAD" >&2
  echo "  THRESHOLD  e.g. Cutoff95 | Cutoff99 | Cutoff9995" >&2
  exit 1
fi

PHENO=$1
THRESHOLD=$2

echo "======================================================"
echo "METAL run started: $(date)"
echo "  Phenotype : ${PHENO}"
echo "  Threshold : ${THRESHOLD}"
echo "======================================================"

ml purge
ml metal/2018-08-28

# paths — edit for your environment
OUTDIR=/path/to/output/${PHENO}
COHORT1_FILE=/path/to/cohort1/regenie/${PHENO}.${THRESHOLD}.regenie.txt
COHORT2_FILE=/path/to/cohort2/regenie/${PHENO}.${THRESHOLD}.regenie.txt

# P-value column name in each cohort's REGENIE output
COHORT1_PVAL_COL=P
COHORT2_PVAL_COL=PVAL

mkdir -p "${OUTDIR}"

METAL_SCRIPT="${OUTDIR}/METAL.${THRESHOLD}.${PHENO}.txt"
OUT_PREFIX="${OUTDIR}/METAANALYSIS.${THRESHOLD}.${PHENO}."

cat > "${METAL_SCRIPT}" <<EOF
MARKER ID
WEIGHT N
ALLELE ALLELE0 ALLELE1
EFFECT BETA
PVAL ${COHORT1_PVAL_COL}
PROCESS ${COHORT1_FILE}

MARKER ID
WEIGHT N
ALLELE ALLELE0 ALLELE1
EFFECT BETA
PVAL ${COHORT2_PVAL_COL}
PROCESS ${COHORT2_FILE}

OUTFILE ${OUT_PREFIX} .txt
ANALYZE
EOF

metal < "${METAL_SCRIPT}"

echo ""
echo "======================================================"
echo "METAL run finished: $(date)"
echo "  Output prefix: ${OUT_PREFIX}"
echo "======================================================"
