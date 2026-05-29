#!/bin/bash
# =============================================================================
# UKBB REGENIE Association Analysis — Neurodegenerative STRs
# =============================================================================
# Gabrielle Altman
#
# Runs REGENIE step 1 + step 2 for TR expansion burden in UKBB WGS TR genotype data.
# Binary trait, Firth approx correction. One submission per pheno x cutoff.
# Similar script used in AoU with different covariates.
#
# Usage: bash runREGENIE.neuroDegen.sh <PHENO> <THRESHOLD>
#   PHENO      e.g. neuroDegen | neuroDegen_noPDorAD
#   THRESHOLD  e.g. Cutoff95 | Cutoff99 | Cutoff9995
#
# HPC: bsub -P acc_PROJECTID -L /bin/bash -q premium -n 18 \
#           -R rusage[mem=10000] -R span[hosts=1] -W 24:00 \
#           bash runREGENIE.neuroDegen.sh neuroDegen Cutoff99
#
# Dependencies: regenie/3.4.1 (loaded via environment modules)
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
TYPE=Binary

echo "======================================================"
echo "REGENIE run started: $(date)"
echo "  Phenotype : ${PHENO}"
echo "  Threshold : ${THRESHOLD}"
echo "  Type      : ${TYPE}"
echo "======================================================"

ml purge
ml regenie/3.4.1

# paths — edit for your environment
INDIR=/path/to/genotype_data
OUTDIR=/path/to/output/${PHENO}_${THRESHOLD}
PGEN=${INDIR}/BinaryTR_Plink/${THRESHOLD}/Autosomes_${THRESHOLD}_BinaryGTForPhewas.BinVar
PHENO_FILE=/path/to/phenotypes/UKBB.phenos.${PHENO}.txt
SAMPLE_FILE=/path/to/phenotypes/UKBB.samples.${PHENO}.txt
COVAR_FILE=${INDIR}/UKB_EUR_Covar.tsv
TR_LIST=/path/to/variant_lists/UKBB.keep.${THRESHOLD}.${PHENO}.txt

COVAR_COLS=Insert_Size,Age,Age_sq,SNP_PC1,SNP_PC2,SNP_PC3,SNP_PC4,SNP_PC5
CAT_COVAR_COLS=Gender,SeqCenter

STEP1_PREFIX=${OUTDIR}/${PHENO}_${TYPE}_${THRESHOLD}_AllTRs.step1
STEP2_PREFIX=${OUTDIR}/${PHENO}.${TYPE}.AllSamples.${THRESHOLD}.AllTRs
LOWMEM_TMP1=${OUTDIR}/tmp_rg1_${PHENO}_${THRESHOLD}
LOWMEM_TMP2=${OUTDIR}/tmp_rg2_${PHENO}_${THRESHOLD}

mkdir -p "${OUTDIR}"

# Flags shared by both REGENIE steps
COMMON_FLAGS=(
    --pgen            "${PGEN}"
    --keep            "${SAMPLE_FILE}"
    --extract         "${TR_LIST}"
    --phenoFile       "${PHENO_FILE}"
    --phenoCol        "${PHENO}"
    --covarFile       "${COVAR_FILE}"
    --covarColList    "${COVAR_COLS}"
    --catCovarList    "${CAT_COVAR_COLS}"
    --bt
    --firth --approx
    --strict
    --bsize           1000
    --loocv
    --lowmem
    --threads         18
    --write-samples
)

# =============================================================================
# Step 1 — Whole-genome regression (null model)
# =============================================================================
echo ""
echo "Running step 1 ... $(date)"

regenie \
    --step 1 \
    "${COMMON_FLAGS[@]}" \
    --lowmem-prefix "${LOWMEM_TMP1}" \
    --out            "${STEP1_PREFIX}"

echo "Step 1 complete: $(date)"

# =============================================================================
# Step 2 — Association testing
# =============================================================================

echo ""
echo "Running step 2 ... $(date)"

regenie \
    --step 2 \
    "${COMMON_FLAGS[@]}" \
    --lowmem-prefix "${LOWMEM_TMP2}" \
    --pred          "${STEP1_PREFIX}_pred.list" \
    --pThresh       0.05 \
    --minMAC        2 \
    --out           "${STEP2_PREFIX}"

echo "Step 2 complete: $(date)"
echo ""
echo "======================================================"
echo "REGENIE run finished: $(date)"
echo "  Output prefix: ${STEP2_PREFIX}"
echo "======================================================"
