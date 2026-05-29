#!/bin/bash
# =============================================================================
# UKBB REGENIE Association Analysis — Neurodegenerative Disease STR Analysis
# =============================================================================
# Description : Runs REGENIE (step 1 + step 2) for tandem repeat (TR)
#               expansion burden association testing in UK Biobank WGS data.
#               Binary trait logistic regression with Firth correction.
#               Designed to be submitted per phenotype and expansion cutoff. 
#               Similar script used to run REGENIE in AoU with different covariates
#
# Usage       : bsub [options] bash UKBB.runREGENIE.neuroDegen.sh \
#                   <PHENO> <THRESHOLD>
#
# Arguments   :
#   PHENO       Phenotype name matching phenotype file and column
#               e.g. neuroDegen | neuroDegen_noPDorAD
#   THRESHOLD   Expansion cutoff label matching input file names
#               e.g. Cutoff95 | Cutoff98 | Cutoff99 | Cutoff995 |
#                    Cutoff998 | Cutoff999 | Cutoff9995
#
# Example HPC submit:
#   bsub -P acc_PROJECTID -L /bin/bash -q premium -n 18 \
#        -R rusage[mem=10000] -R span[hosts=1] -W 24:00 \
#        bash UKBB.runREGENIE.neuroDegen.sh neuroDegen Cutoff99
#
# Dependencies: regenie/3.4.1 (loaded via environment modules)
#
# Author      : Gabrielle Altman
# Date        : July 2025
# =============================================================================

set -euo pipefail

# =============================================================================
# Argument parsing and validation
# =============================================================================

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

# =============================================================================
# Environment
# =============================================================================

ml purge
ml regenie/3.4.1

# =============================================================================
# Configuration — edit these paths to adapt the script to your environment
# =============================================================================

# Root directory containing WGS genotype data (plink2 .pgen/.psam/.pvar files)
INDIR=/path/to/genotype_data

# Output directory for REGENIE results
OUTDIR=/path/to/output/${PHENO}_${THRESHOLD}

# Plink2 binary TR genotype file prefix (no extension)
# Expected: ${INDIR}/.../${THRESHOLD}/...BinVar.{pgen,psam,pvar}
PGEN=${INDIR}/BinaryTR_Plink/${THRESHOLD}/Autosomes_${THRESHOLD}_BinaryGTForPhewas.BinVar

# Phenotype file: tab-delimited, header row, columns FID/IID + phenotype column
PHENO_FILE=/path/to/phenotypes/UKBB.phenos.${PHENO}.txt

# Sample keep file: two-column FID/IID list of samples to include
SAMPLE_FILE=/path/to/phenotypes/UKBB.samples.${PHENO}.txt

# Covariate file: tab-delimited, header row, columns FID/IID + covariate columns
COVAR_FILE=${INDIR}/UKB_EUR_Covar.tsv

# TR variant ID list to extract (one variant ID per line)
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
