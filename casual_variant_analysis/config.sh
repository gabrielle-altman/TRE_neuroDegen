#!/usr/bin/env bash
# config.sh — paths, binaries, and covariates for the TR causal variant pipeline
# source at the top of every %%bash cell: source ./config.sh && source ./common.sh

set -euo pipefail

export WORKSPACE_ROOT="/home/jupyter/workspaces/yourworkspacename"  # EDIT THIS
export PROJECT_ROOT="${WORKSPACE_ROOT}/TR_Causal_Variant_Analysis"
export SNP_DIR="${PROJECT_ROOT}/SNP"
export TR_PLINK_DIR="${PROJECT_ROOT}/TR_plink"
export TRAIT_DIR="${PROJECT_ROOT}/tr_traits"
export REGENIE_DIR="${PROJECT_ROOT}/REGENIE"
export CAVIAR_DIR="${PROJECT_ROOT}/Caviar_FineMapping"
export LD_DIR="${PROJECT_ROOT}/LD_Matrix"
export COND_DIR="${PROJECT_ROOT}/Conditional_Analysis"
export COND_STEP1_DIR="${COND_DIR}/REGENIE_Step1_GlobalModel"
export LOG_DIR="${PROJECT_ROOT}/logs"

export TR_TRAIT_PAIRS_RAW="${PROJECT_ROOT}/TR_Trait_pairs.txt"
export TR_TRAIT_PAIRS="${PROJECT_ROOT}/TR_Trait_pairs_filtered.txt"  # built by step 00
export MISSING_LIST="${PROJECT_ROOT}/missing_regions_list.txt"

export REGIONS_BED="${PROJECT_ROOT}/TRListToUse_250kFlanks.bed"
export SAMPLES_FILE="${PROJECT_ROOT}/samplesToUse.txt"
export TR_BINARY_INDIR="${WORKSPACE_ROOT}/R7_BinaryTR_Plink"
export UPSTREAM_REGENIE_DIR="${WORKSPACE_ROOT}/regenie/AllAncestries"

# forPlink variant has #FID/IID columns; bare variant is REGENIE-formatted
sample_list_for_regenie() { echo "${WORKSPACE_ROOT}/AoU.ALL.sampleList.${1}.tsv"; }
sample_list_for_plink()   { echo "${WORKSPACE_ROOT}/AoU.ALL.sampleList.${1}.forPlink.tsv"; }
phenotype_file()          { echo "${WORKSPACE_ROOT}/AoU_ALL_phenotype.${1}.tsv"; }
export -f sample_list_for_regenie sample_list_for_plink phenotype_file

# plink2 and plink (1.9) are on $PATH in AoU images; regenie and CAVIAR are local builds
export REGENIE_BIN="${WORKSPACE_ROOT}/regenie_3.4.1"
export CAVIAR_BIN="${WORKSPACE_ROOT}/caviar/CAVIAR-C++/CAVIAR"

export COVAR_FILE="${AOU_COVAR_FILE:-/home/jupyter/workspaces/REPLACEME/AoU_ALL_samples_covariates_neuroDegen.tsv.gz}"
[[ "$COVAR_FILE" == *REPLACEME* ]] && { echo "ERROR: set AOU_COVAR_FILE or edit COVAR_FILE in config.sh" >&2; exit 1; }
export QUANT_COVARS="SNP_PC1,SNP_PC2,SNP_PC3,SNP_PC4,SNP_PC5,Age,Insert_Size"
export CAT_COVARS="Gender,Cohort,ancestry_pred_other"

export THREADS=16
export STEP1_PHENOS=(neuroDegen neuroDegen_noPDorAD)
export TYPE="Binary"

mkdir -p "$SNP_DIR" "$TR_PLINK_DIR" "$TRAIT_DIR/$TYPE" \
         "$REGENIE_DIR/$TYPE" "$CAVIAR_DIR" "$LD_DIR" \
         "$COND_DIR" "$COND_STEP1_DIR" "$LOG_DIR"
