#!/bin/bash
# =============================================================================
# ExpansionHunter v5.0.0 — STR Genotyping
# =============================================================================
# Gabrielle N. Altman
#
# Runs ExpansionHunter v5.0.0 on a single sample CRAM file to genotype
# short tandem repeat loci. Run per sample on DNAnexus (UKB) and the
# All of Us Researcher Workbench (AoU).
#
# Usage: bash run_expansionhunter.sh <INPUT_CRAM> <INPUT_RefFasta> <INPUT_STRs> <gender> <sampleName>
#   INPUT_CRAM      Path to input CRAM file
#   INPUT_RefFasta  Path to GRCh38 reference genome FASTA
#   INPUT_STRs      Path to ExpansionHunter variant catalog (JSON)
#                   Catalog used in this study: ExpansionHunter_TRCatalog_NeuroDegen_16871loci_hg38.json (this repository)
#   gender          Sample sex (male/female), used for sex chromosome genotyping
#   sampleName      Sample identifier used for output file naming
#
# Dependencies: ExpansionHunter/5.0.0
# =============================================================================

ExpansionHunter \
    --reads ${INPUT_CRAM} \
    --reference ${INPUT_RefFasta} \
    --variant-catalog ${INPUT_STRs} \
    --sex ${gender} \
    --analysis-mode streaming \
    -n 16 \
    --out ${sampleName}.EHv5.0.0