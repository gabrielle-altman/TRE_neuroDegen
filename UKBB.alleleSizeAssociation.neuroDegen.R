# =============================================================================
# UKBB Allele Size Association — Neurodegenerative Disease STR Analysis
# =============================================================================
# Description : For each tandem repeat (TR) locus that was nominally
#               significant in the AoU+UKBB METAL meta-analysis, tests
#               association between long-allele size (as a binary variable
#               at each repeat-length cutoff >= the 99th percentile) and
#               neurodegenerative disease status using Firth logistic
#               regression (brglm). Run per chromosome.
#
# Usage       : Rscript UKBB.alleleSizeAssociation.neuroDegen.R \
#                   --dir  <main_dir> \
#                   --cohort <phenotype> \
#                   --chrom  <chromosome>
#
# Arguments   :
#   --dir     Main project directory
#   --cohort  Phenotype label: neuroDegen | neuroDegen_noPDorAD
#   --chrom   Chromosome: chr1 ... chr22
#
# Example HPC submit:
#   bsub -P <project> -L /bin/bash -q express -n 4 \
#        -R span[hosts=1] -R rusage[mem=30000] -W 10:00 \
#        Rscript UKBB.alleleSizeAssociation.neuroDegen.R \
#            --dir /path/to/project/dir \
#            --cohort neuroDegen --chrom chr12
#
# Dependencies: tidyverse, data.table, scales, argparser, R.utils,
#               readxl, parallel, brglm
#
# Author      : Gabrielle Altman, adapted from code by Bharati Jadhav
# =============================================================================

cat(paste0("R version: ", getRversion(), "\n"))
cat("Loading libraries...\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(scales)
  library(argparser)
  library(R.utils)
  library(readxl)
  library(parallel)
  library(brglm)
})

# =============================================================================
# Argument parsing and validation
# =============================================================================

p <- arg_parser("Per-allele-size association test for TR loci")
p <- add_argument(p, "--dir",    help = "Main project directory [required]")
p <- add_argument(p, "--cohort", help = "neuroDegen | neuroDegen_noPDorAD [required]")
p <- add_argument(p, "--chrom",  help = "Chromosome, e.g. chr12 [required]")
argv <- parse_args(p)

if (is.na(argv$dir) || is.na(argv$cohort) || is.na(argv$chrom)) {
  print(p)
  stop("ERROR: --dir, --cohort, and --chrom are all required.", call. = FALSE)
}

cat("\nArguments:\n")
cat("  --dir    :", argv$dir,    "\n")
cat("  --cohort :", argv$cohort, "\n")
cat("  --chrom  :", argv$chrom,  "\n\n")


setwd(argv$dir)

# =============================================================================
# Configuration — edit these paths to adapt the script to your environment
# =============================================================================

EH_DIR    <- "/path/to/EH/genotypes/PerChrGT_FilteredEUR_AfterPCA"
COVAR_DIR <- "/path/to/covariates"
META_DIR  <- file.path(argv$dir, "metaAnalysis/neuroDegen/AoU_UKBB")
OUT_DIR   <- file.path(argv$dir, "alleleSizeAssociation")

PHENO_FILE  <- file.path(argv$dir, paste0("UKBB.phenos.", argv$cohort, ".txt"))
COVAR_FILE  <- file.path(COVAR_DIR, "UKB500k_EUR_UnrelSampleListToIncludeAfterClassifierQc_Covar.tsv")
METAL_FILE  <- file.path(META_DIR,  paste0("AoU.UKBB.", argv$cohort,
                           ".METAL.results.nominallySignificant.filtered.xlsx"))

LONGALLELE_SC     <- file.path(EH_DIR, "SC/LongAlleleMatrix",
                                paste0(argv$chrom, "_SC_EUR_LongAlleleMatrix.tsv.gz"))
LONGALLELE_DECODE <- file.path(EH_DIR, "deCODE/LongAlleleMatrix",
                                paste0(argv$chrom, "_deCODE_EUR_LongAlleleMatrix.tsv.gz"))

CLASSIFIER_SC     <- file.path(EH_DIR, "SC/ClassifierQC",
                                paste0(argv$chrom, "_SC_EUR_Cutoff95_ClassifierPredictions.tsv.gz"))
CLASSIFIER_DECODE <- file.path(EH_DIR, "deCODE/ClassifierQC",
                                paste0(argv$chrom, "_deCODE_EUR_Cutoff95_ClassifierPredictions.tsv.gz"))

N_CORES <- min(4L, parallel::detectCores() - 1L)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Load phenotypes, covariates, and target loci
# =============================================================================

cat("Loading phenotypes...\n")
phenotypes <- fread(PHENO_FILE, sep = "\t", check.names = FALSE, header = TRUE) %>%
  select(-FID) %>%
  setNames(c("SampleId", "Phenotype"))

cat("Loading covariates...\n")
covariates <- fread(COVAR_FILE, sep = "\t", check.names = FALSE, header = TRUE) %>%
  select(IID, SNP_PC1, SNP_PC2, SNP_PC3, SNP_PC4, SNP_PC5,
         Age, Age_sq, Insert_Size, predicted_gender, SeqCenter)

cat("Loading METAL target loci...\n")
target_loci <- read_excel(METAL_FILE) %>%
  as.data.frame() %>%
  select(Locus) %>%
  distinct()

# =============================================================================
# Load long-allele matrices and filter to target loci + phenotyped samples
# =============================================================================

cat("Loading long-allele matrices...\n")

load_longallele <- function(filepath) {
  fread(filepath, sep = "\t", check.names = FALSE, header = TRUE) %>%
    rename(SampleId = V1) %>%
    filter(SampleId %in% phenotypes$SampleId) %>%
    data.table::melt(id.vars = "SampleId", variable.name = "Locus", value.name = "LongAllele") %>%
    inner_join(phenotypes, by = "SampleId") %>%
    filter(Locus %in% target_loci$Locus) %>%
    inner_join(covariates, by = c("SampleId" = "IID"))
}

df_SC     <- load_longallele(LONGALLELE_SC)
df_deCODE <- load_longallele(LONGALLELE_DECODE)

# =============================================================================
# Load classifier predictions (Cutoff95) and filter out failed calls ("F")
# =============================================================================

cat("Loading classifier predictions...\n")

classifier_SC     <- fread(CLASSIFIER_SC,     sep = "\t", check.names = FALSE, header = TRUE) %>%
  select(Locus, SampleId, PredictedLabels)
classifier_deCODE <- fread(CLASSIFIER_DECODE, sep = "\t", check.names = FALSE, header = TRUE) %>%
  select(Locus, SampleId, PredictedLabels)

classifier <- rbind(classifier_SC, classifier_deCODE)

# Merge allele data with classifier; retain calls that passed QC ("P") or were not classified (NA)
eh <- rbind(df_SC, df_deCODE) %>%
  left_join(classifier, by = c("SampleId", "Locus")) %>%
  filter(PredictedLabels != "F" | is.na(PredictedLabels))

cat(paste0("Loci to test : ", length(unique(eh$Locus)),   "\n"))
cat(paste0("Samples      : ", length(unique(eh$SampleId)), "\n\n"))

rm(df_SC, df_deCODE, classifier_SC, classifier_deCODE, classifier)
gc()

eh$Locus <- as.character(eh$Locus)

# =============================================================================
# Functions
# =============================================================================

# For a given repeat-length threshold, categorize samples into expanded vs normal
make_allele_bins <- function(repeat_len, df) {
  df %>% mutate(
    RepeatStatus = as.integer(LongAllele >= repeat_len),
    RepeatSize   = repeat_len
  )
}

# Firth logistic regression for one repeat-length cutoff
run_binary_regression <- function(repeat_len, dt) {
  dt$Phenotype    <- relevel(factor(dt$Phenotype),    ref = "0")
  dt$RepeatStatus <- relevel(factor(dt$RepeatStatus), ref = "0")

  if (nlevels(dt$RepeatStatus) < 2 || nlevels(dt$Phenotype) < 2) {
    return(data.frame(repeatLen = repeat_len, LR_StdErr = NA_real_,
                      LR_Zvalue = NA_real_, LR_OddsRatio = NA_real_, LR_Pval = NA_real_))
  }

  fit <- try(
    brglm(Phenotype ~ RepeatStatus + SNP_PC1 + SNP_PC2 + SNP_PC3 + SNP_PC4 + SNP_PC5 +
            Age + Insert_Size + predicted_gender + SeqCenter,
          data = dt, family = binomial, pl = TRUE),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    return(data.frame(repeatLen = repeat_len, LR_StdErr = NA_real_,
                      LR_Zvalue = NA_real_, LR_OddsRatio = NA_real_, LR_Pval = NA_real_))
  }

  coef_row <- summary(fit)$coefficients["RepeatStatus1", ]
  data.frame(
    repeatLen    = repeat_len,
    LR_StdErr    = coef_row["Std. Error"],
    LR_Zvalue    = coef_row["z value"],
    LR_OddsRatio = exp(coef_row["Estimate"]),
    LR_Pval      = coef_row["Pr(>|z|)"]
  )
}

# Test all allele-size cutoffs >= 99th percentile for one locus
test_locus <- function(marker) {
  message("Testing locus: ", marker)
  df      <- eh[eh$Locus == marker, ]
  cutoff  <- quantile(df$LongAllele, probs = 0.99, na.rm = TRUE)
  alleles <- sort(unique(df$LongAllele[df$LongAllele >= cutoff]))

  if (length(alleles) == 0) {
    message("  No alleles above 99th-percentile cutoff — skipping.")
    return(NULL)
  }

  lmdf <- do.call(rbind, lapply(alleles, make_allele_bins, df = df))

  count_summary <- lmdf %>%
    group_by(RepeatSize) %>%
    summarise(
      CaseLargeRepeat    = sum(Phenotype == 1 & RepeatStatus == 1),
      ControlLargeRepeat = sum(Phenotype == 0 & RepeatStatus == 1),
      CaseSmallRepeat    = sum(Phenotype == 1 & RepeatStatus == 0),
      ControlSmallRepeat = sum(Phenotype == 0 & RepeatStatus == 0),
      .groups = "drop"
    )

  regression_results <- parallel::mclapply(alleles, function(len) {
    tryCatch(
      run_binary_regression(len, lmdf[lmdf$RepeatSize == len, ]),
      error = function(e) {
        message("  Error at cutoff ", len, ": ", e$message)
        NULL
      }
    )
  }, mc.cores = N_CORES)

  valid_results <- Filter(Negate(is.null), regression_results)
  if (length(valid_results) == 0) {
    message("  No valid regression results for locus ", marker, " — skipping.")
    return(NULL)
  }

  assoc <- do.call(rbind, valid_results) %>%
    rename(RepeatSize = repeatLen) %>%
    inner_join(count_summary, by = "RepeatSize") %>%
    filter(!is.infinite(LR_OddsRatio), LR_OddsRatio != 0) %>%
    mutate(Locus = marker)

  if (nrow(assoc) == 0) {
    message("  All results filtered out for locus ", marker, " — skipping.")
    return(NULL)
  }

  list(
    all  = assoc,
    best = assoc[which.min(assoc$LR_Pval), ]
  )
}

# =============================================================================
# Run association tests across all loci
# =============================================================================

loci <- unique(eh$Locus)
cat(paste0("Running association tests for ", length(loci), " loci...\n"))

res_list <- lapply(loci, function(loc) {
  tryCatch(
    test_locus(loc),
    error = function(e) {
      message("Error processing locus ", loc, ": ", e$message)
      NULL
    }
  )
})
res_list <- Filter(Negate(is.null), res_list)

all_table  <- do.call(rbind, lapply(res_list, `[[`, "all"))
best_table <- do.call(rbind, lapply(res_list, `[[`, "best"))

cat("Association testing complete.\n\n")

# =============================================================================
# Save results
# =============================================================================

all_file  <- file.path(OUT_DIR, paste0("UKBB.", argv$chrom, ".", argv$cohort,
                         ".alleleSizeAssociation.AllResults.txt.gz"))
best_file <- file.path(OUT_DIR, paste0("UKBB.", argv$chrom, ".", argv$cohort,
                         ".alleleSizeAssociation.BestResults.txt.gz"))

fwrite(all_table,  all_file,  sep = "\t", row.names = FALSE)
fwrite(best_table, best_file, sep = "\t", row.names = FALSE)

cat(paste0("Saved all results  : ", all_file,  "\n"))
cat(paste0("Saved best results : ", best_file, "\n"))
cat("Done.\n")

sessionInfo()
