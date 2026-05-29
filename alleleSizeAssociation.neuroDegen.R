# =============================================================================
# Allele Size Association — Neurodegenerative STRs
# =============================================================================
# Gabrielle Altman, adapted from code by Bharati Jadhav
#
# For each TR locus nominally significant in METAL, tests association between
# long-allele size (binary at each cutoff >= 99th percentile) and disease status
# using Firth logistic regression (brglm). Run per chromosome.
# Similar script used in AoU with different covariates.
#
# Usage:
#   Rscript alleleSizeAssociation.neuroDegen.R \
#     --dir <dir> --cohort neuroDegen --chrom chr12 \
#     --pheno_file phenos.txt --covar_file covars.tsv --metal_file metal.xlsx
#
# Arguments   :
#   --dir        Main project directory
#   --cohort     Phenotype label (e.g. neuroDegen | neuroDegen_noPDorAD)
#   --chrom      Chromosome: chr1 ... chr22
#   --pheno_file Path to phenotype file
#   --covar_file Path to covariates file
#   --metal_file Path to METAL results xlsx
#
# Dependencies: tidyverse, data.table, scales, argparser, R.utils,
#               readxl, parallel, brglm
#
# ============================================================================

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

## args

p <- arg_parser("Per-allele-size association test for TR loci")
p <- add_argument(p, "--dir",        help = "Main project directory [required]")
p <- add_argument(p, "--cohort",     help = "Phenotype label, e.g. neuroDegen [required]")
p <- add_argument(p, "--chrom",      help = "Chromosome, e.g. chr12 [required]")
p <- add_argument(p, "--pheno_file", help = "Path to phenotype file [required]")
p <- add_argument(p, "--covar_file", help = "Path to covariates file [required]")
p <- add_argument(p, "--metal_file", help = "Path to METAL results xlsx [required]")
argv <- parse_args(p)

if (is.na(argv$dir) || is.na(argv$cohort) || is.na(argv$chrom) ||
    is.na(argv$pheno_file) || is.na(argv$covar_file) || is.na(argv$metal_file)) {
  print(p)
  stop("ERROR: --dir, --cohort, --chrom, --pheno_file, --covar_file, and --metal_file are all required.", call. = FALSE)
}

cat("\nArguments:\n")
cat("  --dir        :", argv$dir,        "\n")
cat("  --cohort     :", argv$cohort,     "\n")
cat("  --chrom      :", argv$chrom,      "\n")
cat("  --pheno_file :", argv$pheno_file, "\n")
cat("  --covar_file :", argv$covar_file, "\n")
cat("  --metal_file :", argv$metal_file, "\n\n")


setwd(argv$dir)

# paths — edit for your environment
EH_DIR  <- "/path/to/EH/genotypes/PerChrGT_FilteredEUR_AfterPCA"
OUT_DIR <- file.path(argv$dir, "alleleSizeAssociation")

PHENO_FILE <- argv$pheno_file
COVAR_FILE <- argv$covar_file
METAL_FILE <- argv$metal_file

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

test_locus <- function(marker) {
  message("Testing locus: ", marker)
  df      <- eh[eh$Locus == marker, ]
  cutoff  <- quantile(df$LongAllele, probs = 0.99, na.rm = TRUE)
  alleles <- sort(unique(df$LongAllele[df$LongAllele >= cutoff]))

  if (length(alleles) == 0) {
    message("  No alleles above 99th-percentile cutoff — skipping.")
    return(NULL)
  }

  results_per_allele <- parallel::mclapply(alleles, function(len) {
    dt <- make_allele_bins(len, df)
    counts <- data.frame(
      RepeatSize         = len,
      CaseLargeRepeat    = sum(dt$Phenotype == 1 & dt$RepeatStatus == 1),
      ControlLargeRepeat = sum(dt$Phenotype == 0 & dt$RepeatStatus == 1),
      CaseSmallRepeat    = sum(dt$Phenotype == 1 & dt$RepeatStatus == 0),
      ControlSmallRepeat = sum(dt$Phenotype == 0 & dt$RepeatStatus == 0)
    )
    reg <- tryCatch(
      run_binary_regression(len, dt),
      error = function(e) {
        message("  Error at cutoff ", len, ": ", e$message)
        NULL
      }
    )
    if (is.null(reg)) return(NULL)
    cbind(
      rename(reg, RepeatSize = repeatLen),
      counts[, c("CaseLargeRepeat", "ControlLargeRepeat",
                 "CaseSmallRepeat", "ControlSmallRepeat")]
    )
  }, mc.cores = N_CORES)

  valid_results <- Filter(Negate(is.null), results_per_allele)
  if (length(valid_results) == 0) {
    message("  No valid regression results for locus ", marker, " — skipping.")
    return(NULL)
  }

  assoc <- do.call(rbind, valid_results) %>%
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

all_file  <- file.path(OUT_DIR, paste0(argv$cohort, ".", argv$chrom,
                         ".alleleSizeAssociation.AllResults.txt.gz"))
best_file <- file.path(OUT_DIR, paste0(argv$cohort, ".", argv$chrom,
                         ".alleleSizeAssociation.BestResults.txt.gz"))

fwrite(all_table,  all_file,  sep = "\t", row.names = FALSE)
fwrite(best_table, best_file, sep = "\t", row.names = FALSE)

cat(paste0("Saved all results  : ", all_file,  "\n"))
cat(paste0("Saved best results : ", best_file, "\n"))
cat("Done.\n")

sessionInfo()
