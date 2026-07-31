# =============================================================================
# AMP-PD TR-eQTL Association — Neurodegenerative STRs
# =============================================================================
# Gabrielle N. Altman, adapted from code by Mariya Shadrina, Celine A. Manigbas,
# Alejandro Martin-Trujillo
#
# TR-eQTL association testing in AMP-PD (blood RNA-seq, single tissue):
# expression residualization followed by long-allele association testing at
# host genes and genes within a window of each TR locus. Similar script used
# for GTEx and MESA; see GTEx.eQTLAssociation.neuroDegen.R for the GTEx
# version.
#
# Method:
#   1. Rank-based inverse-normal transform (INT) each needed gene's raw
#      normalized expression, then residualize on cohort, sex, age, median
#      insert size, and genotyping principal components
#      (lm(expr ~ covariates), keep residuals).
#   2. Per locus x gene: expand each TR locus +/- WindowKb, overlap against a
#      gene annotation BED for host-gene + window-gene pairs, then test
#      lm(residual ~ long_allele_length), gated on --MinSamples and
#      --MinDistinctAlleles. A leverage check refits after dropping the
#      single largest-allele sample (beta_drop_top / p_drop_top /
#      delta_beta_drop_top).
#   3. Bonferroni + FDR correction across all tests.
#
# Notes:
#   - gene_id/gene_name are resolved directly from the expression matrix's own
#     gene_id/gene_name columns; ambiguous genes are excluded from the
#     host/window-gene lookup.
#   - --GeneAnnotationFile is a headerless BED (chrom, start, end, gene_id by
#     column position, not header name).
#   - Host genes come from the ANNOVAR closest-gene annotation at each locus
#     and are always tested regardless of distance; window genes are added
#     only within +/-WindowKb.
#   - RNA-seq sample IDs are matched to genotype SampleIds by truncating to
#     the first two hyphen-delimited fields (e.g. "BF-1001-SVM0_5T1" -> "BF-1001").
#
# Usage:
#   Rscript AMPPD.eQTLAssociation.neuroDegen.R \
#     --Dir /path/to/dir --RegionsFile /path/to/regions.txt --FolderName neuroDegen \
#     --GenotypeFile /path/to/genotypes.tsv --ExpressionFile /path/to/expression.tsv \
#     --AnnotationFile /path/to/annotation.tsv
#
# Required: --Dir --RegionsFile --FolderName --GenotypeFile --ExpressionFile
#           --AnnotationFile
#
# Key optional: --WindowKb (100) --MinDistinctAlleles (5) --MinSamples (10)
#   --FDRThreshold (0.1) --GeneAnnotationFile. See --help for the full list.
#
# Dependencies: argparser, data.table, dplyr, tidyr, tibble, readr
#
# HPC example (LSF):
#   bsub -P acc_PROJECTID -q premium -W 4:00 -R "span[hosts=1]" -R "rusage[mem=30000]" \
#     Rscript AMPPD.eQTLAssociation.neuroDegen.R --Dir ... --RegionsFile ... --FolderName neuroDegen \
#       --GenotypeFile ... --ExpressionFile ... --AnnotationFile ...
# =============================================================================

.libPaths(c("~/.Rlib", .libPaths()))

cat(paste0("Loading libraries\n"))
suppressWarnings(suppressMessages({
  library(argparser)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
}))

cat(paste0("\n", Sys.Date(), "\n"))
cat(paste0("R version: ", getRversion(), "\n"))

########################################################################################
# Arguments
########################################################################################

p <- arg_parser("AMP-PD TR-eQTL analysis (long allele, window scan)")
p <- add_argument(p, "--Dir", help = "Main project directory (output written here) [required]")
p <- add_argument(p, "--RegionsFile", help = "File of Locus values (with chr/start/end encoded as chr_start_end) to test [required]")
p <- add_argument(p, "--FolderName", help = "Output subfolder name, or 'no' for none [required]")
p <- add_argument(p, "--GenotypeFile",
                  help = "TR genotype + classifier + covariate file (VARID, SampleId, long_allele, Cohort, Gender, age_at_baseline, MedianInsertSize, SNP_PC1-5) [required]")
p <- add_argument(p, "--ExpressionFile",
                  help = "Normalized gene expression matrix (gene_id, gene_name, plus one sample column per RNA-seq library) [required]")
p <- add_argument(p, "--AnnotationFile",
                  help = "ANNOVAR annotation file with Chr/Start/End/Gene.ensGene/Func.ensGene columns, used to assign each locus its host gene [required]")
p <- add_argument(p, "--GeneAnnotationFile", default = "no",
                  help = "Headerless BED (chrom, start, end, gene_id) for the window scan. 'no' disables the window scan (host gene only).")
p <- add_argument(p, "--WindowKb", default = 100,
                  help = "Window (kb) each side of the TR locus to search for additional genes.")
p <- add_argument(p, "--MinDistinctAlleles", default = 5,
                  help = "Minimum distinct long-allele values to fit the model.")
p <- add_argument(p, "--MinSamples", default = 10,
                  help = "Minimum samples per (locus,gene) to run the test.")
p <- add_argument(p, "--FDRThreshold", default = 0.1,
                  help = "FDR (q-value) threshold for calling a significant TR:gene eQTL.")
argv <- parse_args(p)

if (is.na(argv$Dir) || is.na(argv$RegionsFile) || is.na(argv$FolderName) ||
    is.na(argv$GenotypeFile) || is.na(argv$ExpressionFile) || is.na(argv$AnnotationFile)) {
  print(p)
  stop("ERROR: --Dir, --RegionsFile, --FolderName, --GenotypeFile, --ExpressionFile, and --AnnotationFile are all required.", call. = FALSE)
}

Dir                <- argv$Dir
RegionsFile        <- argv$RegionsFile
FolderName         <- argv$FolderName
GenotypeFile       <- argv$GenotypeFile
ExpressionFile     <- argv$ExpressionFile
AnnotationFile     <- argv$AnnotationFile
GeneAnnotationFile <- argv$GeneAnnotationFile
WindowKb           <- as.numeric(argv$WindowKb)
MinDistinctAlleles <- as.integer(argv$MinDistinctAlleles)
MinSamples         <- as.integer(argv$MinSamples)
FDRThreshold       <- as.numeric(argv$FDRThreshold)

classifierPath     <- GenotypeFile
geneExpressionPath <- ExpressionFile
annotationPath     <- AnnotationFile

if (FolderName != "no") {
  outdir <- file.path(Dir, "geneExpression", FolderName, "eqtlAnalysis")
} else {
  outdir <- file.path(Dir, "geneExpression", "eqtlAnalysis")
}
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

cat("\nArguments:\n")
cat("  Dir:                ", Dir, "\n")
cat("  Genotype file:      ", GenotypeFile, "\n")
cat("  Expression file:    ", ExpressionFile, "\n")
cat("  Annotation file:    ", AnnotationFile, "\n")
cat("  Gene annotation:    ", GeneAnnotationFile, "\n")
cat("  Window (kb):        ", WindowKb, "\n")
cat("  MinDistinctAlleles: ", MinDistinctAlleles, "\n")
cat("  MinSamples:         ", MinSamples, "\n")
cat("  FDR threshold:      ", FDRThreshold, "\n")
cat("  Out dir:            ", outdir, "\n")

########################################################################################
# Helpers
########################################################################################

inverse_normal <- function(x) {
  n <- sum(!is.na(x))
  qnorm((rank(x, na.last = "keep") - 0.5) / n)
}

# Rank-based inverse-normal transform expression, then residualize on covariates
# via lm, keeping the residuals with no further transform.
compute_z <- function(d) {
  d$expr_norm <- inverse_normal(as.numeric(d$expression))

  candidate_covs <- c("Cohort", "Gender", "age_at_baseline",
                      "MedianInsertSize",
                      "SNP_PC1", "SNP_PC2", "SNP_PC3", "SNP_PC4", "SNP_PC5")
  covars <- candidate_covs[candidate_covs %in% names(d)]
  covars <- covars[sapply(covars, function(cv) {
    v <- d[[cv]]
    length(unique(na.omit(v))) > 1 && mean(!is.na(v)) > 0.5
  })]

  for (cv in intersect(covars, c("Cohort", "Gender"))) {
    d[[cv]] <- as.factor(d[[cv]])
  }

  use_cols <- c("expr_norm", covars)
  ok <- complete.cases(d[, use_cols, drop = FALSE])

  if (length(covars) == 0 || sum(ok) < MinSamples) {
    d$resid <- d$expr_norm - mean(d$expr_norm, na.rm = TRUE)
    d$z     <- d$resid
    return(d)
  }
  fmla <- as.formula(paste0("expr_norm ~ ", paste(covars, collapse = " + ")))
  fit  <- tryCatch(lm(fmla, data = d[ok, ]), error = function(e) NULL)
  if (is.null(fit)) {
    d$resid <- d$expr_norm - mean(d$expr_norm, na.rm = TRUE)
  } else {
    d$resid <- d$expr_norm - predict(fit, newdata = d)
  }
  d$z <- d$resid
  d
}

fit_model <- function(d) {
  if (nrow(d) < MinSamples) return(NULL)
  if (length(unique(d$long_allele)) < MinDistinctAlleles) return(NULL)

  fit <- tryCatch(lm(z ~ long_allele, data = d), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  co <- summary(fit)$coefficients
  if (!"long_allele" %in% rownames(co)) return(NULL)

  beta <- co["long_allele", "Estimate"]
  se   <- co["long_allele", "Std. Error"]
  pval <- co["long_allele", "Pr(>|t|)"]
  r_squared <- summary(fit)$r.squared

  # Leverage: refit without top sample
  max_val <- max(d$long_allele, na.rm = TRUE)
  top_idx <- which(d$long_allele == max_val)
  d_drop  <- d[-top_idx[1], , drop = FALSE]
  beta_drop <- NA_real_; p_drop <- NA_real_
  if (nrow(d_drop) >= MinSamples && length(unique(d_drop$long_allele)) >= MinDistinctAlleles) {
    fit2 <- tryCatch(lm(z ~ long_allele, data = d_drop), error = function(e) NULL)
    if (!is.null(fit2)) {
      co2 <- summary(fit2)$coefficients
      if ("long_allele" %in% rownames(co2)) {
        beta_drop <- co2["long_allele", "Estimate"]
        p_drop    <- co2["long_allele", "Pr(>|t|)"]
      }
    }
  }

  list(
    n = nrow(d), max_allele = max_val, n_at_max = length(top_idx),
    n_distinct     = length(unique(d$long_allele)),
    beta           = beta, se = se, p = pval, r_squared = r_squared,
    beta_drop_top  = beta_drop, p_drop_top = p_drop,
    delta_beta_drop_top = beta_drop - beta
  )
}

########################################################################################
# Load data (genotype/classifier, annotation, expression)
########################################################################################

cat("\nLoading regions ...\n")
regionsToKeep <- fread(RegionsFile, sep = "\t", check.names = FALSE, header = TRUE)
loci_tmp <- tempfile(fileext = ".txt")
writeLines(regionsToKeep$Locus, loci_tmp)
cat(nrow(regionsToKeep), "loci\n")

cat("\nLoading AMP-PD classifier data (pre-filtered to test loci) ...\n")
df <- fread(
  cmd = paste0("{ head -1 ", classifierPath,
               "; grep -F -f ", loci_tmp, " ", classifierPath, "; }"),
  sep = "\t", check.names = FALSE, header = TRUE
)
df <- df %>% mutate(long_allele = suppressWarnings(as.numeric(long_allele)))
cat("Filtered to", n_distinct(df$VARID), "VARIDs\n")

cat("\nLoading gene annotation (closest gene per locus) ...\n")
annotation <- fread(annotationPath, sep = "\t", check.names = FALSE, header = TRUE)
annotation$VARID <- paste0(annotation$Chr, "_", annotation$Start, "_", annotation$End)

cat("\nLoading gene expression ...\n")
geneExpression <- fread(geneExpressionPath, sep = "\t", check.names = FALSE, header = TRUE)

expr_ann_cols <- c("gene_id", "#chrom", "start", "end", "strand", "gene_type", "gene_name")
sample_cols   <- setdiff(names(geneExpression), expr_ann_cols)

# RNA-seq IDs look like "BF-1001-SVM0_5T1" -> strip to "BF-1001" to match SampleId
sample_id_map <- data.frame(
  RNAseq_col = sample_cols,
  SampleId   = sub("^([^-]+-[^-]+)-.*", "\\1", sample_cols),
  stringsAsFactors = FALSE
)

covariates_df <- df %>%
  select(SampleId, Cohort, Gender, PDStatus, age_at_baseline,
         MedianInsertSize, SNP_PC1, SNP_PC2, SNP_PC3, SNP_PC4, SNP_PC5) %>%
  distinct()

# gene_id <-> gene_name lookup built directly from the expression matrix -- no
# external crosswalk file needed. Genes appearing more than once (ambiguous) are
# dropped from the lookup so they can't be silently matched to the wrong row.
gene_lookup <- geneExpression %>%
  select(gene_id, gene_name) %>%
  mutate(gene_id = sub("\\..*$", "", gene_id)) %>%
  distinct()
dup_by_name <- gene_lookup %>% count(gene_name) %>% filter(n > 1) %>% pull(gene_name)
dup_by_id   <- gene_lookup %>% count(gene_id)   %>% filter(n > 1) %>% pull(gene_id)
if (length(dup_by_name) > 0) {
  cat("Note:", length(dup_by_name), "gene_name values are ambiguous in the expression matrix; excluded from host/window-gene lookup.\n")
}
gene_lookup <- gene_lookup %>% filter(!(gene_name %in% dup_by_name), !(gene_id %in% dup_by_id))

########################################################################################
# Build the host-gene + window-gene test list
########################################################################################

host_pairs <- annotation %>%
  select(VARID, Gene.ensGene, Func.ensGene) %>%
  rename(Locus = VARID, gene_name = Gene.ensGene, genomicFeature = Func.ensGene) %>%
  filter(Locus %in% regionsToKeep$Locus) %>%
  distinct() %>%
  left_join(gene_lookup, by = "gene_name") %>%
  mutate(is_host_gene = TRUE)

n_no_expr <- sum(is.na(host_pairs$gene_id))
if (n_no_expr > 0) {
  cat("WARNING:", n_no_expr, "host-gene locus/gene pairs have no matching row in the expression matrix (ambiguous or missing gene_name) -- will be skipped:\n")
  print(host_pairs %>% filter(is.na(gene_id)) %>% select(Locus, gene_name))
}
host_pairs <- host_pairs %>% filter(!is.na(gene_id))
cat(nrow(host_pairs), "host-gene (locus, gene) pairs to test\n")

window_gene_map <- NULL
if (GeneAnnotationFile != "no" && file.exists(GeneAnnotationFile)) {
  geneAnnot <- fread(GeneAnnotationFile, sep = "\t", header = FALSE, check.names = FALSE) %>% as.data.frame()
  geneAnnot <- geneAnnot[, 1:4]
  names(geneAnnot) <- c("chr_ann", "start_ann", "end_ann", "gene_ann")
  geneAnnot <- geneAnnot %>%
    mutate(chr_ann = gsub("^chr", "", as.character(chr_ann)),
           start_ann = as.numeric(start_ann), end_ann = as.numeric(end_ann),
           gene_id = sub("\\..*$", "", gene_ann))
  cat("Loaded gene annotation for", nrow(geneAnnot), "genes from", GeneAnnotationFile, "\n")

  locusCoords <- regionsToKeep %>%
    select(Locus) %>%
    distinct() %>%
    separate(Locus, into = c("chr_loc", "start_loc", "end_loc"),
             sep = "_", extra = "drop", remove = FALSE) %>%
    mutate(chr_loc   = gsub("^chr", "", as.character(chr_loc)),
           start_loc = as.numeric(start_loc),
           end_loc   = as.numeric(end_loc)) %>%
    filter(!is.na(chr_loc), !is.na(start_loc), !is.na(end_loc))

  win_bp <- WindowKb * 1000
  window_rows <- vector("list", nrow(locusCoords))
  for (i in seq_len(nrow(locusCoords))) {
    loc <- locusCoords[i, ]
    hits <- geneAnnot %>%
      filter(chr_ann == loc$chr_loc,
             end_ann   >= (loc$start_loc - win_bp),
             start_ann <= (loc$end_loc   + win_bp))
    if (nrow(hits) > 0) {
      window_rows[[i]] <- data.frame(Locus = loc$Locus, gene_id = unique(hits$gene_id), stringsAsFactors = FALSE)
    }
  }
  window_gene_map <- bind_rows(window_rows) %>%
    left_join(gene_lookup, by = "gene_id") %>%
    filter(!is.na(gene_name)) %>%
    mutate(is_host_gene = FALSE, genomicFeature = NA_character_)
  cat("Window scan (+/-", WindowKb, "kb) found ", nrow(window_gene_map),
      " (locus, gene) candidate pairs across ", length(unique(window_gene_map$Locus)),
      " loci.\n", sep = "")
} else {
  cat("GeneAnnotationFile not provided -- window scan DISABLED, host-gene pairs only.\n")
}

test_list <- bind_rows(
  host_pairs      %>% select(Locus, gene_id, gene_name, genomicFeature, is_host_gene),
  window_gene_map %>% select(Locus, gene_id, gene_name, genomicFeature, is_host_gene)
) %>% distinct(Locus, gene_id, .keep_all = TRUE)
cat(nrow(test_list), "total (locus, gene) pairs to test (host + window)\n")

needed_gene_ids <- unique(test_list$gene_id)

########################################################################################
# Residualize every needed gene once (host + window genes only, not the whole
# expression matrix)
########################################################################################

cat("\nResidualizing", length(needed_gene_ids), "needed genes ...\n")
resid_cache <- list()
for (gid in needed_gene_ids) {
  exprRow <- geneExpression %>% filter(sub("\\..*$", "", gene_id) == gid)
  if (nrow(exprRow) != 1) next

  expr_vals <- as.numeric(exprRow[1, sample_cols, with = FALSE])
  exprSub <- data.frame(
    RNAseq_col = sample_cols,
    expression = expr_vals,
    stringsAsFactors = FALSE
  ) %>%
    left_join(sample_id_map,  by = "RNAseq_col") %>%
    inner_join(covariates_df, by = "SampleId")

  if (nrow(exprSub) == 0) next
  exprSub <- compute_z(exprSub)
  resid_cache[[gid]] <- exprSub %>% select(SampleId, expression, z)
}
cat("Residualized", sum(!sapply(resid_cache, is.null)), "of", length(needed_gene_ids), "needed genes.\n")

########################################################################################
# Main loop: test every (locus, gene) pair that maps to a gene we residualized
########################################################################################

per_sample_long <- list()
per_locus_rows  <- list()

pairs_here <- test_list %>% filter(gene_id %in% names(resid_cache))
cat("\nTesting", nrow(pairs_here), "(locus,gene) pairs ...\n")

for (i in seq_len(nrow(pairs_here))) {
  region         <- pairs_here$Locus[i]
  gid            <- pairs_here$gene_id[i]
  gname          <- pairs_here$gene_name[i]
  isHost         <- pairs_here$is_host_gene[i]
  genomicFeature <- pairs_here$genomicFeature[i]

  alleles <- df %>%
    filter(VARID == region) %>%
    select(SampleId, long_allele) %>%
    filter(!is.na(long_allele)) %>%
    distinct()
  if (nrow(alleles) == 0) next

  d <- resid_cache[[gid]] %>% inner_join(alleles, by = "SampleId") %>% filter(!is.na(z))
  if (nrow(d) == 0) next

  per_sample_long[[length(per_sample_long) + 1]] <- d %>%
    mutate(Locus = region, gene_id = gid, gene_name = gname, is_host_gene = isHost) %>%
    select(Locus, gene_id, gene_name, is_host_gene, SampleId, expression, z, long_allele)

  fit_res <- fit_model(d)
  if (is.null(fit_res)) {
    cat(sprintf("[%2d] %-25s gene=%-12s skipped (n=%d, n_distinct=%d)\n",
                i, region, gname, nrow(d), length(unique(d$long_allele))))
    next
  }
  if (isHost) cat("    ", region, " x ", gid, " (", gname, "): association lm used ", fit_res$n, " samples\n", sep = "")

  per_locus_rows[[length(per_locus_rows) + 1]] <- data.frame(
    Locus                = region,
    gene_id              = gid,
    gene_name            = gname,
    is_host_gene         = isHost,
    genomicFeature       = genomicFeature,
    n                    = fit_res$n,
    n_distinct           = fit_res$n_distinct,
    max_allele           = fit_res$max_allele,
    n_at_max             = fit_res$n_at_max,
    beta                 = fit_res$beta,
    se                   = fit_res$se,
    p                    = fit_res$p,
    r_squared            = fit_res$r_squared,
    beta_drop_top        = fit_res$beta_drop_top,
    p_drop_top           = fit_res$p_drop_top,
    delta_beta_drop_top  = fit_res$delta_beta_drop_top,
    stringsAsFactors     = FALSE
  )
}

########################################################################################
# Combine, Bonferroni + FDR-correct, write
########################################################################################

per_sample_df <- bind_rows(per_sample_long)
per_locus_df  <- bind_rows(per_locus_rows)

if (nrow(per_locus_df) > 0) {
  N <- nrow(per_locus_df)
  per_locus_df <- per_locus_df %>%
    mutate(bonferroni      = pmin(p * N, 1),
           q_fdr           = p.adjust(p, method = "fdr"),
           significant_QTL = q_fdr < FDRThreshold) %>%
    arrange(p)
}

suffix <- ".eQTL"

fwrite(per_sample_df,
       file = file.path(outdir, paste0("AMPPD", suffix, ".per_sample_z_longAllele.txt")),
       sep = "\t", row.names = FALSE)
fwrite(per_locus_df %>% rename(`Slope (repeat units)` = beta, `R-squared` = r_squared),
       file = file.path(outdir, paste0("AMPPD", suffix, ".per_locus_lm.txt")),
       sep = "\t", row.names = FALSE)

cat("\n========== Results (sorted by p-value) ==========\n")
if (nrow(per_locus_df) > 0) {
  print(per_locus_df %>%
    select(Locus, gene_id, gene_name, is_host_gene, n, n_distinct, beta, se, p, q_fdr, significant_QTL, delta_beta_drop_top))
}

cat("\n========== Significant TR:gene eQTLs (FDR q<", FDRThreshold, ") ==========\n", sep = "")
if (nrow(per_locus_df) > 0) {
  print(per_locus_df %>% filter(significant_QTL) %>%
    select(Locus, gene_id, gene_name, is_host_gene, n, beta, p, q_fdr))
}

cat("\nOutputs written to:\n")
cat(" ", file.path(outdir, paste0("AMPPD", suffix, ".per_sample_z_longAllele.txt")), "\n")
cat(" ", file.path(outdir, paste0("AMPPD", suffix, ".per_locus_lm.txt")), "\n")
