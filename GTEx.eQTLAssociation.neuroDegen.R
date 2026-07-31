# =============================================================================
# GTEx TR-eQTL Association — Neurodegenerative STRs
# =============================================================================
# Gabrielle N. Altman, adapted from code by Mariya Shadrina, Celine A. Manigbas,
# Alejandro Martin-Trujillo
#
# TR-eQTL association testing in GTEx: per-tissue expression residualization
# followed by long-allele association testing at host genes and genes within a
# window of each TR locus, with cross-tissue meta-analysis. Approach inspired
# by the GTEx TR-eQTL analysis in Manigbas et al. 2024 Nat Commun.
#
# Similar script used for MESA and AMP-PD, differing in expression processing:
#   GTEx:   pre-processed normalized_expression.bed, residualized on sex,
#           insert size, sequencer, genotyping PCs, PEER factors.
#   MESA:   raw TPM, INT-transformed then residualized on sex, insert size,
#           age, analyte isolation batch, predicted ancestry.
#   AMP-PD: raw normalized counts, INT-transformed then residualized on sex,
#           age, cohort (study of origin), insert size, 5 genotyping PCs.
#
# Method:
#   1. Per tissue: load normalized_expression.bed + covariates.txt (every
#      covariate column used), join AGE/SEX from the subject phenotype file
#      and an optional insert-size file, then residualize each needed gene
#      (lm(expression ~ all covariates), keep residuals, no INT).
#   2. Per locus x gene x tissue: expand each TR locus +/- WindowKb, overlap
#      against a gene annotation BED for host-gene + window-gene pairs, then
#      test lm(residual ~ long_allele_length), gated on --MinSamples and
#      --MinDistinctAlleles.
#   3. Per-tissue FDR correction, plus cross-tissue random-effects meta-analysis
#      (metafor::rma; brain/nonbrain/all) for host-gene pairs.
#
# Notes:
#   - gene_id is matched as a version-stripped Ensembl ID; use
#     --GeneSymbolToEnsemblFile if host genes are symbols instead.
#   - --GeneAnnotationFile is a headerless BED (chrom, start, end, gene_id by
#     column position, not header name).
#   - Sample IDs are assumed SUBJID-style ("GTEX-1117F")
#
# Usage:
#   Rscript GTEx.eQTLAssociation.neuroDegen.R \
#     --Dir /path/to/dir --RegionsFile /path/to/regions.txt --FolderName neuroDegen \
#     --GenotypeFile /path/to/genotypes.txt.gz --GenotypeAllFile /path/to/all_alleles.bed.gz \
#     --ExpressionMatrixDir /path/to/GTEx_Analysis_v8_eQTL_expression_matrices \
#     --PEERDir /path/to/GTEx_Analysis_v8_eQTL_covariates
#
# Required: --Dir --RegionsFile --FolderName --GenotypeFile --GenotypeAllFile
#           --ExpressionMatrixDir --PEERDir
#
# Key optional: --WindowKb (100) --MinSamples (10) --FDRThreshold (0.1)
#   --GeneAnnotationFile --GeneSymbolToEnsemblFile --SubjectPhenoFile
#   --InsertSizeFile --Threads (8). See --help for the full list.
#
# Dependencies: argparser, data.table, dplyr, tidyverse, readxl, metafor, parallel
#
# HPC example (LSF):
#   bsub -P acc_PROJECTID -q premium -W 10:00 -n 8 -R "span[hosts=1]" -R "rusage[mem=60000]" \
#     Rscript GTEx.eQTLAssociation.neuroDegen.R --Dir ... --RegionsFile ... --FolderName neuroDegen \
#       --GenotypeFile ... --GenotypeAllFile ... --ExpressionMatrixDir ... --PEERDir ... --Threads 8
# =============================================================================

.libPaths(c("~/.Rlib", .libPaths()))

cat(paste0("Loading libraries\n"))
suppressWarnings(suppressMessages({
  library(argparser)
  library(data.table)
  library(dplyr)
  library(tidyverse)
  library(readxl)
  library(metafor)
  library(parallel)
}))

cat(paste0("\n", Sys.Date(), "\n"))
cat(paste0("R version: ", getRversion(), "\n"))

########################################################################################
# Arguments
########################################################################################

p <- arg_parser("GTEx TR-eQTL analysis (long allele, window scan, host-gene meta-analysis)")
p <- add_argument(p, "--Dir", help = "Main project directory (output + TR genotype input files live here) [required]")
p <- add_argument(p, "--RegionsFile", help = "File path or 'no' -- which Locus values to keep [required]")
p <- add_argument(p, "--FolderName", help = "Output subfolder name, or 'no' for none [required]")
p <- add_argument(p, "--GenotypeFile",
                  help = "TR genotype + annotation file (per-sample allele calls, closestGene/geneName) [required]")
p <- add_argument(p, "--GenotypeAllFile",
                  help = "TR genotype file with the full per-locus allele distribution (Max column) [required]")
p <- add_argument(p, "--ExpressionMatrixDir",
                  help = "Directory with {tissue}.v8.normalized_expression.bed[.gz] files (GTEx_Analysis_v8_eQTL_expression_matrices) [required]")
p <- add_argument(p, "--PEERDir",
                  help = "Directory with {tissue}.v8.covariates.txt files (GTEx_Analysis_v8_eQTL_covariates) [required]")
p <- add_argument(p, "--SubjectPhenoFile", default = "no",
                  help = "GTEx subject phenotypes file, used for AGE (and SEX if not already a covariates-file column). 'no' to skip.")
p <- add_argument(p, "--InsertSizeFile", default = "no",
                  help = "Tab-delimited file: SUBJID + insert size column (e.g. MedianInsertSize.dat). 'no' to skip.")
p <- add_argument(p, "--MinDistinctAlleles", default = 5,
                  help = "Minimum distinct long-allele values in a tissue to fit the per-tissue model.")
p <- add_argument(p, "--MinSamples", default = 10,
                  help = "Minimum samples per (locus,gene,tissue) to run the test.")
p <- add_argument(p, "--Threads", default = 8,
                  help = "Number of tissues to process in parallel (forked worker processes, Linux/Mac only). Match this to the bsub -n core count.")
p <- add_argument(p, "--FDRThreshold", default = 0.1,
                  help = "FDR (q-value) threshold for calling a significant TR:gene eQTL.")
p <- add_argument(p, "--WindowKb", default = 100,
                  help = "Window (kb) each side of the TR locus to search for additional genes.")
p <- add_argument(p, "--GeneAnnotationFile", default = "no",
                  help = "Headerless BED (chrom, start, end, gene_id) for the window scan. Gene IDs should be in the SAME format as your expression bed's gene_id column (Ensembl recommended). 'no' disables the window scan (host gene only).")
p <- add_argument(p, "--GeneSymbolToEnsemblFile", default = "no",
                  help = "Optional symbol,ensembl_id crosswalk to translate host-gene symbols (closestGene/geneName) into your expression bed's gene_id format. 'no' = assume closestGene/geneName already match gene_id directly.")
argv <- parse_args(p)

if (is.na(argv$Dir) || is.na(argv$RegionsFile) || is.na(argv$FolderName) ||
    is.na(argv$GenotypeFile) || is.na(argv$GenotypeAllFile) ||
    is.na(argv$ExpressionMatrixDir) || is.na(argv$PEERDir)) {
  print(p)
  stop("ERROR: --Dir, --RegionsFile, --FolderName, --GenotypeFile, --GenotypeAllFile, --ExpressionMatrixDir, and --PEERDir are all required.", call. = FALSE)
}

Dir                   <- argv$Dir
RegionsFile           <- argv$RegionsFile
FolderName            <- argv$FolderName
GenotypeFile          <- argv$GenotypeFile
GenotypeAllFile       <- argv$GenotypeAllFile
ExpressionMatrixDir   <- argv$ExpressionMatrixDir
PEERDir               <- argv$PEERDir
SubjectPhenoFile      <- argv$SubjectPhenoFile
InsertSizeFile        <- argv$InsertSizeFile
MinDistinctAlleles    <- as.integer(argv$MinDistinctAlleles)
MinSamples            <- as.integer(argv$MinSamples)
Threads               <- max(1L, as.integer(argv$Threads))
FDRThreshold          <- as.numeric(argv$FDRThreshold)
WindowKb              <- as.numeric(argv$WindowKb)
GeneAnnotationFile    <- argv$GeneAnnotationFile
GeneSymbolToEnsemblFile <- argv$GeneSymbolToEnsemblFile

if(FolderName != "no"){
  outdir <- paste0(Dir, "/geneExpression/", FolderName, "/eqtlAnalysis/")
} else {
  outdir <- paste0(Dir, "/geneExpression/eqtlAnalysis/")
}
if(!file.exists(outdir)) dir.create(outdir, recursive = TRUE)

COHORT <- "GTEx"

cat("\nArguments:\n")
cat("  Dir:                    ", Dir, "\n")
cat("  Genotype file:          ", GenotypeFile, "\n")
cat("  Genotype (all) file:    ", GenotypeAllFile, "\n")
cat("  Expression matrix dir:  ", ExpressionMatrixDir, "\n")
cat("  Covariates dir:         ", PEERDir, "\n")
cat("  Subject phenotype file: ", SubjectPhenoFile, "\n")
cat("  Insert size file:       ", InsertSizeFile, "\n")
cat("  MinDistinctAlleles:     ", MinDistinctAlleles, "\n")
cat("  MinSamples:             ", MinSamples, "\n")
cat("  Threads:                ", Threads, "\n")
cat("  FDR threshold:          ", FDRThreshold, "\n")
cat("  Window (kb):            ", WindowKb, "\n")
cat("  Gene annotation:        ", GeneAnnotationFile, "\n")
cat("  Symbol->Ensembl x-walk: ", GeneSymbolToEnsemblFile, "\n")
cat("  Out dir:                ", outdir, "\n")

########################################################################################
# Load TR genotype data
########################################################################################

setwd(Dir)
cat("\nLoading TR genotype file ...\n")

df <- fread(GenotypeFile, sep = "\t", check.names = F, header = T) %>%
  separate(SampleId, into = c("cohort","id"), remove = FALSE, extra = "drop") %>%
  mutate(MergeId = paste0(cohort, "-", id)) %>%
  select(-c(cohort, id))

if(RegionsFile != "no"){
  regionsToKeep <- fread(RegionsFile, sep = "\t", check.names = F, header = T)
  df <- df %>% filter(Locus %in% regionsToKeep$Locus)
}

data_all <- fread(GenotypeAllFile, sep = "\t", check.names = F, header = T) %>%
  separate(SampleId, into = c("cohort","id"), remove = FALSE, extra = "drop") %>%
  mutate(MergeId = paste0(cohort, "-", id)) %>%
  select(-c(cohort, id)) %>%
  mutate(Max_numeric = suppressWarnings(as.numeric(ifelse(Max == ".", NA, Max))))

########################################################################################
# Subject phenotypes (AGE always pulled here; SEX pulled here only if not already a
# column in the per-tissue covariates file 
########################################################################################

has_subj <- FALSE
subjPheno <- NULL
if(SubjectPhenoFile != "no" && file.exists(SubjectPhenoFile)){
  subjPheno <- fread(SubjectPhenoFile, sep = "\t", header = TRUE) %>% as.data.frame()
  has_subj <- TRUE
  if("AGE" %in% names(subjPheno)){
    if(any(grepl("-", subjPheno$AGE))){
      # GTEx AGE sometimes binned ("20-29",...). Convert to bin midpoint if so.
      subjPheno$AGE <- suppressWarnings(as.numeric(sub("-.*","", subjPheno$AGE))) + 5
    } else {
      subjPheno$AGE <- suppressWarnings(as.numeric(subjPheno$AGE))
    }
  }
  cat("Subject phenotype file loaded.\n")
} else {
  cat("SubjectPhenoFile not provided -- proceeding without AGE/SEX fallback.\n")
}

insertSizeDf <- NULL
if(InsertSizeFile != "no" && file.exists(InsertSizeFile)){
  insertSizeDf <- fread(InsertSizeFile, sep = "\t", header = TRUE, check.names = FALSE) %>% as.data.frame()
  cat("Loaded insert size covariate from", InsertSizeFile, "\n")
} else {
  cat("InsertSizeFile not provided -- proceeding WITHOUT insert size covariate.\n")
}

########################################################################################
# Gene annotation + locus coordinates for the +/-WindowKb window scan
########################################################################################

geneAnnot <- NULL
if(GeneAnnotationFile != "no" && file.exists(GeneAnnotationFile)){
  geneAnnot <- fread(GeneAnnotationFile, sep = "\t", header = FALSE, check.names = FALSE) %>% as.data.frame()
  if(ncol(geneAnnot) < 4){
    stop("GeneAnnotationFile must have at least 4 columns (chrom, start, end, gene_id); found ", ncol(geneAnnot))
  }
  # Headerless BED: chrom, start, end, gene_id -- assigned by position, not by name.
  geneAnnot <- geneAnnot[, 1:4]
  names(geneAnnot) <- c("chr_ann", "start_ann", "end_ann", "gene_ann")
  geneAnnot <- geneAnnot %>%
    mutate(chr_ann = gsub("^chr", "", as.character(chr_ann)),
           start_ann = as.numeric(start_ann), end_ann = as.numeric(end_ann),
           gene_ann_stripped = sub("\\..*$", "", gene_ann))
  cat("Loaded gene annotation for", nrow(geneAnnot), "genes from", GeneAnnotationFile, "\n")
} else {
  cat("GeneAnnotationFile not provided -- window scan DISABLED, host-gene pairs only.\n")
}

locusCoords <- NULL
if(!is.null(geneAnnot)){
  if(RegionsFile == "no"){
    stop("Splitting Locus into chr/start/end requires --RegionsFile (need a Locus list to split).")
  }

  locusCoords <- regionsToKeep %>%
    select(Locus) %>%
    distinct() %>%
    separate(Locus, into = c("chr_loc", "start_loc", "end_loc"),
             sep = "_", extra = "drop", remove = FALSE) %>%
    mutate(chr_loc   = gsub("^chr", "", as.character(chr_loc)),
           start_loc = as.numeric(start_loc),
           end_loc   = as.numeric(end_loc))

  n_bad <- sum(is.na(locusCoords$chr_loc) | is.na(locusCoords$start_loc) | is.na(locusCoords$end_loc))
  if(n_bad > 0){
    cat("WARNING:", n_bad, "of", nrow(locusCoords), "Locus values didn't split into 3 usable ",
        "chr/start/end fields. Example Locus values:\n", sep = "")
    print(head(regionsToKeep$Locus, 5))
    locusCoords <- locusCoords %>% filter(!is.na(chr_loc), !is.na(start_loc), !is.na(end_loc))
  }

  if(nrow(locusCoords) == 0){
    stop("Could not parse any chr/start/end from Locus by splitting on '_'. Check the actual Locus format in your RegionsFile.")
  }
  cat("Parsed locus coordinates for", nrow(locusCoords), "loci by splitting Locus on '_'.\n")
}
window_gene_map <- NULL
if(!is.null(geneAnnot) && !is.null(locusCoords)){
  win_bp <- WindowKb * 1000
  window_rows <- vector("list", nrow(locusCoords))

  for(i in seq_len(nrow(locusCoords))){
    loc <- locusCoords[i, ]
    hits <- geneAnnot %>%
      filter(chr_ann == loc$chr_loc,
             end_ann   >= (loc$start_loc - win_bp),
             start_ann <= (loc$end_loc   + win_bp))
    if(nrow(hits) > 0){
      window_rows[[i]] <- data.frame(Locus = loc$Locus, gene_id = hits$gene_ann_stripped,
                                      stringsAsFactors = FALSE)
    }
  }

  window_gene_map <- bind_rows(window_rows)
  cat("Window scan (+/-", WindowKb, "kb) found ", nrow(window_gene_map),
      " (locus, gene) candidate pairs across ", length(unique(window_gene_map$Locus)),
      " loci.\n", sep = "")
}

########################################################################################
# Build the host-gene + window-gene test list, with gene_id resolved to the format
# used in your expression bed files (Ensembl, version-stripped, assumed).
########################################################################################

symbolToEnsembl <- NULL
if(GeneSymbolToEnsemblFile != "no" && file.exists(GeneSymbolToEnsemblFile)){
  symbolToEnsembl <- fread(GeneSymbolToEnsemblFile, sep = "\t", header = TRUE, check.names = FALSE) %>%
    as.data.frame()
  if(!all(c("symbol","ensembl_id") %in% names(symbolToEnsembl))){
    stop("GeneSymbolToEnsemblFile must have columns: symbol, ensembl_id")
  }
  symbolToEnsembl$ensembl_id <- sub("\\..*$", "", symbolToEnsembl$ensembl_id)
  cat("Loaded symbol->Ensembl crosswalk (", nrow(symbolToEnsembl), " genes) from ",
      GeneSymbolToEnsemblFile, "\n", sep = "")
} else {
  cat("No symbol->Ensembl crosswalk provided -- assuming closestGene/geneName already ",
      "match your expression bed's gene_id format directly.\n", sep = "")
}

# Reverse (Ensembl -> symbol) lookup for annotating output files with gene names.
ensemblToSymbol <- NULL
if(!is.null(symbolToEnsembl)){
  ensemblToSymbol <- symbolToEnsembl %>% select(ensembl_id, symbol) %>% distinct(ensembl_id, .keep_all = TRUE)
}
add_gene_name <- function(d){
  if(is.null(ensemblToSymbol) || !("gene_id" %in% names(d)) || nrow(d) == 0) return(d)
  d %>% left_join(ensemblToSymbol, by = c("gene_id" = "ensembl_id")) %>%
    rename(gene_name = symbol) %>%
    mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name)) %>%
    relocate(gene_name, .after = gene_id)
}

host_pairs <- df %>% select(Locus, closestGene, geneName, genomicFeature) %>% unique()
if(RegionsFile != "no"){
  host_pairs <- host_pairs %>% filter(Locus %in% regionsToKeep$Locus)
}
host_pairs <- host_pairs %>%
  mutate(gene_symbol = case_when(
    !is.na(geneName) & geneName == "DM1-AS" ~ "DMPK",
    closestGene == "NAA38,CHD3" ~ "CHD3",
    TRUE ~ closestGene
  ))
cat(nrow(host_pairs), "host-gene (locus, gene) pairs to test\n")

resolve_gene_id <- function(symbol){
  if(is.null(symbolToEnsembl)) return(symbol)
  hit <- symbolToEnsembl$ensembl_id[match(symbol, symbolToEnsembl$symbol)]
  ifelse(is.na(hit), symbol, hit)
}

test_list <- host_pairs %>%
  transmute(Locus, gene_symbol,
            gene_id = resolve_gene_id(gene_symbol),
            is_host_gene = TRUE)

if(!is.null(window_gene_map)){
  window_pairs <- window_gene_map %>%
    filter(Locus %in% host_pairs$Locus) %>%
    rename(gene_id = gene_id) %>%
    anti_join(test_list %>% select(Locus, gene_id), by = c("Locus","gene_id")) %>%
    mutate(gene_symbol = NA_character_, is_host_gene = FALSE) %>%
    select(Locus, gene_symbol, gene_id, is_host_gene)
  cat(nrow(window_pairs), "additional window-gene (locus, gene) pairs (+/-", WindowKb, "kb)\n", sep = "")
  test_list <- bind_rows(test_list, window_pairs)
}

needed_gene_ids <- unique(test_list$gene_id)
cat(length(needed_gene_ids), "unique genes needed across all tests\n")

########################################################################################
# Per-tissue loading helpers
########################################################################################

find_expression_file <- function(tissue){
  cand <- c(
    file.path(ExpressionMatrixDir, paste0(tissue, ".v8.normalized_expression.bed")),
    file.path(ExpressionMatrixDir, paste0(tissue, ".v8.normalized_expression.bed.gz"))
  )
  cand[file.exists(cand)][1]
}
find_covariates_file <- function(tissue){
  file.path(PEERDir, paste0(tissue, ".v8.covariates.txt"))
}

# Load only the rows for genes we actually need, for one tissue's expression bed.
load_tissue_expression <- function(tissue, gene_ids_needed){
  f <- find_expression_file(tissue)
  if(is.na(f) || !file.exists(f)) return(NULL)
  bed <- fread(f, sep = "\t", header = TRUE, check.names = FALSE) %>% as.data.frame()
  if(ncol(bed) < 5) return(NULL)
  bed <- bed[, -(1:3), drop = FALSE]
  names(bed)[1] <- "gene_id_raw"
  bed$gene_id <- sub("\\..*$", "", bed$gene_id_raw)
  bed <- bed %>% filter(gene_id %in% gene_ids_needed)
  if(nrow(bed) == 0) return(NULL)
  bed <- bed %>% select(-gene_id_raw) %>% distinct(gene_id, .keep_all = TRUE)
  sample_cols <- setdiff(names(bed), "gene_id")
  # normalize sample IDs to dash format in case fread mangled dashes to dots
  sample_cols_clean <- gsub("\\.", "-", sample_cols)
  names(bed)[names(bed) %in% sample_cols] <- sample_cols_clean
  list(mat = bed, sample_cols = sample_cols_clean)
}

# Load per-tissue covariates file, transposed to SUBJID x covariate, using EVERY
# available column.
load_tissue_covariates <- function(tissue){
  f <- find_covariates_file(tissue)
  if(!file.exists(f)) return(NULL)
  tab <- tryCatch(fread(f, sep = "\t", header = TRUE, check.names = FALSE),
                  error = function(e) NULL)
  if(is.null(tab) || nrow(tab) == 0) return(NULL)
  rn <- as.character(tab[[1]])
  mat <- as.matrix(tab[, -1, with = FALSE])
  rownames(mat) <- rn
  sub_cov <- as.data.frame(t(mat), stringsAsFactors = FALSE)
  sub_cov[] <- lapply(sub_cov, function(x) suppressWarnings(as.numeric(x)))
  sub_cov$SUBJID <- gsub("\\.", "-", rownames(sub_cov))
  rownames(sub_cov) <- NULL
  sub_cov
}

is_brain <- function(x) grepl("Brain", x)

meta_tissues <- function(tissue_df, label){
  use <- tissue_df %>% filter(!is.na(beta) & !is.na(se) & se > 0)
  if(nrow(use) < 2) return(NULL)
  res <- tryCatch(
    metafor::rma(yi = use$beta, sei = use$se, method = "REML", slab = use$SMTSD),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  if(is.null(res)) return(NULL)
  data.frame(
    tissue_group = label,
    k_tissues    = res$k,
    meta_beta    = as.numeric(res$beta),
    meta_se      = as.numeric(res$se),
    meta_ci_lb   = as.numeric(res$ci.lb),
    meta_ci_ub   = as.numeric(res$ci.ub),
    meta_z       = as.numeric(res$zval),
    meta_p       = as.numeric(res$pval),
    Q            = as.numeric(res$QE),
    Q_p          = as.numeric(res$QEp),
    I2           = as.numeric(res$I2),
    tau2         = as.numeric(res$tau2),
    stringsAsFactors = FALSE
  )
}

########################################################################################
# Determine tissues to process: intersection of ExpressionMatrixDir and PEERDir
########################################################################################

expr_files <- list.files(ExpressionMatrixDir, pattern = "\\.v8\\.normalized_expression\\.bed(\\.gz)?$")
tissues_expr <- unique(sub("\\.v8\\.normalized_expression\\.bed(\\.gz)?$", "", expr_files))
cov_files <- list.files(PEERDir, pattern = "\\.v8\\.covariates\\.txt$")
tissues_cov <- unique(sub("\\.v8\\.covariates\\.txt$", "", cov_files))
tissues <- intersect(tissues_expr, tissues_cov)
cat("\nFound", length(tissues), "tissues with both expression and covariates available\n")

########################################################################################
# Main loop: per tissue, residualize needed genes, then test all (locus,gene) pairs
########################################################################################

per_locus_meta_rows <- list()

run_tissue <- function(tissue){
  cat("\n== Tissue:", tissue, "==\n")

  per_sample_this <- list()
  per_tissue_this  <- list()

  expr <- load_tissue_expression(tissue, needed_gene_ids)
  if(is.null(expr)){
    cat("  No needed genes found in this tissue's expression file -- skipping\n")
    return(list(per_sample = per_sample_this, per_tissue = per_tissue_this))
  }
  cov <- load_tissue_covariates(tissue)
  if(is.null(cov)){
    cat("  No covariates file found for this tissue -- skipping\n")
    return(list(per_sample = per_sample_this, per_tissue = per_tissue_this))
  }

  # Assemble the full per-subject covariate table: covariates file columns (all of
  # them) + AGE (+ SEX if not already present) from subject phenotypes + insert size.
  cov_full <- cov
  if(has_subj){
    add_cols <- c("SUBJID")
    if(!("SEX" %in% names(cov_full)) && "SEX" %in% names(subjPheno)) add_cols <- c(add_cols, "SEX")
    if("AGE" %in% names(subjPheno)) add_cols <- c(add_cols, "AGE")
    cov_full <- cov_full %>% left_join(subjPheno[, intersect(add_cols, names(subjPheno)), drop = FALSE], by = "SUBJID")
  }
  if(!is.null(insertSizeDf)){
    join_key <- if("SUBJID" %in% names(insertSizeDf)) "SUBJID" else names(insertSizeDf)[1]
    cov_full <- cov_full %>% left_join(insertSizeDf, by = setNames(join_key, "SUBJID"))
  }

  covar_cols <- setdiff(names(cov_full), "SUBJID")
  covar_cols <- covar_cols[sapply(covar_cols, function(cv){
    v <- cov_full[[cv]]
    is.numeric(v) && length(unique(na.omit(v))) > 1 && mean(!is.na(v)) > 0.5
  })]
  cat("  Using", length(covar_cols), "covariates:", paste(head(covar_cols, 10), collapse=", "),
      ifelse(length(covar_cols) > 10, "...", ""), "\n")

  sample_cols <- intersect(expr$sample_cols, cov_full$SUBJID)
  cat("  ", length(sample_cols), "subjects with both expression and covariates\n", sep="")
  if(length(sample_cols) < MinSamples){
    return(list(per_sample = per_sample_this, per_tissue = per_tissue_this))
  }

  # Residualize every needed gene present in this tissue, once.
  resid_cache <- list()
  for(gid in expr$mat$gene_id){
    expr_vec <- as.numeric(expr$mat[expr$mat$gene_id == gid, sample_cols, drop = TRUE])
    d <- data.frame(SUBJID = sample_cols, expr = expr_vec, stringsAsFactors = FALSE) %>%
      left_join(cov_full, by = "SUBJID")

    use_cols <- c("expr", covar_cols)
    ok <- complete.cases(d[, use_cols, drop = FALSE])
    if(sum(ok) < MinSamples){
      resid_cache[[gid]] <- NULL
      next
    }
    fmla <- as.formula(paste0("expr ~ ", paste(covar_cols, collapse = " + ")))
    fit <- tryCatch(lm(fmla, data = d[ok, ]), error = function(e) NULL)
    if(is.null(fit)){
      resid_cache[[gid]] <- NULL
      next
    }
    d$resid <- NA_real_
    d$resid[ok] <- residuals(fit)
    resid_cache[[gid]] <- d %>% select(SUBJID, resid)
  }

  # Now test every (locus,gene) pair that maps to a gene we successfully residualized.
  pairs_here <- test_list %>% filter(gene_id %in% names(resid_cache)[!sapply(resid_cache, is.null)])
  if(nrow(pairs_here) == 0){
    cat("  No testable (locus,gene) pairs for this tissue\n")
    return(list(per_sample = per_sample_this, per_tissue = per_tissue_this))
  }

  for(j in seq_len(nrow(pairs_here))){
    region <- pairs_here$Locus[j]
    gid    <- pairs_here$gene_id[j]
    isHost <- pairs_here$is_host_gene[j]

    alleles <- data_all %>% filter(Locus == region) %>%
      select(MergeId, Max_numeric) %>% filter(!is.na(Max_numeric)) %>% unique() %>%
      rename(SUBJID = MergeId)

    d <- resid_cache[[gid]] %>% inner_join(alleles, by = "SUBJID") %>% filter(!is.na(resid))
    if(nrow(d) < MinSamples) next
    if(length(unique(d$Max_numeric)) < MinDistinctAlleles) next

    per_sample_this[[length(per_sample_this) + 1]] <- d %>%
      mutate(Locus = region, gene_id = gid, SMTSD = tissue) %>%
      select(Locus, gene_id, SMTSD, SUBJID, resid, Max_numeric)

    fit <- tryCatch(lm(resid ~ Max_numeric, data = d), error = function(e) NULL)
    if(is.null(fit)) next
    co <- summary(fit)$coefficients
    if(!"Max_numeric" %in% rownames(co)) next
    if(isHost) cat("    ", region, " x ", gid, ": association lm used ", nrow(d), " samples\n", sep = "")

    per_tissue_this[[length(per_tissue_this) + 1]] <- data.frame(
      Locus         = region,
      gene_id       = gid,
      is_host_gene  = isHost,
      SMTSD         = tissue,
      tissue_class  = ifelse(is_brain(tissue), "brain", "nonbrain"),
      n             = nrow(d),
      n_distinct    = length(unique(d$Max_numeric)),
      beta          = co["Max_numeric","Estimate"],
      se            = co["Max_numeric","Std. Error"],
      p             = co["Max_numeric","Pr(>|t|)"],
      r_squared     = summary(fit)$r.squared,
      stringsAsFactors = FALSE
    )
  }
  cat("  Tested", nrow(pairs_here), "(locus,gene) pairs in", tissue, "\n")

  list(per_sample = per_sample_this, per_tissue = per_tissue_this)
}

tissue_results <- if(Threads > 1){
  mclapply(tissues, run_tissue, mc.cores = Threads, mc.preschedule = FALSE)
} else {
  lapply(tissues, run_tissue)
}
names(tissue_results) <- tissues

failed <- sapply(tissue_results, function(x) inherits(x, "try-error"))
if(any(failed)){
  cat("\nWARNING: worker failed for tissue(s):", paste(names(tissue_results)[failed], collapse = ", "), "\n")
  for(t in names(tissue_results)[failed]) cat("  ", t, ": ", attr(tissue_results[[t]], "condition")$message, "\n", sep = "")
  tissue_results <- tissue_results[!failed]
}

per_sample_long <- unlist(lapply(tissue_results, `[[`, "per_sample"), recursive = FALSE)
per_tissue_rows <- unlist(lapply(tissue_results, `[[`, "per_tissue"), recursive = FALSE)

per_sample_df <- bind_rows(per_sample_long) %>% add_gene_name()
per_tissue_df <- bind_rows(per_tissue_rows) %>% add_gene_name()

########################################################################################
# Cross-tissue meta-analysis for HOST-GENE pairs only
########################################################################################

if(nrow(per_tissue_df) > 0){
  host_tissue_df <- per_tissue_df %>% filter(is_host_gene)
  for(region in unique(host_tissue_df$Locus)){
    sub_df <- host_tissue_df %>% filter(Locus == region)
    gid <- sub_df$gene_id[1]
    brain_meta    <- meta_tissues(sub_df %>% filter(tissue_class == "brain"),    "brain")
    nonbrain_meta <- meta_tissues(sub_df %>% filter(tissue_class == "nonbrain"), "nonbrain")
    all_meta      <- meta_tissues(sub_df,                                         "all")
    for(res in list(brain_meta, nonbrain_meta, all_meta)){
      if(is.null(res)) next
      res$Locus   <- region
      res$gene_id <- gid
      per_locus_meta_rows[[length(per_locus_meta_rows) + 1]] <- res
    }
  }
}
per_locus_meta_df <- bind_rows(per_locus_meta_rows)

########################################################################################
# Per-tissue Bonferroni + FDR (primary correction: N = number of tests within that
# single tissue after the MinSamples filter, Bonferroni = p * N, FDR = p.adjust(p,
# "fdr") scoped to that tissue's own test set).
########################################################################################

if(nrow(per_tissue_df) > 0){
  per_tissue_df <- per_tissue_df %>%
    group_by(SMTSD) %>%
    mutate(N_tissue      = n(),
           bonferroni_tissue = pmin(p * N_tissue, 1),
           q_fdr_tissue  = p.adjust(p, method = "fdr")) %>%
    ungroup() %>%
    mutate(significant_QTL = q_fdr_tissue < FDRThreshold)
}

per_locus_summary <- NULL
if(nrow(per_tissue_df) > 0){
  per_locus_summary <- per_tissue_df %>%
    group_by(Locus, gene_id, gene_name, is_host_gene) %>%
    summarise(
      n_tissues_tested      = n(),
      n_tissues_significant  = sum(significant_QTL, na.rm = TRUE),
      any_significant        = n_tissues_significant > 0,
      min_q_tissue           = min(q_fdr_tissue, na.rm = TRUE),
      best_tissue            = SMTSD[which.min(q_fdr_tissue)],
      beta_at_best           = beta[which.min(q_fdr_tissue)],
      p_at_best              = p[which.min(q_fdr_tissue)],
      r_squared_at_best      = r_squared[which.min(q_fdr_tissue)],
      .groups = "drop"
    ) %>%
    arrange(min_q_tissue)
}

if(nrow(per_locus_meta_df) > 0){
  per_locus_meta_df <- per_locus_meta_df %>% add_gene_name() %>% group_by(tissue_group) %>%
    mutate(meta_p_fdr = p.adjust(meta_p, method = "fdr")) %>% ungroup() %>%
    select(Locus, gene_id, gene_name, tissue_group, k_tissues,
           meta_beta, meta_se, meta_ci_lb, meta_ci_ub, meta_z,
           meta_p, meta_p_fdr, Q, Q_p, I2, tau2)
}

########################################################################################
# Output
########################################################################################

suffix <- ".eQTL"

fwrite(per_sample_df,
       file = paste0(outdir, COHORT, suffix, ".per_sample_resid_max.txt"),
       sep = "\t", row.names = FALSE)
fwrite(per_tissue_df %>% rename(`Slope (repeat units)` = beta, `R-squared` = r_squared),
       file = paste0(outdir, COHORT, suffix, ".per_tissue_lm.txt"),
       sep = "\t", row.names = FALSE)
if(!is.null(per_locus_summary)){
  fwrite(per_locus_summary,
         file = paste0(outdir, COHORT, suffix, ".per_locus_anyTissueQTL.txt"),
         sep = "\t", row.names = FALSE)
}
if(nrow(per_locus_meta_df) > 0){
  fwrite(per_locus_meta_df,
         file = paste0(outdir, COHORT, suffix, ".hostGene_crossTissueMeta.txt"),
         sep = "\t", row.names = FALSE)
}

cat("\n========== Significant TR:gene eQTLs (FDR q<", FDRThreshold, ", any tissue) ==========\n", sep="")
if(!is.null(per_locus_summary)){
  print(per_locus_summary %>% filter(any_significant) %>%
        select(Locus, gene_id, gene_name, is_host_gene, n_tissues_tested, n_tissues_significant,
               min_q_tissue, best_tissue, beta_at_best, p_at_best) %>% as.data.frame())
}

cat("\n========== Host-gene cross-tissue meta-analysis (all) ==========\n")
if(nrow(per_locus_meta_df) > 0){
  print(per_locus_meta_df %>% filter(tissue_group == "all") %>%
        arrange(meta_p) %>%
        select(Locus, gene_id, gene_name, k_tissues, meta_beta, meta_se, meta_p, meta_p_fdr, I2, Q_p) %>% as.data.frame())
}

cat("\n========== Host-gene cross-tissue meta-analysis (brain) ==========\n")
if(nrow(per_locus_meta_df) > 0){
  print(per_locus_meta_df %>% filter(tissue_group == "brain") %>%
        arrange(meta_p) %>%
        select(Locus, gene_id, gene_name, k_tissues, meta_beta, meta_se, meta_p, meta_p_fdr, I2, Q_p) %>% as.data.frame())
}

cat("\nOutputs written to:\n")
cat(" ", paste0(outdir, COHORT, suffix, ".per_sample_resid_max.txt"), "\n")
cat(" ", paste0(outdir, COHORT, suffix, ".per_tissue_lm.txt"),        "\n")
if(!is.null(per_locus_summary)){
  cat(" ", paste0(outdir, COHORT, suffix, ".per_locus_anyTissueQTL.txt"), "\n")
}
if(nrow(per_locus_meta_df) > 0){
  cat(" ", paste0(outdir, COHORT, suffix, ".hostGene_crossTissueMeta.txt"), "\n")
}

sessionInfo()
