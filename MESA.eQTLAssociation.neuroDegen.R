# =============================================================================
# MESA TR-eQTL Association — Neurodegenerative STRs
# =============================================================================
# Gabrielle N. Altman, adapted from code by Mariya Shadrina, Celine A. Manigbas,
# Alejandro Martin-Trujillo
#
# TR-eQTL association testing in MESA (Monocytes, PBMCs, T-cells): per-tissue
# expression residualization followed by long-allele association testing at
# host genes and genes within a window of each TR locus. Similar script used
# for GTEx and AMP-PD; see GTEx.eQTLAssociation.neuroDegen.R for the GTEx
# version.
#
# Method:
#   1. Per tissue: load gene TPM matrix + WGS/RNA-seq metadata, rank-based
#      inverse-normal transform (INT) each needed gene's expression, then
#      residualize on sex, analyte isolation batch, and predicted ancestry
#      (lm(expr ~ covariates), keep residuals). Genes are only tested in a
#      tissue if median TPM there clears --MinTPM.
#   2. Per locus x gene x tissue: expand each TR locus +/- WindowKb, overlap
#      against a gene annotation BED for host-gene + window-gene pairs, then
#      test lm(residual ~ long_allele_length), gated on --MinSamples and
#      --MinDistinctAlleles. A leverage check refits after dropping the single
#      largest-allele sample (beta_drop_top / p_drop_top / delta_beta_drop_top).
#   3. Per-tissue Bonferroni + FDR correction.
#
# Notes:
#   - gene_id/gene_name are resolved directly from each tissue's own
#     expression matrix (Name/Description columns); genes ambiguous across
#     tissues are excluded from the host/window-gene lookup.
#   - --GeneAnnotationFile is a headerless BED (chrom, start, end, gene_id by
#     column position, not header name).
#   - Host genes come from the ANNOVAR closest-gene annotation at each locus
#     and are always tested regardless of distance; window genes are added
#     only within +/-WindowKb.
#
# Usage:
#   Rscript MESA.eQTLAssociation.neuroDegen.R \
#     --Dir /path/to/dir --RegionsFile /path/to/regions.txt --FolderName neuroDegen \
#     --GenotypeDir /path/to/EH_output --MetadataDir /path/to/metadata \
#     --ExpressionDir /path/to/RNASeq --AnnotationFile /path/to/annotation.tsv
#
# Required: --Dir --RegionsFile --FolderName --GenotypeDir --MetadataDir
#           --ExpressionDir --AnnotationFile
#
# Key optional: --WindowKb (100) --MinTPM (0.1) --MinDistinctAlleles (5)
#   --MinSamples (10) --FDRThreshold (0.1) --GeneAnnotationFile. See --help
#   for the full list.
#
# Expected input layout:
#   --GenotypeDir/{chr}_MESA_EHv5_GenomewidePolymorphic_ForPhewasGWAS.tsv.gz
#   --MetadataDir/MESA_WGS_meta_genotyped.txt
#   --MetadataDir/MESA_RNASeq_meta_{Mono,PBMC,TCell}.txt
#   --ExpressionDir/MESA.RNASeq.{Monocytes,PBMCs,TCell}.gene_tpm.ExpressedGenes.FINAL.gct.gz
#
# Dependencies: argparser, data.table, dplyr, tidyr, tibble
#
# HPC example (LSF):
#   bsub -P acc_PROJECTID -q premium -W 10:00 -R "span[hosts=1]" -R "rusage[mem=60000]" \
#     Rscript MESA.eQTLAssociation.neuroDegen.R --Dir ... --RegionsFile ... --FolderName neuroDegen \
#       --GenotypeDir ... --MetadataDir ... --ExpressionDir ... --AnnotationFile ...
# =============================================================================

.libPaths(c("~/.Rlib", .libPaths()))

cat("Loading libraries\n")
suppressWarnings(suppressMessages({
  library(argparser)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
}))

cat(paste0("\n", Sys.Date(), "\n"))
cat(paste0("R version: ", getRversion(), "\n"))

########################################################################################
# Arguments
########################################################################################

p <- arg_parser("MESA TR-eQTL analysis (long allele, window scan)")
p <- add_argument(p, "--Dir", help = "Main project directory (output written here) [required]")
p <- add_argument(p, "--RegionsFile", help = "File of Locus values (with chr/start/end encoded as chr_start_end) to test [required]")
p <- add_argument(p, "--FolderName", help = "Output subfolder name, or 'no' for none [required]")
p <- add_argument(p, "--GenotypeDir",
                  help = "Directory containing per-chromosome ExpansionHunter output ({chr}_MESA_EHv5_GenomewidePolymorphic_ForPhewasGWAS.tsv.gz) [required]")
p <- add_argument(p, "--MetadataDir",
                  help = "Directory with MESA_WGS_meta_genotyped.txt and MESA_RNASeq_meta_{Mono,PBMC,TCell}.txt [required]")
p <- add_argument(p, "--ExpressionDir",
                  help = "Directory with MESA.RNASeq.{Monocytes,PBMCs,TCell}.gene_tpm.ExpressedGenes.FINAL.gct.gz [required]")
p <- add_argument(p, "--AnnotationFile",
                  help = "ANNOVAR annotation file with Chr/Start/End/Gene.ensGene/Func.ensGene columns, used to assign each locus its host gene [required]")
p <- add_argument(p, "--GeneAnnotationFile", default = "no",
                  help = "Headerless BED (chrom, start, end, gene_id) for the window scan. 'no' disables the window scan (host gene only).")
p <- add_argument(p, "--WindowKb", default = 100,
                  help = "Window (kb) each side of the TR locus to search for additional genes.")
p <- add_argument(p, "--MinDistinctAlleles", default = 5,
                  help = "Minimum distinct long-allele values in a tissue to fit the model.")
p <- add_argument(p, "--MinTPM", default = 0.1,
                  help = "Minimum median TPM in a tissue to test a gene there.")
p <- add_argument(p, "--MinSamples", default = 10,
                  help = "Minimum samples per (locus,gene,tissue) to run the test.")
p <- add_argument(p, "--FDRThreshold", default = 0.1,
                  help = "FDR (q-value) threshold for calling a significant TR:gene eQTL.")
argv <- parse_args(p)

if (is.na(argv$Dir) || is.na(argv$RegionsFile) || is.na(argv$FolderName) ||
    is.na(argv$GenotypeDir) || is.na(argv$MetadataDir) ||
    is.na(argv$ExpressionDir) || is.na(argv$AnnotationFile)) {
  print(p)
  stop("ERROR: --Dir, --RegionsFile, --FolderName, --GenotypeDir, --MetadataDir, --ExpressionDir, and --AnnotationFile are all required.", call. = FALSE)
}

Dir                <- argv$Dir
RegionsFile        <- argv$RegionsFile
FolderName         <- argv$FolderName
GenotypeDir        <- argv$GenotypeDir
MetadataDir        <- argv$MetadataDir
ExpressionDir      <- argv$ExpressionDir
AnnotationFile     <- argv$AnnotationFile
GeneAnnotationFile <- argv$GeneAnnotationFile
WindowKb           <- as.numeric(argv$WindowKb)
MinDistinctAlleles <- as.integer(argv$MinDistinctAlleles)
MinTPM             <- as.numeric(argv$MinTPM)
MinSamples         <- as.integer(argv$MinSamples)
FDRThreshold       <- as.numeric(argv$FDRThreshold)

ehBaseDir  <- GenotypeDir
metaDir    <- MetadataDir
rnaseqDir  <- ExpressionDir
annotationPath <- AnnotationFile

tissues <- list(
  Monocytes = list(
    expr = file.path(rnaseqDir, "MESA.RNASeq.Monocytes.gene_tpm.ExpressedGenes.FINAL.gct.gz"),
    meta = file.path(metaDir,   "MESA_RNASeq_meta_Mono.txt")
  ),
  PBMCs = list(
    expr = file.path(rnaseqDir, "MESA.RNASeq.PBMCs.gene_tpm.ExpressedGenes.FINAL.gct.gz"),
    meta = file.path(metaDir,   "MESA_RNASeq_meta_PBMC.txt")
  ),
  TCell = list(
    expr = file.path(rnaseqDir, "MESA.RNASeq.TCell.gene_tpm.ExpressedGenes.FINAL.gct.gz"),
    meta = file.path(metaDir,   "MESA_RNASeq_meta_TCell.txt")
  )
)

if (FolderName != "no") {
  outdir <- file.path(Dir, "geneExpression", FolderName, "eqtlAnalysis")
} else {
  outdir <- file.path(Dir, "geneExpression", "eqtlAnalysis")
}
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

cat("\nArguments:\n")
cat("  Dir:                ", Dir, "\n")
cat("  Genotype dir:       ", GenotypeDir, "\n")
cat("  Metadata dir:       ", MetadataDir, "\n")
cat("  Expression dir:     ", ExpressionDir, "\n")
cat("  Annotation file:    ", AnnotationFile, "\n")
cat("  Gene annotation:    ", GeneAnnotationFile, "\n")
cat("  Window (kb):        ", WindowKb, "\n")
cat("  MinDistinctAlleles: ", MinDistinctAlleles, "\n")
cat("  MinTPM:             ", MinTPM, "\n")
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

  candidate_covs <- c("AGE_AT_COLLECTION", "SEX",
                      "ANALYTE_ISOLATION_BATCH_ID", "PredictedAncestry")
  covars <- candidate_covs[candidate_covs %in% names(d)]
  covars <- covars[sapply(covars, function(cv) {
    v <- d[[cv]]
    length(unique(na.omit(v))) > 1 && mean(!is.na(v)) > 0.5
  })]

  for (cv in intersect(covars, c("SEX", "ANALYTE_ISOLATION_BATCH_ID", "PredictedAncestry")))
    d[[cv]] <- as.factor(d[[cv]])

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
  if (length(unique(d$LongAllele)) < MinDistinctAlleles) return(NULL)

  fit <- tryCatch(lm(z ~ LongAllele, data = d), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  co <- summary(fit)$coefficients
  if (!"LongAllele" %in% rownames(co)) return(NULL)

  beta <- co["LongAllele", "Estimate"]
  se   <- co["LongAllele", "Std. Error"]
  pval <- co["LongAllele", "Pr(>|t|)"]
  r_squared <- summary(fit)$r.squared

  max_val <- max(d$LongAllele, na.rm = TRUE)
  top_idx <- which(d$LongAllele == max_val)
  d_drop  <- d[-top_idx[1], , drop = FALSE]
  beta_drop <- NA_real_; p_drop <- NA_real_
  if (nrow(d_drop) >= MinSamples && length(unique(d_drop$LongAllele)) >= MinDistinctAlleles) {
    fit2 <- tryCatch(lm(z ~ LongAllele, data = d_drop), error = function(e) NULL)
    if (!is.null(fit2)) {
      co2 <- summary(fit2)$coefficients
      if ("LongAllele" %in% rownames(co2)) {
        beta_drop <- co2["LongAllele", "Estimate"]
        p_drop    <- co2["LongAllele", "Pr(>|t|)"]
      }
    }
  }

  list(
    n = nrow(d), max_allele = max_val, n_at_max = length(top_idx),
    n_distinct          = length(unique(d$LongAllele)),
    beta                = beta, se = se, p = pval, r_squared = r_squared,
    beta_drop_top       = beta_drop, p_drop_top = p_drop,
    delta_beta_drop_top = beta_drop - beta
  )
}

########################################################################################
# Load regions + EH allele data
########################################################################################

cat("\nLoading regions ...\n")
regionsToKeep <- fread(RegionsFile, sep = "\t", check.names = FALSE, header = TRUE)
loci_vec      <- regionsToKeep$Locus
chroms_needed <- unique(sub("_.*", "", loci_vec))
cat(length(loci_vec), "loci across", length(chroms_needed), "chromosomes\n")

cat("\nLoading MESA EH allele lengths ...\n")
data_EH <- rbindlist(lapply(chroms_needed, function(chr) {
  f <- file.path(ehBaseDir,
                 paste0(chr, "_MESA_EHv5_GenomewidePolymorphic_ForPhewasGWAS.tsv.gz"))
  if (!file.exists(f)) { cat("  Missing:", f, "\n"); return(NULL) }
  dt <- fread(f, select = c("SampleId", "chrom", "start", "end", "LongAllele"),
              showProgress = FALSE)
  dt[, VARID := paste0(chrom, "_", start, "_", end)]
  dt[VARID %in% loci_vec, .(SampleId, VARID, LongAllele)]
}))
data_EH <- data_EH %>% mutate(LongAllele = suppressWarnings(as.numeric(LongAllele)))
cat("EH data:", nrow(data_EH), "rows,", n_distinct(data_EH$VARID), "VARIDs,",
    n_distinct(data_EH$SampleId), "WGS samples\n")

cat("\nLoading WGS metadata (SampleId -> SubjectID + SEX + PredictedAncestry) ...\n")
meta_WGS <- fread(file.path(metaDir, "MESA_WGS_meta_genotyped.txt"),
                  sep = "\t", header = TRUE, check.names = FALSE) %>%
  select(SampleID, SubjectID, SEX, PredictedAncestry) %>%
  rename(SampleId = SampleID) %>%
  mutate(SubjectID = as.character(SubjectID))

cat("\nLoading gene annotation (closest gene per locus) ...\n")
annotation <- fread(annotationPath, sep = "\t", check.names = FALSE, header = TRUE)
annotation$VARID <- paste0(annotation$Chr, "_", annotation$Start, "_", annotation$End)

########################################################################################
# Load RNA-seq + metadata per tissue (once, outside main loop)
########################################################################################

cat("\nLoading RNA-seq data for all tissues ...\n")
tissue_data <- list()
gene_lookup_rows <- list()

for (tname in names(tissues)) {
  cat("  Loading", tname, "...\n")

  expr_mat  <- fread(tissues[[tname]]$expr, sep = "\t", check.names = FALSE,
                     header = TRUE, showProgress = FALSE)
  ann_cols  <- c("Name", "Description")
  samp_cols <- setdiff(names(expr_mat), ann_cols)

  tmeta <- fread(tissues[[tname]]$meta, sep = "\t", header = TRUE,
                 check.names = FALSE, showProgress = FALSE) %>%
    select(SampleID, SubjectID, AGE_AT_COLLECTION, ANALYTE_ISOLATION_BATCH_ID) %>%
    rename(RNAseq_col = SampleID) %>%
    mutate(SubjectID = as.character(SubjectID))

  samp_cols <- intersect(samp_cols, tmeta$RNAseq_col)

  gene_lookup_rows[[tname]] <- expr_mat %>%
    select(Name, Description) %>%
    rename(gene_id = Name, gene_name = Description) %>%
    mutate(gene_id = sub("\\..*$", "", gene_id))

  tissue_data[[tname]] <- list(
    expr      = expr_mat,
    samp_cols = samp_cols,
    meta      = tmeta
  )
  cat("   ", length(samp_cols), "samples with expression and metadata\n")
}

# gene_id <-> gene_name lookup built directly from MESA's own expression matrices 
gene_lookup <- bind_rows(gene_lookup_rows) %>% distinct()
dup_by_name <- gene_lookup %>% count(gene_name) %>% filter(n > 1) %>% pull(gene_name)
dup_by_id   <- gene_lookup %>% count(gene_id)   %>% filter(n > 1) %>% pull(gene_id)
if (length(dup_by_name) > 0) {
  cat("Note:", length(dup_by_name), "gene_name values are ambiguous across tissues' expression matrices; excluded from host/window-gene lookup.\n")
}
gene_lookup <- gene_lookup %>% filter(!(gene_name %in% dup_by_name), !(gene_id %in% dup_by_id))

########################################################################################
# Build the host-gene + window-gene test list
########################################################################################

host_pairs <- annotation %>%
  select(VARID, Gene.ensGene, Func.ensGene) %>%
  rename(Locus = VARID, gene_name = Gene.ensGene, genomicFeature = Func.ensGene) %>%
  filter(Locus %in% loci_vec) %>%
  distinct() %>%
  left_join(gene_lookup, by = "gene_name") %>%
  mutate(is_host_gene = TRUE)

n_no_expr <- sum(is.na(host_pairs$gene_id))
if (n_no_expr > 0) {
  cat("WARNING:", n_no_expr, "host-gene locus/gene pairs have no matching row in any tissue's expression matrix (ambiguous or missing gene_name) -- will be skipped:\n")
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
# Residualize every needed gene once per tissue (host + window genes only)
########################################################################################

cat("\nResidualizing", length(needed_gene_ids), "needed genes across", length(tissues), "tissues ...\n")
resid_cache <- list()
for (tname in names(tissues)) {
  td <- tissue_data[[tname]]
  resid_cache[[tname]] <- list()

  for (gid in needed_gene_ids) {
    exprRow <- td$expr %>% filter(sub("\\..*$", "", Name) == gid)
    if (nrow(exprRow) != 1) next

    expr_vals <- as.numeric(exprRow[1, td$samp_cols, with = FALSE])
    if (median(expr_vals, na.rm = TRUE) < MinTPM) next

    exprSub <- data.frame(
      RNAseq_col = td$samp_cols,
      expression = expr_vals,
      stringsAsFactors = FALSE
    ) %>% inner_join(td$meta, by = "RNAseq_col")

    if (nrow(exprSub) == 0) next
    exprSub <- compute_z(exprSub)
    resid_cache[[tname]][[gid]] <- exprSub %>% select(RNAseq_col, SubjectID, expression, z)
  }
  cat("  ", tname, ": residualized", sum(!sapply(resid_cache[[tname]], is.null)),
      "of", length(needed_gene_ids), "needed genes.\n")
}

########################################################################################
# Main loop: test every (locus, gene, tissue) triple that maps to a gene 
# residualized in that tissue
########################################################################################

per_sample_long <- list()
per_tissue_rows <- list()

for (i in seq_len(nrow(test_list))) {
  region         <- test_list$Locus[i]
  gid            <- test_list$gene_id[i]
  gname          <- test_list$gene_name[i]
  isHost         <- test_list$is_host_gene[i]
  genomicFeature <- test_list$genomicFeature[i]

  alleles <- data_EH %>%
    filter(VARID == region, !is.na(LongAllele)) %>%
    select(SampleId, LongAllele) %>%
    inner_join(meta_WGS, by = "SampleId") %>%
    select(SubjectID, LongAllele) %>%
    distinct()
  if (nrow(alleles) == 0) next

  tissue_results <- list()

  for (tname in names(tissues)) {
    if (is.null(resid_cache[[tname]][[gid]])) next

    d <- resid_cache[[tname]][[gid]] %>% inner_join(alleles, by = "SubjectID") %>% filter(!is.na(z))
    if (nrow(d) == 0) next

    per_sample_long[[length(per_sample_long) + 1]] <- d %>%
      mutate(Locus = region, gene_id = gid, gene_name = gname, is_host_gene = isHost, tissue = tname) %>%
      select(Locus, gene_id, gene_name, is_host_gene, tissue, RNAseq_col, SubjectID, expression, z, LongAllele)

    fit_res <- fit_model(d)
    if (is.null(fit_res)) {
      cat(sprintf("[%2d] %-25s | %-10s gene=%-12s skipped (n=%d, n_distinct=%d)\n",
                  i, region, tname, gname, nrow(d), length(unique(d$LongAllele))))
      next
    }
    if (isHost) cat("    ", region, " x ", gid, " (", gname, ") | ", tname, ": association lm used ", fit_res$n, " samples\n", sep = "")

    tissue_results[[length(tissue_results) + 1]] <- data.frame(
      Locus                = region,
      gene_id              = gid,
      gene_name            = gname,
      is_host_gene         = isHost,
      genomicFeature       = genomicFeature,
      tissue               = tname,
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

  if (length(tissue_results) == 0) next
  per_tissue_rows[[length(per_tissue_rows) + 1]] <- bind_rows(tissue_results)
}

per_sample_df <- bind_rows(per_sample_long)
per_tissue_df <- bind_rows(per_tissue_rows)

########################################################################################
# Per-tissue Bonferroni + FDR (primary correction: N = tests within that tissue)
########################################################################################

if (nrow(per_tissue_df) > 0) {
  per_tissue_df <- per_tissue_df %>%
    group_by(tissue) %>%
    mutate(N_tissue          = n(),
           bonferroni_tissue = pmin(p * N_tissue, 1),
           q_fdr_tissue      = p.adjust(p, method = "fdr")) %>%
    ungroup() %>%
    mutate(significant_QTL = q_fdr_tissue < FDRThreshold) %>%
    arrange(p)
}

########################################################################################
# Output
########################################################################################

suffix <- ".eQTL"

fwrite(per_sample_df,
       file = file.path(outdir, paste0("MESA", suffix, ".per_sample_z_longAllele.txt")),
       sep = "\t", row.names = FALSE)
fwrite(per_tissue_df %>% rename(`Slope (repeat units)` = beta, `R-squared` = r_squared),
       file = file.path(outdir, paste0("MESA", suffix, ".per_tissue_lm.txt")),
       sep = "\t", row.names = FALSE)

cat("\n========== Significant TR:gene eQTLs (per-tissue FDR q<", FDRThreshold, ") ==========\n", sep = "")
if (nrow(per_tissue_df) > 0) {
  print(per_tissue_df %>% filter(significant_QTL) %>%
    select(Locus, gene_id, gene_name, is_host_gene, tissue, n, beta, p, q_fdr_tissue))
}

cat("\nOutputs written to:\n")
cat(" ", file.path(outdir, paste0("MESA", suffix, ".per_sample_z_longAllele.txt")), "\n")
cat(" ", file.path(outdir, paste0("MESA", suffix, ".per_tissue_lm.txt")), "\n")
