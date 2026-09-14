# Neurodegenerative disease TRE association analysis

This code is part of the following preprint:

> Altman GN, Jadhav B, Garg P, Shadrina M, Manigbas CA, Lee W, Kandoi S, Martin-Trujillo A, Sharp AJ. **Tandem repeat expansions in *DAPK1*, *ANK3*, and *RPL14* are associated with diverse neurodegenerative diseases.** medRxiv 2026.08.06.26358503; doi: [10.64898/2026.08.06.26358503](https://doi.org/10.64898/2026.08.06.26358503) *(currently under review)*

Scripts included, in the order the corresponding analyses appear in the Methods:

1. **TR expansion association via REGENIE** (`runREGENIE.neuroDegen.sh`) — association test using REGENIE's two-step Firth approx logistic regression, applied to binary TR genotype data in UKB. Similar script used in AoU.
2. **Cross-cohort meta-analysis via METAL** (`runMETAL.neuroDegen.sh`) — fixed-effects meta-analysis of REGENIE association results across cohorts (e.g. AoU + UKB).
3. **Allele-size association** (`alleleSizeAssociation.neuroDegen.R`) — per-locus Firth logistic regression testing association of varying expansion size cutoff (≥ 99th-percentile repeat length) with disease status, run on UKB WGS TR genotype data. Similar script used in AoU.
4. **Causal variant analysis** (`causal_variant_analysis/`) — fine-mapping (CAVIAR) and conditional analysis (REGENIE) pipeline for identifying whether TR signals are independent of nearby SNPs, run in the *All of Us* Researcher Workbench. Similar pipeline used in UKB.
5. **GTEx TR-eQTL analysis** (`GTEx.eQTLAssociation.neuroDegen.R`) — per-tissue expression residualization and long-allele eQTL testing (host gene + window scan) in GTEx v8, with cross-tissue meta-analysis.
6. **MESA TR-eQTL analysis** (`MESA.eQTLAssociation.neuroDegen.R`) — per-tissue (Monocytes, PBMCs, T-cells) expression residualization and long-allele eQTL testing (host gene + window scan) in MESA.
7. **AMP-PD TR-eQTL analysis** (`AMPPD.eQTLAssociation.neuroDegen.R`) — single-tissue (blood) expression residualization and long-allele eQTL testing (host gene + window scan) in AMP-PD.

---

## Repository layout

```
TRE_neuroDegen/
├── runREGENIE.neuroDegen.sh                    # TR expansion association via REGENIE (bash)
├── runMETAL.neuroDegen.sh                      # Cross-cohort METAL meta-analysis (bash)
├── alleleSizeAssociation.neuroDegen.R          # Per-locus allele-size association (R)
├── causal_variant_analysis/
│   ├── TR_Causal_Variant_Analysis.ipynb        # End-to-end causal variant pipeline (Jupyter)
│   ├── config.sh                               # All paths, binaries, covariate names
│   ├── common.sh                               # Shared bash functions
│   └── README.md                               # Causal variant pipeline details
├── GTEx.eQTLAssociation.neuroDegen.R           # GTEx TR-eQTL analysis (R)
├── MESA.eQTLAssociation.neuroDegen.R           # MESA TR-eQTL analysis (R)
└── AMPPD.eQTLAssociation.neuroDegen.R          # AMP-PD TR-eQTL analysis (R)
```

---

## 1 · TR expansion association via REGENIE (`runREGENIE.neuroDegen.sh`)

### What it does

Runs REGENIE step 1 (whole-genome regression null model) and step 2 (association testing) for TR expansion in UKB WGS data. Binary trait logistic regression with Firth approx correction. Designed to be run once per phenotype × expansion cutoff combination.

### Usage

```bash
bash runREGENIE.neuroDegen.sh <PHENO> <THRESHOLD>
```

| Argument | Description |
|---|---|
| `PHENO` | Phenotype name matching the phenotype file column, e.g. `neuroDegen` |
| `THRESHOLD` | Expansion cutoff label matching input file names, e.g. `Cutoff99` |

### Paths to configure

Edit the `Configuration` block in the script:

| Variable | Description |
|---|---|
| `INDIR` | Root directory of WGS genotype data (plink2 `.pgen/.psam/.pvar` files) |
| `OUTDIR` | Output directory for REGENIE results |
| `PGEN` | Plink2 binary TR genotype file prefix |
| `PHENO_FILE` | Tab-delimited phenotype file (`FID IID <PHENO>`) |
| `SAMPLE_FILE` | Two-column `FID IID` sample keep list |
| `COVAR_FILE` | Tab-delimited covariate file |
| `TR_LIST` | One-column variant ID extract list |

Covariate columns used: `Insert_Size, Age, Age_sq, SNP_PC1–SNP_PC5` (quantitative); `Gender, SeqCenter` (categorical).

Run on DNAnexus.

### Dependencies

`regenie/3.4.1` (loaded via environment modules: `ml regenie/3.4.1`)

---

## 2 · Cross-cohort meta-analysis via METAL (`runMETAL.neuroDegen.sh`)

### What it does

Runs METAL to combine TR expansion association results (REGENIE step 2 output, e.g. from `runREGENIE.neuroDegen.sh`) across two cohorts using fixed-effects (inverse-variance-weighted) meta-analysis. One submission per phenotype × expansion cutoff combination.

### Usage

```bash
bash runMETAL.neuroDegen.sh <PHENO> <THRESHOLD>
```

| Argument | Description |
|---|---|
| `PHENO` | Phenotype name matching the REGENIE output file naming, e.g. `neuroDegen` |
| `THRESHOLD` | Expansion cutoff label matching input file names, e.g. `Cutoff99` |

### Paths to configure

Edit the `Configuration` block in the script:

| Variable | Description |
|---|---|
| `OUTDIR` | Output directory for METAL results and the generated METAL script |
| `COHORT1_FILE` | REGENIE step 2 output for the first cohort |
| `COHORT2_FILE` | REGENIE step 2 output for the second cohort |
| `COHORT1_PVAL_COL` | P-value column name in `COHORT1_FILE` (default `P`) |
| `COHORT2_PVAL_COL` | P-value column name in `COHORT2_FILE` (default `PVAL`) |

Both input files are expected to have `ID`, `N`, `ALLELE0`, `ALLELE1`, and `BETA` columns (standard REGENIE step 2 output); only the p-value column name is configurable, since it can differ if one file was post-processed to rename `LOG10P` to a linear p-value.

### Output

Written to `<OUTDIR>/METAANALYSIS.<THRESHOLD>.<PHENO>..txt`, plus the generated METAL input script at `<OUTDIR>/METAL.<THRESHOLD>.<PHENO>.txt`.

### HPC example (LSF)

```bash
bsub -P acc_PROJECTID -L /bin/bash -q express -n 1 \
     -R rusage[mem=8000] -W 4:00 \
     bash runMETAL.neuroDegen.sh neuroDegen Cutoff99
```

### Dependencies

`metal` (loaded via environment modules: `ml metal/2018-08-28`)

---

## 3 · Allele-size association (`alleleSizeAssociation.neuroDegen.R`)

### What it does

For each tandem repeat locus that was nominally significant in a METAL meta-analysis, tests association between long-allele size (binary variable at each repeat-length cutoff ≥ the 99th percentile) and neurodegenerative disease status using Firth logistic regression (`brglm`). Run per chromosome.

### Usage

```bash
Rscript alleleSizeAssociation.neuroDegen.R \
    --dir        /path/to/project/dir \
    --cohort     neuroDegen \
    --chrom      chr12 \
    --pheno_file /path/to/phenos.txt \
    --covar_file /path/to/covariates.tsv \
    --metal_file /path/to/metal_results.xlsx
```

### Required arguments

| Argument | Description |
|---|---|
| `--dir` | Main project directory (output written here) |
| `--cohort` | Phenotype label, e.g. `neuroDegen` or `neuroDegen_noPDorAD` |
| `--chrom` | Chromosome, e.g. `chr12` |
| `--pheno_file` | Tab-delimited phenotype file (`FID IID Phenotype`, binary 0/1) |
| `--covar_file` | Tab-delimited covariate file (`IID SNP_PC1..5 Age Age_sq Insert_Size predicted_gender SeqCenter`) |
| `--metal_file` | METAL results `.xlsx` with a `Locus` column |

### Input file paths to configure

Edit the `Configuration` block near the top of the script and set `EH_DIR` to the directory containing per-chromosome long-allele matrices and classifier QC files.

Expected input structure under `EH_DIR`:

```
EH_DIR/
├── SC/LongAlleleMatrix/      chrN_SC_EUR_LongAlleleMatrix.tsv.gz
├── SC/ClassifierQC/          chrN_SC_EUR_Cutoff95_ClassifierPredictions.tsv.gz
├── deCODE/LongAlleleMatrix/  chrN_deCODE_EUR_LongAlleleMatrix.tsv.gz
└── deCODE/ClassifierQC/      chrN_deCODE_EUR_Cutoff95_ClassifierPredictions.tsv.gz
```

### Output

Written to `<dir>/alleleSizeAssociation/`:

| File | Contents |
|---|---|
| `<cohort>.<chrom>.alleleSizeAssociation.AllResults.txt.gz` | All allele-size cutoffs tested per locus |
| `<cohort>.<chrom>.alleleSizeAssociation.BestResults.txt.gz` | Best (lowest p-value) cutoff per locus |

Run on DNAnexus.

### Dependencies

`tidyverse`, `data.table`, `scales`, `argparser`, `R.utils`, `readxl`, `parallel`, `brglm`

---

## 4 · Causal variant analysis (`causal_variant_analysis/`)

See [`causal_variant_analysis/README.md`](causal_variant_analysis/README.md) for the full pipeline description.

Fine-mapping (CAVIAR) and conditional analysis (REGENIE) of TR/trait pairs to test whether TR associations are independent of nearby SNP signals. Designed to run inside a Jupyter notebook on the *All of Us* Researcher Workbench.

---

## 5 · GTEx TR-eQTL analysis (`GTEx.eQTLAssociation.neuroDegen.R`)

### What it does

Per-tissue expression residualization followed by long-allele association testing at host genes and genes within a window of each TR locus, with cross-tissue meta-analysis. Approach inspired by the GTEx TR-eQTL analysis in Manigbas et al. 2024 Nat Commun. A similar script was used for the MESA and AMP-PD gene-expression association analyses (§6, §7), with the following differences in how expression was processed: in GTEx, pre-processed `normalized_expression.bed` files were residualized on sex, median insert size, sequencer, genotyping principal components, and PEER factors; in MESA, raw TPM values were INT-transformed within the analysis and residualized on sex, median insert size, age, analyte isolation batch, and predicted ancestry; in AMP-PD, raw normalized counts were likewise INT-transformed within the analysis and residualized on sex, age, cohort (study of origin), median insert size, and five genotyping principal components.

- **Per tissue**: loads GTEx's normalized expression matrix and covariates file (every covariate column used, no manual selection), joins AGE/SEX from the subject phenotype file and an optional insert-size covariate, then residualizes each needed gene's expression against all covariates.
- **Per (locus, gene, tissue)**: expands each TR locus ± a window and overlaps it against a gene annotation BED to build host-gene and window-gene candidate pairs, then tests `lm(residual ~ long_allele_length)` for each pair.
- **Across tissues**: applies per-tissue FDR correction, plus an optional random-effects meta-analysis (`metafor::rma`, brain/nonbrain/all) for host-gene pairs.

### Usage

```bash
Rscript GTEx.eQTLAssociation.neuroDegen.R \
    --Dir                      /path/to/project/dir \
    --RegionsFile              /path/to/regions.txt \
    --FolderName                neuroDegen \
    --GenotypeFile               /path/to/GTEx.filtered.alleles.anno.exps.classifierPredictions.txt.gz \
    --GenotypeAllFile            /path/to/GTEx.QC.filtered.alleles.bed.gz \
    --ExpressionMatrixDir       /path/to/GTEx_Analysis_v8_eQTL_expression_matrices \
    --PEERDir                   /path/to/GTEx_Analysis_v8_eQTL_covariates \
    --GeneAnnotationFile         /path/to/all_genes.bed \
    --GeneSymbolToEnsemblFile    /path/to/symbol_to_ensembl.txt \
    --SubjectPhenoFile           /path/to/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt
```

### Required arguments

| Argument | Description |
|---|---|
| `--Dir` | Main project directory (output written here) |
| `--RegionsFile` | File of `Locus` values to keep, or `no` to keep all |
| `--FolderName` | Output subfolder name, or `no` for none |
| `--GenotypeFile` | TR genotype + annotation file (per-sample allele calls, closestGene/geneName) |
| `--GenotypeAllFile` | TR genotype file with the full per-locus allele distribution (`Max` column) |
| `--ExpressionMatrixDir` | Directory with `{tissue}.v8.normalized_expression.bed[.gz]` files |
| `--PEERDir` | Directory with `{tissue}.v8.covariates.txt` files |

### Key optional arguments

| Argument | Description |
|---|---|
| `--WindowKb` | ± window (kb) for the window-gene scan (default `100`) |
| `--MinSamples` | Minimum samples per (locus, gene, tissue) test (default `10`) |
| `--MinDistinctAlleles` | Minimum distinct long-allele values per tissue (default `5`) |
| `--FDRThreshold` | FDR (q) threshold for calling a significant eQTL (default `0.1`) |
| `--GeneAnnotationFile` | Headerless BED (`chrom, start, end, gene_id`); `no` disables the window scan |
| `--GeneSymbolToEnsemblFile` | `symbol, ensembl_id` crosswalk for host-gene symbols |
| `--SubjectPhenoFile` | GTEx subject phenotypes file for AGE/SEX; `no` to skip |
| `--InsertSizeFile` | `SUBJID` + insert-size covariate file; `no` to skip |
| `--Threads` | Tissues to process in parallel (forked workers; Linux/Mac only) |

Run with `--help` for the full argument list.

### Notes / assumptions

- The expression bed's `gene_id` column is expected to be a versioned Ensembl ID (`.version` suffix stripped before matching). If your host genes are gene symbols rather than Ensembl IDs, supply `--GeneSymbolToEnsemblFile`.
- `--GeneAnnotationFile` is read as a headerless BED with columns assigned by position (`chrom, start, end, gene_id`), not by header name.
- Sample IDs are assumed to be `SUBJID`-style (`GTEX-1117F`); dots are converted back to dashes in case `fread` mangled them.

### Output

Written to `<Dir>/geneExpression/<FolderName>/eqtlAnalysis/`:

| File | Contents |
|---|---|
| `GTEx.eQTL.per_sample_resid_max.txt` | Per-sample residual expression + long-allele length for every tested pair |
| `GTEx.eQTL.per_tissue_lm.txt` | Per-(locus, gene, tissue) association stats, with per-tissue FDR |
| `GTEx.eQTL.per_locus_anyTissueQTL.txt` | Per-(locus, gene) summary across tissues |
| `GTEx.eQTL.hostGene_crossTissueMeta.txt` | Cross-tissue meta-analysis for host-gene pairs (brain/nonbrain/all) |

### HPC example (LSF)

```bash
bsub -P acc_PROJECTID -q premium -W 10:00 -n 8 \
     -R "span[hosts=1]" -R "rusage[mem=60000]" \
     Rscript GTEx.eQTLAssociation.neuroDegen.R \
         --Dir /path/to/project/dir --RegionsFile /path/to/regions.txt \
         --FolderName neuroDegen \
         --GenotypeFile /path/to/GTEx.filtered.alleles.anno.exps.classifierPredictions.txt.gz \
         --GenotypeAllFile /path/to/GTEx.QC.filtered.alleles.bed.gz \
         --ExpressionMatrixDir /path/to/GTEx_Analysis_v8_eQTL_expression_matrices \
         --PEERDir /path/to/GTEx_Analysis_v8_eQTL_covariates \
         --Threads 8
```

### Dependencies

`argparser`, `data.table`, `dplyr`, `tidyverse`, `readxl`, `metafor`, `parallel`

---

## 6 · MESA TR-eQTL analysis (`MESA.eQTLAssociation.neuroDegen.R`)

### What it does

Per-tissue (Monocytes, PBMCs, T-cells) expression residualization followed by long-allele association testing at host genes and genes within a window of each TR locus. Part of the same cis-eQTL analysis as the GTEx and AMP-PD scripts (§5, §7); see §5 for how expression processing differs across the three cohorts.

- **Per tissue**: rank-based inverse-normal transform (INT) each needed gene's TPM, then residualizes on sex, analyte isolation batch, and predicted ancestry. Genes are only tested in a tissue if median TPM there clears `--MinTPM`.
- **Per (locus, gene, tissue)**: expands each TR locus ± a window and overlaps it against a gene annotation BED to build host-gene and window-gene candidate pairs, then tests `lm(residual ~ long_allele_length)` for each pair. A leverage check refits after dropping the single largest-allele sample (`beta_drop_top` / `p_drop_top` / `delta_beta_drop_top`).
- **Across tissues**: applies per-tissue FDR correction.

### Usage

```bash
Rscript MESA.eQTLAssociation.neuroDegen.R \
    --Dir            /path/to/project/dir \
    --RegionsFile    /path/to/regions.txt \
    --FolderName     neuroDegen \
    --GenotypeDir    /path/to/EH_output \
    --MetadataDir    /path/to/metadata \
    --ExpressionDir  /path/to/RNASeq \
    --AnnotationFile /path/to/annotation.tsv
```

### Required arguments

| Argument | Description |
|---|---|
| `--Dir` | Main project directory (output written here) |
| `--RegionsFile` | File of `Locus` values (chr/start/end encoded as `chr_start_end`) to test |
| `--FolderName` | Output subfolder name, or `no` for none |
| `--GenotypeDir` | Directory with per-chromosome ExpansionHunter output |
| `--MetadataDir` | Directory with WGS + per-tissue RNA-seq metadata |
| `--ExpressionDir` | Directory with per-tissue gene TPM matrices |
| `--AnnotationFile` | ANNOVAR annotation file (host gene per locus) |

### Key optional arguments

| Argument | Description |
|---|---|
| `--WindowKb` | ± window (kb) for the window-gene scan (default `100`) |
| `--MinTPM` | Minimum median TPM in a tissue to test a gene there (default `0.1`) |
| `--MinDistinctAlleles` | Minimum distinct long-allele values per tissue (default `5`) |
| `--MinSamples` | Minimum samples per (locus, gene, tissue) test (default `10`) |
| `--FDRThreshold` | FDR (q) threshold for calling a significant eQTL (default `0.1`) |
| `--GeneAnnotationFile` | Headerless BED (`chrom, start, end, gene_id`); `no` disables the window scan |

Run with `--help` for the full argument list.

### Expected input layout

```
--GenotypeDir/{chr}_MESA_EHv5_GenomewidePolymorphic_ForPhewasGWAS.tsv.gz
--MetadataDir/MESA_WGS_meta_genotyped.txt
--MetadataDir/MESA_RNASeq_meta_{Mono,PBMC,TCell}.txt
--ExpressionDir/MESA.RNASeq.{Monocytes,PBMCs,TCell}.gene_tpm.ExpressedGenes.FINAL.gct.gz
```

### Notes / assumptions

- `gene_id`/`gene_name` are resolved directly from each tissue's own expression matrix (`Name`/`Description` columns); genes ambiguous across tissues are excluded from the host/window-gene lookup.
- `--GeneAnnotationFile` is read as a headerless BED with columns assigned by position, not by header name.
- Host genes come from the ANNOVAR closest-gene annotation at each locus and are always tested regardless of distance; window genes are added only within `--WindowKb`.

### Output

Written to `<Dir>/geneExpression/<FolderName>/eqtlAnalysis/`:

| File | Contents |
|---|---|
| `MESA.eQTL.per_sample_z_longAllele.txt` | Per-sample residual expression + long-allele length for every tested pair |
| `MESA.eQTL.per_tissue_lm.txt` | Per-(locus, gene, tissue) association stats, with per-tissue FDR |

### HPC example (LSF)

```bash
bsub -P acc_PROJECTID -q premium -W 10:00 -R "span[hosts=1]" -R "rusage[mem=60000]" \
    Rscript MESA.eQTLAssociation.neuroDegen.R \
        --Dir /path/to/project/dir --RegionsFile /path/to/regions.txt --FolderName neuroDegen \
        --GenotypeDir /path/to/EH_output --MetadataDir /path/to/metadata \
        --ExpressionDir /path/to/RNASeq --AnnotationFile /path/to/annotation.tsv
```

### Dependencies

`argparser`, `data.table`, `dplyr`, `tidyr`, `tibble`

---

## 7 · AMP-PD TR-eQTL analysis (`AMPPD.eQTLAssociation.neuroDegen.R`)

### What it does

Expression residualization (blood RNA-seq, single tissue) followed by long-allele association testing at host genes and genes within a window of each TR locus. Part of the same cis-eQTL analysis as the GTEx and MESA scripts (§5, §6).

- Rank-based inverse-normal transform (INT) each needed gene's raw normalized expression, then residualizes on cohort, sex, age, median insert size, and genotyping principal components.
- **Per (locus, gene)**: expands each TR locus ± a window and overlaps it against a gene annotation BED to build host-gene and window-gene candidate pairs, then tests `lm(residual ~ long_allele_length)` for each pair. A leverage check refits after dropping the single largest-allele sample (`beta_drop_top` / `p_drop_top` / `delta_beta_drop_top`).
- Bonferroni + FDR correction across all tests.

### Usage

```bash
Rscript AMPPD.eQTLAssociation.neuroDegen.R \
    --Dir            /path/to/project/dir \
    --RegionsFile    /path/to/regions.txt \
    --FolderName     neuroDegen \
    --GenotypeFile   /path/to/genotypes.tsv \
    --ExpressionFile /path/to/expression.tsv \
    --AnnotationFile /path/to/annotation.tsv
```

### Required arguments

| Argument | Description |
|---|---|
| `--Dir` | Main project directory (output written here) |
| `--RegionsFile` | File of `Locus` values (chr/start/end encoded as `chr_start_end`) to test |
| `--FolderName` | Output subfolder name, or `no` for none |
| `--GenotypeFile` | TR genotype + classifier + covariate file |
| `--ExpressionFile` | Normalized gene expression matrix |
| `--AnnotationFile` | ANNOVAR annotation file (host gene per locus) |

### Key optional arguments

| Argument | Description |
|---|---|
| `--WindowKb` | ± window (kb) for the window-gene scan (default `100`) |
| `--MinDistinctAlleles` | Minimum distinct long-allele values (default `5`) |
| `--MinSamples` | Minimum samples per (locus, gene) test (default `10`) |
| `--FDRThreshold` | FDR (q) threshold for calling a significant eQTL (default `0.1`) |
| `--GeneAnnotationFile` | Headerless BED (`chrom, start, end, gene_id`); `no` disables the window scan |

Run with `--help` for the full argument list.

### Notes / assumptions

- `gene_id`/`gene_name` are resolved directly from the expression matrix's own columns; ambiguous genes are excluded from the host/window-gene lookup.
- Host genes come from the ANNOVAR closest-gene annotation at each locus and are always tested regardless of distance; window genes are added only within `--WindowKb`.
- RNA-seq sample IDs are matched to genotype SampleIds by truncating to the first two hyphen-delimited fields (e.g. `BF-1001-SVM0_5T1` -> `BF-1001`).

### Output

Written to `<Dir>/geneExpression/<FolderName>/eqtlAnalysis/`:

| File | Contents |
|---|---|
| `AMPPD.eQTL.per_sample_z_longAllele.txt` | Per-sample residual expression + long-allele length for every tested pair |
| `AMPPD.eQTL.per_locus_lm.txt` | Per-(locus, gene) association stats, with Bonferroni + FDR |

### HPC example (LSF)

```bash
bsub -P acc_PROJECTID -q premium -W 4:00 -R "span[hosts=1]" -R "rusage[mem=30000]" \
    Rscript AMPPD.eQTLAssociation.neuroDegen.R \
        --Dir /path/to/project/dir --RegionsFile /path/to/regions.txt --FolderName neuroDegen \
        --GenotypeFile /path/to/genotypes.tsv --ExpressionFile /path/to/expression.tsv \
        --AnnotationFile /path/to/annotation.tsv
```

### Dependencies

`argparser`, `data.table`, `dplyr`, `tidyr`, `tibble`, `readr`

---

## Data availability

Input genotype and phenotype data are derived from UK Biobank (UKB), the *All of Us* Research Program, GTEx, MESA, and AMP-PD, and are not publicly distributable. Access to the UKB, *All of Us*, MESA, and AMP-PD data requires approval from the respective data access committees; GTEx v8 expression matrices and covariates used by the eQTL script are available from the GTEx Portal (open-access tier), while individual-level genotypes require dbGaP access:

- **UK Biobank**: [ukbiobank.ac.uk](https://www.ukbiobank.ac.uk/)
- **All of Us**: [researchallofus.org](https://www.researchallofus.org/)
- **GTEx**: [gtexportal.org](https://www.gtexportal.org/)
- **MESA**: genotype/RNA-seq data via [dbGaP](https://www.ncbi.nlm.nih.gov/gap/)
- **AMP-PD**: application required via the AMP-PD Knowledge Platform

---

## Authors

Gabrielle N. Altman, Bharati Jadhav, Paras Garg, Mariya Shadrina, Celine A. Manigbas, William Lee, Shrishtee Kandoi, Alejandro Martin-Trujillo, Andrew J. Sharp — Icahn School of Medicine at Mount Sinai.


---

## License

See [LICENSE](LICENSE).
