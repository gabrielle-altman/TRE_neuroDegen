# TRE_neuroDegen

Scripts for tandem repeat expansion (TRE) association analysis in neurodegenerative disease, as described in:

> *[manuscript citation to be added]*

Two complementary analyses are included:

1. **Allele-size association** (`alleleSizeAssociation.neuroDegen.R`) — per-locus Firth logistic regression testing association of varying expansion size cutoff (≥ 99th-percentile repeat length) with disease status, run on UKBB WGS TR genotype data. Similar script used in AoU.
2. **Burden association via REGENIE** (`runREGENIE.neuroDegen.sh`) — genome-wide association test using REGENIE's two-step Firth approx logistic regression, applied to binary TR genotype data in UKBB. Similar script used in AoU.
3. **Causal variant analysis** (`casual_variant_analysis/`) — fine-mapping (CAVIAR) and conditional analysis (REGENIE) pipeline for identifying whether TR signals are independent of nearby SNPs, run in the *All of Us* Researcher Workbench. Similar pipeline used in UKBB.

---

## Repository layout

```
TRE_neuroDegen/
├── alleleSizeAssociation.neuroDegen.R          # Per-locus allele-size association (R)
├── runREGENIE.neuroDegen.sh                    # REGENIE burden association (bash)
└── casual_variant_analysis/
    ├── TR_Causal_Variant_Analysis.ipynb        # End-to-end causal variant pipeline (Jupyter)
    ├── config.sh                               # All paths, binaries, covariate names
    ├── common.sh                               # Shared bash functions
    └── README.md                               # Causal variant pipeline details
```

---

## 1 · Allele-size association (`alleleSizeAssociation.neuroDegen.R`)

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

### HPC example (LSF)

```bash
bsub -P acc_PROJECTID -L /bin/bash -q express -n 4 \
     -R span[hosts=1] -R rusage[mem=30000] -W 10:00 \
     Rscript alleleSizeAssociation.neuroDegen.R \
         --dir        /path/to/project/dir \
         --cohort     neuroDegen \
         --chrom      chr12 \
         --pheno_file /path/to/phenos.txt \
         --covar_file /path/to/covariates.tsv \
         --metal_file /path/to/metal_results.xlsx
```

### Dependencies

`tidyverse`, `data.table`, `scales`, `argparser`, `R.utils`, `readxl`, `parallel`, `brglm`

---

## 2 · REGENIE burden association (`runREGENIE.neuroDegen.sh`)

### What it does

Runs REGENIE step 1 (whole-genome regression null model) and step 2 (association testing) for TR expansion burden in UKBB WGS data. Binary trait logistic regression with Firth approx correction. Designed to be run once per phenotype × expansion cutoff combination.

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

### HPC example (LSF)

```bash
bsub -P acc_PROJECTID -L /bin/bash -q premium -n 18 \
     -R rusage[mem=10000] -R span[hosts=1] -W 24:00 \
     bash runREGENIE.neuroDegen.sh neuroDegen Cutoff99
```

### Dependencies

`regenie/3.4.1` (loaded via environment modules: `ml regenie/3.4.1`)

---

## 3 · Causal variant analysis (`casual_variant_analysis/`)

See [`casual_variant_analysis/README.md`](casual_variant_analysis/README.md) for the full pipeline description.

Fine-mapping (CAVIAR) and conditional analysis (REGENIE) of TR/trait pairs to test whether TR associations are independent of nearby SNP signals. Designed to run inside a Jupyter notebook on the *All of Us* Researcher Workbench.

---

## Data availability

Input genotype and phenotype data are derived from UK Biobank (UKBB) and the *All of Us* Research Program and are not publicly distributable. Access to these datasets requires approval from the respective data access committees:

- **UK Biobank**: [ukbiobank.ac.uk](https://www.ukbiobank.ac.uk/)
- **All of Us**: [researchallofus.org](https://www.researchallofus.org/)

---

## Authors

Gabrielle Altman, Bharati Jadhav, Paras Garg, Alejandro Martin-Trujillo, Celine Manigbas, Mariya Shadrina, William Lee — Icahn School of Medicine at Mount Sinai.


---

## License

See [LICENSE](LICENSE).
