# TR Causal Variant Analysis (AoU)

End-to-end pipeline for fine-mapping (CAVIAR) and conditional analysis (REGENIE) of tandem-repeat / trait pairs in *All of Us*. Designed to run inside a Jupyter notebook on the AoU Researcher Workbench. Similar pipeline used in UKBB.

## Layout

```
.
├── TR_Causal_Variant_Analysis.ipynb   # main notebook — all pipeline steps
├── config.sh                          # paths, binaries, covariate names
├── common.sh                          # shared bash functions
└── README.md
```

Every `%%bash` cell in the notebook starts with
```bash
source ./config.sh
source ./common.sh
```
so configuration lives in one place and shared logic (TR + SNP merge, REGENIE column parsing, psam header fix) lives in another. Variables don't persist between Jupyter cells, so sourcing on every cell is the workaround.

## Setup

1. Place this directory in your AoU workspace.
2. Edit `config.sh`:
   - Set filenames
3. Confirm binaries:
   - `plink2` and `plink` (1.9) — on `$PATH` in AoU images.
   - `regenie_3.4.1` — local build, path set by `REGENIE_BIN`.
   - `CAVIAR` — local build, path set by `CAVIAR_BIN`.
4. Open `TR_Causal_Variant_Analysis.ipynb` and run the sanity-check cell first.

## Input file schemas

- `TR_Trait_pairs.txt` — TSV with header. Columns 1–3 used: `TR`, `THRESH`, `PHENO`. Additional columns ignored.
- `TRListToUse_250kFlanks.bed` — BED with TR ID in column 4. Chromosome may have `chr` prefix; the pipeline strips it where needed. Genomic regions are TR locus +/- 250k bp. Genomic window can be adjusted when you create this file.
- `*_Binary_PhenoList` — one per TR; one phenotype per line in column 1.
- `samplesToUse.txt` — `#FID IID` (header + 2 columns)
- Per-phenotype sample lists and phenotype files in the parent directory (`../AoU.ALL.sampleList.${PHENO}.tsv` etc.); see helper functions in `config.sh`.

## Running

The cells are numbered 00–14 and can be re-run independently as long as their upstream outputs exist. Logs for the long-running REGENIE/CAVIAR calls land in `./TR_Causal_Variant_Analysis/logs/`.

## Notes / known constraints

- The pipeline is binary-only (`TYPE=Binary`)
- `--minMAC` differs between the two conditional step-2 cells (2 for SNP, 5 for TR) — see comments in cell 14.
- The `%%bash` cells use `set -euo pipefail` (via `config.sh`) so any plink2/REGENIE failure aborts the cell instead of silently producing empty output.
- Similar pipeline was used in UKBB

## Authors

Gabrielle Altman, Bharati Jadhav, Paras Garg, Celine Manigbas — Icahn School of Medicine at Mount Sinai.