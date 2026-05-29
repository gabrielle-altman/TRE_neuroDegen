#!/usr/bin/env bash
# common.sh — shared functions for the TR causal variant pipeline
# source after config.sh

# emit TR THRESH PHENO triples from TR_TRAIT_PAIRS
iterate_pairs() {
    awk -F'\t' 'NR>1 && NF>=3 {print $1"\t"$2"\t"$3}' "$TR_TRAIT_PAIRS" | sort -u
}

unique_trs() {
    awk -F'\t' 'NR>1 {print $1}' "$TR_TRAIT_PAIRS" | sort -u
}

chroms_from_bed() {
    local bed="$1"
    awk '{sub(/^chr/,"",$1); print $1}' "$bed" | sort -u
}

# parse REGENIE output by column name rather than index — handles column reordering
regenie_top_lines() {
    local infile="$1"
    awk '
        NR==1 {
            for (i=1; i<=NF; i++) col[$i]=i;
            for (k in needed) if (!(needed[k] in col)) {
                print "ERROR: column "needed[k]" missing in "FILENAME > "/dev/stderr";
                exit 1
            }
            next
        }
        BEGIN { split("LOG10P ID BETA SE", needed, " ") }
        $col["LOG10P"] != "NA" && $col["LOG10P"] != "" {
            print $col["LOG10P"], $col["ID"], $col["BETA"], $col["SE"]
        }
    ' "$infile"
}

# plink2 writes "#IID" only when no family info; REGENIE needs "#FID IID"
fix_psam_header() {
    local psam="$1"
    local header
    header=$(head -n1 "$psam")
    if [[ "$header" == *"#FID"* && "$header" == *"IID"* ]]; then
        echo "psam $psam already has #FID/IID — skipping" >&2
        return 0
    fi
    awk -v OFS='\t' '
        NR==1 { print "#FID","IID", substr($0, index($0,$2)); next }
        { print $1, $1, substr($0, index($0,$2)) }
    ' "$psam" > "${psam}.tmp" && mv "${psam}.tmp" "$psam"
}

# merge TR pgen + locus SNP pgen, restricted to shared samples
# optional 5th arg: SNP extract list (fine-mapping uses top-100; conditional uses all)
build_tr_snp_merged_pgen() {
    local tr="$1" pheno="$2" thresh="$3" out_prefix="$4"
    local snp_extract_list="${5:-}"

    local tr_pgen="${TR_PLINK_DIR}/TR_${thresh}"
    local snp_pgen="${TRAIT_DIR}/${TYPE}/${tr}_${pheno}_snps"
    local sample_keep="${out_prefix}_SampleKeep"

    # Build the intersection of TR and SNP samples (uniq -d keeps lines
    # appearing twice, i.e. in both psams).
    {
        printf "#FID\tIID\n"
        cat "${tr_pgen}.psam" "${snp_pgen}.psam" \
            | grep -v '^#FID' \
            | awk -v OFS='\t' '{print $1, $2}' \
            | sort \
            | uniq -d
    } > "$sample_keep"

    # TR ID must exist in the TR pvar.
    if ! grep -qw "$tr" "${tr_pgen}.pvar"; then
        echo "ERROR: TR '$tr' not found in ${tr_pgen}.pvar" >&2
        return 1
    fi

    plink2 --pfile "$tr_pgen" \
           --snps "$tr" \
           --keep "$sample_keep" \
           --make-bed \
           --out "${out_prefix}_TR"

    local extract_args=()
    if [[ -n "$snp_extract_list" ]]; then
        extract_args+=(--extract "$snp_extract_list")
    fi
    plink2 --pfile "$snp_pgen" \
           "${extract_args[@]}" \
           --keep "$sample_keep" \
           --make-bed \
           --out "${out_prefix}_snps"

    plink --bfile "${out_prefix}_snps" \
          --bmerge "${out_prefix}_TR" \
          --make-bed \
          --out "${out_prefix}_snps_TR"

    rm -f "${out_prefix}_TR".{bed,bim,fam} "${out_prefix}_snps".{bed,bim,fam}
}

run_regenie_step2_conditional() {
    local tr="$1" pheno="$2" thresh="$3"
    local condition_file="$4" out_suffix="$5" min_mac="$6"

    local pgen="${COND_DIR}/${tr}_${pheno}_${TYPE}_${thresh}_snps_TR"
    local sample_list; sample_list=$(sample_list_for_regenie "$pheno")
    local pheno_file; pheno_file=$(phenotype_file "$pheno")
    local exclude_snp="${TRAIT_DIR}/${TYPE}/${tr}_${pheno}_drop_snp.txt"
    local pred_list="${COND_STEP1_DIR}/${pheno}_${TYPE}_SNPs.step1_pred.list"

    "$REGENIE_BIN" --step 2 \
        --pgen "$pgen" \
        --keep "$sample_list" \
        --exclude "$exclude_snp" \
        --phenoFile "$pheno_file" \
        --phenoCol "$pheno" \
        --strict \
        --condition-list "$condition_file" \
        --covarFile "$COVAR_FILE" \
        --catCovarList "$CAT_COVARS" \
        --covarColList "$QUANT_COVARS" \
        --out "${COND_STEP1_DIR}/${tr}_${pheno}_${TYPE}_${thresh}_${out_suffix}" \
        --bt --bsize 1000 \
        --firth --approx \
        --loocv \
        --lowmem \
        --lowmem-prefix "${COND_STEP1_DIR}/tmp_rg_${pheno}_${tr}" \
        --threads "$THREADS" \
        --pred "$pred_list" \
        --pThresh 0.05 \
        --minMAC "$min_mac"
}
