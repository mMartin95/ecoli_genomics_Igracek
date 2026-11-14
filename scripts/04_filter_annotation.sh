#!/usr/bin/env bash
# =============================================================================
# Script: 04_annotation.sh
# Description: Quality Hard-filtering and SnpEff Annotation using custom local DB
# =============================================================================

set -e 

# 1. CONFIGURATION & VARIABLES 
VAR_DIR="results/variants2"
SNPEFF_DB="ATCC25922_whole"
DATA_DIR="./data"

RAW_VCF="${VAR_DIR}/all_samples_joint_raw_variants.vcf"
FILTERED_VCF="${VAR_DIR}/all_samples_joint_filtered.vcf"
ANNOTATED_VCF="${VAR_DIR}/all_samples_joint_annotated.vcf"
CONFIG_FILE="./snpEff.config"

# 2. INPUT CHECK 
if [ ! -f "${RAW_VCF}" ]; then
    echo "[ERROR] Raw joint VCF file was not found: ${RAW_VCF}"
    echo "Make sure script 03 finished successfully."
    exit 1
fi

mkdir -p "${VAR_DIR}"

# 3. STEP 1: QUALITY HARD-FILTERING
echo ""
echo "[1/2] Quality Hard-Filtering (QUAL > 30, FMT/DP >= 10, MQ > 30)"
bcftools filter -i 'QUAL > 30 && DP >= 10 && MQ > 30' "${RAW_VCF}" -O v -o "${FILTERED_VCF}"

PASSED_COUNT=$(bcftools view -H "${FILTERED_VCF}" | wc -l)
echo "📊 High-quality variants remaining after filtering: ${PASSED_COUNT}"

# 4. STEP 2: FUNCTIONAL ANNOTATION (SnpEff)
echo ""
echo "[2/2] Functional Variant Annotation (SnpEff)"
snpEff ann -v -c "${CONFIG_FILE}" -dataDir "${DATA_DIR}" "${SNPEFF_DB}" "${FILTERED_VCF}" > "${ANNOTATED_VCF}"

echo ""
echo "[SUCCESS] Pipeline finished! Final annotated output: ${ANNOTATED_VCF}"