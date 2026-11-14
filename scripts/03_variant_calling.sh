#!/usr/bin/env bash
# =============================================================================
# Script: 03_variant_calling.sh
# Description: Multi-sample variant calling (bcftools mpileup + call), BCF indexing, and VCF conversion.
# =============================================================================

set -e

# 1. CONFIGURATION & VARIABLES
SAMPLES=("C40" "ZnO40" "ZnONPs40" "ZnO20+20" "ZnONPs20+20")
THREADS=4

# Relative pathways (synchronized with previous steps)
REF_DIR="data/reference2"
REF_FASTA="${REF_DIR}/ATCC25922.fasta"  # Ujisti se, že název FASTA souboru odpovídá tomu, co máš ve složce reference
ALN_DIR="results/mapped2"
VAR_DIR="results/variants2"

# 2. GLOBAL INPUT CHECK FOR REFERENCE
if [ ! -f "${REF_FASTA}" ]; then
    echo "[ERROR] Reference genome was not found in the directory: ${REF_FASTA}"
    exit 1
fi

# Create output directory
mkdir -p "${VAR_DIR}"

echo "Starting Variant Calling Pipeline for all samples..."

# 3. MAIN LOOP OVER ALL SAMPLES
for SAMPLE in "${SAMPLES[@]}"; do
    echo ""
    echo "=========================================================="
    echo "Processing sample: ${SAMPLE}"
    echo "=========================================================="

    # Input check per sample
    if [ ! -f "${ALN_DIR}/${SAMPLES}_final.bam" ]; then
        echo "[ERROR] Input final BAM was not found: ${ALN_DIR}/${SAMPLE}_final.bam"
        echo "Skipping this sample. Make sure script 02 ran successfully."
        continue
    fi

    # STEP A: VARIANT CALLING TO BCF FORMAT (-Ob)
    echo "[1/3] Computing Genotype Likelihoods & Calling Variants to BCF (${SAMPLE})"
    bcftools mpileup --threads ${THREADS} -f "${REF_FASTA}" "${ALN_DIR}/${SAMPLE}_final.bam" | \
    bcftools call --threads ${THREADS} -m -v -Ob -o "${VAR_DIR}/${SAMPLE}_raw_variants.bcf"

    # STEP B: INDEXING BCF FILE
    echo "[2/3] Indexing BCF File (${SAMPLE})"
    bcftools index "${VAR_DIR}/${SAMPLE}_raw_variants.bcf"

    # STEP C: CONVERTING BCF TO VCF FOR INSPECTION (-Ov)
    echo "[3/3] Converting BCF to readable VCF (${SAMPLE})"
    bcftools view -Ov -o "${VAR_DIR}/${SAMPLE}_raw_variants.vcf" "${VAR_DIR}/${SAMPLE}_raw_variants.bcf"

    # Count total variants
    TOTAL_VARIANTS=$(bcftools view -H "${VAR_DIR}/${SAMPLE}_raw_variants.vcf" | wc -l)
    echo "📊 Total raw variants detected for ${SAMPLE}: ${TOTAL_VARIANTS}"

done

echo ""
echo "=========================================================="
echo "[SUCCESS] Variant calling and VCF conversion complete for all samples!"
echo "=========================================================="

# 4. JOINT VARIANT CALLING (All samples combined)
echo ""
echo "=========================================================="
echo "Running Joint Variant Calling for ALL samples together..."
echo "=========================================================="

JOINT_OUT="${VAR_DIR}/all_samples_joint_raw_variants"

# Spuštění mpileup nad všemi BAM soubory současně
bcftools mpileup --threads ${THREADS} -f "${REF_FASTA}" \
    "${ALN_DIR}/C40_final.bam" \
    "${ALN_DIR}/ZnO40_final.bam" \
    "${ALN_DIR}/ZnONPs40_final.bam" \
    "${ALN_DIR}/ZnO20+20_final.bam" \
    "${ALN_DIR}/ZnONPs20+20_final.bam" | \
bcftools call --threads ${THREADS} -m -v -Ob -o "${JOINT_OUT}.bcf"

# Indexace a převod joint BCF na VCF
bcftools index "${JOINT_OUT}.bcf"
bcftools view -Ov -o "${JOINT_OUT}.vcf" "${JOINT_OUT}.bcf"

TOTAL_JOINT_VARIANTS=$(bcftools view -H "${JOINT_OUT}.vcf" | wc -l)
echo "📊 Total joint variants detected across all samples: ${TOTAL_JOINT_VARIANTS}"

echo "[SUCCESS] Both individual and joint variant calling pipelines completed!"


