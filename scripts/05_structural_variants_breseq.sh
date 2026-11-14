#!/usr/bin/env bash

# ==============================================================================
# Pipeline: Breseq Variant Calling, Structural Variant Filtration & gdtools
# Reference: E. coli ATCC 25922 (GenBank: NZ_CP009072 / .gb format)
# ==============================================================================

set -euo pipefail

# 1. Configuration of paths and variables
REF="data/reference2/ATCC25922.gb"
FASTQ_DIR="results/trimmed"
OUT_DIR="results/breseq"
THREADS=4

# Check if reference file exists
if [ ! -f "$REF" ]; then
    echo "Error: Reference file $REF not found!" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

# Sample list (adjust according to actual filenames without _R1/_R2.fastq.gz)
SAMPLES=("C40" "ZnO40" "ZnONPs40" "ZnO20_20" "ZnONPs20_20")

echo "=== STEP 1: Running breseq for all samples ==="
for SAMPLE in "${SAMPLES[@]}"; do
    R1="${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
    R2="${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"
    
    if [ -f "$R1" ] && [ -f "$R2" ]; then
        echo "--> Processing sample: $SAMPLE"
        breseq -r "$REF" -o "${OUT_DIR}/${SAMPLE}_output" -j "$THREADS" "$R1" "$R2"
    else
        echo "Warning: FASTQ files for $SAMPLE not found, skipping." >&2
    fi
done

echo "=== STEP 2: Extraction and filtering of structural variants (MC/JC) ==="
mkdir -p "${OUT_DIR}/structural_variants"

# Example: Filtering ZnONPs40 (or ZnO40) against the C40 control
TREATMENT="ZnONPs40"
CONTROL="C40"

TREAT_GD="${OUT_DIR}/${TREATMENT}_output/output/output.gd"
CTRL_GD="${OUT_DIR}/${CONTROL}_output/output/output.gd"
SV_OUT="${OUT_DIR}/structural_variants"

if [ -f "$TREAT_GD" ] && [ -f "$CTRL_GD" ]; then
    echo "--> Extracting Missing Coverage (MC) and Junctions (JC) for $TREATMENT"
    
    # Treatment extraction (Note: breseq outputs annotated.gd inside output/ directory)
    awk 'BEGIN{OFS="\t"} /^MC/ {print $4, $5, $6}' "$TREAT_GD" > "${SV_OUT}/${TREATMENT}_annotated_deletions.bed"
    awk 'BEGIN{OFS="\t"} /^JC/ {print $4, $5, $5+1}' "$TREAT_GD" > "${SV_OUT}/${TREATMENT}_annotated_junctions.bed"
    
    # Control extraction
    awk 'BEGIN{OFS="\t"} /^MC/ {print $4, $5, $6}' "$CTRL_GD" > "${SV_OUT}/${CONTROL}_annotated_deletions.bed"
    
    echo "--> Running bedtools subtract (removing deletions present in control)"
    bedtools subtract -a "${SV_OUT}/${TREATMENT}_annotated_deletions.bed" \
                      -b "${SV_OUT}/${CONTROL}_annotated_deletions.bed" \
                      > "${SV_OUT}/${TREATMENT}_cleaned.bed"
                      
    echo "--> Checking read depth for cleaned deletions via bedtools coverage"
    BAM_FILE="${OUT_DIR}/${TREATMENT}_output/data/reference.bam"
    if [ -f "$BAM_FILE" ]; then
        bedtools coverage -a "${SV_OUT}/${TREATMENT}_cleaned.bed" \
                          -b "$BAM_FILE" \
                          > "${SV_OUT}/${TREATMENT}_MC_coverage.tsv"
    else
        echo "Warning: BAM file for coverage not found at $BAM_FILE" >&2
    fi
else
    echo "Warning: Output .gd files for $TREATMENT or $CONTROL are missing." >&2
fi

echo "=== STEP 3: gdtools SUBTRACT & COMPARE ==="
mkdir -p "${OUT_DIR}/gdtools_comparisons"

# Remove background mutations using gdtools SUBTRACT
if [ -f "${OUT_DIR}/ZnO40_output/output/output.gd" ] && [ -f "${OUT_DIR}/C40_output/output/output.gd" ]; then
    gdtools SUBTRACT -o "${OUT_DIR}/gdtools_comparisons/ZnO40_clean.gd" \
                     "${OUT_DIR}/ZnO40_output/output/output.gd" \
                     "${OUT_DIR}/C40_output/output/output.gd"
    echo "--> Created cleaned .gd file for ZnO40 without control background."
fi

# Compare all treated lineages in a single HTML summary
gdtools COMPARE -r "$REF" \
                -o "${OUT_DIR}/gdtools_comparisons/all_treatments_comparison.html" \
                "${OUT_DIR}/ZnO40_output/output/output.gd" \
                "${OUT_DIR}/ZnONPs40_output/output/output.gd" \
                "${OUT_DIR}/ZnO20_20_output/output/output.gd" \
                "${OUT_DIR}/ZnONPs20_20_output/output/output.gd"

echo "=== PIPELINE COMPLETED SUCCESSFULLY ==="