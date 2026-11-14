#!/usr/bin/env bash
# =============================================================================
# Script: 02_mapping_and_coverage.sh
# Description: Multi-sample reads alignment (BWA-MEM), BAM sorting, duplicate removal, indexing, and depth/coverage metrics evaluation.
# =============================================================================

set -e

# 1. CONFIGURATION AND VARIABLES
SAMPLES=("C40" "ZnO40" "ZnONPs40" "ZnO20+20" "ZnONPs20+20")
THREADS=4

# Paths synchronized with Script 01
REF_DIR="data/reference2"
REF_FASTA="${REF_DIR}/ATCC25922.fasta"
TRIM_DIR="results/trimmed"
ALN_DIR="results/mapped2"
QC_DIR="qc"

MIN_COVERAGE=20

# 2. GLOBAL INPUT CHECK
if [ ! -f "${REF_FASTA}" ]; then
    echo "[ERROR] Reference genome was not found in the directory : ${REF_FASTA}"
    echo "Download or prepare reference genome of E. coli ATCC 25922."
    exit 1
fi

# Preparation of directories
mkdir -p "${ALN_DIR}" "${QC_DIR}"

# 3. REFERENCE INDEXING (Executed ONCE globally before the loop)
echo "=========================================================="
echo "Checking Reference Genome Index"
echo "=========================================================="
if [ ! -f "${REF_FASTA}.bwt" ]; then
    echo "Building BWA index for reference genome..."
    bwa index "${REF_FASTA}"
else
    echo "BWA index already exists. Skipping."
fi

if [ ! -f "${REF_FASTA}.fai" ]; then
    echo "Building Samtools FASTA index (.fai)..."
    samtools faidx "${REF_FASTA}"
else
    echo "Samtools FASTA index already exists. Skipping."
fi

# 4. MAIN LOOP OVER ALL SAMPLES
for SAMPLE in "${SAMPLES[@]}"; do
    echo ""
    echo "=========================================================="
    echo "Processing sample: ${SAMPLE}"
    echo "=========================================================="

    # Check if trimmed files exist in results/trimmed (fixed path alignment)
    if [ ! -f "${TRIM_DIR}/${SAMPLE}_R1_trimmed.fastq.gz" ] || [ ! -f "${TRIM_DIR}/${SAMPLE}_R2_trimmed.fastq.gz" ]; then
        echo "[ERROR] Trimmed FASTQ files for ${SAMPLE} were not found in ${TRIM_DIR}/!"
        echo "Skipping this sample. Make sure script 01 ran successfully."
        continue
    fi

    # MAPPING, CONVERSION & COORDINATE SORTING 
    echo "[1/5] Aligning reads with BWA-MEM and sorting to BAM (${SAMPLE})"
    bwa mem -t ${THREADS} -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
        "${REF_FASTA}" \
        "${TRIM_DIR}/${SAMPLE}_R1_trimmed.fastq.gz" \
        "${TRIM_DIR}/${SAMPLE}_R2_trimmed.fastq.gz" | \
        samtools view -u - | \
        samtools sort -@ ${THREADS} -o "${ALN_DIR}/${SAMPLE}_sorted.bam" -

    # MARK & REMOVE PCR DUPLICATES
    echo "[2/5] Marking and Removing PCR Duplicates (${SAMPLE})"
    samtools collate -o "${ALN_DIR}/${SAMPLE}_collate.bam" "${ALN_DIR}/${SAMPLE}_sorted.bam"
    samtools fixmate -m "${ALN_DIR}/${SAMPLE}_collate.bam" "${ALN_DIR}/${SAMPLE}_fixmate.bam"
    samtools sort -@ ${THREADS} -o "${ALN_DIR}/${SAMPLE}_fixmate_sorted.bam" "${ALN_DIR}/${SAMPLE}_fixmate.bam"
    samtools markdup -r -@ ${THREADS} "${ALN_DIR}/${SAMPLE}_fixmate_sorted.bam" "${ALN_DIR}/${SAMPLE}_final.bam"

    # Removing generated intermediate BAM files to save disk space
    rm "${ALN_DIR}/${SAMPLE}_sorted.bam" "${ALN_DIR}/${SAMPLE}_collate.bam" "${ALN_DIR}/${SAMPLE}_fixmate.bam" "${ALN_DIR}/${SAMPLE}_fixmate_sorted.bam"

    # INDEXING FINAL BAM
    echo "[3/5] Indexing Final BAM File (${SAMPLE})"
    samtools index "${ALN_DIR}/${SAMPLE}_final.bam"

    # COVERAGE EVALUATION & QUALITY CHECKPOINT
    echo "[4/5] Calculating Coverage and Depth Metrics (${SAMPLE})"
    COVERAGE_REPORT="${QC_DIR}/${SAMPLE}_coverage_stats.txt"
    samtools coverage "${ALN_DIR}/${SAMPLE}_final.bam" > "${COVERAGE_REPORT}"

    # ALIGNMENT STATISTICS
    echo "[5/5] Generating Alignment Statistics (samtools flagstat) (${SAMPLE})"
    FLAGSTAT_REPORT="${QC_DIR}/${SAMPLE}_flagstat_stats.txt"
    samtools flagstat "${ALN_DIR}/${SAMPLE}_final.bam" > "${FLAGSTAT_REPORT}"

    # Extraction of average mean depth
    MEAN_DEPTH=$(awk 'NR==2 {print $7}' "${COVERAGE_REPORT}")
    echo "Mean Coverage Depth for ${SAMPLE}: ${MEAN_DEPTH}x"

    # Comparison of calculated coverage with the minimal coverage (MIN_COVERAGE)
    IS_SUFFICIENT=$(awk -v depth="${MEAN_DEPTH}" -v min="${MIN_COVERAGE}" 'BEGIN {if (depth >= min) print "YES"; else print "NO"}')

    if [ "${IS_SUFFICIENT}" == "YES" ]; then
        echo "[SUCCESS] Sample ${SAMPLE} passed coverage checkpoint (${MEAN_DEPTH}x >= ${MIN_COVERAGE}x)."
    else
        echo "[WARNING] Low coverage detected for ${SAMPLE} (${MEAN_DEPTH}x < ${MIN_COVERAGE}x)!"
        echo "Variant calling in Step 03 may yield false positives or missing calls for this lineage."
    fi
done

echo ""
echo "=========================================================="
echo "Pipeline execution completed for all samples!"
echo "=========================================================="