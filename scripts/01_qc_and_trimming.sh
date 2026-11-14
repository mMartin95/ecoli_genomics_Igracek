#!/usr/bin/env bash
# =============================================================================
# Script: 01_qc_and_trimming.sh (Multi-sample version)
# Description: Quality control and adapter/quality trimming (fastp) for multiple E. coli experimental lineages
# =============================================================================

set -e

## 1. CONFIGURATION & VARIABLES 
# Definice pole se všemi vzorky najednou
SAMPLES=("C40" "ZnO40" "ZnONPs40" "ZnO20+20" "ZnONPs20+20")
THREADS=4

RAW_DIR="data/raw"
TRIM_DIR="results/trimmed"
QC_DIR="qc"

# Vytvoření globálních adresářů
mkdir -p "${QC_DIR}/before_trimming" "${QC_DIR}/after_trimming" "${TRIM_DIR}"

## 2. MAIN PROCESSING LOOP
for SAMPLE in "${SAMPLES[@]}"; do
    echo "=========================================================="
    echo "Processing sample: ${SAMPLE}"
    echo "=========================================================="

    ## Input check per sample
    if [ ! -f "${RAW_DIR}/${SAMPLE}_R1.fastq.gz" ] || [ ! -f "${RAW_DIR}/${SAMPLE}_R2.fastq.gz" ]; then
        echo "[ERROR] Input files for sample ${SAMPLE} were not detected in ${RAW_DIR}/!"
        echo "Make sure ${SAMPLE}_R1.fastq.gz and ${SAMPLE}_R2.fastq.gz exist."
        exit 1
    fi

    ## 3. RAW DATA QC (FastQC) 
    echo "--- Running FastQC on Raw Reads (${SAMPLE}) ---"
    fastqc "${RAW_DIR}/${SAMPLE}_R1.fastq.gz" "${RAW_DIR}/${SAMPLE}_R2.fastq.gz" \
        --outdir "${QC_DIR}/before_trimming" \
        -t ${THREADS}

    ## 4. TRIMMING & FILTERING (fastp)
    echo "--- Trimming adapters and low-quality bases (fastp) ---"
    fastp \
        -i "${RAW_DIR}/${SAMPLE}_R1.fastq.gz" \
        -I "${RAW_DIR}/${SAMPLE}_R2.fastq.gz" \
        -o "${TRIM_DIR}/${SAMPLE}_R1_trimmed.fastq.gz" \
        -O "${TRIM_DIR}/${SAMPLE}_R2_trimmed.fastq.gz" \
        -3 -W 4 -M 20 -l 50 \
        -j "${QC_DIR}/after_trimming/fastp_${SAMPLE}.json" \
        -h "${QC_DIR}/after_trimming/fastp_${SAMPLE}.html" \
        -w ${THREADS}

    ## 5. TRIMMED DATA QC (FastQC)
    echo "--- Running FastQC on Trimmed Reads (${SAMPLE}) ---"
    fastqc "${TRIM_DIR}/${SAMPLE}_R1_trimmed.fastq.gz" "${TRIM_DIR}/${SAMPLE}_R2_trimmed.fastq.gz" \
        --outdir "${QC_DIR}/after_trimming" \
        -t ${THREADS}

    echo "[SUCCESS] Sample ${SAMPLE} finished successfully."
    echo ""
done

## 6. GLOBAL MULTIQC (Aggregating all samples at once)
echo "=== Generating Global MultiQC Reports ==="
multiqc "${QC_DIR}/before_trimming" \
    -o "${QC_DIR}" \
    -n raw_data_multiqc_report.html \
    --force

multiqc "${QC_DIR}/after_trimming" \
    -o "${QC_DIR}" \
    -n trimmed_data_multiqc_report.html \
    --force

echo "[SUCCESS] All samples processed and MultiQC reports generated!"