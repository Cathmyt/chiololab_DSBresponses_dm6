#! /bin/bash

set -euo pipefail
eval "$(conda shell.bash hook)"

# ============================
# Default parameters
# ============================
THREADS=20
REF_GENOME="/Users/chiololab/wdir/JacobMac/ReferenceGenomes/dm6.fa"

# ============================
# Validate Conda environment
# ============================
validate_env() {
    local env_name=$1
    conda activate "$env_name" || { echo "Error: Conda environment '$env_name' not found!"; exit 1; }
}

# ============================
# Main
# ============================
cd "$(pwd)"

for FILE in *.fastq; do
    BASE="${FILE%.fastq}" 
    BASE_TR="${BASE}_trimmed"
    timestamp "Processing: $BASE"
    
    # --- Trimming ---
    trim_galore --illumina --fastqc "${BASE}.fastq"

    # --- Alignment ---
    validate_env bwaenv
    bwa mem -t $THREADS "$REF_GENOME" "${BASE_TR}.fq" > "${BASE_TR}_bwa_aligned_se.sam"

    # --- Convert SAM to BAM + Sort ---
    validate_env samtoolsenv
    samtools view -b -@ $THREADS "${BASE_TR}_bwa_aligned_se.sam" \
      | samtools sort -@ $THREADS -o "${BASE_TR}_bwa_aligned_se.sam.sorted.bam"
    rm -f "${BASE_TR}_bwa_aligned_se.sam"
    samtools index "${BASE_TR}_bwa_aligned_se.sam.sorted.bam"

    # --- Remove duplicates ---
    validate_env picardenv
    picard MarkDuplicates \
        INPUT="${BASE_TR}_bwa_aligned_se.sam.sorted.bam" \
        METRICS_FILE="${BASE_TR}_bwa_aligned_se.sam.sorted.bam.markup.metrics" \
        OUTPUT="${BASE_TR}_bwa_aligned_se.sam.sorted.bam.nodup.bam" \
        REMOVE_DUPLICATES=true ASSUME_SORTED=true
    rm -f "${BASE_TR}_bwa_aligned_se.sam.sorted.bam.markup.metrics"

    validate_env samtoolsenv
    samtools index "${BASE_TR}_bwa_aligned_se.sam.sorted.bam.nodup.bam"

    # --- Coverage Tracks (Parallelized) ---
    validate_env deeptoolsenv
    bamCoverage -b "${BASE_TR}_bwa_aligned_se.sam.sorted.bam.nodup.bam" \
        -o "${BASE_TR}_bwa_aligned_se_nodup.bw" -p $THREADS --extendReads 200 &
    bamCoverage -b "${BASE_TR}_bwa_aligned_se.sam.sorted.bam" \
        -o "${BASE_TR}_bwa_aligned_se_wdup.bw" -p $THREADS --extendReads 200 &
    wait

    timestamp "Processing DONE for $BASE"
done

timestamp "Alignment completed."
