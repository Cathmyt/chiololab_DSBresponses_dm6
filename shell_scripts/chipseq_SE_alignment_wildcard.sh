#!/usr/bin/env bash
set -euo pipefail
eval "$(conda shell.bash hook)"

# ===========================================================
# Script: Single-end FASTQ → Trim → Align → Dedup → Coverage

# Environment setup (Conda):
#   YMalignmentenv : bwa, samtools
#   picardenv      : picard
#   samtoolsenv    : samtools
#   deeptoolsenv   : deepTools (bamCoverage, etc.)
#   cutadapt       : trim_galore (requires cutadapt + fastqc)
#
# Description:
#   - Accepts .fastq / .fq / .fastq.gz / .fq.gz
#   - Trims adapters (Trim Galore, always outputs gzip)
#   - Aligns with BWA MEM (using compressed FASTQ directly)
#   - Sorts, indexes, removes duplicates (Picard)
#   - Generates normalized coverage tracks (bigWig) with/without duplicates
#

# ============================
# Parameters
# ============================
THREADS=20
REF_GENOME="/Users/chiololab/wdir/JacobMac/ReferenceGenomes/dm6.fa"
EXTEND=200

# ============================
# Utilities
# ============================
timestamp() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

validate_env() {
    local env_name=$1
    conda activate "$env_name" || { echo "Error: Conda environment '$env_name' not found!"; exit 1; }
}

# ============================
# Main
# ============================
shopt -s nullglob
for FILE in *.fastq *.fq *.fastq.gz *.fq.gz; do
    [ -e "$FILE" ] || continue

    BASE="${FILE%.fastq.gz}"
    BASE="${BASE%.fq.gz}"
    BASE="${BASE%.fastq}"
    BASE="${BASE%.fq}"
    BASE_TR="${BASE}_trimmed"

    timestamp "Processing: $BASE"

    # --- Trimming ---
    trim_galore --illumina --fastqc --gzip "$FILE"

    TRIMMED="${BASE_TR}.fq.gz"
    if [[ ! -f "$TRIMMED" ]]; then
        echo "ERROR: Trimmed file $TRIMMED not found."
        exit 1
    fi

    # --- Alignment ---
    validate_env YMalignmentenv
    bwa mem -t $THREADS "$REF_GENOME" "$TRIMMED" \
        | samtools view -b -@ $THREADS - \
        | samtools sort -@ $THREADS -o "${BASE_TR}_bwa_aligned_se.sorted.bam"
    samtools index "${BASE_TR}_bwa_aligned_se.sorted.bam"

    # --- Remove duplicates ---
    validate_env picardenv
    picard MarkDuplicates \
        INPUT="${BASE_TR}_bwa_aligned_se.sorted.bam" \
        METRICS_FILE="${BASE_TR}_bwa_aligned_se.markdup.metrics" \
        OUTPUT="${BASE_TR}_bwa_aligned_se.nodup.bam" \
        REMOVE_DUPLICATES=true ASSUME_SORTED=true
    rm -f "${BASE_TR}_bwa_aligned_se.markdup.metrics"

    validate_env samtoolsenv
    samtools index "${BASE_TR}_bwa_aligned_se.nodup.bam"

    # --- Coverage Tracks ---
    validate_env deeptoolsenv
    bamCoverage -b "${BASE_TR}_bwa_aligned_se.nodup.bam" \
        -o "${BASE_TR}_nodup.bw" -p $THREADS --extendReads $EXTEND &
    bamCoverage -b "${BASE_TR}_bwa_aligned_se.nodup.bam" \
        -o "${BASE_TR}_CPM_nodup.bw" -p $THREADS --extendReads $EXTEND --normalizeUsing CPM &
    bamCoverage -b "${BASE_TR}_bwa_aligned_se.sorted.bam" \
        -o "${BASE_TR}_wdup.bw" -p $THREADS --extendReads $EXTEND &
    wait

    timestamp "Processing DONE for $BASE"
done

timestamp "Alignment completed."
