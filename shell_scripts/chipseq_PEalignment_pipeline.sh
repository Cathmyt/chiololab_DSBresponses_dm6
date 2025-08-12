#! /bin/bash

set -euo pipefail
eval "$(conda shell.bash hook)"

# ===========================================================
# Script: Paired-end FASTQ → Trim → Align → Dedup → Coverage → Compare
#
# Environment setup (Conda):
#   bwaenv        : bwa
#   samtoolsenv   : samtools
#   picardenv     : picard
#   deeptoolsenv  : deepTools (bamCoverage, bamCompare)
#   cutadaptenv   : trim_galore (requires cutadapt + fastqc)
#
# Description:
#   - Accepts .fastq / .fq / .fastq.gz / .fq.gz
#   - Trims adapters (Trim Galore, paired-end, gzip output)
#   - Aligns with BWA MEM
#   - Sorts, indexes, removes duplicates (Picard)
#   - Generates normalized coverage tracks (bigWig) with/without duplicates
#   - Compares treatment vs untreated (bamCompare, SES normalization)
#
# Usage:
#   ./script.sh --trt <treatment_keyword> --unt <untreated_keyword> \
#               --output <output_name> [--threads N] [--binSize B]
#
# ===========================================================



# ============================
# Default parameters
# ============================
THREADS=20
BINSIZE=30
TRT_KEYWORD=""
UNT_KEYWORD=""
FILENAME=""

# ============================
# Parse inputs
# ============================
while [[ $# -gt 0 ]]; do
    case $1 in
        --trt) TRT_KEYWORD="$2"; shift 2 ;;
        --unt) UNT_KEYWORD="$2"; shift 2 ;;
        --output) FILENAME="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --binSize) BINSIZE="$2"; shift 2 ;;
        *) echo "Unknown option $1"; exit 1 ;;
    esac
done

if [[ -z $TRT_KEYWORD || -z $UNT_KEYWORD || -z $FILENAME ]]; then
    echo "Usage: $0 --trt <treatment_keyword> --unt <untreated_keyword> --output <output_name> [--threads N] [--binSize B]"
    exit 1
fi

timestamp "Using THREADS=$THREADS, BINSIZE=$BINSIZE"

# ============================
# Locate input files
# ============================
TRT_FILE1=$(ls | grep "$TRT_KEYWORD" | grep "_1\.f\(ast\)\?q\.gz\?$" | head -n1)
UNT_FILE1=$(ls | grep "$UNT_KEYWORD" | grep "_1\.f\(ast\)\?q\.gz\?$" | head -n1)

if [[ -z $TRT_FILE1 || -z $UNT_FILE1 ]]; then
    echo "Error: Matching paired-end files not found for keywords."
    exit 1
fi

TRT_RAW="${TRT_FILE1%_*1.f*}"
UNT_RAW="${UNT_FILE1%_*1.f*}"

timestamp "Treatment sample: $TRT_RAW"
timestamp "Untreated sample: $UNT_RAW"

# ============================
# Process function
# ============================
process_sample() {
    local RAW=$1
    local BASE=$2

    timestamp "Processing: $RAW as $BASE"

    mkdir -p trimming_reports
    validate_env cutadaptenv
    trim_galore --illumina --fastqc --paired "${RAW}_1.fq"* "${RAW}_2.fq"* -o trimming_reports

    validate_env bwaenv
    bwa mem -t $THREADS "$REF_GENOME" \
        trimming_reports/${RAW}_1_val_1.fq.gz trimming_reports/${RAW}_2_val_2.fq.gz \
        > ${BASE}_bwa_aligned_pe.sam

    validate_env samtoolsenv
    samtools view -b -@ $THREADS ${BASE}_bwa_aligned_pe.sam \
        | samtools sort -@ $THREADS -o ${BASE}_bwa_aligned_se.sam.sorted.bam
    rm ${BASE}_bwa_aligned_pe.sam
    samtools index ${BASE}_bwa_aligned_se.sam.sorted.bam

    validate_env picardenv
    picard MarkDuplicates \
        INPUT=${BASE}_bwa_aligned_se.sam.sorted.bam \
        METRICS_FILE=${BASE}_markdup.metrics \
        OUTPUT=${BASE}_bwa_aligned_se.sam.sorted.bam.nodup.bam \
        REMOVE_DUPLICATES=true ASSUME_SORTED=true
    rm -f ${BASE}_markdup.metrics

    validate_env samtoolsenv
    samtools index ${BASE}_bwa_aligned_se.sam.sorted.bam.nodup.bam

    validate_env deeptoolsenv
    bamCoverage -b ${BASE}_bwa_aligned_se.sam.sorted.bam \
        -o ${BASE}_wdup.bw -p $THREADS --extendReads $EXTEND &
    bamCoverage -b ${BASE}_bwa_aligned_se.sam.sorted.bam.nodup.bam \
        -o ${BASE}_nodup.bw -p $THREADS --extendReads $EXTEND &
    wait
}

# ============================
# Run processing
# ============================
process_sample "$TRT_RAW" "${TRT_KEYWORD}_${FILENAME}_TRT"
process_sample "$UNT_RAW" "${UNT_KEYWORD}_${FILENAME}_UNT"

# ============================
# Generate bamCompare
# ============================
TRT_BASE="${TRT_KEYWORD}_${FILENAME}_TRT"
UNT_BASE="${UNT_KEYWORD}_${FILENAME}_UNT"

timestamp "Generating bamCompare outputs"
validate_env deeptoolsenv

bamCompare -b1 ${TRT_BASE}_bwa_aligned_se.sam.sorted.bam \
           -b2 ${UNT_BASE}_bwa_aligned_se.sam.sorted.bam \
           --scaleFactorsMethod SES --binSize $BINSIZE --smoothLength 300 \
           -p $THREADS --extendReads $EXTEND \
           -o ${FILENAME}_${TRT_KEYWORD}${UNT_KEYWORD}_wdup.bw

bamCompare -b1 ${TRT_BASE}_bwa_aligned_se.sam.sorted.bam.nodup.bam \
           -b2 ${UNT_BASE}_bwa_aligned_se.sam.sorted.bam.nodup.bam \
           --scaleFactorsMethod SES --binSize $BINSIZE --smoothLength 300 \
           -p $THREADS --extendReads $EXTEND --ignoreDuplicates \
           -o ${FILENAME}_${TRT_KEYWORD}${UNT_KEYWORD}_nodup.bw

timestamp "Alignment completed."