#! /bin/bash

set -euo pipefail
eval "$(conda shell.bash hook)"

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
        --trt)
            TRT_KEYWORD="$2"
            shift 2
            ;;
        --unt)
            UNT_KEYWORD="$2"
            shift 2
            ;;
        --output)
            FILENAME="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --binSize)
            BINSIZE="$2"
            shift 2
            ;;
        *)
            echo "Unknown option $1"
            exit 1
            ;;
    esac
done

if [[ -z $TRT_KEYWORD || -z $UNT_KEYWORD || -z $FILENAME ]]; then
    echo "Usage: $0 --trt <treatment_keyword> --unt <untreated_keyword> --output <formatted_FILENAME> [--threads N] [--binSize B]"
    exit 1
fi

# ============================
# Validate Conda environments
# ============================
validate_env() {
    local env_name=$1
    conda activate "$env_name" || { echo "Error: Conda environment '$env_name' not found!"; exit 1; }
}
echo "Using THREADS=$THREADS, BINSIZE=$BINSIZE"



cd "$(pwd)"
REF_GENOME="/Users/chiololab/wdir/JacobMac/ReferenceGenomes/dm6.fa"

# ============================
# Find input files
# ============================
TRT_FILE=$(ls | grep "$TRT_KEYWORD" | grep "_1.fq.gz")
UNT_FILE=$(ls | grep "$UNT_KEYWORD" | grep "_1.fq.gz")

if [[ -z $TRT_FILE || -z $UNT_FILE ]]; then
    echo "Error: Matching files not found for the keywords"
    exit 1
fi

# Extract base names
TRT_RAW="${TRT_FILE%_*1.fq.gz}"
UNT_RAW="${UNT_FILE%_*1.fq.gz}"

echo "Processing Treatment: $TRT_RAW"
echo "Processing Untreated: $UNT_RAW"


# ============================
# Process Treatment and Untreated
# ============================
for RAW in $TRT_RAW $UNT_RAW; do
    if [[ "$RAW" == "$TRT_RAW" ]]; then
        BASE="${TRT_KEYWORD}_${FILENAME}_TRT"
    else
        BASE="${UNT_KEYWORD}_${FILENAME}_UNT"
    fi
    echo "Processing: $RAW as $BASE"

    mkdir -p trimming_reports
    cd trimming_reports
    trim_galore --illumina --fastqc --paired ../${RAW}_1.fq.gz ../${RAW}_2.fq.gz
    cd ..

    validate_env bwaenv
    bwa mem -t $THREADS $REF_GENOME trimming_reports/${RAW}_1_val_1.fq.gz trimming_reports/${RAW}_2_val_2.fq.gz > ${BASE}_bwa_aligned_pe.sam

    validate_env samtoolsenv
    samtools view -b -@ 20 ${BASE}_bwa_aligned_pe.sam | samtools sort -@ 20 -o ${BASE}_bwa_aligned_se.sam.sorted.bam
    rm ${BASE}_bwa_aligned_pe.sam
    samtools index ${BASE}_bwa_aligned_se.sam.sorted.bam

    validate_env picardenv
    picard MarkDuplicates INPUT=${BASE}_bwa_aligned_se.sam.sorted.bam \
    METRICS_FILE=${BASE}_bwa_aligned_se.sam.sorted.bam.markup.metrics \
    OUTPUT=${BASE}_bwa_aligned_se.sam.sorted.bam.nodup.bam \
    REMOVE_DUPLICATES=true ASSUME_SORTED=true
    rm -f ${BASE}_bwa_aligned_se.sam.sorted.bam.markup.metrics

    validate_env samtoolsenv
    samtools index ${BASE}_bwa_aligned_se.sam.sorted.bam.nodup.bam

    validate_env deeptoolsenv
    bamCoverage -b ${BASE}_bwa_aligned_se.sam.sorted.bam \
        -o ${BASE}_bwa_aligned_se_wdup.bw -p $THREADS --extendReads 200 &
    bamCoverage -b ${BASE}_bwa_aligned_se.sam.sorted.bam.nodup.bam \
        -o ${BASE}_bwa_aligned_se_nodup.bw -p $THREADS --extendReads 200 &
    wait

    echo "Processing DONE for $BASE"
done

# ============================
# Generate bamCompare
# ============================
TRT_BASE="${TRT_KEYWORD}_${FILENAME}_TRT"
UNT_BASE="${UNT_KEYWORD}_${FILENAME}_UNT"

echo "Generating bamCompare for $TRT_BASE and $UNT_BASE"
bamCompare -b1 ${TRT_BASE}_bwa_aligned_se.sam.sorted.bam \
           -b2 ${UNT_BASE}_bwa_aligned_se.sam.sorted.bam \
           --scaleFactorsMethod SES --binSize $BINSIZE --smoothLength 300 \
           -p $THREADS --extendReads 200 \
           -o ${FILENAME}_${TRT_KEYWORD}${UNT_KEYWORD}_5hrVsUNT_bwa_aligned_PE_SES_wdup.bw

bamCompare -b1 ${TRT_BASE}_bwa_aligned_se.sam.sorted.bam.nodup.bam \
           -b2 ${UNT_BASE}_bwa_aligned_se.sam.sorted.bam.nodup.bam \
           --scaleFactorsMethod SES --binSize $BINSIZE --smoothLength 300 \
           -p $THREADS --extendReads 200 --ignoreDuplicates \
           -o ${FILENAME}_${TRT_KEYWORD}${UNT_KEYWORD}_5hrVsUNT_bwa_aligned_PE_SES_nodup.bw

echo "Alignment completed."