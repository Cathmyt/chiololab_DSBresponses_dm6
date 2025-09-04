SAMPLES=("Blanks_60_trim3" "Blanks_0_beta_plus_trim3" "Blanks_0_beta_minus_trim3" "WT_0_trim3")

conda activate YMalignmentenv
for s in "${SAMPLES[@]}"; do
  echo "Processing sample: $s"
  # trim to 21 nt
  cutadapt -m 21 -M 21 -o fastq/${s}.L21.fastq.gz fastq/${s}.fastq.gz > align/${s}.L21.cutadapt.txt

  # map to TE consensus (no mismatches)
  bwa aln -n 0 -t 8 /Users/chiololab/wdir/Cathy/reference_genomes/transposon_sequence_set.fa fastq/${s}.L21.fastq.gz > align/${s}.TE.sai
  bwa samse /Users/chiololab/wdir/Cathy/reference_genomes/transposon_sequence_set.fa align/${s}.TE.sai fastq/${s}.L21.fastq.gz \
    | samtools view -bS - > align/${s}.TE.bam

  # keep unmapped-to-TE reads → FASTQ for genome step
  samtools view -b -f 4 align/${s}.TE.bam > align/${s}.nonTE.bam
  samtools fastq align/${s}.nonTE.bam | gzip > fastq/${s}.L21.nonTE.fastq.gz

  # perfect match only
  bwa aln -n 0 -t 8 /Users/chiololab/wdir/Cathy/reference_genomes/dm6.fa fastq/${s}.L21.nonTE.fastq.gz > align/${s}.dm6.sai
  bwa samse /Users/chiololab/wdir/Cathy/reference_genomes/dm6.fa align/${s}.dm6.sai fastq/${s}.L21.nonTE.fastq.gz \
    | samtools view -bS - > align/${s}.dm6.L21.nonTE.bam

  # keep unique mappers (bwa-aln unique ≈ MAPQ >= 37)
  samtools view -b -q 37 align/${s}.dm6.L21.nonTE.bam \
    | samtools sort -o align/${s}.dm6.L21.nonTE.unique.sorted.bam
  samtools index align/${s}.dm6.L21.nonTE.unique.sorted.bam

  # make CPM-normalized bigWig (1 bp bins)
bamCoverage -b align/${s}.dm6.L21.nonTE.unique.sorted.bam \
  -o tracks/${s}.plus.cpm.bw \
  --binSize 1 --normalizeUsing CPM --filterRNAstrand forward
bamCoverage -b align/${s}.dm6.L21.nonTE.unique.sorted.bam \
  -o tracks/${s}.minus.cpm.bw \
  --binSize 1 --normalizeUsing CPM --filterRNAstrand reverse
bamCoverage -b align/${s}.dm6.L21.nonTE.unique.sorted.bam \
  -o tracks/${s}.cpm.bw \
  --binSize 1 --normalizeUsing CPM
done
