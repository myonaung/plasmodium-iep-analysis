#!/bin/bash
#PBS -l walltime=08:00:00
#PBS -l mem=96gb
#PBS -l ncpus=48
#PBS -q normal
#PBS -N trim_align
#PBS -o trim_align.out
#PBS -e trim_align.err

# =============================================================================
# 01_trim_align.sh
#
# Step 1 of 4: trims paired-end reads with Trimmomatic, aligns to the
# alignment reference with BWA-MEM, sorts, marks duplicates, and subsets to
# the core-genome regions. Produces one core-region BAM per sample, ready
# for 02_bqsr.sh.
#
# Configure paths in pipeline_config.sh or with environment variables.
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/pipeline_config.sh"

# Required tools must be available on PATH: bwa, samtools, java and gatk.

mkdir -p "$TRIMMED_DIR" "$BAM_DIR" "$LOCK_DIR"

# ===================== Reference Index =====================
if [ ! -f "${ALIGN_REF_DIR}/${ALIGN_REFERENCE}.bwt" ]; then
    echo "Indexing alignment reference..."
    bwa index "${ALIGN_REF_DIR}/${ALIGN_REFERENCE}"
fi

# =================== Pipeline Start ========================
find "$RAW_FASTQ_DIR" -maxdepth 1 -type f -name "*_R1.fastq.gz" | sort -r | while read -r R1_FILE; do
  BASE=$(basename "$R1_FILE" "_R1.fastq.gz")
  DIR=$(dirname "$R1_FILE")
  R2_FILE="${DIR}/${BASE}_R2.fastq.gz"

  FWD_PAIRED="${TRIMMED_DIR}/${BASE}_forward_paired.fq.gz"
  FWD_UNPAIRED="${TRIMMED_DIR}/${BASE}_forward_unpaired.fq.gz"
  REV_PAIRED="${TRIMMED_DIR}/${BASE}_reverse_paired.fq.gz"
  REV_UNPAIRED="${TRIMMED_DIR}/${BASE}_reverse_unpaired.fq.gz"

  FINAL_BAM="${BAM_DIR}/${BASE}.sorted.dup.core.bam"
  LOCKFILE="${LOCK_DIR}/${BASE}.lock"

  # Final BAM exists -> skip
  if [ -f "$FINAL_BAM" ]; then
      echo "Skipping $BASE: final BAM already exists."
      continue
  fi

  # Intermediate files exist -> skip
  if ls "${WORK_DIR}/${BASE}"*.bam "${WORK_DIR}/${BASE}"*.sam "${TRIMMED_DIR}/${BASE}"*.fq.gz 2>/dev/null | grep -q .; then
      echo "Skipping $BASE: intermediate files suggest in-progress or failed run."
      continue
  fi

  # Lockfile already exists -> skip
  if [ -f "$LOCKFILE" ]; then
      echo "Skipping $BASE: lockfile exists (already running elsewhere)."
      continue
  fi

  # Create lock
  touch "$LOCKFILE"

  echo "==> Trimming $BASE..."
  java -jar "$TRIMMOMATIC_JAR" PE -threads 48 -phred33 \
    "$R1_FILE" "$R2_FILE" \
    "$FWD_PAIRED" "$FWD_UNPAIRED" \
    "$REV_PAIRED" "$REV_UNPAIRED" \
    HEADCROP:20 ILLUMINACLIP:"$ADAPTER_FILE":2:30:10 \
    LEADING:3 TRAILING:3 MINLEN:36

  echo "==> Aligning $BASE..."
  bwa mem -t 48 -M -R "@RG\tID:${BASE}\tLB:${BASE}\tPL:ILLUMINA\tSM:${BASE}\tPU:${BASE}" \
    "${ALIGN_REF_DIR}/${ALIGN_REFERENCE}" "$FWD_PAIRED" "$REV_PAIRED" > "${WORK_DIR}/${BASE}_aligned.sam"

  samtools view -@ 4 -bS "${WORK_DIR}/${BASE}_aligned.sam" > "${WORK_DIR}/${BASE}_aligned.bam"
  rm "${WORK_DIR}/${BASE}_aligned.sam"

  gatk --java-options "-Xmx32g" CleanSam \
    -R "${ALIGN_REF_DIR}/${ALIGN_REFERENCE}" \
    -I "${WORK_DIR}/${BASE}_aligned.bam" \
    -O "${WORK_DIR}/${BASE}.clean.bam"
  rm "${WORK_DIR}/${BASE}_aligned.bam"

  gatk --java-options "-Xmx32g" SortSam \
    -R "${ALIGN_REF_DIR}/${ALIGN_REFERENCE}" \
    -I "${WORK_DIR}/${BASE}.clean.bam" \
    -O "${WORK_DIR}/${BASE}.sorted.bam" \
    -SO coordinate --CREATE_INDEX true
  rm "${WORK_DIR}/${BASE}.clean.bam"
  rm -f "${WORK_DIR}/${BASE}.sorted.bam.bai"

  gatk --java-options "-Xmx32g" MarkDuplicates \
    -R "${ALIGN_REF_DIR}/${ALIGN_REFERENCE}" \
    -I "${WORK_DIR}/${BASE}.sorted.bam" \
    -O "${WORK_DIR}/${BASE}.sorted.dup.bam" \
    -M "${WORK_DIR}/${BASE}_dup_metrics.tsv"
  rm "${WORK_DIR}/${BASE}.sorted.bam"

  samtools view -@ 4 -b -h "${WORK_DIR}/${BASE}.sorted.dup.bam" \
    -T "${ALIGN_REF_DIR}/${ALIGN_REFERENCE}" \
    -L "$CORE_REGIONS" \
    > "$FINAL_BAM"
  rm "${WORK_DIR}/${BASE}.sorted.dup.bam"

  samtools index "$FINAL_BAM"

  echo "==> Cleaning up $BASE..."
  rm -f "$FWD_PAIRED" "$FWD_UNPAIRED" "$REV_PAIRED" "$REV_UNPAIRED"
  rm -f "${WORK_DIR}/${BASE}_dup_metrics.tsv"
  rm -f "$LOCKFILE"

  echo "==> Completed $BASE"
done

echo "All samples processed and intermediate files cleaned."
