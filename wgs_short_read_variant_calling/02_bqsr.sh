#!/bin/bash
#PBS -l walltime=48:00:00
#PBS -l mem=100gb
#PBS -l ncpus=1
#PBS -l jobfs=1gb
#PBS -q normal
#PBS -N BQSR
#PBS -o BQSR.out
#PBS -e BQSR.err

# =============================================================================
# 02_bqsr.sh
#
# Step 2 of 4: Base Quality Score Recalibration. Builds a recalibration
# table from known variant sites, applies it to each core-region BAM
# produced by 01_trim_align.sh, and generates before/after covariate plots
# for QC. Produces one recalibrated BAM per sample, ready for
# 03_haplotype_caller.sh.
#
# Configure paths in pipeline_config.sh or with environment variables.
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/pipeline_config.sh"

mkdir -p "$RECAL_BAM_DIR"

# Required tools must be available on PATH: R, Java and GATK.

# ========================== Base Recalibration =============================
cd "$BAM_DIR" || exit

# Build the --known-sites argument list; KNOWN_SITES2 is optional.
known_sites_args=(--known-sites "$KNOWN_SITES1")
if [[ -n "$KNOWN_SITES2" ]]; then
    known_sites_args+=(--known-sites "$KNOWN_SITES2")
fi

find . -name "*.sorted.dup.core.bam" | while read -r sample_path; do
    sample=$(basename "$sample_path")
    prefix=$(basename "$sample" .sorted.dup.core.bam)
    lockfile="${RECAL_BAM_DIR}/${prefix}.bqsr.lock"
    output_bam="${RECAL_BAM_DIR}/${prefix}.recal_reads.bam"

    if [[ -f "$output_bam" ]]; then
        echo "Skipping $prefix (already processed)"
        continue
    fi

    if [[ -f "$lockfile" ]]; then
        echo "Skipping $prefix (lock file exists)"
        continue
    fi

    touch "$lockfile"

    echo "Processing sample: $prefix"

    # Step 1: Generate BQSR table
    gatk BaseRecalibrator \
        -R "${REF_DIR}/${REFERENCE}" \
        -I "$sample" \
        "${known_sites_args[@]}" \
        --use-original-qualities \
        -O "${RECAL_BAM_DIR}/${prefix}.recal_data.table"

    # Step 2: Apply BQSR
    gatk ApplyBQSR \
        -R "${REF_DIR}/${REFERENCE}" \
        -I "$sample" \
        --bqsr-recal-file "${RECAL_BAM_DIR}/${prefix}.recal_data.table" \
        --use-original-qualities \
        -O "$output_bam"

    # Step 3: Post-BQSR table (QC)
    gatk BaseRecalibrator \
        -R "${REF_DIR}/${REFERENCE}" \
        -I "$output_bam" \
        "${known_sites_args[@]}" \
        --use-original-qualities \
        -O "${RECAL_BAM_DIR}/${prefix}.post_recal_data.table"

    # Step 4: Analyze covariates
    gatk AnalyzeCovariates \
        -before "${RECAL_BAM_DIR}/${prefix}.recal_data.table" \
        -after "${RECAL_BAM_DIR}/${prefix}.post_recal_data.table" \
        -csv "${RECAL_BAM_DIR}/${prefix}.recalibration.csv"

    echo "Completed: $prefix"
    rm -f "$lockfile"
done

echo "All samples processed with base recalibration!"
