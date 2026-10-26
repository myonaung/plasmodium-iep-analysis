#!/bin/bash
#PBS -l walltime=48:00:00
#PBS -l mem=100gb
#PBS -l jobfs=5gb
#PBS -l ncpus=20
#PBS -q normal
#PBS -N haplotype_caller
#PBS -o haplotype_caller.out
#PBS -e haplotype_caller.err

# =============================================================================
# 03_haplotype_caller.sh
#
# Step 3 of 4: runs GATK HaplotypeCaller in GVCF mode per sample, per
# chromosome, chunked into ~10kb intervals for parallel-friendly progress
# tracking and safe resuming. Produces one GVCF per sample per chromosome,
# ready for 04_vqsr.sh.
#
# The `intervals` array below defines the callable chromosome ranges to
# call over — these must match the reference genome set in
# pipeline_config.sh (REF_DIR/REFERENCE). Replace the example
# chromosome/start/end triples with the ones appropriate for your reference
# and organism.
#
# Configure paths in pipeline_config.sh or with environment variables.
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/pipeline_config.sh"

mkdir -p "$INTERVAL_DIR" "$GVCF_DIR" "$LOG_DIR"

# Required tools must be available on PATH: Java and GATK.

export OMP_NUM_THREADS=20  # GATK native HMM threads

# ========================== Define Callable Genome Intervals ==================
# Example intervals below are for the Plasmodium falciparum 3D7 v3 reference
# core genome. Replace with "<chrom> <start> <end>" triples matching your
# own reference.
intervals=(
"Pf3D7_01_v3 92901 457931"
"Pf3D7_01_v3 460312 575900"
"Pf3D7_02_v3 105801 447300"
"Pf3D7_02_v3 450451 862500"
"Pf3D7_03_v3 70631 597816"
"Pf3D7_03_v3 600276 1003060"
"Pf3D7_04_v3 91421 545800"
"Pf3D7_04_v3 614901 642003"
"Pf3D7_04_v3 644530 935030"
"Pf3D7_04_v3 983081 1143990"
"Pf3D7_05_v3 37901 455740"
"Pf3D7_05_v3 457253 1321390"
"Pf3D7_06_v3 72351 478652"
"Pf3D7_06_v3 480972 723117"
"Pf3D7_06_v3 742801 1294830"
"Pf3D7_07_v3 77101 508360"
"Pf3D7_07_v3 605651 809245"
"Pf3D7_07_v3 811717 1381600"
"Pf3D7_08_v3 73561 299079"
"Pf3D7_08_v3 301404 427430"
"Pf3D7_08_v3 467341 1365730"
"Pf3D7_09_v3 79101 1242137"
"Pf3D7_09_v3 1244484 1473560"
"Pf3D7_10_v3 68971 1571815"
"Pf3D7_11_v3 110001 831968"
"Pf3D7_11_v3 834246 2003320"
"Pf3D7_12_v3 60301 766654"
"Pf3D7_12_v3 780451 1282773"
"Pf3D7_12_v3 1285068 1688600"
"Pf3D7_12_v3 1745531 2163700"
"Pf3D7_13_v3 74414 1168127"
"Pf3D7_13_v3 1170426 2791900"
"Pf3D7_14_v3 35775 1071523"
"Pf3D7_14_v3 1075090 3255710"
)

# ======================= Generate Interval Files ==============================
for entry in "${intervals[@]}"; do
    read -r chrom start end <<< "$entry"
    interval_file="$INTERVAL_DIR/${chrom}.intervals"
    [[ -f "$interval_file" ]] && continue
    > "$interval_file"
    chunk_start=$start
    while [[ "$chunk_start" -lt "$end" ]]; do
        chunk_end=$((chunk_start + 9999))
        [[ "$chunk_end" -gt "$end" ]] && chunk_end=$end
        echo "${chrom}:${chunk_start}-${chunk_end}" >> "$interval_file"
        chunk_start=$((chunk_end + 1))
    done
    echo "Created: $interval_file"
done

# ===================== Run HaplotypeCaller per Sample ========================
cd "$RECAL_BAM_DIR" || { echo "Cannot access $RECAL_BAM_DIR"; exit 1; }

shopt -s nullglob
bam_files=(*.recal_reads.bam)
if (( ${#bam_files[@]} == 0 )); then
    echo "No BAM files found in $RECAL_BAM_DIR"
    exit 1
fi

for recal_bam in *.recal_reads.bam; do
    prefix=$(basename "$recal_bam" .recal_reads.bam)
    sample_log="${LOG_DIR}/${prefix}.log"
    echo "Processing $prefix" | tee -a "$sample_log"

    for entry in "${intervals[@]}"; do
        read -r chrom _ _ <<< "$entry"
        interval_file="$INTERVAL_DIR/${chrom}.intervals"
        chr_dir="$GVCF_DIR/${chrom}"
        mkdir -p "$chr_dir"
        output_gvcf="${chr_dir}/${prefix}_${chrom}.g.vcf.gz"
        lockfile="${output_gvcf}.lock"

        if [[ -f "$output_gvcf" ]]; then
            echo "  Skipping $chrom: already exists." | tee -a "$sample_log"
            continue
        fi

        if [[ -f "$lockfile" ]]; then
            echo "  Skipping $chrom: lock exists." | tee -a "$sample_log"
            continue
        fi

        touch "$lockfile"
        echo "  Calling $prefix @ $chrom" | tee -a "$sample_log"

        gatk --java-options "-Xmx16G" HaplotypeCaller \
            -R "${REF_DIR}/${REFERENCE}" \
            -I "$recal_bam" \
            -L "$interval_file" \
            -ERC GVCF \
            -ploidy 2 \
            --native-pair-hmm-threads 4 \
            -O "$output_gvcf" \
            -G StandardAnnotation \
            -G AS_StandardAnnotation \
            --assembly-region-padding 100 \
            --max-num-haplotypes-in-population 128 \
            --kmer-size 10 \
            --kmer-size 25 \
            --min-dangling-branch-length 4 \
            --heterozygosity 0.001 \
            --indel-heterozygosity 0.000125 \
            --min-assembly-region-size 50 \
            --contamination-fraction-to-filter 0 \
            -mbq 5 \
            --base-quality-score-threshold 12 \
            2>> "$sample_log"

        rm -f "$lockfile"
        echo "  Done $chrom" | tee -a "$sample_log"
    done

    echo "Finished: $prefix" | tee -a "$sample_log"
done

echo "All GVCF generation completed."
