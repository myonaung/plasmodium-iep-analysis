#!/usr/bin/env bash
set -euo pipefail

# Standalone long-read amplicon variant-calling workflow used by sm-SNIPER.
# This script is provided as a transparent command-line representation of the
# analysis steps and does not require Snakemake.
#
# Expected input: one or more FASTQ files named SAMPLE.fastq or SAMPLE.fastq.gz.

# -------------------------------------------------------------------------
# Command-line settings
# -------------------------------------------------------------------------
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FASTQ_DIR=""
REFERENCE=""
BED="${PIPELINE_DIR}/resources/pf_target_panel.bed"
OUTDIR=""
THREADS=8

MIN_COVERAGE=5
MAX_COVERAGE=5000
MIN_ALT_FRAC=0.01
ONT_CHEMISTRY="ont_r9_guppy5_sup"

PEPPER_CONTAINER=""
CONTAINER_ENGINE="singularity"  # use "apptainer" if preferred

SCRIPTS="${PIPELINE_DIR}/scripts"
RESOURCES="${PIPELINE_DIR}/resources"

usage() {
  cat <<'EOF'
Usage: run_pf_long_read_variant_calling.sh --fastq-dir DIR --reference FASTA \
  --pepper-container IMAGE --outdir DIR [options]

Options:
  --bed FILE             Target BED (default: resources/pf_target_panel.bed)
  --threads N            Threads (default: 8)
  --container-engine CMD Container engine (default: singularity)
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --fastq-dir) FASTQ_DIR="$2"; shift 2 ;;
    --reference) REFERENCE="$2"; shift 2 ;;
    --pepper-container) PEPPER_CONTAINER="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --bed) BED="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --container-engine) CONTAINER_ENGINE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done


for value in FASTQ_DIR REFERENCE PEPPER_CONTAINER OUTDIR; do
  [[ -n "${!value}" ]] || { echo "ERROR: required argument missing: $value" >&2; usage >&2; exit 2; }
done

# -------------------------------------------------------------------------
# Dependency checks
# -------------------------------------------------------------------------
for cmd in minimap2 samtools longshot intersectBed bcftools bgzip tabix vcftools whatshap Rscript python; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "ERROR: required command not found on PATH: ${cmd}" >&2
    exit 1
  fi
done

if ! command -v "${CONTAINER_ENGINE}" >/dev/null 2>&1; then
  echo "ERROR: container engine not found on PATH: ${CONTAINER_ENGINE}" >&2
  exit 1
fi

if [ ! -f "${REFERENCE}" ]; then
  echo "ERROR: reference FASTA not found: ${REFERENCE}" >&2
  exit 1
fi

if [ ! -f "${BED}" ]; then
  echo "ERROR: BED file not found: ${BED}" >&2
  exit 1
fi

if [ ! -f "${PEPPER_CONTAINER}" ]; then
  echo "ERROR: PEPPER container not found: ${PEPPER_CONTAINER}" >&2
  exit 1
fi

if [ ! -d "${FASTQ_DIR}" ]; then
  echo "ERROR: FASTQ directory not found: ${FASTQ_DIR}" >&2
  exit 1
fi

mkdir -p "${OUTDIR}"

# -------------------------------------------------------------------------
# Per-sample analysis
# -------------------------------------------------------------------------
for FASTQ in "${FASTQ_DIR}"/*.fastq "${FASTQ_DIR}"/*.fastq.gz; do
  [ -e "${FASTQ}" ] || continue

  SAMPLE="$(basename "${FASTQ}")"
  SAMPLE="${SAMPLE%.gz}"
  SAMPLE="${SAMPLE%.fastq}"

  echo "Processing sample: ${SAMPLE}"

  SAMPLE_OUT="${OUTDIR}/${SAMPLE}/out"
  BAM_DIR="${OUTDIR}/bam"
  PEPPER_DIR="${OUTDIR}/PEPPER/${SAMPLE}"

  mkdir -p \
    "${BAM_DIR}" \
    "${SAMPLE_OUT}/QC" \
    "${SAMPLE_OUT}/confident_snv" \
    "${SAMPLE_OUT}/text_file" \
    "${SAMPLE_OUT}/snps_to_evaluate/longshot_pass" \
    "${SAMPLE_OUT}/snps_to_evaluate/pepper_pass" \
    "${SAMPLE_OUT}/snps_to_evaluate/svm_train_data" \
    "${SAMPLE_OUT}/final_vcf_SNIPER" \
    "${PEPPER_DIR}"

  # 1. Map Nanopore reads to the reference genome.
  minimap2 -B 2 -t "${THREADS}" -m 20 -ax map-ont \
    "${REFERENCE}" "${FASTQ}" > "${BAM_DIR}/${SAMPLE}_minimap2.sam"

  samtools view -b "${BAM_DIR}/${SAMPLE}_minimap2.sam" \
    -o "${BAM_DIR}/${SAMPLE}_minimap2.bam"

  samtools sort \
    -o "${BAM_DIR}/${SAMPLE}.fastq_sorted.minimap2.bam" \
    "${BAM_DIR}/${SAMPLE}_minimap2.bam"

  samtools index "${BAM_DIR}/${SAMPLE}.fastq_sorted.minimap2.bam"

  samtools view -q 40 -F 2048 -F 256 -b \
    "${BAM_DIR}/${SAMPLE}.fastq_sorted.minimap2.bam" \
    > "${BAM_DIR}/${SAMPLE}_filtered.bam"

  samtools index "${BAM_DIR}/${SAMPLE}_filtered.bam"

  # 2. Alignment quality-control summaries.
  samtools bedcov "${BED}" "${BAM_DIR}/${SAMPLE}.fastq_sorted.minimap2.bam" \
    > "${SAMPLE_OUT}/QC/${SAMPLE}_amplicon_total_depth.txt"

  samtools depth -b "${BED}" "${BAM_DIR}/${SAMPLE}.fastq_sorted.minimap2.bam" \
    > "${SAMPLE_OUT}/QC/${SAMPLE}_amplicon_depth.txt"

  samtools flagstat "${BAM_DIR}/${SAMPLE}.fastq_sorted.minimap2.bam" \
    > "${SAMPLE_OUT}/QC/${SAMPLE}_mapping_summary.txt"

  awk -f "${SCRIPTS}/average_depth.awk" \
    "${SAMPLE_OUT}/QC/${SAMPLE}_amplicon_depth.txt" \
    | cut -d" " -f3 \
    > "${SAMPLE_OUT}/QC/${SAMPLE}_mean_depth.txt"

  # 3. Longshot variant calling and restriction to target intervals.
  longshot -S \
    -b "${BAM_DIR}/${SAMPLE}.fastq_sorted.minimap2.bam" \
    -f "${REFERENCE}" \
    -o "${SAMPLE_OUT}/${SAMPLE}.vcf" \
    -c "${MIN_COVERAGE}" \
    -F \
    -E "${MIN_ALT_FRAC}" \
    -C "${MAX_COVERAGE}" \
    -I 3 \
    --sample_id "${SAMPLE}"

  intersectBed -a "${SAMPLE_OUT}/${SAMPLE}.vcf" -b "${BED}" -header \
    > "${SAMPLE_OUT}/${SAMPLE}_amplicon.vcf"

  # 4. PEPPER variant calling.
  samtools addreplacerg \
    -r "ID:${SAMPLE}_cohort" \
    -r "LB:${SAMPLE}_LSK109" \
    -r "PL:nanopore" \
    -r "PU:${SAMPLE}_probe" \
    -r "SM:${SAMPLE}" \
    -o "${SAMPLE_OUT}/${SAMPLE}_RG.bam" \
    "${BAM_DIR}/${SAMPLE}.fastq_sorted.minimap2.bam"

  samtools sort -o "${SAMPLE_OUT}/${SAMPLE}_RG_sorted.bam" \
    "${SAMPLE_OUT}/${SAMPLE}_RG.bam"

  samtools index "${SAMPLE_OUT}/${SAMPLE}_RG_sorted.bam"

  "${CONTAINER_ENGINE}" exec --bind "${PWD}" "${PEPPER_CONTAINER}" \
    run_pepper_margin_deepvariant call_variant \
    --only_pepper \
    -b "${SAMPLE_OUT}/${SAMPLE}_RG_sorted.bam" \
    -f "${REFERENCE}" \
    -o "${PEPPER_DIR}/pepper_work" \
    -p "${SAMPLE}" \
    -s "${SAMPLE}" \
    -t "${THREADS}" \
    --"${ONT_CHEMISTRY}"

  PEPPER_VCF="$(find "${PEPPER_DIR}/pepper_work" -type f -name "*FULL.vcf.gz" | sort | head -n 1)"
  if [ -z "${PEPPER_VCF}" ]; then
    echo "ERROR: no PEPPER VCF found for ${SAMPLE}" >&2
    exit 1
  fi

  cp "${PEPPER_VCF}" "${PEPPER_DIR}/PEPPER_VARIANT_FULL.vcf.gz"

  intersectBed -a "${PEPPER_DIR}/PEPPER_VARIANT_FULL.vcf.gz" -b "${BED}" -header \
    > "${SAMPLE_OUT}/${SAMPLE}_PEPPER_VARIANT_FULL.vcf"

  # 5. High-confidence variants supported by both Longshot and PEPPER.
  bgzip -c "${SAMPLE_OUT}/${SAMPLE}_amplicon.vcf" \
    > "${SAMPLE_OUT}/${SAMPLE}_amplicon.vcf.gz"
  tabix "${SAMPLE_OUT}/${SAMPLE}_amplicon.vcf.gz"

  bgzip -c "${SAMPLE_OUT}/${SAMPLE}_PEPPER_VARIANT_FULL.vcf" \
    > "${SAMPLE_OUT}/${SAMPLE}_PEPPER_VARIANT_FULL.vcf.gz"
  tabix "${SAMPLE_OUT}/${SAMPLE}_PEPPER_VARIANT_FULL.vcf.gz"

  bcftools isec -c all \
    "${SAMPLE_OUT}/${SAMPLE}_amplicon.vcf.gz" \
    "${SAMPLE_OUT}/${SAMPLE}_PEPPER_VARIANT_FULL.vcf.gz" \
    -p "${SAMPLE_OUT}/confident_snv/${SAMPLE}" \
    -n=2 -w 1

  bgzip -c "${SAMPLE_OUT}/confident_snv/${SAMPLE}/0000.vcf" \
    > "${SAMPLE_OUT}/confident_snv/${SAMPLE}_confident_snps.vcf.gz"
  tabix "${SAMPLE_OUT}/confident_snv/${SAMPLE}_confident_snps.vcf.gz"

  # 6. Convert single-caller variants to SVM feature tables.
  python "${SCRIPTS}/vcf_to_svm_tables.py" \
    --sample "${SAMPLE}" \
    --longshot-vcf "${SAMPLE_OUT}/${SAMPLE}_amplicon.vcf" \
    --pepper-vcf "${SAMPLE_OUT}/${SAMPLE}_PEPPER_VARIANT_FULL.vcf" \
    --outdir "${SAMPLE_OUT}/text_file"

  Rscript "${SCRIPTS}/unique_variants.R" \
    "${SAMPLE_OUT}/text_file" \
    "${SAMPLE_OUT}/snps_to_evaluate"

  cp "${RESOURCES}/svm_training_longshot.txt" "${SAMPLE_OUT}/snps_to_evaluate/svm_train_data/"
  cp "${RESOURCES}/svm_training_pepper.txt" "${SAMPLE_OUT}/snps_to_evaluate/svm_train_data/"

  Rscript "${SCRIPTS}/svm_longshot.R" "${SAMPLE_OUT}/snps_to_evaluate/"
  Rscript "${SCRIPTS}/svm_pepper.R" "${SAMPLE_OUT}/snps_to_evaluate/"

  "${SCRIPTS}/post_SVM.sh" \
    "${SAMPLE_OUT}/snps_to_evaluate/longshot_pass/${SAMPLE}_longshot_pass.txt" \
    > "${SAMPLE_OUT}/snps_to_evaluate/longshot_pass/${SAMPLE}_longshot_svm_pass_sites"

  "${SCRIPTS}/post_SVM.sh" \
    "${SAMPLE_OUT}/snps_to_evaluate/pepper_pass/${SAMPLE}_pepper_pass.txt" \
    > "${SAMPLE_OUT}/snps_to_evaluate/pepper_pass/${SAMPLE}_pepper_svm_pass_sites"

  # 7. Recover SVM-passed variants and merge with confident variants.
  vcftools --gzvcf "${SAMPLE_OUT}/${SAMPLE}_amplicon.vcf.gz" \
    --positions "${SAMPLE_OUT}/snps_to_evaluate/longshot_pass/${SAMPLE}_longshot_svm_pass_sites" \
    --recode \
    --out "${SAMPLE_OUT}/longshot_final"

  mv "${SAMPLE_OUT}/longshot_final.recode.vcf" \
    "${SAMPLE_OUT}/${SAMPLE}_longshot_SVM_pass.vcf"

  bgzip -c "${SAMPLE_OUT}/${SAMPLE}_longshot_SVM_pass.vcf" \
    > "${SAMPLE_OUT}/${SAMPLE}_longshot_SVM_pass.vcf.gz"
  tabix "${SAMPLE_OUT}/${SAMPLE}_longshot_SVM_pass.vcf.gz"

  bcftools concat -a \
    "${SAMPLE_OUT}/${SAMPLE}_longshot_SVM_pass.vcf.gz" \
    "${SAMPLE_OUT}/confident_snv/${SAMPLE}_confident_snps.vcf.gz" \
    -o "${SAMPLE_OUT}/${SAMPLE}_final_longshot.vcf"

  bgzip -c "${SAMPLE_OUT}/${SAMPLE}_final_longshot.vcf" \
    > "${SAMPLE_OUT}/${SAMPLE}_final_longshot.vcf.gz"
  tabix "${SAMPLE_OUT}/${SAMPLE}_final_longshot.vcf.gz"

  vcftools --gzvcf "${SAMPLE_OUT}/${SAMPLE}_PEPPER_VARIANT_FULL.vcf.gz" \
    --positions "${SAMPLE_OUT}/snps_to_evaluate/pepper_pass/${SAMPLE}_pepper_svm_pass_sites" \
    --recode \
    --out "${SAMPLE_OUT}/pepper_final"

  mv "${SAMPLE_OUT}/pepper_final.recode.vcf" \
    "${SAMPLE_OUT}/${SAMPLE}_PEPPER_SVM_pass.vcf"

  bgzip -c "${SAMPLE_OUT}/${SAMPLE}_PEPPER_SVM_pass.vcf" \
    > "${SAMPLE_OUT}/${SAMPLE}_PEPPER_SVM_pass.vcf.gz"
  tabix "${SAMPLE_OUT}/${SAMPLE}_PEPPER_SVM_pass.vcf.gz"

  longshot \
    -b "${BAM_DIR}/${SAMPLE}.fastq_sorted.minimap2.bam" \
    -f "${REFERENCE}" \
    -v "${SAMPLE_OUT}/${SAMPLE}_PEPPER_SVM_pass.vcf.gz" \
    -o "${SAMPLE_OUT}/${SAMPLE}_pepper_to_longshot.vcf" \
    --sample_id "${SAMPLE}" \
    -E 0.01 -e 2 -C 5000 -c 1 -F -P 0.00001 -n \
    --output-ref -q 5 -B 50 -l 20 -y 10 -a 10

  bgzip -c "${SAMPLE_OUT}/${SAMPLE}_pepper_to_longshot.vcf" \
    > "${SAMPLE_OUT}/${SAMPLE}_pepper_to_longshot.vcf.gz"
  tabix "${SAMPLE_OUT}/${SAMPLE}_pepper_to_longshot.vcf.gz"

  # 8. Final sm-SNIPER callset and haplotype phasing.
  bcftools concat -a \
    "${SAMPLE_OUT}/${SAMPLE}_final_longshot.vcf.gz" \
    "${SAMPLE_OUT}/${SAMPLE}_pepper_to_longshot.vcf.gz" \
    -o "${SAMPLE_OUT}/final_vcf_SNIPER/${SAMPLE}_final_SNIPER.vcf"

  whatshap phase \
    -o "${SAMPLE_OUT}/final_vcf_SNIPER/${SAMPLE}_final_SNIPER_phased.vcf" \
    --reference "${REFERENCE}" \
    "${SAMPLE_OUT}/final_vcf_SNIPER/${SAMPLE}_final_SNIPER.vcf" \
    "${BAM_DIR}/${SAMPLE}_filtered.bam" \
    --ignore-read-groups \
    --mapq 40

  echo "Finished sample: ${SAMPLE}"
done
