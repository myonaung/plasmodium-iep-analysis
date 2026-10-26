#!/bin/bash

##############################################################################
# run_clair3_variant_calling.sh
#
# Clair3 v2.0.0 variant calling for Oxford Nanopore (ONT) reads, run inside
# a Singularity container. Loops over every sorted, indexed BAM in
# INPUT_DIR and writes one Clair3 output folder per sample under
# OUTPUT_DIR.
#
# Generic standalone script — no job-scheduler directives. If you're
# running this on an HPC cluster via a scheduler (PBS/Slurm/etc.), add
# your own resource-request directives/wrapper around this script, or
# submit it directly if your scheduler accepts plain shell scripts.
#
# Supply input paths on the command line; no scheduler or module layout is assumed.
##############################################################################

set -euo pipefail

INPUT_DIR=""
OUTPUT_DIR=""
REF=""
BED_FILE=""
SINGULARITY_IMAGE=""
CONTAINER_ENGINE="singularity"

# Clair3 basecaller model, as named inside the container's /opt/models/.
MODEL_NAME="r1041_e82_400bps_sup_v500"

THREADS=4

usage() {
  cat <<'EOF'
Usage: run_pv_long_read_variant_calling.sh --bam-dir DIR --reference FASTA \
  --clair3-container IMAGE --outdir DIR [options]

Options:
  --bed FILE             Optional target BED
  --model NAME           Clair3 model name
  --threads N            Threads (default: 4)
  --container-engine CMD Container engine (default: singularity)
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam-dir) INPUT_DIR="$2"; shift 2 ;;
    --reference) REF="$2"; shift 2 ;;
    --clair3-container) SINGULARITY_IMAGE="$2"; shift 2 ;;
    --outdir) OUTPUT_DIR="$2"; shift 2 ;;
    --bed) BED_FILE="$2"; shift 2 ;;
    --model) MODEL_NAME="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --container-engine) CONTAINER_ENGINE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

for value in INPUT_DIR REF SINGULARITY_IMAGE OUTPUT_DIR; do
  [[ -n "${!value}" ]] || { echo "ERROR: required argument missing: $value" >&2; usage >&2; exit 2; }
done

# ---------------- Sanity checks ----------------

if ! command -v "${CONTAINER_ENGINE}" >/dev/null 2>&1; then
  echo "ERROR: container engine not found on PATH: ${CONTAINER_ENGINE}" >&2
  exit 1
fi

if [ ! -d "${INPUT_DIR}" ]; then
  echo "ERROR: INPUT_DIR not found: ${INPUT_DIR}" >&2
  exit 1
fi

if [ ! -f "${REF}" ]; then
  echo "ERROR: reference FASTA not found: ${REF}" >&2
  exit 1
fi

if [ ! -f "${SINGULARITY_IMAGE}" ]; then
  echo "ERROR: Clair3 Singularity image not found: ${SINGULARITY_IMAGE}" >&2
  exit 1
fi

# ---------------- Start ----------------

mkdir -p "${OUTPUT_DIR}"

echo "==========================================="
echo "Clair3 v2.0.0 STARTED"
echo "Date: $(date)"
echo "==========================================="

cd "${INPUT_DIR}" || exit 1

# ---------------- Loop BAM files ----------------

for file in *.minimap2.sorted.bam
do
    [ -e "$file" ] || { echo "No BAM files found"; exit 1; }

    BASE=$(basename "${file}" .minimap2.sorted.bam)

    echo "-------------------------------------------"
    echo "Processing: ${BASE}"

    SAMPLE_OUT="${OUTPUT_DIR}/${BASE}"
    mkdir -p "${SAMPLE_OUT}"

    # Check BAM index
    if [ ! -f "${file}.bai" ]; then
        echo "Missing index for ${file} -> skipping"
        continue
    fi

    # ---------------- Run Clair3 ----------------

    bed_args=()
    if [[ -n "${BED_FILE}" ]]; then
        [[ -f "${BED_FILE}" ]] || { echo "ERROR: BED file not found: ${BED_FILE}" >&2; exit 1; }
        bed_args+=(--bed_fn="${BED_FILE}")
    fi

    "${CONTAINER_ENGINE}" exec "${SINGULARITY_IMAGE}" \
        /opt/bin/run_clair3.sh \
        --bam_fn="${INPUT_DIR}/${file}" \
        --ref_fn="${REF}" \
        --threads="${THREADS}" \
        --platform=ont \
        --model_path="/opt/models/${MODEL_NAME}" \
        --output="${SAMPLE_OUT}" \
        --sample_name="${BASE}" \
        "${bed_args[@]}" \
        --include_all_ctgs \
        --no_phasing_for_fa \
        --haploid_sensitive \
        --print_ref_calls

    echo "Finished: ${BASE}"
done

echo "==========================================="
echo "ALL SAMPLES COMPLETED"
echo "==========================================="

# For more stringent variant calling, consider adding:
#        --snp_min_af=0.0 \
#        --indel_min_af=0.0 \
