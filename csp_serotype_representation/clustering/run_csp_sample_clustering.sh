#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: run_csp_sample_clustering.sh --metadata FILE --vcf FILE --outdir DIR [options]

Options:
  --threads N                 bcftools threads (default: 4)
  --distance-threshold X      Maximum within-cluster distance (default: 0.001)
  --year-threshold N          Maximum sampling-year difference (default: 0)
EOF
}

METADATA=""
VCF=""
OUTDIR=""
THREADS=4
DISTANCE_THRESHOLD=0.001
YEAR_THRESHOLD=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --metadata) METADATA="$2"; shift 2 ;;
    --vcf) VCF="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --distance-threshold) DISTANCE_THRESHOLD="$2"; shift 2 ;;
    --year-threshold) YEAR_THRESHOLD="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

for value in METADATA VCF OUTDIR; do
  [[ -n "${!value}" ]] || { echo "ERROR: missing --${value,,}" >&2; usage >&2; exit 2; }
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HELPERS="${SCRIPT_DIR}/scripts"
mkdir -p "${OUTDIR}"

for cmd in Rscript python3 bcftools; do
  command -v "${cmd}" >/dev/null 2>&1 || { echo "ERROR: ${cmd} is required" >&2; exit 1; }
done
for file in "${METADATA}" "${VCF}"; do
  [[ -f "${file}" ]] || { echo "ERROR: input not found: ${file}" >&2; exit 1; }
done

Rscript "${HELPERS}/prepare_pf7_metadata.R" \
  --pf7-metadata "${METADATA}" \
  --out "${OUTDIR}/csp_sample_metadata.tsv" \
  --sample-col Sample --study-col Study --country-col Country \
  --admin1-col "Admin level 1" --region-col Population --year-col Year \
  --qc-col "QC pass" --qc-keep TRUE

python3 "${HELPERS}/distance_matrix.py" \
  --vcf "${VCF}" \
  --out-prefix "${OUTDIR}/csp_sample_distance" \
  --bcftools bcftools --threads "${THREADS}" \
  --site-chunk 256 --sample-block 1024 \
  --min-total-depth 5 --min-allele-depth 2 --min-allele-frac 0.05

python3 "${HELPERS}/cluster_similar_samples.py" \
  --matrix "${OUTDIR}/csp_sample_distance.npy" \
  --sample-ids "${OUTDIR}/csp_sample_distance.sample_ids.txt" \
  --metadata "${OUTDIR}/csp_sample_metadata.tsv" \
  --out-clusters "${OUTDIR}/csp_sample_clusters.tsv" \
  --out-pairs "${OUTDIR}/csp_linked_sample_pairs.tsv" \
  --distance-threshold "${DISTANCE_THRESHOLD}" \
  --site-col site --year-col sample_year --year-threshold "${YEAR_THRESHOLD}"

echo "Cluster table: ${OUTDIR}/csp_sample_clusters.tsv"
