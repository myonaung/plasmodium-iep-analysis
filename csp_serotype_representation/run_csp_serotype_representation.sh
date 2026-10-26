#!/usr/bin/env bash
set -euo pipefail

# CSP putative-serotype representation analysis.
# This wrapper calls immune-escape profiles and estimates global and
# country-level representation. Plotting can be performed from the output TSVs.

WORKFLOW_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="${WORKFLOW_DIR}/scripts"

# -------------------------------------------------------------------------
# User inputs
# -------------------------------------------------------------------------
usage() {
  cat <<'EOF'
Usage: run_csp_serotype_representation.sh --annotated-summary FILE --vcf FILE \
  --metadata FILE --outdir DIR [options]

Options:
  --clusters FILE       Optional CSP sample-cluster table
  --gene ID             Gene identifier (default: PF3D7_0304600)
  --antigen NAME        Antigen label (default: CSP)
  --sites LIST          Comma-separated amino-acid positions (default: 352,357)
  --cores N             Parallel chains/cores (default: 4)
EOF
}

ANNOTATED_SUMMARY=""
VCF=""
METADATA=""
CLUSTERS=""
OUTDIR=""

# Example shown for PfCSP C-terminal LD haplotype sites used in vaccine
# coverage/sieve analyses. For data-driven candidate IEPs, leave SITE empty
# and use the posterior/SNP/DAF filters below.
GENE="PF3D7_0304600"
ANTIGEN="CSP"
SITE="352,357"

# Candidate IEP filters used when SITE is empty.
PROB_SYM_GT_ASYM_THRESHOLD="0.80"
PROB_SYM_GT_EXPECTED_THRESHOLD="0"
REQUIRE_SYM_GT_EXPECTED_OVER_ASYM="true"
MAX_DOMINANT_ALLELE_FREQ="0.99"

TOP_N="10"
COUNTRY_TOP_N="5"
COUNTRY_ESTIMATION="dirichlet"

DRAWS="2000"
ITER="4000"
WARMUP="2000"
CHAINS="4"
CORES="4"
ADAPT_DELTA="0.995"
MAX_TREEDEPTH="15"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --annotated-summary) ANNOTATED_SUMMARY="$2"; shift 2 ;;
    --vcf) VCF="$2"; shift 2 ;;
    --metadata) METADATA="$2"; shift 2 ;;
    --clusters) CLUSTERS="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --gene) GENE="$2"; shift 2 ;;
    --antigen) ANTIGEN="$2"; shift 2 ;;
    --sites) SITE="$2"; shift 2 ;;
    --cores) CORES="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

for value in ANNOTATED_SUMMARY VCF METADATA OUTDIR; do
  [[ -n "${!value}" ]] || { echo "ERROR: required argument missing: $value" >&2; usage >&2; exit 2; }
done

RESULT_DIR="${OUTDIR}/${ANTIGEN}_serotype_representation"
mkdir -p "${OUTDIR}"

for cmd in python3 Rscript awk sort paste; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "ERROR: required command not found on PATH: ${cmd}" >&2
    exit 1
  fi
done

for f in "${ANNOTATED_SUMMARY}" "${VCF}" "${METADATA}"; do
  if [ ! -f "${f}" ]; then
    echo "ERROR: input file not found: ${f}" >&2
    exit 1
  fi
done

representation_cmd=(
  bash "${SCRIPTS}/run_single_antigen_serotype_representation.sh"
  --annotated-summary "${ANNOTATED_SUMMARY}"
  --vcf "${VCF}"
  --gene "${GENE}"
  --antigen "${ANTIGEN}"
  --metadata "${METADATA}"
  --outdir "${RESULT_DIR}"
  --prob-sym-gt-asym-threshold "${PROB_SYM_GT_ASYM_THRESHOLD}"
  --prob-sym-gt-expected-threshold "${PROB_SYM_GT_EXPECTED_THRESHOLD}"
  --require-pf7-ns-snp true
  --require-sym-gt-expected-over-asym "${REQUIRE_SYM_GT_EXPECTED_OVER_ASYM}"
  --max-dominant-allele-freq "${MAX_DOMINANT_ALLELE_FREQ}"
  --top-n "${TOP_N}"
  --country-top-n "${COUNTRY_TOP_N}"
  --country-estimation "${COUNTRY_ESTIMATION}"
  --draws "${DRAWS}"
  --iter "${ITER}"
  --warmup "${WARMUP}"
  --chains "${CHAINS}"
  --cores "${CORES}"
  --adapt-delta "${ADAPT_DELTA}"
  --max-treedepth "${MAX_TREEDEPTH}"
)

if [ -n "${SITE}" ]; then
  representation_cmd+=(--site "${SITE}")
fi

if [ -n "${CLUSTERS}" ]; then
  [[ -f "${CLUSTERS}" ]] || { echo "ERROR: cluster file not found: ${CLUSTERS}" >&2; exit 1; }
  representation_cmd+=(--clusters "${CLUSTERS}")
fi

"${representation_cmd[@]}"

echo "CSP serotype-representation outputs written to: ${RESULT_DIR}"
