#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash scripts/run_single_antigen_serotype_representation.sh \
    --annotated-summary /path/to/all_bayesian_site_summaries_annotated.csv \
    --vcf PF3D7_1133400_AMA1.vcf.gz \
    --gene PF3D7_1133400 \
    --antigen AMA1 \
    [--site 197] or [--site 352,357] \
    --metadata pf7_clonal_sample_metadata.tsv \
    [--clusters csp_sample_clusters.tsv] \
    --outdir AMA1_serotype_representation \
    [--prob-sym-gt-asym-threshold 0.80] \
    [--prob-sym-gt-expected-threshold 0.2] \
    [--require-high-confidence true] \
    [--require-pf7-ns-snp true] \
    [--require-pv4-snps true] \
    [--require-snp-col Pv4_snps] \
    [--max-dominant-allele-freq 0.95] \
    [--require-sym-gt-expected-over-asym true] \
    [--top-n 10] \
    [--country-top-n 5] \
    [--draws 2000] \
    [--iter 6000] \
    [--warmup 3000] \
    [--chains 4] \
    [--cores 4] \
    [--adapt-delta 0.995] \
    [--max-treedepth 15] \
    [--country-estimation brms|dirichlet]

This analysis wrapper:
  1. selects immune-escape amino-acid positions from the annotated CSV
  2. calls aa-position serotypes from the gene VCF
  3. estimates global Bayesian serotype representation
  4. estimates country-level Bayesian serotype representation

When --site is supplied, the listed amino-acid positions are treated as a
manual serotype-defining panel and bypass posterior/SNP/DAF filters.
EOF
  exit 1
}

ANNOTATED_SUMMARY=""
VCF=""
GENE=""
ANTIGEN=""
SITE=""
METADATA=""
CLUSTERS=""
OUTDIR=""
PROB_SYM_GT_ASYM_THRESHOLD="0.80"
PROB_SYM_GT_EXPECTED_THRESHOLD="0.2"
REQUIRE_HIGH_CONFIDENCE="ignore"
REQUIRE_PF7_NS_SNP="true"
REQUIRE_PV4_SNPS="ignore"
REQUIRE_SNP_COL=""
MAX_DOMINANT_ALLELE_FREQ=""
REQUIRE_SYM_GT_EXPECTED_OVER_ASYM="true"
TOP_N="10"
COUNTRY_TOP_N="5"
DRAWS="2000"
ITER="6000"
WARMUP="3000"
CHAINS="4"
CORES="4"
ADAPT_DELTA="0.995"
MAX_TREEDEPTH="15"
COUNTRY_ESTIMATION="brms"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --annotated-summary) ANNOTATED_SUMMARY="$2"; shift 2 ;;
    --vcf) VCF="$2"; shift 2 ;;
    --gene) GENE="$2"; shift 2 ;;
    --antigen) ANTIGEN="$2"; shift 2 ;;
    --site) SITE="$2"; shift 2 ;;
    --metadata) METADATA="$2"; shift 2 ;;
    --clusters) CLUSTERS="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --prob-sym-gt-asym-threshold) PROB_SYM_GT_ASYM_THRESHOLD="$2"; shift 2 ;;
    --prob_sym_gt_asym_threshold) PROB_SYM_GT_ASYM_THRESHOLD="$2"; shift 2 ;;
    --prob-sym-gt-expected-threshold) PROB_SYM_GT_EXPECTED_THRESHOLD="$2"; shift 2 ;;
    --prob_sym_gt_expected_threshold) PROB_SYM_GT_EXPECTED_THRESHOLD="$2"; shift 2 ;;
    --require-high-confidence) REQUIRE_HIGH_CONFIDENCE="$2"; shift 2 ;;
    --require_high_confidence) REQUIRE_HIGH_CONFIDENCE="$2"; shift 2 ;;
    --require-pf7-ns-snp) REQUIRE_PF7_NS_SNP="$2"; shift 2 ;;
    --require_pf7_ns_snp) REQUIRE_PF7_NS_SNP="$2"; shift 2 ;;
    --require-pv4-snps) REQUIRE_PV4_SNPS="$2"; shift 2 ;;
    --require_pv4_snps) REQUIRE_PV4_SNPS="$2"; shift 2 ;;
    --require-snp-col) REQUIRE_SNP_COL="$2"; shift 2 ;;
    --require_snp_col) REQUIRE_SNP_COL="$2"; shift 2 ;;
    --max-dominant-allele-freq) MAX_DOMINANT_ALLELE_FREQ="$2"; shift 2 ;;
    --max_dominant_allele_freq) MAX_DOMINANT_ALLELE_FREQ="$2"; shift 2 ;;
    --require-sym-gt-expected-over-asym) REQUIRE_SYM_GT_EXPECTED_OVER_ASYM="$2"; shift 2 ;;
    --require_prob_sym_gt_expected_over_asym) REQUIRE_SYM_GT_EXPECTED_OVER_ASYM="$2"; shift 2 ;;
    --top-n) TOP_N="$2"; shift 2 ;;
    --country-top-n) COUNTRY_TOP_N="$2"; shift 2 ;;
    --draws) DRAWS="$2"; shift 2 ;;
    --iter) ITER="$2"; shift 2 ;;
    --warmup) WARMUP="$2"; shift 2 ;;
    --chains) CHAINS="$2"; shift 2 ;;
    --cores) CORES="$2"; shift 2 ;;
    --adapt-delta) ADAPT_DELTA="$2"; shift 2 ;;
    --max-treedepth) MAX_TREEDEPTH="$2"; shift 2 ;;
    --country-estimation) COUNTRY_ESTIMATION="$2"; shift 2 ;;
    --country_estimation) COUNTRY_ESTIMATION="$2"; shift 2 ;;
    *) echo "Unknown option: $1" >&2; usage ;;
  esac
done

if [[ -z "$ANNOTATED_SUMMARY" || -z "$VCF" || -z "$GENE" || -z "$ANTIGEN" || -z "$METADATA" || -z "$OUTDIR" ]]; then
  usage
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HELPER_DIR="${SEROTYPE_REPRESENTATION_HELPER_DIR:-$SCRIPT_DIR}"
if [[ ! -f "$HELPER_DIR/extract_candidate_iep_positions.py" && -f "$(dirname "$SCRIPT_DIR")/extract_candidate_iep_positions.py" ]]; then
  HELPER_DIR="$(dirname "$SCRIPT_DIR")"
fi
required_helpers=(
  "$HELPER_DIR/extract_candidate_iep_positions.py"
  "$HELPER_DIR/call_antigen_serotypes_from_annotated_vcf.py"
  "$HELPER_DIR/estimate_global_serotype_representation_bayesian.R"
  "$HELPER_DIR/estimate_country_serotype_representation_bayesian.R"
)

missing_helpers=()
for helper in "${required_helpers[@]}"; do
  if [[ ! -f "$helper" ]]; then
    missing_helpers+=("$helper")
  fi
done

if [[ ${#missing_helpers[@]} -gt 0 ]]; then
  echo "Missing required helper script(s):" >&2
  printf "  %s\n" "${missing_helpers[@]}" >&2
  exit 1
fi

mkdir -p "$OUTDIR"/{derived,global_bayesian,country_bayesian}

POSITIONS_TSV="$OUTDIR/derived/immune_escape_positions.tsv"
CALLS_TSV="$OUTDIR/derived/aa_positions_calls.tsv"
VARIANT_TSV="$OUTDIR/derived/aa_positions_variants.tsv"
REQUESTED_POSITIONS_TXT="$OUTDIR/derived/requested_aa_positions.txt"
EFFECTIVE_POSITIONS_TXT="$OUTDIR/derived/effective_aa_positions.txt"

echo "[1/4] Extracting immune-escape positions for ${ANTIGEN}"
extract_cmd=(
  python3 "$HELPER_DIR/extract_candidate_iep_positions.py"
  --input "$ANNOTATED_SUMMARY"
  --antigen "$ANTIGEN"
  --gene-id "$GENE"
  --prob-sym-gt-asym-threshold "$PROB_SYM_GT_ASYM_THRESHOLD"
  --prob-sym-gt-expected-threshold "$PROB_SYM_GT_EXPECTED_THRESHOLD"
  --require-high-confidence "$REQUIRE_HIGH_CONFIDENCE"
  --require-pf7-ns-snp "$REQUIRE_PF7_NS_SNP"
  --require-pv4-snps "$REQUIRE_PV4_SNPS"
  --require-snp-col "$REQUIRE_SNP_COL"
  --require-sym-gt-expected-over-asym "$REQUIRE_SYM_GT_EXPECTED_OVER_ASYM"
  --out "$POSITIONS_TSV"
)

if [[ -n "$SITE" ]]; then
  extract_cmd+=(--site "$SITE")
fi

if [[ -n "$MAX_DOMINANT_ALLELE_FREQ" ]]; then
  extract_cmd+=(--max-dominant-allele-freq "$MAX_DOMINANT_ALLELE_FREQ")
fi

AA_POSITIONS="$("${extract_cmd[@]}")"

if [[ -z "$AA_POSITIONS" ]]; then
  echo "No immune-escape positions passed the thresholds for ${ANTIGEN}." >&2
  exit 2
fi

echo "Selected positions: $AA_POSITIONS"
printf "%s\n" "$AA_POSITIONS" > "$REQUESTED_POSITIONS_TXT"

echo "[2/4] Calling AA-position serotypes"
python3 "$HELPER_DIR/call_antigen_serotypes_from_annotated_vcf.py" \
  --vcf "$VCF" \
  --gene "$GENE" \
  --mode aa_positions \
  --aa-positions "$AA_POSITIONS" \
  --out "$CALLS_TSV" \
  --variant-summary "$VARIANT_TSV"

AVAILABLE_POSITIONS="$(
  awk -F'\t' 'NR > 1 && $6 != "" { print $6 }' "$VARIANT_TSV" | sort -n -u | paste -sd, -
)"

if [[ -z "$AVAILABLE_POSITIONS" ]]; then
  echo "No amino-acid variant positions from the selected panel were present in the VCF for ${ANTIGEN}." >&2
  exit 2
fi

if [[ "$AVAILABLE_POSITIONS" != "$AA_POSITIONS" ]]; then
  echo "Reducing requested positions to VCF-available positions: $AVAILABLE_POSITIONS"
  python3 "$HELPER_DIR/call_antigen_serotypes_from_annotated_vcf.py" \
    --vcf "$VCF" \
    --gene "$GENE" \
    --mode aa_positions \
    --aa-positions "$AVAILABLE_POSITIONS" \
    --out "$CALLS_TSV" \
    --variant-summary "$VARIANT_TSV"
fi

printf "%s\n" "$AVAILABLE_POSITIONS" > "$EFFECTIVE_POSITIONS_TXT"

echo "[3/4] Estimating global Bayesian serotype representation"
global_cmd=(
Rscript "$HELPER_DIR/estimate_global_serotype_representation_bayesian.R"
  --calls "$CALLS_TSV" \
  --metadata "$METADATA" \
  --outdir "$OUTDIR/global_bayesian" \
  --top-n "$TOP_N" \
  --draws "$DRAWS" \
  --iter "$ITER" \
  --warmup "$WARMUP" \
  --chains "$CHAINS" \
  --cores "$CORES" \
  --adapt-delta "$ADAPT_DELTA" \
  --max-treedepth "$MAX_TREEDEPTH"
)

echo "[4/4] Estimating country-level Bayesian serotype representation"
country_cmd=(
Rscript "$HELPER_DIR/estimate_country_serotype_representation_bayesian.R"
  --calls "$CALLS_TSV" \
  --metadata "$METADATA" \
  --out-prefix "$OUTDIR/country_bayesian/representation" \
  --top-n "$COUNTRY_TOP_N" \
  --draws "$DRAWS" \
  --iter "$ITER" \
  --warmup "$WARMUP" \
  --chains "$CHAINS" \
  --cores "$CORES" \
  --adapt-delta "$ADAPT_DELTA" \
  --max-treedepth "$MAX_TREEDEPTH" \
  --estimation "$COUNTRY_ESTIMATION"
)

if [[ -n "$CLUSTERS" ]]; then
  global_cmd+=(--clusters "$CLUSTERS")
  country_cmd+=(--clusters "$CLUSTERS")
fi

"${global_cmd[@]}"
"${country_cmd[@]}"

echo "Done. Analysis outputs written to: $OUTDIR"
echo "Use the generated TSV files for downstream serotype-representation plots."
