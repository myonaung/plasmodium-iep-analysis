#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  run_genomewide_iep_prediction.sh --features FILE --training-sites FILE --vcf FILE \
    --outdir DIR [--threads N] [--prefix NAME] [--python PYTHON] [--bcftools PATH]

Required inputs:
  --features        Residue-level feature table (TSV)
  --training-sites  Candidate-IEP training-site table (TSV)
  --vcf             SnpEff-annotated SNP VCF containing INFO/ANN
  --outdir          Output directory

Optional:
  --threads         Worker processes (default: 4)
  --prefix          Output prefix basename (default: pf7_purf_iep_like_sites)
  --python          Python executable (default: python3)
  --bcftools        bcftools executable (default: bcftools)
EOF
}

FEATURES=""
TRAINING_SITES=""
VCF=""
OUTDIR=""
THREADS=4
PREFIX_NAME="pf7_purf_iep_like_sites"
PYTHON_BIN="python3"
BCFTOOLS_BIN="bcftools"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --features) FEATURES="$2"; shift 2 ;;
    --training-sites) TRAINING_SITES="$2"; shift 2 ;;
    --vcf) VCF="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --prefix) PREFIX_NAME="$2"; shift 2 ;;
    --python) PYTHON_BIN="$2"; shift 2 ;;
    --bcftools) BCFTOOLS_BIN="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

for value in FEATURES TRAINING_SITES VCF OUTDIR; do
  if [[ -z "${!value}" ]]; then
    echo "ERROR: --$(tr '[:upper:]_' '[:lower:]-' <<<"$value") is required" >&2
    usage >&2
    exit 2
  fi
done

for file in "$FEATURES" "$TRAINING_SITES" "$VCF"; do
  [[ -s "$file" ]] || { echo "ERROR: input not found or empty: $file" >&2; exit 1; }
done

command -v "$PYTHON_BIN" >/dev/null || { echo "ERROR: Python not found: $PYTHON_BIN" >&2; exit 1; }
command -v "$BCFTOOLS_BIN" >/dev/null || { echo "ERROR: bcftools not found: $BCFTOOLS_BIN" >&2; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"
PREFIX="${OUTDIR}/${PREFIX_NAME}"

export TMPDIR="${TMPDIR:-${OUTDIR}/tmp}"
export JOBLIB_TEMP_FOLDER="${JOBLIB_TEMP_FOLDER:-${TMPDIR}/joblib}"
mkdir -p "$TMPDIR" "$JOBLIB_TEMP_FOLDER"

echo "1/3 Fitting residue-level PURF model"
"$PYTHON_BIN" "${SCRIPT_DIR}/run_purf_site_model.py" \
  --positive-unit site \
  --input "$FEATURES" \
  --out-prefix "$PREFIX" \
  --gene-col gene_id \
  --label-col gene_id \
  --position-col position_std \
  --training-sites "$TRAINING_SITES" \
  --training-gene-col gene_id \
  --training-positive-col train_positive \
  --training-site-col position \
  --min-positive-sites 1 \
  --exclude-regex '(^gene_|_key$|Gene_Name_annotated|^position$|position_std|aa_position|site_index|positive|known_iep|candidate_iep|potential_iep|positive_or_potential|training|training_source|iep_evidence|prob_|sym_|asym_|expected|dominant_allele|background_|purf_|reference|tree_filter|^antigen_label$|^sp_prob$|^gpi_specificity$|signal|gpi|t[_\. -]?cell|tcell|t_epitope|ctl_epitope|cd4|cd8|mhc|hla|netmhc|netctl|net_charge|tm|c?nonsyn_syn_snp_ratio)' \
  --cor-cutoff 0.70 \
  --prefer-feature-regex 'TajimaD|rsa|B_epitope|IUPRED2_score|discotope|ss_is_H|ss_is_E|ss_is_C' \
  --positive-levels 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 \
  --select-positive-level-metric reference_auc \
  --n-trees-tuning 500 \
  --n-trees-final 2000 \
  --max-features sqrt \
  --min-samples-leaf 1 \
  --variable-space-weighting TRUE \
  --max-positive-replicates 8 \
  --positive-replicate-power 1 \
  --tree-filter FALSE \
  --min-trees-after-filter 100 \
  --tree-filter-min-reference-recall 0.80 \
  --tree-filter-fallback best \
  --threads "$THREADS" \
  --top-n 250

echo "2/3 Computing permutation importance"
"$PYTHON_BIN" "${SCRIPT_DIR}/compute_purf_permutation_importance.py" \
  --scores "${PREFIX}.all_gene_scores.tsv" \
  --feature-importance "${PREFIX}.feature_importance.tsv" \
  --out-prefix "${PREFIX}.site_permutation_accuracy" \
  --known-col purf_reference_positive \
  --score-col purf_score \
  --top-n-features 0 \
  --negative-score-max 0.20 \
  --negative-ratio 20 \
  --max-negatives 50000 \
  --n-estimators 500 \
  --n-repeats 10 \
  --scoring accuracy \
  --threads "$THREADS"

echo "3/3 Adding genomic coordinates"
"$PYTHON_BIN" "${SCRIPT_DIR}/add_vcf_coordinates.py" \
  --site-scores "${PREFIX}.all_gene_scores.tsv" \
  --vcf "$VCF" \
  --out "${PREFIX}.predicted_iep_like_snps.tsv" \
  --simple-out "${PREFIX}.predicted_iep_like_snps.simple.tsv" \
  --score-threshold 0.80 \
  --top-n 250 \
  --gene-col source_gene_id_normalised \
  --position-col position_std \
  --score-col purf_score \
  --rank-col purf_rank \
  --known-col purf_reference_positive \
  --bcftools "$BCFTOOLS_BIN"

echo "Completed. Outputs use prefix: $PREFIX"
