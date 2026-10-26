# Haplotype and putative-serotype profiles

This workflow compares full sequence haplotypes with profiles defined only at
candidate IEP positions and produces the flow tables used for serotype analyses.

## Inputs

- **Candidate-site summaries:** one Pf and one Pv CSV containing antigen/gene
  identifiers, amino-acid positions, posterior support values and SNP status.
- **Antigen sequence data:** per-antigen VCF files, or aligned FASTA files, for
  each species.
- **Sample metadata:** sample identifier and geographic grouping columns when
  restricting the analysis to a population; metadata are unnecessary for a
  global sample scope.

## Example

```bash
Rscript haplotype_serotype_profiles/build_haplotype_serotype_profiles.R \
  --pf-summary data/pf_candidate_sites.csv \
  --pv-summary data/pv_candidate_sites.csv \
  --pf-vcf-dir data/pf_antigen_vcfs \
  --pv-vcf-dir data/pv_antigen_vcfs \
  --pf-metadata data/pf_sample_metadata.tsv \
  --pv-metadata data/pv_sample_metadata.tsv \
  --pf-sample-scope png \
  --pv-sample-scope png \
  --max-iep-hamming-distance 2 \
  --output results/haplotype_serotype_profiles.pdf
```

The script also writes summary, flow, raw-flow and diagnostic TSV files. It
does not require a related-sample clustering table.

## Accumulation curves

`scripts/plot_serotype_accumulation.R` uses the raw and displayed flow tables
from the profile-building step to estimate how rapidly distinct putative
serotypes accumulate as sequence haplotypes are sampled.

```bash
Rscript haplotype_serotype_profiles/scripts/plot_serotype_accumulation.R \
  --raw-flows results/profiles_raw_flows.tsv \
  --alluvial-flows results/profiles_flows.tsv \
  --output results/serotype_accumulation.pdf \
  --output-svg results/serotype_accumulation.svg \
  --summary-output results/serotype_accumulation.tsv \
  --min-serotypes-exclusive 1 \
  --max-serotypes-exclusive 11 \
  --max-rare-original-haplotype-percent 91 \
  --permutations 10000 \
  --ribbon-type sd \
  --max-iep-hamming-distance 2 \
  --seed 1
```

## Serotype-discriminating IEPs

The companion scripts prepare an IEP-position map and train antigen-specific
random-forest classifiers to identify positions that distinguish the inferred
putative-serotype profiles.

```bash
Rscript haplotype_serotype_profiles/scripts/make_candidate_iep_position_map.R \
  --pf-summary data/pf_candidate_sites.csv \
  --pv-summary data/pv_candidate_sites.csv \
  --out results/candidate_iep_position_map.tsv

python3 haplotype_serotype_profiles/scripts/prepare_serotype_rf_inputs.py \
  --raw-flows results/profiles_raw_flows.tsv \
  --alluvial-flows results/profiles_flows.tsv \
  --position-map results/candidate_iep_position_map.tsv \
  --out-prefix results/serotype_rf_input

Rscript haplotype_serotype_profiles/scripts/train_serotype_random_forest.R \
  --raw-flows results/serotype_rf_input.raw_flows.tsv \
  --alluvial-flows results/serotype_rf_input.alluvial_flows.tsv \
  --iep-position-map results/serotype_rf_input.iep_position_map.tsv \
  --out-prefix results/serotype_rf \
  --max-iep-hamming-distance 2 \
  --trees 1000 \
  --seed 1
```

`scripts/summarise_haplotype_serotype_profiles.R` optionally combines the
profile summary and flow tables into a compact publication table.
