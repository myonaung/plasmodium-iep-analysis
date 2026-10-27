# Candidate IEP inference

This analysis compares allele turnover across within-participant infection
pairs preceding symptomatic and asymptomatic infections using a Bayesian
mixed-effects model.

## Inputs

- **Aligned antigen FASTA:** one aligned sequence per sampled infection. The
  sequence header begins with the sample identifier used in the metadata.
- **Sample metadata CSV:** one row per sample, with sample ID, participant ID,
  sampling time and clinical status columns.
- **Exposure covariate CSV:** one row per participant containing molecular
  force of infection or an equivalent numeric exposure covariate.

Column names are configurable, so the public repository does not require the
original cohort naming scheme.

## Example

```r
source("candidate_iep_inference/infer_candidate_ieps_bayesian.R")

result <- allele_specific_immunity_bayesian(
  data_dir = "data/example_antigen",
  output_dir = "results/candidate_ieps",
  fasta_file = "aligned_sequences.fasta",
  metadata_file = "sample_metadata.csv",
  molfoi_file = "exposure_covariate.csv",
  Antigen = "ANTIGEN_NAME",
  sample_col = "sample_id",
  studyid_col = "participant_id",
  timepoint_col = "timepoint",
  status_col = "clinical_status",
  molfoi_col = "molFOI",
  min_prob_sym_gt_asym = 0.80
)
```

Outputs include infection-pair comparisons, model-ready site records,
posterior site summaries, candidate sites, the fitted model and a summary plot.
