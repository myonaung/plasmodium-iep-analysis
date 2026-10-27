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

The `example_data/` directory contains a reduced demonstration dataset. All
participant and sample identifiers have been replaced, sampling times are
represented only by within-participant order, unnecessary cohort variables
have been removed, and the exposure values are synthetic. The aligned AMA1
sequences contain only residues 150-329 and are provided solely to demonstrate
the input structure and workflow.

## Example

```r
source("candidate_iep_inference/infer_candidate_ieps_bayesian.R")

result <- allele_specific_immunity_bayesian(
  data_dir = "candidate_iep_inference/example_data",
  output_dir = "results/candidate_ieps",
  fasta_file = "ama1_example_alignment.fasta",
  metadata_file = "sample_metadata.csv",
  molfoi_file = "exposure_covariate.csv",
  Antigen = "AMA1_EXAMPLE",
  sample_col = "sample_id",
  studyid_col = "participant_id",
  timepoint_col = "timepoint",
  status_col = "clinical_status",
  molfoi_col = "molFOI",
  min_prob_sym_gt_asym = 0.80
)
```

The public example can be regenerated from authorised source files with
`prepare_public_example_inputs.py`. No identifier crosswalk is written:

```bash
python candidate_iep_inference/prepare_public_example_inputs.py \
  --metadata /path/to/sample_metadata.csv \
  --exposure /path/to/exposure_covariate.csv \
  --fasta /path/to/aligned_sequences.fasta \
  --outdir candidate_iep_inference/example_data
```

Before releasing any derived sequence data, confirm that sharing is permitted
by the relevant participant consent, ethics approval and data-use agreement.

Outputs include infection-pair comparisons, model-ready site records,
posterior site summaries, candidate sites, the fitted model and a summary plot.
