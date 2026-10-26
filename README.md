# Analysis code

This directory contains publication-facing analysis workflows for the accompanying manuscript. Paths are supplied at run time; no institutional filesystem layout, project allocation or user-specific Conda installation is assumed.

## Contents

- `candidate_iep_inference/`: Bayesian mixed-effects inference of candidate immune-escape polymorphisms.
- `haplotype_serotype_profiles/`: profile construction, serotype accumulation and random-forest position importance.
- `genomewide_iep_prediction/`: positive-unlabelled random-forest prediction of IEP-like residues.
- `csp_serotype_representation/`: CSP putative-serotype calling, optional sample clustering, population-representation analysis and geographic maps.
- `targeted_long_read_variant_calling/`: targeted long-read SNV calling and caller-specific SVM filtering.
- `wgs_short_read_variant_calling/`: GATK whole-genome short-read variant calling.

## General requirements

The scripts use command-line arguments or editable example configuration files. Large input datasets, reference genomes, VCF files and intermediate results are not included. Obtain MalariaGEN data under the applicable data-access terms and download reference genomes from the original provider.

The principal Python dependencies are listed in `requirements.txt`. PURF itself
must also be installed for the residue-level model. Genomic workflows require
the command-line programs named in their README files or dependency checks. R
package requirements are declared near the beginning of each R script.

## PURF residue-level analysis

Run locally or within an allocated compute job:

```bash
bash genomewide_iep_prediction/run_genomewide_iep_prediction.sh \
  --features data/residue_features.tsv \
  --training-sites data/training_sites.tsv \
  --vcf data/annotated_snps.vcf.gz \
  --outdir results/genomewide_iep_prediction \
  --threads 4
```

The `genomewide_iep_prediction/submit_genomewide_iep_prediction.pbs.example`
file is an optional PBS template. Scheduler allocation, storage and environment
activation are intentionally left to the user.

## Candidate immune-escape analysis

`candidate_iep_inference/infer_candidate_ieps_bayesian.R` is a sourceable R module. It constructs all eligible
within-participant infection pairs, fits the Bayesian mixed-effects model and
writes pairwise, model-input, posterior-summary and candidate-site outputs.

```r
source("candidate_iep_inference/infer_candidate_ieps_bayesian.R")

result <- allele_specific_immunity_bayesian(
  data_dir = "data/antigen",
  output_dir = "results/immune_escape",
  fasta_file = "aligned_sequences.fasta",
  metadata_file = "sample_metadata.csv",
  molfoi_file = "moi_estimates.csv",
  Antigen = "AMA1",
  min_prob_sym_gt_asym = 0.80
)
```

The function default for its high-confidence export is `0.95`. The downstream
cross-antigen classification used for the manuscript is configured explicitly
in `haplotype_serotype_profiles/build_haplotype_serotype_profiles.R` at `0.80`; pass the intended threshold explicitly
when recreating manuscript outputs.

## Haplotype and putative-serotype inputs

`haplotype_serotype_profiles/build_haplotype_serotype_profiles.R` combines filtered Pf7 and Pv4 antigen data with
candidate-IEP summaries, collapses sequence haplotypes by their IEP-site
profiles and writes plotting, flow and diagnostic tables. The required
arguments and a generic command are documented in the script header; all input
paths are supplied at run time.

## Reproducibility notes

- Random seeds and analytical thresholds remain explicit in the scripts.
- SVM training tables and the target-panel BED used by the long-read workflow are included as analysis resources. Confirm that these resource files may be redistributed under the repository's chosen licence before public release.
- The reference genome is intentionally not bundled. Record its release and checksum in the repository release notes.
- Generated outputs, scheduler logs, caches and large sequence files are excluded through `.gitignore`.

## Before depositing a release

1. Add the repository licence agreed by all contributors.
2. Add the manuscript citation and permanent archive DOI when available.
3. Record tested software versions and input-data accession or release identifiers.
4. Confirm redistribution permission for the two SVM training tables and target-panel BED.
5. Verify that threshold values and reference releases match the final Methods.
