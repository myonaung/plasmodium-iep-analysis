# Plasmodium immune-escape analysis

Publication-facing workflows for identifying candidate immune-escape
polymorphisms (IEPs), reducing parasite haplotypes to putative serotypes and
prioritising immune-escape-like residues across the *Plasmodium* genome.

The repository accompanies the manuscript **Genomic signatures of immune
escape reveal antigenic redundancy within Plasmodium sequence diversity**.

## Workflows

| Workflow | Purpose |
|---|---|
| [`iep_inference/`](iep_inference/) | Bayesian mixed-effects inference of  IEPs from longitudinal infection pairs |
| [`haplotype_serotype_profiles/`](haplotype_serotype_profiles/) | Haplotype reduction, serotype accumulation and random-forest position importance |
| [`genomewide_iep_prediction/`](genomewide_iep_prediction/) | Positive-unlabelled random-forest prediction of IEP-like residues |
| [`csp_serotype_representation/`](csp_serotype_representation/) | CSP putative-serotype calling, optional sample clustering and geographic representation |
| [`targeted_long_read_variant_calling/`](targeted_long_read_variant_calling/) | Targeted Nanopore SNV calling and caller-specific SVM filtering |
| [`wgs_short_read_variant_calling/`](wgs_short_read_variant_calling/) | GATK whole-genome short-read variant calling |



## Requirements

Paths are supplied at run time; no institutional filesystem layout, scheduler
allocation or user-specific Conda installation is assumed. The principal
Python dependencies are listed in [`requirements.txt`](requirements.txt).
PURF must be installed separately for the residue-level model. Genomic
workflows also require the command-line tools named in their README files.
R-package requirements are declared at the beginning of each R script.

Large MalariaGEN datasets and reference genomes are not duplicated here.
Obtain them under the applicable source-data access terms and record the
release and reference checksums used for each analysis.

## IEP inference

`iep_inference/infer_candidate_ieps_bayesian.R` constructs eligible
within-participant infection pairs, fits the Bayesian mixed-effects model and
writes pair-level, model-input, posterior-summary and candidate-site outputs.

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

The function default for its high-confidence export is `0.95`. The manuscript
cross-antigen classification is configured explicitly at `0.80` in
`haplotype_serotype_profiles/build_haplotype_serotype_profiles.R`; pass the
intended threshold explicitly when recreating manuscript outputs.

## Haplotype and putative-serotype profiling

`haplotype_serotype_profiles/build_haplotype_serotype_profiles.R` combines
filtered Pf7 and Pv4 antigen data with IEP summaries, collapses
sequence haplotypes by their IEP-site profiles and writes plotting, flow and
diagnostic tables. Its required arguments and an example command are provided
in the script header and workflow README.

## PURF residue-level analysis

The required PURF inputs are stored in
[`PURF_residue_level_analysis_input/`](PURF_residue_level_analysis_input/):

- `purf_residue_features.tsv`
- `purf_training_sites.tsv`
- `purf_annotated_ns_snps.vcf.gz` and its index

Run from the repository root:

```bash
bash genomewide_iep_prediction/run_genomewide_iep_prediction.sh \
  --features PURF_residue_level_analysis_input/purf_residue_features.tsv \
  --training-sites PURF_residue_level_analysis_input/purf_training_sites.tsv \
  --vcf PURF_residue_level_analysis_input/purf_annotated_ns_snps.vcf.gz \
  --outdir results/genomewide_iep_prediction \
  --threads 4
```

The PBS example uses these repository-relative inputs by default. The feature
and training tables are analysis-ready derivatives; provenance and column
requirements are documented in the input-directory README.

## Documentation

- [`SCRIPT_MANIFEST.md`](SCRIPT_MANIFEST.md): released scripts and renamed files
- [`CONTRIBUTING.md`](CONTRIBUTING.md): contributions and reproducibility checks


## Reproducibility

- Random seeds and analytical thresholds remain explicit in the scripts.
- Generated results, scheduler logs, caches and large general sequence files
  are excluded through `.gitignore`.
- The PURF input VCF is tracked specifically through Git LFS; run
  `git lfs install` before adding it.
- The reference genome is not bundled. Record its provider, release and
  checksum in the archived repository release.

## Support

Use GitHub Issues for reproducible code problems and workflow questions. Do
not post controlled-access data, participant information or credentials in an
issue.
