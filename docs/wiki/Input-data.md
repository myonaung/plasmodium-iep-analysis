# Input data

## PURF inputs distributed with the repository

The `PURF_residue_level_analysis_input/` directory contains:

- `purf_residue_features.tsv`: residue-level predictor matrix;
- `purf_training_sites.tsv`: candidate-IEP training labels;
- `purf_annotated_ns_snps.vcf.gz`: annotated SNP coordinates;
- the corresponding VCF index and, when supplied, `input_manifest.tsv`.

The VCF is required only to restore genomic coordinates to predicted residues.
It does not contribute predictors or labels during model fitting.

## External inputs

The remaining workflows use MalariaGEN releases, longitudinal cohort sequence
data, reference genomes and accompanying metadata. These files are not
duplicated in the repository. Obtain each dataset from its original provider
under the applicable access terms.

Do not upload participant identifiers or controlled-access source data to
GitHub Issues or pull requests.
