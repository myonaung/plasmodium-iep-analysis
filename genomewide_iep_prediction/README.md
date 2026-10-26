# Residue-level PURF analysis

The workflow trains a positive-unlabelled random forest on known candidate IEP sites, ranks unlabelled nonsynonymous sites, estimates permutation importance with a supervised surrogate model, and links high-scoring amino-acid positions to genomic coordinates in an annotated VCF.

## Inputs

- **Residue-feature TSV:** one nonsynonymous residue per row with gene ID,
  amino-acid position and numeric immunological, structural or
  population-genetic predictors.
- **Training-site TSV:** gene ID, amino-acid position and a binary indicator for
  candidate IEP training sites.
- **Annotated VCF:** SNP records with SnpEff `ANN` annotations, used to restore
  genomic coordinates to ranked residues.

## Example

```bash
bash run_genomewide_iep_prediction.sh \
  --features data/residue_features.tsv \
  --training-sites data/candidate_iep_training_sites.tsv \
  --vcf data/annotated_snps.vcf.gz \
  --outdir results/genomewide_iep_prediction \
  --threads 4
```

`run_genomewide_iep_prediction.sh` exposes the required paths as command-line arguments. It does not activate a Conda environment or assume a scheduler.

Required Python packages include NumPy, pandas, scikit-learn and the package providing `purf.pu_ensemble.PURandomForestClassifier`. `bcftools` is required to link residue scores to SnpEff `ANN` annotations.

The PBS file is an example submission wrapper only. Set scheduler allocation options at submission time and pass `REPOSITORY_DIR`, `FEATURES`, `TRAINING_SITES`, `VCF` and `OUTDIR` as environment variables.
