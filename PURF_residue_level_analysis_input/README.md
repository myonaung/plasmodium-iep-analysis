# PURF residue-level analysis inputs

This directory contains the analysis-ready inputs used by the genome-wide
residue-level PURF workflow.

| File | Role in the analysis |
|---|---|
| `purf_residue_features.tsv` | Nonsynonymous residue-level predictors used to train and apply PURF |
| `purf_training_sites.tsv` | Gene, amino-acid position and binary candidate-IEP training label |
| `purf_annotated_ns_snps.vcf.gz` | SnpEff-annotated nonsynonymous SNPs used to restore genomic coordinates to ranked residues |
| `purf_annotated_ns_snps.vcf.gz.csi` or `.tbi` | Index for the compressed VCF |
| `input_manifest.tsv` | File sizes and SHA-256 checksums for release verification |

The VCF is used only during coordinate annotation; it is not used to train the
PURF model. The input files should not contain participant identifiers or other
directly identifying metadata.

Because GitHub applies a per-file size limit, the compressed VCF and its index
are configured for Git LFS in the repository `.gitattributes` file. Install Git
LFS before adding or cloning these files:

```bash
git lfs install
git lfs pull
```

Run the workflow from the repository root as documented in
[`../genomewide_iep_prediction/README.md`](../genomewide_iep_prediction/README.md).
