# Targeted long-read variant calling

## Inputs

- **Pf workflow:** Oxford Nanopore FASTQ files, a matching reference FASTA,
  target-region BED and PEPPER container image.
- **Pv workflow:** sorted and indexed Oxford Nanopore BAM files, a matching
  reference FASTA, optional target BED and Clair3 container image.

Reference sequences and containers are external resources and are not bundled.

## *P. falciparum*

`run_pf_long_read_variant_calling.sh` maps Oxford Nanopore reads, calls target-region SNVs with Longshot and PEPPER, retains concordant calls, applies caller-specific SVM models to discordant calls and phases the final VCF.

```bash
bash run_pf_long_read_variant_calling.sh \
  --fastq-dir data/fastq \
  --reference references/Pfalciparum3D7.fasta \
  --pepper-container containers/pepper_deepvariant_r0.7.sif \
  --outdir results/pf_long_reads
```

The workflow expects sample files named `SAMPLE.fastq` or `SAMPLE.fastq.gz`. Longshot calling uses 5--5,000x depth and a minimum alternate-allele fraction of 1%. PEPPER uses its caller-specific model and the discordant calls from each caller are evaluated with separate radial-kernel SVMs.

## *P. vivax*

`run_pv_long_read_variant_calling.sh` runs Clair3 on sorted, indexed BAM files.

```bash
bash run_pv_long_read_variant_calling.sh \
  --bam-dir data/bam \
  --reference references/Pvivax.fasta \
  --clair3-container containers/clair3_v2.0.0.sif \
  --outdir results/pv_long_reads \
  --bed resources/pv_target_panel.bed
```

## External resources

Reference FASTA files and caller containers are not distributed here. The included SVM training tables and target BED should be reviewed for redistribution permission before creating a public release.
