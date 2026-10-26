# CSP sample clustering

This optional workflow groups highly similar infections sampled at the same site
and time. The resulting `sample_id`/`cluster_id` table can be supplied to the CSP
representation workflow to prevent closely related infections from dominating
serotype ranking. Raw prevalence is still estimated from all callable samples.

## Inputs

- Pf7 sample metadata with sample, site, country, region, year and QC fields.
- A PASS-filtered multisample VCF containing `FORMAT/AD` values.

## Example

```bash
bash clustering/run_csp_sample_clustering.sh \
  --metadata data/Pf7_samples.txt \
  --vcf data/pf7_pass_snps.vcf.gz \
  --outdir results/csp_sample_clustering \
  --threads 4
```

The main output is `csp_sample_clusters.tsv`, containing `sample_id` and
`cluster_id` columns required by `--clusters`.
