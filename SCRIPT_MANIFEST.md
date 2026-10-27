# Script manifest

This release contains one maintained copy of each publication-facing workflow.
Names were shortened where helpful, but analytical code was otherwise retained.

| Release path | Purpose |
|---|---|
| `candidate_iep_inference/` | Bayesian within-host candidate-IEP inference |
| `haplotype_serotype_profiles/` | Haplotype-to-IEP-profile reduction, accumulation curves and serotype-discriminating IEP analysis |
| `genomewide_iep_prediction/` | Residue-level positive-unlabelled random-forest analysis |
| `csp_serotype_representation/` | CSP profile calling, CSP-specific sample clustering, population representation and maps |
| `targeted_long_read_variant_calling/` | Targeted Nanopore SNV calling and SVM filtering |
| `wgs_short_read_variant_calling/` | Four-stage GATK WGS workflow |

Repository support files include:

| Path | Purpose |
|---|---|
| `PURF_residue_level_analysis_input/` | Analysis-ready inputs for the residue-level PURF workflow |
| `docs/wiki/` | Wiki-ready workflow and reproducibility documentation |
| `.github/` | Issue forms and pull-request template |
| `CONTRIBUTING.md` | Contribution and validation guidance |

The PURF scripts were renamed as follows:

| Original name | Release name |
|---|---|
| `run_pf7_purf_gene_model_rtchou.py` | `genomewide_iep_prediction/run_purf_site_model.py` |
| `compute_purf_site_permutation_auc_importance.py` | `genomewide_iep_prediction/compute_purf_permutation_importance.py` |
| `merge_rtchou_site_scores_with_vcf_positions.py` | `genomewide_iep_prediction/add_vcf_coordinates.py` |

The following development artefacts were deliberately excluded:

- `.snakemake` state, scheduler logs, caches and compiled Python files;
- duplicated files under the historical long-read `workflow/` directory;
- generated outputs and temporary files;
- the bundled reference FASTA, which should be obtained from its source and
  identified by release and checksum.

The retained `resources/` files are analysis inputs rather than executable code.
Confirm their redistribution terms before applying a public repository licence.
