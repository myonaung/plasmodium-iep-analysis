# Workflow overview

## Candidate IEP inference

Longitudinal sequences and infection metadata are converted into eligible
within-participant infection pairs. A Bayesian mixed-effects model estimates
site-level support for greater allele turnover before symptomatic than
asymptomatic infections and relative to antigen-specific background turnover.

## Haplotype and putative-serotype profiling

Full sequence haplotypes are collapsed by their amino-acid profiles at
candidate IEP positions. Accumulation analyses and random-forest importance
summarise serotype diversity and the residues that distinguish profiles.

## Genome-wide PURF prediction

Known candidate IEP sites are treated as labelled positives and other
nonsynonymous residues as unlabelled. PURF ranks residues with similar
immunological, structural and population-genetic features, after which genomic
coordinates are restored from the annotated VCF.

## CSP representation and mapping

CSP profiles are called from annotated VCF data and summarised globally and by
country. Optional sample clustering reduces over-representation of closely
related infections during serotype ranking. Raw prevalence remains based on
all callable samples.

## Variant-processing workflows

The long-read workflow combines concordant caller output with caller-specific
SVM filtering. The short-read workflow provides the GATK processing stages used
for whole-genome variant calling.
