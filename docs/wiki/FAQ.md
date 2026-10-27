# Frequently asked questions

## Are non-IEP polymorphisms true negatives?

No. In the PURF analysis they are treated as unlabelled because they may arise
through immune selection, other selective pressures or neutral processes.

## Are putative serotypes experimentally validated serotypes?

No. They are computational groupings defined by shared amino-acid profiles at
candidate IEP positions and should be interpreted as predicted immune-relevant
antigenic groups.

## Is CSP sample clustering required for every workflow?

No. It is specific to CSP population-representation and mapping analyses. It
is not required for candidate-IEP inference or general haplotype profiling.

## Why is the PURF VCF stored with Git LFS?

Compressed genomic files can exceed standard GitHub file limits. Git LFS keeps
the repository history manageable while preserving a versioned pointer to the
file.

## Where should outputs be written?

Use a directory under `results/` or an external analysis location. Generated
outputs are intentionally excluded from source control.
