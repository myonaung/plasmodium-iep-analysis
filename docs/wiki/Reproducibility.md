# Reproducibility

Record the following for each analysis run:

- repository commit;
- input release, accession and checksum;
- reference-genome release and checksum;
- software and package versions;
- command-line arguments;
- random seed;
- scheduler resources when relevant.

The released scripts retain manuscript thresholds and seeds. Validate bundled
PURF inputs against `input_manifest.tsv` when that file is present. Generated
results should be written under `results/`, which is excluded from version
control.

Before creating a repository release, run syntax checks, inspect the repository
for private paths and archive the release with a permanent DOI where possible.
