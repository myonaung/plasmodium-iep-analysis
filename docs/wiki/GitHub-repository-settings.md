# Recommended GitHub repository settings

These settings make effective use of GitHub Pro while keeping the scientific
workflow straightforward.

## Repository features

Under **Settings > General > Features**, enable:

- Issues;
- Discussions, if broader methodological questions are expected;
- Wikis;
- Projects only if the team will actively maintain a task board.

Add a concise repository description and topics such as `malaria`,
`population-genomics`, `immune-escape`, `vaccine-design`, `r` and `python`.

## Default-branch ruleset

Under **Settings > Rules > Rulesets**, create an active branch ruleset targeting
the default branch. Recommended rules are:

- require a pull request before merging;
- require at least one approving review for collaborative changes;
- dismiss stale approvals after new commits;
- require conversation resolution;
- require the `Syntax and privacy checks` and `PURF input schema` status checks;
- block force pushes and branch deletion;
- require linear history if squash or rebase merging is used consistently.

Run the validation workflow once before making its checks mandatory so their
names become selectable in the ruleset interface.

## Code ownership

GitHub Pro supports code-owner review in private repositories. Add
`.github/CODEOWNERS` after confirming contributor GitHub usernames. A suitable
starting structure is:

```text
*                                      @REPOSITORY_MAINTAINER
/candidate_iep_inference/              @IEP_METHOD_MAINTAINER
/genomewide_iep_prediction/            @PURF_MAINTAINER
/targeted_long_read_variant_calling/   @SEQUENCING_MAINTAINER
```

Do not commit placeholder usernames because an invalid CODEOWNERS file gives a
false impression that review is enforced.

## Automated maintenance

The repository includes:

- a validation workflow in `.github/workflows/validate.yml`;
- monthly Dependabot checks for Python and GitHub Actions;
- issue and pull-request templates;
- automatic release-note categories in `.github/release.yml`.

Enable Dependabot alerts and security updates under **Settings > Code security
and analysis**. Review automated dependency pull requests rather than merging
them without rerunning the scientific tests.

## Releases and archiving

Create versioned releases for manuscript submission and publication. Attach
small source-data tables where licensing permits, record checksums, and archive
the release with Zenodo or another DOI-issuing repository. Avoid attaching
controlled-access or participant-level data to GitHub releases.

## Git LFS and Actions usage

The PURF VCF is configured for Git LFS. Confirm that LFS objects are present
after cloning with `git lfs pull`. The validation workflow deliberately avoids
downloading large LFS objects, which reduces Actions storage and bandwidth use.
