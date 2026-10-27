# Publishing the GitHub Wiki

First enable **Wikis** under the GitHub repository's **Settings > General >
Features** section. Create the first wiki page in GitHub, then clone the wiki's
separate Git repository:

```bash
git clone git@github.com:OWNER/REPOSITORY.wiki.git
```

Copy the prepared pages and push them:

```bash
cp docs/wiki/*.md REPOSITORY.wiki/
cd REPOSITORY.wiki
git add .
git commit -m "Add analysis documentation"
git push
```

Replace `OWNER` and `REPOSITORY` with the final GitHub organisation or account
and repository name. GitHub uses `Home.md` as the landing page and
`_Sidebar.md` as wiki navigation.

For the main repository, also enable Issues, add a short repository
description and topics, protect the default branch, and create a DOI-backed
release when the manuscript is published.
