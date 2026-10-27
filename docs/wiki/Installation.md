# Installation

Clone the repository and enter its root directory:

```bash
git clone REPOSITORY_URL
cd REPOSITORY_NAME
```

Install Git LFS before retrieving the bundled PURF VCF:

```bash
git lfs install
git lfs pull
```

Create a Python environment and install the listed dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

PURF must be installed separately from its source repository. Individual
workflows may additionally require R, `bcftools`, GATK, alignment software or a
job scheduler. Consult the README within each workflow directory before
running it.
