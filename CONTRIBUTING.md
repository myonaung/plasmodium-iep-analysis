# Contributing

Contributions that improve reproducibility, portability or documentation are
welcome.

## Before opening a pull request

1. Create a focused branch and keep unrelated changes separate.
2. Do not commit credentials, participant identifiers, controlled-access data,
   scheduler logs or institutional filesystem paths.
3. Preserve documented defaults unless the scientific rationale for changing
   them is stated in the pull request.
4. Run syntax checks for every language touched.
5. Update the relevant README when inputs, outputs or command-line arguments
   change.

## Basic checks

```bash
bash -n path/to/script.sh
python3 -m py_compile path/to/script.py
Rscript -e 'parse(file="path/to/script.R")'
```

For analytical changes, include the command used, software versions, random
seed, a small non-sensitive test input and the expected result. Large or
controlled datasets should be referenced by accession or release rather than
attached to an issue.
