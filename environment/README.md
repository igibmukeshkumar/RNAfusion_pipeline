# Environment files

- `workflow.yml` declares the Conda packages needed by the Snakemake workflow.
- `python-requirements.txt` declares the direct Python dependencies for literature ingestion.
- `r-packages.txt` lists R packages used by the R scripts.
- `reported-versions.yml` records the user-supplied Ubuntu, R, and RStudio versions.

The supplied OS/R/RStudio versions do not specify the R package, Snakemake,
STAR-Fusion, FastQC, MultiQC, Python package, or CTAT versions. The dependency
declarations are not lockfiles. The chosen package versions can affect analysis
outputs. Before a publication release, recover the missing versions if possible,
test the intended environment against representative inputs, and record the
resolved versions, CTAT/genome build, external database releases, and checksums.
Do not treat installing the latest available versions as evidence of reproducibility.
