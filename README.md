# RNA Fusion Pipeline

This repository contains three related but independently run components:

1. A Snakemake workflow that integrates paired-end FASTQ quality control, STAR-Fusion/FusionInspector, fusion-table combination, filtering, and other project-specific analysis steps.
2. R analysis scripts for combining calls, filtering, HGNC mapping, and protein-domain data preparation.
3. PubMed/PMC literature-mining scripts, including an optional NVIDIA Nemotron curation step.

The [CTAT Genome Resource Library (CTAT Lib)](https://toolshed.g2.bx.psu.edu/view/trinity_ctat) is a public resource developed and maintained by the Trinity CTAT team at the Broad Institute, led primarily by [Brian J. Haas](https://www.google.com/search?kgmid=%2Fg%2F11ftc28wtm&q=ctat+lib+by+whome+brian+j.+haas). This repository does not develop or maintain CTAT Lib; it configures STAR-Fusion and FusionInspector to use an externally supplied CTAT genome library. Our pipeline integrates these public tools with project-specific FusionInspector parsing, fusion-table combination, filtering, HGNC mapping, protein-domain preparation, and literature mining.

The repository does not include sequencing reads, a CTAT genome library, the external case/control tables used by the R analysis, or a validated software lockfile. The example configuration is a template, not a runnable dataset. Do not interpret a successful syntax check as validation of scientific reproducibility.

## Repository map

| Path | Purpose |
|---|---|
| `workflow/fusion.smk` | FastQC, MultiQC, and STAR-Fusion workflow |
| `workflow/config.yaml` | Workflow paths and thread settings |
| `workflow/scripts/` | Fusion table combination and R analysis scripts |
| `reference_data/` | Curated AML fusion list and Pfam domain table |
| `literature_mining/fusin_text_mining.R` | PubMed/PMC retrieval and text mining |
| `llm_mining/ingest.py` | XML ingestion and model-assisted curation |
| `llm_mining/run_ingest.sh` | Environment-configured runner for `ingest.py` |
| `examples/samples.tsv` | Sample-sheet format example; replace its paths before running |
| `tests/` | Lightweight repository and script tests |
| `environment/` | Dependency declarations; these are not exact environment locks |

## Installation

### Workflow tools

Create a Conda environment with the core workflow executables:

```bash
mamba env create -f environment/workflow.yml
mamba activate rnafusion-workflow
```

The supplied analysis environment is Ubuntu 24.04, R 4.6.2 (2026-10-29), and RStudio 2025.9.2.418; these values are recorded in `environment/reported-versions.yml`. The environment specifies Snakemake, STAR-Fusion, FastQC, and MultiQC, but does not pin their versions. R package and Python package versions were not supplied. For publication, record and validate the actual versions used; do not treat these environment declarations as a lockfile.

Capture the versions and checksums available in a configured environment with:

```bash
bash environment/capture_versions.sh
```

Save that report alongside the run configuration. It records the current environment; it does not recover the environment used for earlier analyses.

### Python literature curation

```bash
python3 -m pip install -r environment/python-requirements.txt
```

The R scripts use `dplyr`, `tidyr`, `readr`, `pubmed.mineR`, `rentrez`, `xml2`, `biomaRt`, `chimeraviz`, and `EnsDb.Hsapiens.v86`. The package names are listed in `environment/r-packages.txt`. A starting installation command is:

```r
install.packages(c("dplyr", "tidyr", "readr", "pubmed.mineR", "rentrez", "xml2"))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("biomaRt", "chimeraviz", "EnsDb.Hsapiens.v86"))
```

These are unpinned installation instructions, not a validated lock. Record and test the resolved package versions before publication.

## FASTQ workflow

1. Copy `examples/samples.tsv` and replace both FASTQ paths with existing paired-end reads. The tab-separated columns are `sampleID`, `pathR1`, and `pathR2`.
2. Edit `workflow/config.yaml`; set `samplesheet` and `ref_ctat` to the sample sheet and the STAR-Fusion CTAT genome-library directory. Paths are interpreted from the repository root.
3. From the repository root, run:

```bash
snakemake --snakefile workflow/fusion.smk \
  --configfile workflow/config.yaml \
  --cores 8 \
  --resources star_fusion_slot=2
```

The workflow expects FastQC, MultiQC, STAR-Fusion, and the configured CTAT library to be available. It does not perform read trimming. The STAR-Fusion rule currently uses:

| STAR-Fusion option | Current value |
|---|---:|
| `--min_junction_reads` | 2 |
| `--min_sum_frags` | 5 |
| `--min_FFPM` | 0.1 |
| `--min_novel_junction_support` | 3 |

These are the existing workflow settings, not recommendations or values inferred from the downstream R script.

Expected outputs include per-sample STAR-Fusion predictions under `results/fusion_result/`, per-sample FastQC results under `results/post_trim_fastqc/`, `results/multiqc/multiqc_report.html`, and logs under `results/logs/`.

## Fusion table combination and downstream R analysis

The combiner accepts one or more study directories and finds matching FusionInspector abridged TSV files. It appends the source path and sample name:

```bash
bash workflow/scripts/combine_fusions.sh \
  -i /path/to/study_directory \
  -o results/all_samples_abridged.tsv \
  -p '*.FusionInspector.fusions.abridged.tsv'
```

The R analysis in `workflow/scripts/filter_data.R` is a separate downstream analysis. It expects case/control combined fusion TSVs and FusionGDB lists that are not included in this repository. Its first `setwd` path is host-specific and must be configured before running. It writes FT0–FT4 and known-fusion tables. Its existing criteria include:

- FT1: junction reads at least 5, spanning fragments at least 2, FFPM at least 0.1, and `LargeAnchorSupport == "YES"`, with additional annotation/gene exclusions.
- FT2: FFPM at least 0.1, maximum left/right FAR at least 0.2, the existing microhomology condition, and `PROT_FUSION_TYPE == "INFRAME"`.
- Later stages exclude fusions using the external FusionGDB files.

These R criteria are distinct from the STAR-Fusion settings above. They are documented as currently written; no filtering criteria have been harmonized or changed.

`hgnc_map_symbols.R` retrieves the current HGNC table. `prepare_protein_bed.R` retrieves Ensembl data and writes Pfam annotations; the repository also contains a separate pre-generated Pfam table. These external sources can change over time, so record their release/version and provenance for each analysis. The protein-preparation script also includes a plotting example that refers to an externally defined `fusion1` object; it is not a standalone plotting workflow.

## Literature mining

`literature_mining/fusin_text_mining.R` queries PubMed and PMC, downloads XML, and writes a fusion-pair frequency table. It needs network access and local HGNC mapping data; edit its absolute paths before use. It currently computes an AML/RNA/fusion flag but does not use that flag in the final selection, which is based on extracted gene pairs.

For model-assisted XML curation, supply `NVIDIA_API_KEY` through a secret manager or
the process environment. Set the other required variables and run
`llm_mining/run_ingest.sh`:

```bash
export INPUT_DIR=/path/to/xml_directory
export GENE_MAP=/path/to/HGNC_symbol_mapping.tsv
export OUTPUT=/path/to/output.tsv
# Optional:
export NVIDIA_MODEL='nvidia/nemotron-3.5-lightning-30b-a3b'
export NVIDIA_BASE_URL='https://integrate.api.nvidia.com/v1'
export NVIDIA_THRESHOLD=0.5
bash llm_mining/run_ingest.sh
```

Do not put API keys in scripts, configuration files, or sample data. The ingestion program reads non-recursive `*.xml` files from `INPUT_DIR`; it writes a TSV and requires a gene-map TSV with a recognized symbol column.

The output TSV contains article identifiers, gene-pair fields, model confidence/confirmation, HGNC status, evidence, and source path. `NVIDIA_THRESHOLD` (default `0.5`) sets the confidence threshold used for the `nvidia_confirmed` output flag; it does not remove other parsed pairs from the TSV.

## Known limitations requiring scientific review

- The FusionGDB file in `filter_data.R` is read with different header settings in different stages.
- The literature-mining final selection does not apply its calculated AML/RNA/fusion flag.
- The HGNC mapping script selects `alias_symbol` but does not currently use it in the emitted mapping.

These are documented rather than changed here because corrections could change fusion inclusion or gene normalization. Resolve them only with explicit scientific decisions and validation.

## Tests and CI

Run the repository tests with:

```bash
python3 -m unittest discover -s tests -v
bash -n llm_mining/run_ingest.sh workflow/scripts/combine_fusions.sh
```

The tests cover repository-level behavior only; they do not run STAR-Fusion, validate biological findings, or establish end-to-end reproducibility. The CI workflow likewise performs static checks and these lightweight tests.

## Citation, license, and release notes

If you use this repository or the associated pipeline, cite the manuscript as follows:

```text
Citing RNAfusion_pipeline: Kumar M, Seth S, Sasi A, Ganguly S, Bakhshi S. Uncovering Novel RNA Fusion Events and the Prevalence of RNA Fusions in Acute Myeloid Leukemia Through Integrative Analysis of Public Transcriptomic Datasets. Under Communication. NA. DOI: NA.
```

```bibtex
@article{kumar_etal_rnafusion_NA,
  author = {Kumar, Mukesh and Seth, Surabhi and Sasi, Archana and Ganguly, Shuvadeep and Bakhshi, Sameer},
  title = {Uncovering Novel RNA Fusion Events and the Prevalence of RNA Fusions in Acute Myeloid Leukemia Through Integrative Analysis of Public Transcriptomic Datasets},
  year = {NA},
  doi = {NA}
}
```

The manuscript was described as “Under Communication”; no venue, publication year, or DOI was supplied. `CITATION.cff` records the software and manuscript title, version `1.0`, authors, and affiliation supplied by the authors. The author symbols `$` and `#` were not defined, so they are not represented in the structured citation. Confirm those details before release.

No license is declared for the code or bundled reference data, and this repository does not grant reuse permission. The code rights holder would choose a software license; terms for bundled reference datasets must be checked with their original providers. Check institutional and publisher requirements before distribution. No license file has been added.

A credential-like literal was found in prior Git history and removed from the current runner. It remains in history; treat it as compromised if valid, rotate/revoke it, and arrange coordinated history cleanup before public release. The credential value is intentionally not reproduced here.

Before release, record exact software versions, CTAT/genome build, sample and reference checksums, HGNC/Ensembl releases, model identifier and endpoint, and the commands/configuration used. A release must not imply that excluded input data or unavailable external dependencies were tested.
