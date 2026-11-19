#!/usr/bin/env bash
set -uo pipefail

if [[ -r /etc/os-release ]]; then
    . /etc/os-release
    printf '%-18s%s %s\n' "Operating system" "${NAME:-unknown}" "${VERSION_ID:-unknown}"
fi

report() {
    local label="$1"
    shift
    printf '%-18s' "$label"
    if ! command -v "$1" >/dev/null 2>&1; then
        printf '%s\n' 'not installed'
        return
    fi
    "$@" 2>&1 || printf '%s\n' '[version command failed]'
}

report "Python" python3 --version
report "Snakemake" snakemake --version
report "STAR-Fusion" STAR-Fusion --version
report "FastQC" fastqc --version
report "MultiQC" multiqc --version
report "R" R --version

if command -v python3 >/dev/null 2>&1; then
    python3 - <<'PY'
from importlib.metadata import PackageNotFoundError, version

for package in ("lxml", "openai"):
    try:
        found = version(package)
    except PackageNotFoundError:
        found = "not installed"
    print(f"{package:18}{found}")
PY
fi

if command -v Rscript >/dev/null 2>&1; then
    Rscript - <<'RS'
cat(sprintf("%-18s%s\n", "Rscript", getRversion()))
rstudio_version <- tryCatch({
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    as.character(rstudioapi::versionInfo()$version)
  } else {
    "not installed"
  }
}, error = function(e) "not running in RStudio")
cat(sprintf("%-18s%s\n", "RStudio", rstudio_version))
packages <- c(
  "dplyr", "tidyr", "readr", "pubmed.mineR", "rentrez", "xml2",
  "biomaRt", "chimeraviz", "EnsDb.Hsapiens.v86"
)
for (package in packages) {
  found <- if (requireNamespace(package, quietly = TRUE)) {
    as.character(packageVersion(package))
  } else {
    "not installed"
  }
  cat(sprintf("%-18s%s\n", package, found))
}
RS
fi

printf '\nReference file checksums:\n'
for path in \
    reference_data/known_AML_ref_fusion_list.txt \
    reference_data/pfam_domains.bed.gz
do
    if [[ -f "$path" ]]; then
        sha256sum "$path"
    fi
done
