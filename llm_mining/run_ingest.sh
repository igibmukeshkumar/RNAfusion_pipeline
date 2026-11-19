#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

: "${NVIDIA_API_KEY:?Set NVIDIA_API_KEY through your environment or secret manager}"
: "${INPUT_DIR:?Set INPUT_DIR to a directory containing XML files}"
: "${GENE_MAP:?Set GENE_MAP to the HGNC gene-map TSV}"
: "${OUTPUT:?Set OUTPUT to the output TSV path}"

NVIDIA_BASE_URL="${NVIDIA_BASE_URL:-https://integrate.api.nvidia.com/v1}"
NVIDIA_MODEL="${NVIDIA_MODEL:-nvidia/nemotron-3.5-lightning-30b-a3b}"
NVIDIA_THRESHOLD="${NVIDIA_THRESHOLD:-0.5}"

python3 "$SCRIPT_DIR/ingest.py" \
    -i "$INPUT_DIR" \
    --gene-map "$GENE_MAP" \
    --nvidia-threshold "$NVIDIA_THRESHOLD" \
    --nvidia-model "$NVIDIA_MODEL" \
    --nvidia-base-url "$NVIDIA_BASE_URL" \
    --skip-processed \
    -o "$OUTPUT"

echo "Output: $OUTPUT"
