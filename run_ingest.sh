#!/bin/bash
set -e

# ============================================================
# Configuration
# ============================================================
echo "=============================================="
echo " NVIDIA Nemotron RNA Fusion Curation"
echo "=============================================="

cd /home/adminsb/disk2/meta-analysis/rnafusion/nvidia_nemotron/code

export NVIDIA_API_KEY="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"
export NVIDIA_BASE_URL="https://integrate.api.nvidia.com/v1"
export NVIDIA_BASE_URL="https://integrate.api.nvidia.com/v1"
export NVIDIA_MODEL="nvidia/nemotron-3.5-lightning-30b-a3b"

echo "Model: $NVIDIA_MODEL"
echo "Endpoint: $NVIDIA_BASE_URL"

# ============================================================
# Run Nemotron
# ============================================================
python ingest.py \
    -i /home/adminsb/disk2/meta-analysis/rnafusion/nvidia_nemotron/litData_9oct/test_mat/v1_68 \
    --gene-map /home/adminsb/disk2/meta-analysis/rnafusion/pubmed_miner/code_pubmedminer/HGNC_symbol_mapping.tsv \
    --nvidia-threshold 0.5 \
    --nvidia-model "$NVIDIA_MODEL" \
    --nvidia-base-url "$NVIDIA_BASE_URL" --skip-processed \
    -o /home/adminsb/disk2/meta-analysis/rnafusion/nvidia_nemotron/results/test_68ids_matrix_nemotron_only.txt

# ============================================================
# Done
# ============================================================
echo
echo "=============================================="
echo " Pipeline completed"
echo "=============================================="
echo "Output:"
echo "/home/adminsb/disk2/meta-analysis/rnafusion/nvidia_nemotron/results/fusion_results_nemotron_only.txt"
echo
echo "Checkpoint:"
echo "/home/adminsb/disk2/meta-analysis/rnafusion/nvidia_nemotron/results/fusion_results_nemotron_only.txt.partial"
