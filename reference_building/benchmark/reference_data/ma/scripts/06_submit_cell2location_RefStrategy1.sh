#!/usr/bin/env bash
# Submits cell2location_deconv.py for Ma et al. 2025 Visium SN samples (05_format_visium_for_cell2location.py output).
# Stage 1 trains the reference to produce inf_aver.csv.gz
# Stages 2+3 loop per sample, reusing the inf_aver.csv.gz from Stage 1.

set -euo pipefail

REPO_ROOT="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
cd "$REPO_ROOT"

PYTHON="python3"
DECONVOLUTE_SCRIPT="docker/cell2location_tools/scripts/cell2location_deconv.py"
SPATIAL_DIR="reference_building/benchmark/reference_data/ma/data/reformatted/spatial_h5ad"
DATA_DIR="reference_building/benchmark/reference_data/ma/data"
LOG_FILE="$DATA_DIR/06_submit_cell2location_RefStrategy1.log"

SC_RNASEQ_REFERENCE_ADATA="reference_building/reference_h5ad/strategy_1/deconvolution_reference_300k_c2lUda_Kamath200Rand.h5ad"
OUTPUT_DIR="reference_building/benchmark/reference_data/ma/data/reformatted/cell2location_Strategy1"

# Hyperparameters match this repo's prior cell2location_deconv.py runs.
BATCH_KEY="batch_id"
CELL_TYPE_KEY="cell_type_pred"
MAX_EPOCHS_REF=250
MAX_EPOCHS_SPATIAL=30000

mkdir -p "$OUTPUT_DIR" "$DATA_DIR"
printf '[%s] Invoked: %s %s\n\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$0" "$*" > "$LOG_FILE"

# Records each Stage 1 / Stages 2+3 command to LOG_FILE -- not their own
# output, which already goes to OUTPUT_DIR's per-stage/sample *.stdout.log.
log_command() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >> "$LOG_FILE"
}

REF_STEM="$(basename "$SC_RNASEQ_REFERENCE_ADATA" .h5ad)"
INF_AVER="$OUTPUT_DIR/${REF_STEM}_inf_aver.csv.gz"

# --- Stage 1: train reference signatures once -------------------------------
if [[ -f "$INF_AVER" ]]; then
    echo "=== Stage 1: $INF_AVER already exists, skipping reference training ==="
else
    echo "=== Stage 1: training reference signatures ==="
    STAGE1_CMD=("$PYTHON" "$DECONVOLUTE_SCRIPT" \
        --run-stages 1 \
        --sc-rnaseq-reference-adata "$SC_RNASEQ_REFERENCE_ADATA" \
        --batch-key "$BATCH_KEY" \
        --cell-type-key "$CELL_TYPE_KEY" \
        --max-epochs-ref "$MAX_EPOCHS_REF" \
        --output-dir "$OUTPUT_DIR" \
        --promote)
    log_command "${STAGE1_CMD[@]}"
    "${STAGE1_CMD[@]}" 2>&1 | tee "$OUTPUT_DIR/stage1_${REF_STEM}.stdout.log"
fi

# --- Stages 2+3: spatial mapping + clustering, per sample -------------------
INPUT_H5ADS=("$SPATIAL_DIR"/*.spatial.h5ad)
number_of_samples="${#INPUT_H5ADS[@]}"
i=0
    for spatial_h5ad in "${INPUT_H5ADS[@]}"; do
    sample_id="$(basename "$spatial_h5ad" .spatial.h5ad)"
    echo "=== Stages 2+3: $sample_id ($((++i))/$number_of_samples) ==="
    STAGE23_CMD=("$PYTHON" "$DECONVOLUTE_SCRIPT" \
        --run-stages 2,3 \
        --sample-id "$sample_id" \
        --spatial-h5ad "$spatial_h5ad" \
        --sc-rnaseq-reference-cell2location-model "$INF_AVER" \
        --max-epochs-spatial "$MAX_EPOCHS_SPATIAL" \
        --output-dir "$OUTPUT_DIR" \
        --promote)
    log_command "${STAGE23_CMD[@]}"
    "${STAGE23_CMD[@]}" 2>&1 | tee "$OUTPUT_DIR/${sample_id}.stdout.log"
done
