#!/usr/bin/env bash
# Submits cell2location_deconv.py for asap-raw-team-scherzer-pmdbs-spatial-visium-mtg Visium samples.
# Stage 1 is NOT run, instead, it uses Strategy 3 reference signatures to produce the inf_aver.csv.gz produced for the benchmark analysis
# Stages 2+3 loop per sample, reusing the inf_aver.csv.gz from Stage 1.

set -euo pipefail

REPO_ROOT="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
cd "$REPO_ROOT"

PYTHON="python3"
DECONVOLUTE_SCRIPT="docker/cell2location_tools/scripts/cell2location_deconv.py"
SPATIAL_DIR="spatial_samples/asap-raw-team-scherzer-pmdbs-spatial-visium-mtg/workflow_execution/spatial_visium/preprocess/counts_to_adata/1.0.0"
SC_RNASEQ_REFERENCE_SIGNATURES="reference_building/benchmark/reference_data/ma/data/reformatted/cell2location_Strategy3/deconvolution_reference_strategy3_inf_aver.csv.gz"
OUTPUT_DIR="spatial_samples/asap-raw-team-scherzer-pmdbs-spatial-visium-mtg/workflow_execution/spatial_visium/cell2location/Strategy3"
LOG_FILE="$OUTPUT_DIR/06_submit_cell2location_ASAP_RefStrategy3.log"

# Hyperparameters match this repo's prior cell2location_deconv.py runs.
MAX_EPOCHS_SPATIAL=30000

mkdir -p "$OUTPUT_DIR"
printf '[%s] Invoked: %s %s\n\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$0" "$*" > "$LOG_FILE"

# Records each Stages 2+3 command to LOG_FILE -- not their own
# output, which already goes to OUTPUT_DIR's per-stage/sample *.stdout.log.
log_command() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >> "$LOG_FILE"
}

# --- Stages 2+3: spatial mapping + clustering, per sample -------------------
INPUT_H5ADS=("$SPATIAL_DIR"/*.cleaned_unfiltered.h5ad)
number_of_samples="${#INPUT_H5ADS[@]}"
i=0
for spatial_h5ad in "${INPUT_H5ADS[@]}"; do
    sample_id="$(basename "$spatial_h5ad" .cleaned_unfiltered.h5ad)"
    echo "=== Stages 2+3: $sample_id ($((++i))/$number_of_samples) ==="
    STAGE23_CMD=("$PYTHON" "$DECONVOLUTE_SCRIPT" \
        --run-stages 2,3 \
        --sample-id "$sample_id" \
        --spatial-h5ad "$spatial_h5ad" \
        --sc-rnaseq-reference-cell2location-model "$SC_RNASEQ_REFERENCE_SIGNATURES" \
        --max-epochs-spatial "$MAX_EPOCHS_SPATIAL" \
        --output-dir "$OUTPUT_DIR" \
        --promote)
    log_command "${STAGE23_CMD[@]}"
    "${STAGE23_CMD[@]}" 2>&1 | tee "$OUTPUT_DIR/${sample_id}.stdout.log"
done
