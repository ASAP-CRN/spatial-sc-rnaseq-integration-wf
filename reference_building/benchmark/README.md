# reference_building/benchmark/

Benchmarks this repo's Cell2Location deconvolution pipeline against Ma et al.
2025 (*Nat Commun*) — 10 human substantia nigra Visium sections with
published RCTD/BayesSpace cell-type proportions.

**Caveat:** Ma et al.'s RCTD proportions are aggregated per-BayesSpace-cluster means, not true per-spot values — every spot sharing a cluster gets an identical proportion vector. See `03_build_benchmark.py`'s docstring.

## Workflow

Fetches the paper's public data, reformats it alongside this repo's own
Visium samples, runs Cell2Location, then compares the two sets of per-spot
cell-type calls.

This is a one-time benchmarking analysis, not part of the spatial
deconvolution WDL pipeline.

![Workflow diagram](docs/benchmark.svg "Workflow diagram")

See slides for full analysis context: [Cell2Location vs Ma (RCTD) overview](https://docs.google.com/presentation/d/1NrLuYe5S31BgBTzDqgeqwHes8uX3g8db5uUCJzg6Jec/edit).

## Directory contents, inputs, and outputs

| Path | What | Inputs | Key outputs |
|---|---|---|---|
| `scripts/01_fetch_al_dalahmah_zenodo_data.py` | Step 1. Fetch the two Zenodo artifacts needed for the benchmark | Zenodo records 19454950 / 19457732 (network) | `data/original/{processed_st_obj.rds, 41467_2025_62478_MOESM3_ESM.xlsx}` |
| `scripts/02_export_spot_metadata.R` | Step 2. Export per-spot Region/Condition/BayesSpace-cluster metadata from the Seurat object list | `data/original/processed_st_obj.rds` (step 1) | `data/reformatted/spot_metadata_all_samples.csv.gz` |
| `scripts/03_build_benchmark.py` | Step 3. Join spot metadata to Supplementary Data 1's per-cluster RCTD proportions | Spot metadata (step 2) + Supplementary Data 1 xlsx (step 1) | `data/reformatted/per_spot_benchmark_al_dalahmah_sn.csv.gz` |
| `scripts/04_plot_sn_dopa_gaba_meg3_pd_vs_control.py` | Step 4. (Optional) Standalone exploratory PD-vs-control stats/plots on the RCTD benchmark (not required for steps 5-7) | Benchmark table (step 3) | `data/analysis/*.html` |
| `scripts/05_format_visium_for_cell2location.py` | Step 5. Convert each sample's Space Ranger `outs/` zip into a Cell2Location-ready h5ad | `data/original/*_result-*.zip` | `data/reformatted/spatial_h5ad/<sample>.spatial.h5ad` |
| `scripts/06_submit_cell2location.sh` | Step 6. Run this repo's `cell2location_deconv.py` — reference training once, then spatial mapping + clustering per sample | Spatial h5ads (step 5) + a scRNA-seq reference h5ad | `data/reformatted/cell2location/<sample>.deconvoluted.h5ad` |
| `scripts/07_compare_uber_kamath_pred_vs_ma_labels.py` | Step 7. Pool both sides into a shared cell-type vocabulary and compare as continuous proportions (PCC/RMSE/Jensen-Shannon) | Cell2Location outputs (step 6) + benchmark table (step 3) | `data/reformatted/collapsed_agreement_*.csv[.gz]` |
| `<repo_root>/notebooks/visualize_deconvolution.ipynb` | Interactive Jupyter notebook to explore benchmark analysis | Outputs from steps 06 and 07 | Confusion matrix heatmap, spatial-plots |
| `docs/` | Workflow diagram and other supporting documentation assets. | — | — |