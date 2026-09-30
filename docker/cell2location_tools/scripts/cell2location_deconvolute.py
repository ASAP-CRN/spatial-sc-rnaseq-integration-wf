#!/usr/bin/env python3
"""
Deconvolute one spatial transcriptomics sample using Cell2location (Visium)
against an scRNA-seq cohort reference.

Steps
-----
1. Build NB regression model on scRNA-seq reference to estimate cell type signatures.
2. Map cell types onto spatial spots using Cell2location.
3. (Optional) Cluster spots by cell type abundance using Leiden.

Outputs
-------
{sample_id}.integrated.h5ad
    Spatial AnnData with Cell2location results in .obsm and per-spot fraction in .obs.
{sample_id}.integration_metadata.csv
    Per-spot table: cell type fractions + obs columns for downstream reporting.

Usage (dry-run)
---------------
cell2location_deconvolute \\
    --sample-id SAMPLE001 \\
    --spatial-h5ad sample.qc.h5ad \\
    --sc-rnaseq-reference-adata final_all_genes.h5ad \\
    --sc-rnaseq-cell-types-parquet scanvi_cell_types.parquet \\
    --spatial-technology visium \\
    --output-dir results/

Usage (write outputs)
---------------------
cell2location_deconvolute ... --promote
"""

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.sparse import issparse

repo_root = next(p for p in Path(__file__).resolve().parents if (p / ".git").exists())
sys.path.insert(0, str(repo_root))

from utils.accelerator import get_accelerator
from utils.genomics_constants import MT_GENE_ENSEMBL_IDS
from utils.logging_extra import log_run_command


# ---------------------------------------------------------------------------
# Cell2location training defaults
# ---------------------------------------------------------------------------

# Source https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_short_demo.html
_CELLS_PER_SPOT = 8 # Used for mouse brain tissue. Note that cell2location_tutorial.html uses N_cells_per_location=30. But that's for human lymph node, which is denser than brain.

# Source https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html
_DETECTION_ALPHA = 100      # detection rate prior. Tutorial recommends detection_alpha between 20 (High technical variability in RNA detection) and 200 (low technical variability in RNA detection)
_REF_N_SAMPLES   = 1_000    # posterior samples for reference model. Tutorial recommends num_samples=1000
_EPOCHS_REF     = 250       # NB regression. Tutorial recommends max_epochs=250 for reference signature estimation.
_EPOCHS_SPATIAL = 30_000    # Cell2location spatial mapping. Tutorial recommends max_epochs=30000.
_REF_BATCH_SIZE  = 2_500    # posterior sampling. Tutorial recommends batch_size=2500, which is higher than scVI.
_TRAIN_SIZE      = 1        # training batch size for spatial model. Tutorial recommends train_size=1.

_ACCELERATOR     = get_accelerator()  # auto-selects: "gpu" (CUDA) > "mps" (Apple Silicon) > "cpu"


def _load_reference(ref_h5ad: str, cell_type_key: str) -> ad.AnnData:
    """
    Load scRNA-seq reference and validate that cell type labels are present in .obs.

    Parameters
    ----------
    ref_h5ad : str
        Path to all-gene reference h5ad (QC-filtered cells × all genes, raw counts).
        Must contain `cell_type_key` in .obs.
    cell_type_key : str
        Column in .obs to use as the cell type label.

    Returns
    -------
    anndata.AnnData
        Reference AnnData with a `cell_type` column in .obs.
    """
    print(f"[INFO] Loading reference: {ref_h5ad}")
    adata = ad.read_h5ad(ref_h5ad)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    if cell_type_key not in adata.obs.columns:
        available = adata.obs.columns.tolist()
        raise ValueError(
            f"Column '{cell_type_key}' not found in reference .obs. "
            f"Available columns: {available}"
        )

    adata.obs["cell_type"] = adata.obs[cell_type_key]
    n_missing = adata.obs["cell_type"].isna().sum()
    if n_missing > 0:
        print(f"  [WARN] Dropping {n_missing:,} cells with NaN cell type label.")
        adata = adata[adata.obs["cell_type"].notna()].copy()

    n_types = adata.obs["cell_type"].nunique()
    print(f"  {adata.n_obs:,} cells retained; {n_types} cell types")
    return adata


def _correct_sample_id(sample_id: str) -> str:
    """Return `sample_id` with the known ASAP typo PMBDS → PMDBS corrected."""
    if sample_id.startswith("ASAP_PMBDS_"):
        return "ASAP_PMDBS_" + sample_id[len("ASAP_PMBDS_"):]
    return sample_id


def _patch_sample_id(adata: ad.AnnData, original: str, corrected: str) -> None:
    """
    Replace `original` with `corrected` in all string .obs columns, in place.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData to patch.
    original : str
        The erroneous sample ID string to replace (e.g. "ASAP_PMBDS_000002_s007_Rep1").
    corrected : str
        The canonical sample ID string to substitute (e.g. "ASAP_PMDBS_000002_s007_Rep1").

    Returns
    -------
    None
    """
    for col in adata.obs.select_dtypes("object").columns:
        if adata.obs[col].str.contains(original, regex=False).any():
            adata.obs[col] = adata.obs[col].str.replace(original, corrected, regex=False)
            print(f"  [PATCH] obs['{col}']: '{original}' → '{corrected}'")


def _filter_mt_genes(adata: ad.AnnData, label: str) -> ad.AnnData:
    """Remove mitochondrial genes from `adata` (symbol- or Ensembl-indexed)."""
    mt_mask = (
        adata.var_names.str.upper().str.startswith("MT-")
        | adata.var_names.isin(MT_GENE_ENSEMBL_IDS)
    )
    n_mt = mt_mask.sum()
    if n_mt:
        print(f"  [{label}] Removing {n_mt} mitochondrial genes.")
        adata = adata[:, ~mt_mask].copy()
    return adata


def _load_spatial(
    spatial_h5ad: str,
    sample_id: str,
) -> tuple:
    """
    Load a Visium spatial h5ad and prepare it for Cell2location.

    If `var["gene_ids"]` is present, var_names are replaced with Ensembl IDs
    so they intersect with the reference (which uses Ensembl IDs as var_names).
    Auto-corrects the known ASAP typo PMBDS → PMDBS in .obs string columns.
    Removes mitochondrial genes and validates raw integer counts.

    Parameters
    ----------
    spatial_h5ad : str
        Path to per-sample Visium h5ad (raw UMI counts, post-QC).
    sample_id : str
        Sample identifier; if it starts with "ASAP_PMBDS_" it is corrected to
        "ASAP_PMDBS_" and obs columns are patched accordingly.

    Returns
    -------
    tuple
        (anndata.AnnData, str) — prepared spatial AnnData and the
        (possibly corrected) canonical sample_id.
    """
    corrected_id = _correct_sample_id(sample_id)

    print(f"[INFO] Loading spatial adata: {spatial_h5ad}")
    adata = ad.read_h5ad(spatial_h5ad)
    print(f"  {adata.n_obs:,} spots × {adata.n_vars:,} genes")

    if "gene_ids" in adata.var.columns:
        adata.var_names = adata.var["gene_ids"].astype(str)
        adata.var_names.name = "gene_id"
        print(f"  var_names set to Ensembl IDs from var['gene_ids']")

    if corrected_id != sample_id:
        print(f"  [PATCH] Correcting sample ID: '{sample_id}' → '{corrected_id}'")
        _patch_sample_id(adata, sample_id, corrected_id)

    adata = _filter_mt_genes(adata, "spatial")
    _ensure_raw_int_counts(adata, "spatial")
    return adata, corrected_id


def _ensure_raw_int_counts(adata: ad.AnnData, label: str) -> None:
    """Raise if X does not contain raw integer counts."""
    x = adata.X
    mat = x.data if issparse(x) else x.ravel()
    if not np.issubdtype(mat.dtype, np.integer):
        non_int = mat[mat % 1 != 0]
        if non_int.size > 0:
            raise ValueError(
                f"[{label}] X contains non-integer values (e.g. {non_int[:5].tolist()}). "
                "Cell2location requires raw integer counts."
            )


def _build_reference_signatures(
    adata_ref: ad.AnnData,
    batch_key: str,
    max_epochs: int,
    accelerator: str,
) -> pd.DataFrame:
    """
    Train NB regression model to estimate per-cell-type gene expression signatures.

    Parameters
    ----------
    adata_ref : anndata.AnnData
        Reference with raw counts and `cell_type` in .obs.
    batch_key : str
        Obs column used as batch covariate in the regression model.
    max_epochs : int
        Training epochs for the NB regression model.
    accelerator : str
        PyTorch Lightning accelerator ("mps", "gpu", or "cpu").

    Returns
    -------
    pandas.DataFrame
        Estimated signatures (genes × cell types), the `inf_aver` matrix
        passed to Cell2location.
    """
    import cell2location
    from cell2location.models import RegressionModel

    print(f"[INFO] Setting up NB regression model ({adata_ref.n_obs:,} cells)...")
    RegressionModel.setup_anndata(
        adata_ref,
        batch_key=batch_key,
        labels_key="cell_type",
    )
    ref_model = RegressionModel(adata_ref)
    ref_model.train(max_epochs=max_epochs, accelerator=accelerator)

    print("[INFO] Exporting reference posterior...")
    adata_ref = ref_model.export_posterior(
        adata_ref,
        sample_kwargs={"num_samples": _REF_N_SAMPLES, "batch_size": _REF_BATCH_SIZE},
    )

    # Extract inf_aver: genes × cell types
    factor_names = adata_ref.uns["mod"]["factor_names"]
    inf_aver = pd.DataFrame(
        adata_ref.varm["means_per_cluster_mu_fg"],
        index=adata_ref.var_names,
        columns=factor_names,
    )
    print(f"  Signatures: {inf_aver.shape[0]:,} genes × {inf_aver.shape[1]} cell types")
    return inf_aver


def _run_cell2location_visium(
    adata_spatial: ad.AnnData,
    inf_aver: pd.DataFrame,
    n_cells_per_spot: int,
    detection_alpha: int,
    max_epochs: int,
    accelerator: str,
) -> ad.AnnData:
    """
    Map cell types onto Visium spots using Cell2location.

    Parameters
    ----------
    adata_spatial : anndata.AnnData
        Per-sample Visium adata (spots × genes, raw UMI counts).
    inf_aver : pandas.DataFrame
        Reference signatures from `_build_reference_signatures` (genes × cell types).
    n_cells_per_spot : int
        Prior on cells per spot (N_cells_per_location).
    detection_alpha : int
        Prior on detection rate consistency across spots.
    max_epochs : int
        Training epochs for the Cell2location model.
    accelerator : str
        PyTorch Lightning accelerator ("mps", "gpu", or "cpu").

    Returns
    -------
    anndata.AnnData
        Spatial adata with posterior estimates in .obsm.
    """
    from cell2location.models import Cell2location

    # Filter to genes present in both spatial and reference signatures
    shared_genes = adata_spatial.var_names.intersection(inf_aver.index)
    if len(shared_genes) == 0:
        raise ValueError("No shared genes between spatial adata and reference signatures.")
    print(f"  Shared genes: {len(shared_genes):,} / {adata_spatial.n_vars:,} spatial, "
          f"{len(inf_aver):,} reference")

    adata_spatial = adata_spatial[:, shared_genes].copy()
    inf_aver_filtered = inf_aver.loc[shared_genes]

    print(f"[INFO] Setting up Cell2location model ({adata_spatial.n_obs:,} spots)...")
    Cell2location.setup_anndata(adata_spatial, batch_key=None)
    model = Cell2location(
        adata_spatial,
        cell_state_df=inf_aver_filtered,
        N_cells_per_location=n_cells_per_spot,
        detection_alpha=detection_alpha,
    )
    model.train(
        max_epochs=max_epochs,
        batch_size=None,
        train_size=_TRAIN_SIZE,
        accelerator=accelerator,
    )

    print("[INFO] Exporting spatial posterior...")
    adata_spatial = model.export_posterior(
        adata_spatial,
        sample_kwargs={"num_samples": _REF_N_SAMPLES, "batch_size": adata_spatial.n_obs},
    )
    return adata_spatial


def _add_cell_type_fractions(adata_spatial: ad.AnnData) -> None:
    """
    Normalise q05 abundance estimates to per-spot cell type fractions and add to .obs.

    Modifies `adata_spatial` in place.

    Parameters
    ----------
    adata_spatial : anndata.AnnData
        Spatial adata after Cell2location posterior export.

    Returns
    -------
    None
    """
    abundances = adata_spatial.obsm["q05_cell_abundance_w_sf"]
    fractions  = abundances.div(abundances.sum(axis=1), axis=0)
    fractions.columns = [f"frac_{c}" for c in fractions.columns]
    adata_spatial.obsm["cell_type_fractions"] = fractions
    # Flatten fractions into .obs for easy CSV export
    for col in fractions.columns:
        adata_spatial.obs[col] = fractions[col].values


def _downstream_analysis(
    adata_spatial: ad.AnnData,
    n_comps: int,
    n_neighbors: int,
    leiden_res: list,
) -> None:
    """
    PCA → neighbours → UMAP → Leiden clustering on cell type abundance vectors.

    Uses `.obsm["cell_type_fractions"]` as the feature matrix so that spatial
    domains are defined by cell type composition rather than raw gene expression.
    Modifies `adata_spatial` in place.

    Parameters
    ----------
    adata_spatial : anndata.AnnData
        Spatial adata with `cell_type_fractions` in .obsm.
    n_comps : int
        PCA components.
    n_neighbors : int
        Neighbourhood size for graph construction.
    leiden_res : list
        Leiden resolutions; one clustering per resolution is added to .obs.

    Returns
    -------
    None
    """
    print("[INFO] Running downstream analysis on cell type fractions...")
    fractions = adata_spatial.obsm["cell_type_fractions"]
    # Embed fractions as a temporary X for scanpy's PCA / neighbours / UMAP
    adata_tmp = ad.AnnData(X=fractions.values, obs=adata_spatial.obs)
    sc.tl.pca(adata_tmp, n_comps=min(n_comps, fractions.shape[1] - 1))
    sc.pp.neighbors(adata_tmp, n_neighbors=n_neighbors, n_pcs=adata_tmp.obsm["X_pca"].shape[1])
    sc.tl.umap(adata_tmp)
    adata_spatial.obsm["X_pca_fractions"]  = adata_tmp.obsm["X_pca"]
    adata_spatial.obsm["X_umap_fractions"] = adata_tmp.obsm["X_umap"]
    adata_spatial.obsp["connectivities"]   = adata_tmp.obsp["connectivities"]
    adata_spatial.obsp["distances"]        = adata_tmp.obsp["distances"]
    adata_spatial.uns["neighbors"]         = adata_tmp.uns["neighbors"]
    for res in leiden_res:
        key = f"leiden_res_{str(res).replace('.', '_')}"
        sc.tl.leiden(adata_tmp, resolution=res, key_added=key)
        adata_spatial.obs[key] = adata_tmp.obs[key].values
        print(f"  Leiden {res}: {adata_tmp.obs[key].nunique()} clusters")


def integrate_visium(
    sample_id: str,
    spatial_h5ad: str,
    ref_h5ad: str | None,
    reference_model: str | None,
    batch_key: str,
    cell_type_key: str,
    n_comps: int,
    n_neighbors: int,
    leiden_res: list,
    n_cells_per_spot: int,
    detection_alpha: int,
    epochs_ref: int,
    epochs_spatial: int,
    accelerator: str,
    output_dir: str,
    promote: bool,
) -> None:
    """
    Full Cell2location pipeline for one Visium sample.

    Parameters
    ----------
    sample_id : str
        Sample identifier; "ASAP_PMBDS_" is auto-corrected to "ASAP_PMDBS_".
    spatial_h5ad : str
        Path to per-sample Visium adata (raw UMI counts, post-QC).
    ref_h5ad : str or None
        Path to scRNA-seq reference adata (QC-filtered cells × all genes, raw counts).
        Must contain `cell_type_key` in .obs. Mutually exclusive with `reference_model`.
    reference_model : str or None
        Path to a pre-computed `inf_aver.csv.gz` from a previous run. When provided,
        Stage 1 NB regression is skipped. Mutually exclusive with `ref_h5ad`.
    batch_key : str
        Obs column used as batch covariate in the NB regression model.
    cell_type_key : str
        Parquet column to use as the cell type label.
    n_comps : int
        PCA components for downstream clustering.
    n_neighbors : int
        Neighbourhood size for downstream clustering.
    leiden_res : list
        Leiden resolutions for downstream clustering.
    n_cells_per_spot : int
        N_cells_per_location prior for Cell2location.
    detection_alpha : int
        Detection rate prior for Cell2location.
    epochs_ref : int
        NB regression training epochs.
    epochs_spatial : int
        Cell2location training epochs.
    accelerator : str
        PyTorch Lightning accelerator ("mps", "gpu", or "cpu").
    output_dir : str
        Directory to write output files.
    promote : bool
        If False, print plan without writing any files.

    Returns
    -------
    None
    """
    sample_id = _correct_sample_id(sample_id)
    out_h5ad      = os.path.join(output_dir, f"{sample_id}.deconvoluted.h5ad")
    out_csv       = os.path.join(output_dir, f"{sample_id}.deconvolution_metadata.csv")
    out_inf_aver  = os.path.join(output_dir, "inf_aver.csv.gz")

    if not promote:
        print(f"[DRY-RUN] Would write:")
        print(f"  {out_h5ad}")
        print(f"  {out_csv}")
        if ref_h5ad:
            print(f"  {out_inf_aver}")
        print("Re-run with --promote to write outputs.")
        return

    os.makedirs(output_dir, exist_ok=True)
    out_timing = os.path.join(output_dir, f"{sample_id}.timing.tsv")
    timing = []
    t_total = time.time()

    # --- Step 1: reference signatures ---
    t0 = time.time()
    if reference_model:
        print(f"[INFO] Loading pre-computed signatures: {reference_model}")
        inf_aver = pd.read_csv(reference_model, index_col=0, compression="gzip")
        print(f"  Signatures: {inf_aver.shape[0]:,} genes × {inf_aver.shape[1]} cell types")
        timing.append(("stage1_nb_regression_skipped", time.time() - t0))
    else:
        adata_ref = _load_reference(ref_h5ad, cell_type_key)
        adata_ref = _filter_mt_genes(adata_ref, "reference")
        _ensure_raw_int_counts(adata_ref, "reference")
        inf_aver = _build_reference_signatures(adata_ref, batch_key, epochs_ref, accelerator)
        del adata_ref
        print(f"[INFO] Writing {out_inf_aver} ...")
        inf_aver.to_csv(out_inf_aver, compression="gzip")
        timing.append(("stage1_nb_regression", time.time() - t0))

    # --- Step 2: spatial mapping ---
    t0 = time.time()
    adata_spatial, sample_id = _load_spatial(spatial_h5ad, sample_id)
    adata_spatial = _run_cell2location_visium(
        adata_spatial, inf_aver, n_cells_per_spot, detection_alpha, epochs_spatial, accelerator
    )
    timing.append(("stage2_spatial_mapping", time.time() - t0))

    # --- Step 3: downstream analysis ---
    t0 = time.time()
    _add_cell_type_fractions(adata_spatial)
    _downstream_analysis(adata_spatial, n_comps, n_neighbors, leiden_res)
    timing.append(("stage3_downstream", time.time() - t0))

    # --- Write outputs ---
    t0 = time.time()
    print(f"[INFO] Writing {out_h5ad} ...")
    adata_spatial.write_h5ad(out_h5ad, compression="gzip")
    meta_cols = (
        [c for c in adata_spatial.obs.columns if c.startswith("frac_")]
        + [c for c in adata_spatial.obs.columns if c.startswith("leiden_")]
        + [c for c in adata_spatial.obs.columns if not c.startswith("frac_") and not c.startswith("leiden_")]
    )
    adata_spatial.obs[meta_cols].to_csv(out_csv)
    print(f"[INFO] Writing {out_csv} ...")
    timing.append(("write_outputs", time.time() - t0))
    timing.append(("total", time.time() - t_total))

    with open(out_timing, "w") as f:
        f.write("step\tduration_seconds\n")
        for step, secs in timing:
            f.write(f"{step}\t{secs:.1f}\n")
    print(f"[INFO] Timing written to {out_timing}")
    print(f"[DONE] {sample_id}")


def main() -> None:
    """Parse arguments and dispatch to the appropriate integration function."""
    args = _parse_args()

    ref_adata = args.sc_rnaseq_reference_adata
    ref_model = args.sc_rnaseq_reference_cell2location_model
    if ref_adata and ref_model:
        raise SystemExit(
            "error: --sc-rnaseq-reference-adata and --sc-rnaseq-reference-cell2location-model "
            "are mutually exclusive."
        )
    if not ref_adata and not ref_model:
        raise SystemExit(
            "error: one of --sc-rnaseq-reference-adata or "
            "--sc-rnaseq-reference-cell2location-model is required."
        )

    log_run_command(Path(args.output_dir) / f"{Path(sys.argv[0]).stem}.log")

    integrate_visium(
        sample_id          = args.sample_id,
        spatial_h5ad       = args.spatial_h5ad,
        ref_h5ad           = ref_adata,
        reference_model    = ref_model,
        batch_key          = args.batch_key,
        cell_type_key      = args.cell_type_key,
        n_comps            = args.n_comps,
        n_neighbors        = args.n_neighbors,
        leiden_res         = args.leiden_res,
        n_cells_per_spot   = args.n_cells_per_spot,
        detection_alpha    = args.detection_alpha,
        epochs_ref         = args.max_epochs_ref,
        epochs_spatial     = args.max_epochs_spatial,
        accelerator        = args.accelerator,
        output_dir         = args.output_dir,
        promote            = args.promote,
    )


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--sample-id",                    required=True,
        help="Spatial sample identifier. 'ASAP_PMBDS_' is auto-corrected to 'ASAP_PMDBS_'.")
    parser.add_argument("--spatial-h5ad",                required=True,
        help="Per-sample spatial h5ad (raw UMI counts, post-QC).")
    parser.add_argument("--sc-rnaseq-reference-adata",    default=None,
        help="scRNA-seq reference h5ad (QC-filtered cells × all genes, raw counts).\n"
             "Trains Stage 1 NB regression and saves inf_aver.csv.gz to --output-dir.\n"
             "Mutually exclusive with --sc-rnaseq-reference-cell2location-model.")
    parser.add_argument("--sc-rnaseq-reference-cell2location-model", default=None,
        help="Path to a pre-computed inf_aver.csv.gz produced by a previous run.\n"
             "Skips Stage 1 NB regression entirely.\n"
             "Mutually exclusive with --sc-rnaseq-reference-adata.")
    parser.add_argument("--spatial-technology",           default="visium",
        choices=["visium"],
        help="Spatial technology. Currently only visium is supported.")
    parser.add_argument("--batch-key",                    default="batch_id",
        help="Obs column used as batch covariate in the NB regression model.\n"
             "Default: batch_id")
    parser.add_argument("--cell-type-key",                default="cell_type_pred",
        help="Column in the reference .obs to use as the cell type label.\n"
             "Default: cell_type_pred")
    parser.add_argument("--n-comps",                      type=int, default=30,
        help="PCA components for downstream spot clustering. Default: 30")
    parser.add_argument("--n-neighbors",                  type=int, default=15,
        help="Neighbourhood size for downstream spot clustering. Default: 15")
    parser.add_argument("--leiden-res",                   type=float, nargs="+",
        default=[0.05, 0.1, 0.2, 0.4],
        help="Leiden resolutions for downstream spot clustering. Default: 0.05 0.1 0.2 0.4")
    parser.add_argument("--n-cells-per-spot",             type=int, default=_CELLS_PER_SPOT,
        help=f"Cell2location N_cells_per_location prior. Default: {_CELLS_PER_SPOT}")
    parser.add_argument("--detection-alpha",              type=int, default=_DETECTION_ALPHA,
        help=f"Cell2location detection rate prior (tutorial range 20–200).\n"
             f"Higher values assume low technical variability in RNA detection.\n"
             f"Default: {_DETECTION_ALPHA}")
    parser.add_argument("--max-epochs-ref",               type=int, default=_EPOCHS_REF,
        help=f"NB regression training epochs. Default: {_EPOCHS_REF}")
    parser.add_argument("--max-epochs-spatial",           type=int, default=_EPOCHS_SPATIAL,
        help=f"Cell2location training epochs. Default: {_EPOCHS_SPATIAL}")
    parser.add_argument("--accelerator",                  default=_ACCELERATOR,
        help=f"PyTorch Lightning accelerator: 'mps' (Apple Silicon), 'gpu' (CUDA), 'cpu'.\n"
             f"Default: {_ACCELERATOR}")
    parser.add_argument("--output-dir",                   default=".",
        help="Directory to write output files. Default: current directory")
    parser.add_argument("--promote",                      action="store_true", default=False,
        help="Write output files. Default is dry-run (prints what would happen).")
    return parser.parse_args()


if __name__ == "__main__":
    main()
