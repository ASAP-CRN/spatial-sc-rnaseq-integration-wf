#!/usr/bin/env python3
"""
Deconvolute one spatial transcriptomics sample using Cell2location against an scRNA-seq reference.
Only Visium is supported in this version.

Steps
-----
1. Build NB regression model on scRNA-seq reference to estimate cell type signatures.
2. Map cell types onto spatial spots using Cell2location.
3. (Optional) Cluster spots by cell type abundance using Leiden.

Inputs
------
--sample-id                               Spatial sample identifier.
--spatial-h5ad                            Per-sample Visium adata path (raw UMI counts, post-QC from `spatial-transcriptomics-wf`).
--sc-rnaseq-reference-adata               scRNA-seq reference h5ad; trains Stage 1 NB regression.
--sc-rnaseq-reference-cell2location-model Pre-computed inf_aver.csv.gz; skips Stage 1.
--build-sc-rnaseq-reference-only          Train Stage 1 and write inf_aver.csv.gz, then exit.
                                          --sample-id and --spatial-h5ad are not required.
--skip-clustering                        Skip Stage 3 (Leiden clustering); write Stage 1/2 outputs only.

Outputs
-------
{sample_id}.deconvoluted.h5ad
    Spatial AnnData with Cell2location results in .obsm and per-spot fraction in .obs.
{sample_id}.deconvolution_metadata.csv.gz
    Per-spot table: cell type fractions + obs columns for downstream reporting.
{ref_stem}_inf_aver.csv.gz
    Reference signatures (genes × cell types); only written when --sc-rnaseq-reference-adata
    is provided. {ref_stem} is the stem of the reference h5ad filename.
stage1_ref_training_history.png / stage2_spatial_training_history.png
    ELBO training-loss curves for each stage, for convergence inspection.
    Stage 1 only written when --sc-rnaseq-reference-adata is provided.

Usage (dry-run)
---------------
cell2location_deconvolute \\
    --sample-id SAMPLE001 \\
    --spatial-h5ad sample.qc.h5ad \\
    --sc-rnaseq-reference-adata final_all_genes.h5ad \\
    --spatial-technology visium \\
    --output-dir results/

Usage (write outputs)
---------------------
cell2location_deconvolute ... --promote
"""

import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import cell2location
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from cell2location.models import Cell2location, RegressionModel
from scipy.sparse import issparse

# Container builds (docker/cell2location_tools/Dockerfile) copy utils/ next to this
# script and set PYTHONPATH accordingly, with no .git present. Only add repo_root to
# sys.path when running from a git clone (local dev); otherwise rely on PYTHONPATH.
repo_root = next((p for p in Path(__file__).resolve().parents if (p / ".git").exists()), None)
if repo_root is not None:
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
_REF_N_SAMPLES   = 1_000    # posterior samples (Monte Carlo draws) for reference model. Tutorial recommends num_samples=1000
_EPOCHS_REF     = 250       # NB regression. Tutorial recommends max_epochs=250 for reference signature estimation.
_EPOCHS_SPATIAL = 30_000    # Cell2location spatial mapping. Tutorial recommends max_epochs=30000.
_REF_BATCH_SIZE  = 2_500    # posterior sampling (number of cells per sampling batch). Tutorial recommends batch_size=2500, which is higher than scVI.
_TRAIN_SIZE      = 1        # training batch size for spatial model. Tutorial recommends train_size=1.

_ACCELERATOR     = get_accelerator()  # auto-selects: "gpu" (CUDA) > "mps" (Apple Silicon) > "cpu"

_EXPORT_POSTERIOR_SUBPROCESS_SCRIPT = Path(__file__).resolve().parent / "_export_posterior_subprocess.py"


def _run_export_posterior_subprocess(
    model,
    model_class: str,
    ckpt_dir: str,
    export_h5ad: str,
    num_samples: int,
    batch_size: int,
    keep_checkpoints: bool = False,
) -> ad.AnnData:
    """
    Save `model` checkpoint, run `_export_posterior_subprocess.py` as a child
    process with `accelerator="cpu"`, and return the resulting AnnData.

    Export runs in a separate process because Pyro's global parameter store
    retains GPU/MPS tensors from training; reloading the model in the same
    process does not clear that state and causes intermittent all-NaN posteriors.

    Parameters
    ----------
    model : RegressionModel or Cell2location
        Trained model to checkpoint and export the posterior for.
    model_class : str
        "RegressionModel" or "Cell2location" — tells the subprocess which class to load.
    ckpt_dir : str
        Directory to write the model checkpoint; removed on success unless `keep_checkpoints`.
    export_h5ad : str
        Path where the subprocess writes the posterior-exported AnnData; removed on success unless `keep_checkpoints`.
    keep_checkpoints : bool
        If True, skip removal of `ckpt_dir` and `export_h5ad` after reading the result.
    num_samples : int
        Posterior samples to draw.
    batch_size : int
        Batch size for posterior sampling.

    Returns
    -------
    anndata.AnnData
        AnnData with posterior export results in .varm/.obsm/.uns.
    """
    print(f"[INFO] Saving model checkpoint to {ckpt_dir} ...")
    model.save(ckpt_dir, overwrite=True, save_anndata=True)

    print("[INFO] Exporting posterior in a fresh subprocess ...")
    subprocess.run(
        [
            sys.executable,
            str(_EXPORT_POSTERIOR_SUBPROCESS_SCRIPT),
            "--checkpoint-dir", ckpt_dir,
            "--model-class", model_class,
            "--output-h5ad", export_h5ad,
            "--posterior-num-samples", str(num_samples),
            "--batch-size", str(batch_size),
        ],
        check=True,
    )
    result = ad.read_h5ad(export_h5ad)
    if not keep_checkpoints:
        os.remove(export_h5ad)
        shutil.rmtree(ckpt_dir)
    return result


def _save_training_history_plot(model, output_path: str, title: str) -> None:
    """Save ELBO training-loss curve to `output_path` as a PNG."""
    fig, ax = plt.subplots(figsize=(6, 4))
    model.plot_history(ax=ax)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"[INFO] Wrote training history plot to {output_path}")


def _load_reference(ref_h5ad: str, cell_type_key: str, min_cells_per_type: int) -> ad.AnnData:
    """
    Load scRNA-seq reference and validate that cell type labels are present in .obs.

    Parameters
    ----------
    ref_h5ad : str
        Path to all-gene reference h5ad (QC-filtered, raw counts). Must contain `cell_type_key` in .obs.
    cell_type_key : str
        Column in .obs to use as the cell type label.
    min_cells_per_type : int
        Drop cell types with fewer than this many cells (scvi-tools warns below 3).

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

    counts = adata.obs["cell_type"].value_counts()
    sparse_types = counts[counts < min_cells_per_type]
    if len(sparse_types) > 0:
        print(f"  [WARN] Dropping {len(sparse_types)} cell type(s) with fewer than "
              f"{min_cells_per_type} cells: "
              + ", ".join(f"{t} ({n})" for t, n in sparse_types.items()))
        adata = adata[~adata.obs["cell_type"].isin(sparse_types.index)].copy()

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
        Erroneous sample ID to replace.
    corrected : str
        Canonical sample ID to substitute.

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
    Load a Visium spatial h5ad, swap var_names to Ensembl IDs if available,
    correct PMBDS → PMDBS typo, remove MT genes, and validate integer counts.

    Parameters
    ----------
    spatial_h5ad : str
        Path to per-sample Visium h5ad (raw UMI counts, post-QC).
    sample_id : str
        Sample identifier; "ASAP_PMBDS_" is auto-corrected to "ASAP_PMDBS_".

    Returns
    -------
    tuple
        (anndata.AnnData, str) — prepared spatial AnnData and corrected sample_id.
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
    output_dir: str,
    num_samples: int,
    posterior_batch_size: int,
    keep_checkpoints: bool = False,
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
    output_dir : str
        Directory to write Stage 1 checkpoint and posterior-export AnnData.

    Returns
    -------
    pandas.DataFrame
        Estimated signatures (genes × cell types) — the `inf_aver` matrix.
    """
    print(f"[INFO] Setting up NB regression model ({adata_ref.n_obs:,} cells)...")
    RegressionModel.setup_anndata(
        adata_ref,
        batch_key=batch_key,
        labels_key="cell_type",
    )
    ref_model = RegressionModel(adata_ref)
    ref_model.train(max_epochs=max_epochs,
                    accelerator=accelerator,
                    early_stopping=True,
                    early_stopping_monitor="elbo_train")
    _save_training_history_plot(
        ref_model,
        os.path.join(output_dir, "stage1_ref_training_history.png"),
        "Stage 1: NB regression (reference signatures)",
    )

    adata_ref = _run_export_posterior_subprocess(
        ref_model,
        model_class="RegressionModel",
        ckpt_dir=os.path.join(output_dir, "stage1_ref_model_checkpoint"),
        export_h5ad=os.path.join(output_dir, "stage1_ref_posterior_export.h5ad"),
        num_samples=num_samples,
        batch_size=posterior_batch_size,
        keep_checkpoints=keep_checkpoints,
    )

    # np.asarray() needed: varm returns a prefixed-column DataFrame; columns= reindexes not relabels.
    factor_names = adata_ref.uns["mod"]["factor_names"]
    inf_aver = pd.DataFrame(
        np.asarray(adata_ref.varm["means_per_cluster_mu_fg"]),
        index=adata_ref.var_names,
        columns=factor_names,
    )
    if inf_aver.isna().values.all():
        raise RuntimeError(
            "Reference signatures (inf_aver) are entirely NaN after posterior export. "
            "This indicates a numerical issue in posterior sampling, not a data problem."
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
    output_dir: str,
    num_samples: int,
    posterior_batch_size: int,
    keep_checkpoints: bool = False,
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
    output_dir : str
        Directory to write Stage 2 checkpoint and posterior-export AnnData.

    Returns
    -------
    anndata.AnnData
        Spatial adata with posterior estimates in .obsm.
    """
    # Filter to genes present in both spatial and reference signatures
    shared_genes = adata_spatial.var_names.intersection(inf_aver.index)
    if len(shared_genes) == 0:
        raise ValueError("No shared genes between spatial adata and reference signatures.")
    print(f"  Shared genes: {len(shared_genes):,} / {adata_spatial.n_vars:,} spatial, "
          f"{len(inf_aver):,} reference")

    adata_spatial = adata_spatial[:, shared_genes].copy()
    # MPS only supports float32; cast before Cell2location registers cell_state_df as a buffer.
    inf_aver_filtered = inf_aver.loc[shared_genes].astype(np.float32)

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
        early_stopping=True,
        early_stopping_monitor="elbo_train",
    )
    _save_training_history_plot(
        model,
        os.path.join(output_dir, "stage2_spatial_training_history.png"),
        "Stage 2: Cell2location (spatial mapping)",
    )

    adata_spatial = _run_export_posterior_subprocess(
        model,
        model_class="Cell2location",
        ckpt_dir=os.path.join(output_dir, "stage2_cell2location_model_checkpoint"),
        export_h5ad=os.path.join(output_dir, "stage2_spatial_posterior_export.h5ad"),
        num_samples=num_samples,
        batch_size=posterior_batch_size,
        keep_checkpoints=keep_checkpoints,
    )
    return adata_spatial


def _add_cell_type_fractions(adata_spatial: ad.AnnData) -> None:
    """Normalise q05 abundances to per-spot cell type fractions; add to .obsm and .obs. In place."""
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
    PCA → neighbours → UMAP → Leiden on `cell_type_fractions`. Modifies `adata_spatial` in place.

    Parameters
    ----------
    adata_spatial : anndata.AnnData
        Spatial adata with `cell_type_fractions` in .obsm.
    n_comps : int
        PCA components.
    n_neighbors : int
        Neighbourhood size for graph construction.
    leiden_res : list
        Leiden resolutions; one clustering per resolution added to .obs.

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
        sc.tl.leiden(adata_tmp, resolution=res, key_added=key,
                     flavor="igraph", n_iterations=2, directed=False)
        adata_spatial.obs[key] = adata_tmp.obs[key].values
        print(f"  Leiden {res}: {adata_tmp.obs[key].nunique()} clusters")


def deconvolute_visium(
    sample_id: str | None,
    spatial_h5ad: str | None,
    ref_h5ad: str | None,
    reference_model: str | None,
    batch_key: str,
    cell_type_key: str,
    min_cells_per_type: int,
    n_comps: int,
    n_neighbors: int,
    leiden_res: list,
    n_cells_per_spot: int,
    detection_alpha: int,
    epochs_ref: int,
    epochs_spatial: int,
    num_samples: int,
    posterior_batch_size: int,
    accelerator: str,
    output_dir: str,
    promote: bool,
    keep_checkpoints: bool = False,
    build_reference_only: bool = False,
    skip_clustering: bool = False,
) -> None:
    """
    Full Cell2location deconvolution pipeline for one Visium sample.

    Parameters
    ----------
    sample_id : str
        Sample identifier; "ASAP_PMBDS_" is auto-corrected to "ASAP_PMDBS_".
    spatial_h5ad : str
        Per-sample Visium adata path (raw UMI counts, post-QC from `spatial-transcriptomics-wf`).
    ref_h5ad : str or None
        scRNA-seq reference h5ad. Mutually exclusive with `reference_model`.
    reference_model : str or None
        Pre-computed `<ref_stem>_inf_aver.csv.gz`; skips Stage 1. Mutually exclusive with `ref_h5ad`.
    batch_key : str
        Obs column used as batch covariate in NB regression.
    cell_type_key : str
        Reference .obs column to use as cell type label.
    min_cells_per_type : int
        Drop reference cell types with fewer than this many cells before training.
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
    keep_checkpoints : bool
        If True, retain model checkpoint directories and intermediate posterior h5ad files.
    build_reference_only : bool
        If True, run Stage 1 and write `inf_aver.csv.gz`, then exit without Stage 2/3.
        Requires `ref_h5ad`; `sample_id` and `spatial_h5ad` are ignored.
    skip_clustering : bool
        If True, skip Stage 3 downstream clustering (PCA/neighbours/UMAP/Leiden) and
        write Stage 1/2 outputs only. Useful to isolate Stage 3-only failures (e.g.
        missing igraph/leidenalg) without losing a completed Stage 1/2 run.

    Returns
    -------
    None
    """
    out_inf_aver = (
        os.path.join(output_dir, f"{Path(ref_h5ad).stem}_inf_aver.csv.gz")
        if ref_h5ad else None
    )
    if not build_reference_only:
        sample_id = _correct_sample_id(sample_id)
        out_h5ad  = os.path.join(output_dir, f"{sample_id}.deconvoluted.h5ad")
        out_csv   = os.path.join(output_dir, f"{sample_id}.deconvolution_metadata.csv.gz")

    if not promote:
        print(f"[DRY-RUN] Would write:")
        if build_reference_only:
            print(f"  {out_inf_aver}")
        else:
            print(f"  {out_h5ad}")
            print(f"  {out_csv}")
            if ref_h5ad:
                print(f"  {out_inf_aver}")
        print("Re-run with --promote to write outputs.")
        return

    os.makedirs(output_dir, exist_ok=True)
    ref_stem   = Path(ref_h5ad).stem if ref_h5ad else (sample_id or "run")
    out_timing = os.path.join(
        output_dir,
        f"{ref_stem}_ref_build.timing.tsv" if build_reference_only else f"{sample_id}.timing.tsv",
    )
    timing = []
    t_total = time.time()

    # --- Step 1: reference signatures ---
    t0 = time.time()
    if reference_model:
        print(f"[INFO] Loading pre-computed scRNA-seq reference signatures: {reference_model}")
        inf_aver = pd.read_csv(reference_model, index_col=0, compression="gzip")
        print(f"  Signatures: {inf_aver.shape[0]:,} genes × {inf_aver.shape[1]} cell types")
        timing.append(("stage1_nb_regression_skipped", time.time() - t0))
    else:
        adata_ref = _load_reference(ref_h5ad, cell_type_key, min_cells_per_type)
        adata_ref = _filter_mt_genes(adata_ref, "reference")
        _ensure_raw_int_counts(adata_ref, "reference")
        inf_aver = _build_reference_signatures(adata_ref, batch_key, epochs_ref, accelerator, output_dir, num_samples, posterior_batch_size, keep_checkpoints)
        del adata_ref
        print(f"[INFO] Writing {out_inf_aver} ...")
        inf_aver.to_csv(out_inf_aver, compression="gzip")
        timing.append(("stage1_nb_regression", time.time() - t0))
        if build_reference_only:
            timing.append(("total", time.time() - t_total))
            with open(out_timing, "w") as f:
                f.write("step\tduration_seconds\n")
                for step, secs in timing:
                    f.write(f"{step}\t{secs:.1f}\n")
            print(f"[INFO] Timing written to {out_timing}")
            print(f"[DONE] Reference built: {out_inf_aver}")
            return

    # --- Step 2: spatial mapping ---
    t0 = time.time()
    adata_spatial, sample_id = _load_spatial(spatial_h5ad, sample_id)
    adata_spatial = _run_cell2location_visium(
        adata_spatial, inf_aver, n_cells_per_spot, detection_alpha, epochs_spatial, accelerator, output_dir, num_samples, posterior_batch_size, keep_checkpoints
    )
    timing.append(("stage2_spatial_mapping", time.time() - t0))

    # --- Step 3: downstream analysis (optional clustering) ---
    t0 = time.time()
    _add_cell_type_fractions(adata_spatial)
    if skip_clustering:
        print("[INFO] Skipping Stage 3 downstream clustering (--skip-clustering).")
        timing.append(("stage3_downstream_skipped", time.time() - t0))
    else:
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
    adata_spatial.obs[meta_cols].to_csv(out_csv, compression="gzip")
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
    """Parse arguments and run deconvolute_visium."""
    args = _parse_args()

    ref_adata = args.sc_rnaseq_reference_adata
    ref_model = args.sc_rnaseq_reference_cell2location_model
    build_reference_only = args.build_sc_rnaseq_reference_only
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
    if build_reference_only and not ref_adata:
        raise SystemExit(
            "error: --build-sc-rnaseq-reference-only requires --sc-rnaseq-reference-adata."
        )
    if not build_reference_only:
        if not args.sample_id:
            raise SystemExit("error: --sample-id is required unless --build-sc-rnaseq-reference-only is set.")
        if not args.spatial_h5ad:
            raise SystemExit("error: --spatial-h5ad is required unless --build-sc-rnaseq-reference-only is set.")

    log_run_command(Path(args.output_dir) / f"{Path(sys.argv[0]).stem}.log")

    deconvolute_visium(
        sample_id          = args.sample_id,
        spatial_h5ad       = args.spatial_h5ad,
        ref_h5ad           = ref_adata,
        reference_model    = ref_model,
        batch_key          = args.batch_key,
        cell_type_key      = args.cell_type_key,
        min_cells_per_type = args.min_cells_per_type,
        n_comps            = args.n_comps,
        n_neighbors        = args.n_neighbors,
        leiden_res         = args.leiden_res,
        n_cells_per_spot   = args.n_cells_per_spot,
        detection_alpha    = args.detection_alpha,
        epochs_ref         = args.max_epochs_ref,
        epochs_spatial     = args.max_epochs_spatial,
        num_samples        = args.posterior_num_samples,
        posterior_batch_size = args.posterior_batch_size,
        accelerator        = args.accelerator,
        output_dir         = args.output_dir,
        promote              = args.promote,
        keep_checkpoints     = args.keep_checkpoints,
        build_reference_only = build_reference_only,
        skip_clustering       = args.skip_clustering,
    )


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--sample-id",                    default=None,
        help="Spatial sample identifier. 'ASAP_PMBDS_' is auto-corrected to 'ASAP_PMDBS_'.\n"
             "Required unless --build-sc-rnaseq-reference-only is set.")
    parser.add_argument("--spatial-h5ad",                default=None,
        help="Per-sample spatial h5ad (raw UMI counts, post-QC).\n"
             "Required unless --build-sc-rnaseq-reference-only is set.")
    parser.add_argument("--sc-rnaseq-reference-adata",    default=None,
        help="scRNA-seq reference h5ad (QC-filtered cells × all genes, raw counts).\n"
             "Trains Stage 1 NB regression and saves <ref_stem>_inf_aver.csv.gz to --output-dir.\n"
             "Mutually exclusive with --sc-rnaseq-reference-cell2location-model.")
    parser.add_argument("--sc-rnaseq-reference-cell2location-model", default=None,
        help="Path to a pre-computed scRNA-seq reference signatures file <ref_stem>_inf_aver.csv.gz produced by a previous run.\n"
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
    parser.add_argument("--min-cells-per-type",           type=int, default=10,
        help="Drop reference cell types with fewer than this many cells before\n"
             "training. scvi-tools warns NB regression may not train properly\n"
             "for categories with fewer than 3 cells. Default: 10")
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
    parser.add_argument("--posterior-num-samples",                  type=int, default=_REF_N_SAMPLES,
        help=f"Posterior samples drawn during export_posterior. Reduce to lower peak RAM\n"
             f"(e.g. 200–300) at the cost of noisier estimates. Default: {_REF_N_SAMPLES}")
    parser.add_argument("--posterior-batch-size",         type=int, default=_REF_BATCH_SIZE,
        help=f"Batch size for posterior sampling. Reduce to lower peak RAM\n"
             f"(e.g. 500–1000). Default: {_REF_BATCH_SIZE}")
    parser.add_argument("--keep-checkpoints",             action="store_true", default=False,
        help="Retain model checkpoint directories and intermediate posterior h5ad files\n"
             "after each stage. Useful for debugging. Default: remove on success.")
    parser.add_argument("--build-sc-rnaseq-reference-only", action="store_true", default=False,
        help="Train Stage 1 NB regression and write <ref_stem>_inf_aver.csv.gz, then exit.\n"
             "Requires --sc-rnaseq-reference-adata. --sample-id and --spatial-h5ad are not required.\n"
             "Use --sc-rnaseq-reference-cell2location-model <ref_stem>_inf_aver.csv.gz in subsequent runs.")
    parser.add_argument("--skip-clustering",              action="store_true", default=False,
        help="Skip Stage 3 downstream clustering (PCA/neighbours/UMAP/Leiden) on cell type\n"
             "fractions. Requires the igraph/leidenalg packages when not skipped. Still\n"
             "writes Stage 1/2 outputs. Default: False")
    return parser.parse_args()


if __name__ == "__main__":
    main()
