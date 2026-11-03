#!/usr/bin/env python3
"""Deconvolute one Visium sample using Cell2location against an scRNA-seq reference.

--run-stages examples
---------------------
1,2,3 (default)  Full pipeline.
1                Build reference only.
1,2              Train, deconvolute.
2,3              Load a pre-computed reference, deconvolute, cluster.
2                Load a pre-computed reference, deconvolute only.
3                Load a pre-computed deconvolution, cluster only.

Outputs
-------
{sample_id}.deconvoluted.h5ad
    Spatial AnnData with Cell2location results in .obsm and per-spot fraction in .obs.
    Written when Stage 2 runs without Stage 3.
{sample_id}.deconvolution_metadata.csv.gz
    Per-spot table: cell type fractions + obs columns for downstream reporting.
    Written when Stage 2 runs without Stage 3.
{sample_id}.deconvoluted_clustered.h5ad
    Same as above, plus downstream Leiden clustering.
    Written when Stage 3 runs.
{sample_id}.deconvolution_clustering_metadata.csv.gz
    Same as above, plus Leiden cluster assignments.
    Written when Stage 3 runs.
{ref_stem}_inf_aver.csv.gz
    Reference signatures (genes × cell types); only written when 1 is in --run-stages.
    {ref_stem} is the stem of the reference h5ad filename.
stage1_ref_training_history.png / {sample_id}.stage2_spatial_training_history.png
    ELBO training-loss curves per stage; only written when that stage trains.
    Stage 1's plot is {ref_stem}-prefixed (not bare) when Stage 2 also runs.
{sample_id}.timing.tsv / {ref_stem}_ref_build.timing.tsv
    Per-step duration (seconds); skipped stages recorded as 0.0 under a "_skipped" name.
    {ref_stem} replaces {sample_id} in Stage-1-only mode.
cell2location_deconv.log
    The exact invocation, for provenance.
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
import scvi
import torch
from cell2location.models import Cell2location, RegressionModel
from scipy.sparse import issparse

# Container builds copy utils/ next to this script (no .git present); only insert repo_root
# when running from a git clone (local dev) -- otherwise rely on PYTHONPATH.
repo_root = next((p for p in Path(__file__).resolve().parents if (p / ".git").exists()), None)
if repo_root is not None:
    sys.path.insert(0, str(repo_root))

from utils.accelerator import get_accelerator
from utils.genomics_constants import MT_GENE_ENSEMBL_IDS
from utils.logging_extra import log_run_command


# Source https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_short_demo.html
_CELLS_PER_SPOT = 8  # mouse brain (tutorial's N=30 is for denser human lymph node)

# Source https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html
_DETECTION_ALPHA = 100      # tutorial range 20 (high tech. variability) - 200 (low)
_REF_N_SAMPLES   = 1_000    # Monte Carlo posterior samples (tutorial default)
_EPOCHS_REF     = 250       # NB regression epochs (tutorial default, early_stopping is included)
_EPOCHS_SPATIAL = 30_000    # spatial mapping epochs (tutorial default, early_stopping is included)
_REF_BATCH_SIZE  = 2_500    # posterior sampling batch size (tutorial default)
_TRAIN_SIZE      = 1        # tutorial default

_ACCELERATOR     = get_accelerator()  # auto-selects: "gpu" (CUDA) > "mps" (Apple Silicon) > "cpu"
_SEED            = 0  # seed for same machine system. Note that "gpu" vs. "mps" still introduces non-determinism.

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
    Checkpoint `model`, run `_export_posterior_subprocess.py` as a CPU child process, and
    return the resulting AnnData.

    A separate process is needed because Pyro's global parameter store retains GPU/MPS
    tensors from training, causing intermittent all-NaN posteriors if reloaded in-process.

    Parameters
    ----------
    model : RegressionModel or Cell2location
        Trained model to checkpoint and export the posterior for.
    model_class : str
        "RegressionModel" or "Cell2location"; tells the subprocess which class to load.
    ckpt_dir : str
        Directory to write the model checkpoint; removed on success unless `keep_checkpoints`.
    export_h5ad : str
        Path for the subprocess's posterior-exported AnnData; removed on success unless `keep_checkpoints`.
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
        Reference h5ad (raw UMI counts, all genes--not only HVGs); must contain `cell_type_key` in .obs.
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
    """Return `sample_id` with the known ASAP typo PMBDS → PMDBS corrected.
    This is a temporary fix while spatial-transcriptomics-wf is re-run to
    fix underlying typos. """
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
    Load a Visium h5ad: swap var_names to Ensembl IDs if available, correct the
    PMBDS/PMDBS typo, remove MT genes, and validate integer counts.

    Parameters
    ----------
    spatial_h5ad : str
        Path to per-sample Visium h5ad (raw UMI counts, all genes--not only HVGs).
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
    stem: str | None = None,
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
    stem : str or None
        Prefixes the training-history plot filename if given; bare filename otherwise.

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
    plot_filename = (
        "stage1_ref_training_history.png" if stem is None
        else f"{stem}.stage1_ref_training_history.png"
    )
    _save_training_history_plot(
        ref_model,
        os.path.join(output_dir, plot_filename),
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
    sample_id: str = "sample",
) -> ad.AnnData:
    """
    Map cell types onto Visium spots using Cell2location.

    Parameters
    ----------
    adata_spatial : anndata.AnnData
        Per-sample Visium adata (spots × genes, raw UMI counts, all genes--not only HVGs).
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
    sample_id : str
        Prefixes the training-history plot filename.

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
        os.path.join(output_dir, f"{sample_id}.stage2_spatial_training_history.png"),
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
    adata_spatial.obsm["cell2location_X_pca_fractions"]  = adata_tmp.obsm["X_pca"]
    adata_spatial.obsm["cell2location_X_umap_fractions"] = adata_tmp.obsm["X_umap"]
    adata_spatial.obsp["cell2location_connectivities"]   = adata_tmp.obsp["connectivities"]
    adata_spatial.obsp["cell2location_distances"]        = adata_tmp.obsp["distances"]
    adata_spatial.uns["cell2location_neighbors"]         = adata_tmp.uns["neighbors"]
    for res in leiden_res:
        key = f"cell2location_leiden_res_{str(res).replace('.', '_')}"
        sc.tl.leiden(adata_tmp, resolution=res, key_added=key,
                     flavor="igraph", n_iterations=2, directed=False)
        adata_spatial.obs[key] = adata_tmp.obs[key].values
        print(f"  Leiden {res}: {adata_tmp.obs[key].nunique()} clusters")


def _write_timing_tsv(path: str, timing: list) -> None:
    """Write a (step, duration_seconds) list to `path` as a 2-column TSV."""
    with open(path, "w") as f:
        f.write("step\tduration_seconds\n")
        for step, secs in timing:
            f.write(f"{step}\t{secs:.1f}\n")


def deconvolute_visium(
    sample_id: str | None,
    spatial_h5ad: str | None,
    ref_h5ad: str | None,
    reference_model: str | None,
    spatial_deconvolution_adata: str | None,
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
    run_stages: set,
    write_h5ad: bool = True,
    keep_checkpoints: bool = False,
) -> None:
    """
    Full Cell2location deconvolution pipeline for one Visium sample.

    Parameters
    ----------
    sample_id : str
        Sample identifier; "ASAP_PMBDS_" is auto-corrected to "ASAP_PMDBS_".
        Required when Stage 2 and/or 3 is in `run_stages`.
    spatial_h5ad : str
        Per-sample Visium adata path. Required when Stage 2 is in `run_stages`.
    ref_h5ad : str or None
        scRNA-seq reference h5ad. Required when Stage 1 is in `run_stages`.
    reference_model : str or None
        Pre-computed `<ref_stem>_inf_aver.csv.gz`. Required when Stage 2 is in `run_stages`
        without Stage 1.
    spatial_deconvolution_adata : str or None
        Pre-computed `.deconvoluted.h5ad`. Required when Stage 3 is in `run_stages` without Stage 2.
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
    run_stages : set of int
        Non-empty subset of {1, 2, 3} selecting which stages run; see module docstring for
        full semantics and examples.
    write_h5ad : bool
        If False, skip writing the `.h5ad` output and only write the metadata CSV.
        Has no effect when neither Stage 2 nor Stage 3 is in `run_stages`.
    keep_checkpoints : bool
        If True, retain model checkpoint directories and intermediate posterior h5ad files.

    Returns
    -------
    None
    """
    run_stage1 = 1 in run_stages
    run_stage2 = 2 in run_stages
    run_stage3 = 3 in run_stages
    write_sample_outputs = run_stage2 or run_stage3

    out_inf_aver = (
        os.path.join(output_dir, f"{Path(ref_h5ad).stem}_inf_aver.csv.gz")
        if ref_h5ad else None
    )
    if write_sample_outputs:
        sample_id = _correct_sample_id(sample_id)
        if run_stage3:
            out_h5ad = os.path.join(output_dir, f"{sample_id}.deconvoluted_clustered.h5ad")
            out_csv  = os.path.join(output_dir, f"{sample_id}.deconvolution_clustering_metadata.csv.gz")
        else:
            out_h5ad = os.path.join(output_dir, f"{sample_id}.deconvoluted.h5ad")
            out_csv  = os.path.join(output_dir, f"{sample_id}.deconvolution_metadata.csv.gz")

    if not promote:
        print(f"[DRY-RUN] Would write:")
        if run_stage1:
            print(f"  {out_inf_aver}")
        if write_sample_outputs:
            if write_h5ad:
                print(f"  {out_h5ad}")
            print(f"  {out_csv}")
        print("Re-run with --promote to write outputs.")
        return

    os.makedirs(output_dir, exist_ok=True)
    ref_stem   = Path(ref_h5ad).stem if ref_h5ad else (sample_id or "run")
    out_timing = os.path.join(
        output_dir,
        f"{sample_id}.timing.tsv" if write_sample_outputs else f"{ref_stem}_ref_build.timing.tsv",
    )
    timing = []
    t_total = time.time()

    # --- Stage 1: reference signatures ---
    t0 = time.time()
    if run_stage1:
        adata_ref = _load_reference(ref_h5ad, cell_type_key, min_cells_per_type)
        adata_ref = _filter_mt_genes(adata_ref, "reference")
        _ensure_raw_int_counts(adata_ref, "reference")
        inf_aver = _build_reference_signatures(
            adata_ref, batch_key, epochs_ref, accelerator, output_dir, num_samples, posterior_batch_size,
            keep_checkpoints, stem=(ref_stem if write_sample_outputs else None),
        )
        del adata_ref
        print(f"[INFO] Writing {out_inf_aver} ...")
        inf_aver.to_csv(out_inf_aver, compression="gzip")
        timing.append(("stage1_nb_regression", time.time() - t0))
    elif run_stage2:
        print(f"[INFO] Loading pre-computed scRNA-seq reference signatures: {reference_model}")
        inf_aver = pd.read_csv(reference_model, index_col=0, compression="gzip")
        print(f"  Signatures: {inf_aver.shape[0]:,} genes × {inf_aver.shape[1]} cell types")
        timing.append(("stage1_nb_regression_skipped", time.time() - t0))
    else:
        timing.append(("stage1_nb_regression_skipped", time.time() - t0))

    if run_stage1 and not write_sample_outputs:
        timing.append(("total", time.time() - t_total))
        _write_timing_tsv(out_timing, timing)
        print(f"[INFO] Timing written to {out_timing}")
        print(f"[DONE] Reference built: {out_inf_aver}")
        return

    # --- Stage 2: spatial mapping ---
    t0 = time.time()
    if run_stage2:
        adata_spatial, sample_id = _load_spatial(spatial_h5ad, sample_id)
        adata_spatial = _run_cell2location_visium(
            adata_spatial, inf_aver, n_cells_per_spot, detection_alpha, epochs_spatial, accelerator, output_dir, num_samples, posterior_batch_size,
            keep_checkpoints, sample_id=sample_id,
        )
        timing.append(("stage2_spatial_mapping", time.time() - t0))
    else:
        print(f"[INFO] Loading pre-computed Stage 2 output: {spatial_deconvolution_adata}")
        adata_spatial = ad.read_h5ad(spatial_deconvolution_adata)
        timing.append(("stage2_spatial_mapping_skipped", time.time() - t0))

    # --- Stage 3: downstream analysis (clustering) ---
    t0 = time.time()
    _add_cell_type_fractions(adata_spatial)
    if run_stage3:
        _downstream_analysis(adata_spatial, n_comps, n_neighbors, leiden_res)
        timing.append(("stage3_downstream", time.time() - t0))
    else:
        timing.append(("stage3_downstream_skipped", time.time() - t0))

    # --- Write outputs ---
    t0 = time.time()
    if write_h5ad:
        print(f"[INFO] Writing {out_h5ad} ...")
        adata_spatial.write_h5ad(out_h5ad, compression="lzf")
    meta_cols = (
        [c for c in adata_spatial.obs.columns if c.startswith("frac_")]
        + [c for c in adata_spatial.obs.columns if c.startswith("cell2location_leiden_res_")]
        + [c for c in adata_spatial.obs.columns if not c.startswith("frac_") and not c.startswith("cell2location_leiden_res_")]
    )
    adata_spatial.obs[meta_cols].to_csv(out_csv, compression="gzip")
    print(f"[INFO] Writing {out_csv} ...")
    timing.append(("write_outputs", time.time() - t0))
    timing.append(("total", time.time() - t_total))

    _write_timing_tsv(out_timing, timing)
    print(f"[INFO] Timing written to {out_timing}")
    print(f"[DONE] {sample_id}")


def _validate_run_stages(
    run_stages: str,
    sc_rnaseq_reference_adata: str | None,
    sc_rnaseq_reference_cell2location_model: str | None,
    spatial_h5ad: str | None,
    spatial_deconvolution_adata: str | None,
    sample_id: str | None,
) -> set[int]:
    """
    Parse `--run-stages` and validate every stage-dependent CLI argument
    to ensure that all required and no extra inputs for each stage are correctly specified.

    Parameters
    ----------
    run_stages : str
        Raw `--run-stages` value, e.g. "2,3".
    sc_rnaseq_reference_adata : str or None
        `--sc-rnaseq-reference-adata` value; H5AD to train Stage 1's reference signature.
    sc_rnaseq_reference_cell2location_model : str or None
        `--sc-rnaseq-reference-cell2location-model` value; pre-computed signature (inf_aver file)
        for Stage 2 without Stage 1.
    spatial_h5ad : str or None
        `--spatial-h5ad` value; raw spatial input consumed by Stage 2.
    spatial_deconvolution_adata : str or None
        `--spatial-deconvolution-adata` value; pre-computed Stage 2 output for Stage 3
        without Stage 2.
    sample_id : str or None
        `--sample-id` value; names per-sample output files for Stage 2 and/or 3.

    Returns
    -------
    set[int]
        Parsed run_stages, e.g. `{2, 3}`.

    Raises
    ------
    SystemExit
        If any stage-dependent validation check fails.
    """
    try:
        parsed_run_stages = {int(s) for s in run_stages.split(",")}
    except ValueError:
        raise SystemExit(
            f"error: --run-stages must be a comma-separated list of integers from {{1,2,3}}, "
            f"got '{run_stages}'."
        )
    if not parsed_run_stages or not parsed_run_stages.issubset({1, 2, 3}):
        raise SystemExit(
            f"error: --run-stages must be a non-empty comma-separated subset of {{1,2,3}}, "
            f"got '{run_stages}'."
        )
    run_stage1 = 1 in parsed_run_stages
    run_stage2 = 2 in parsed_run_stages
    run_stage3 = 3 in parsed_run_stages

    # Gate 1 - scRNA-seq reference: either H5AD or precomputed inf_aver.
    if run_stage1:
        if not sc_rnaseq_reference_adata:
            raise SystemExit("error: --sc-rnaseq-reference-adata is required when Stage 1 is in --run-stages.")
        if sc_rnaseq_reference_cell2location_model:
            raise SystemExit("error: --sc-rnaseq-reference-cell2location-model is ignored when Stage 1 is in --run-stages.")
    else:
        if sc_rnaseq_reference_adata:
            raise SystemExit("error: --sc-rnaseq-reference-adata is ignored when Stage 1 is not in --run-stages.")
        if run_stage2:
            if not sc_rnaseq_reference_cell2location_model:
                raise SystemExit(
                    "error: --sc-rnaseq-reference-cell2location-model is required when Stage 2 "
                    "is in --run-stages without Stage 1."
                )
        elif sc_rnaseq_reference_cell2location_model:
            raise SystemExit(
                "error: --sc-rnaseq-reference-cell2location-model is only used when Stage 2 "
                "is in --run-stages without Stage 1."
            )

    # Gate 2 - Raw spatial H5AD
    if run_stage2:
        if not spatial_h5ad:
            raise SystemExit("error: --spatial-h5ad is required when Stage 2 is in --run-stages.")
    elif spatial_h5ad:
        raise SystemExit("error: --spatial-h5ad is only used when Stage 2 is in --run-stages.")

    # Gate 3 - Deconvoluted spatial H5AD
    if run_stage2:
        if spatial_deconvolution_adata:
            raise SystemExit(
                "error: --spatial-deconvolution-adata is only used when Stage 3 is in "
                "--run-stages without Stage 2."
            )
    elif run_stage3:
        if not spatial_deconvolution_adata:
            raise SystemExit(
                "error: --spatial-deconvolution-adata is required when Stage 3 is in "
                "--run-stages without Stage 2."
            )
    elif spatial_deconvolution_adata:
        raise SystemExit(
            "error: --spatial-deconvolution-adata is only used when Stage 3 is in "
            "--run-stages without Stage 2."
        )

    # sample_id
    if (run_stage2 or run_stage3) and not sample_id:
        raise SystemExit("error: --sample-id is required when Stage 2 and/or 3 is in --run-stages.")

    return parsed_run_stages

def main() -> None:
    args = _parse_args()

    # Fix seed for same system. VM cuda vs. Mac mps still can have differences.
    scvi.settings.seed = args.seed
    if torch.backends.mps.is_available():
        torch.mps.manual_seed(args.seed)

    run_stages = _validate_run_stages(
        run_stages=args.run_stages,
        sc_rnaseq_reference_adata=args.sc_rnaseq_reference_adata,
        sc_rnaseq_reference_cell2location_model=args.sc_rnaseq_reference_cell2location_model,
        spatial_h5ad=args.spatial_h5ad,
        spatial_deconvolution_adata=args.spatial_deconvolution_adata,
        sample_id=args.sample_id,
    )

    log_run_command(Path(args.output_dir) / f"{Path(sys.argv[0]).stem}.log")

    deconvolute_visium(
        sample_id          = args.sample_id,
        spatial_h5ad       = args.spatial_h5ad,
        ref_h5ad           = args.sc_rnaseq_reference_adata,
        reference_model    = args.sc_rnaseq_reference_cell2location_model,
        spatial_deconvolution_adata = args.spatial_deconvolution_adata,
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
        promote            = args.promote,
        run_stages         = run_stages,
        write_h5ad         = args.write_h5ad,
        keep_checkpoints   = args.keep_checkpoints,
    )


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    general = parser.add_argument_group("General (always applicable)")
    general.add_argument("--run-stages",                   type=str, default="1,2,3",
        help="Comma-separated subset of {1,2,3} — which stages to run this invocation:\n"
             "1 = NB regression (reference signatures)\n"
             "2 = Cell2location spatial mapping,\n"
             "3 = downstream Leiden clustering.\n"
             "Default: 1,2,3 (full pipeline).\n")
    general.add_argument("--output-dir",                  default=".",
        help="Directory to write output files. Default: current directory")
    general.add_argument("--promote",                     action="store_true", default=False,
        help="Write output files. Default is dry-run (prints what would happen).")
    general.add_argument("--write-h5ad",                  action=argparse.BooleanOptionalAction, default=True,
        help="Write the per-sample `.h5ad` output in addition to the metadata CSV.\n"
             "Use --no-write-h5ad to write only the metadata CSV. Has no effect in Stage-1-only mode.\n"
             "Default: True"),
    general.add_argument("--seed",                        type=int, default=_SEED,
        help="Random seed for same system torch/numpy/scvi-tools training.\n"
             f"Default: {_SEED}")

    stage1 = parser.add_argument_group("Stage 1 — NB regression (reference signatures)")
    stage1.add_argument("--sc-rnaseq-reference-adata",    default=None,
        help="scRNA-seq reference h5ad (QC-filtered cells × all genes--not only HVGs, raw UMI counts).\n"
             "Trains Stage 1 NB regression and saves <ref_stem>_inf_aver.csv.gz to --output-dir.\n"
             "Required when Stage 1 is in --run-stages.")
    stage1.add_argument("--batch-key",                    default="batch_id",
        help="Obs column used as batch covariate in the NB regression model.\n"
             "Default: batch_id")
    stage1.add_argument("--cell-type-key",                default="cell_type_pred",
        help="Column in the reference .obs to use as the cell type label.\n"
             "Default: cell_type_pred")
    stage1.add_argument("--min-cells-per-type",           type=int, default=10,
        help="Drop reference cell types with fewer than this many cells before\n"
             "training. scvi-tools warns NB regression may not train properly\n"
             "for categories with fewer than 3 cells. Default: 10")
    stage1.add_argument("--max-epochs-ref",               type=int, default=_EPOCHS_REF,
        help=f"NB regression training epochs. Default: {_EPOCHS_REF}")

    stage2 = parser.add_argument_group("Stage 2 — Cell2location spatial mapping")
    stage2.add_argument("--spatial-h5ad",                 default=None,
        help="Per-sample spatial h5ad (raw UMI counts, all genes--not only HVGs).\n"
             "Required when Stage 2 is in --run-stages.")
    stage2.add_argument("--sc-rnaseq-reference-cell2location-model", default=None,
        help="Path to a pre-computed scRNA-seq reference signatures file <ref_stem>_inf_aver.csv.gz produced by a previous run.\n"
             "Required when Stage 2 is in --run-stages without Stage 1.")
    stage2.add_argument("--n-cells-per-spot",             type=int, default=_CELLS_PER_SPOT,
        help=f"Cell2location N_cells_per_location prior. Default: {_CELLS_PER_SPOT}")
    stage2.add_argument("--detection-alpha",              type=int, default=_DETECTION_ALPHA,
        help=f"Cell2location detection rate prior (tutorial range 20–200).\n"
             f"Higher values assume low technical variability in RNA detection.\n"
             f"Default: {_DETECTION_ALPHA}")
    stage2.add_argument("--max-epochs-spatial",           type=int, default=_EPOCHS_SPATIAL,
        help=f"Cell2location training epochs. Default: {_EPOCHS_SPATIAL}")

    stage3 = parser.add_argument_group("Stage 3 — downstream Leiden clustering")
    stage3.add_argument("--n-comps",                      type=int, default=30,
        help="PCA components for downstream spot clustering. Default: 30")
    stage3.add_argument("--n-neighbors",                  type=int, default=15,
        help="Neighbourhood size for downstream spot clustering. Default: 15")
    stage3.add_argument("--leiden-res",                   type=float, nargs="+",
        default=[0.2, 0.5, 1.0],
        help="Leiden resolutions for downstream spot clustering. Default: 0.2 0.5 1.0")
    stage3.add_argument("--spatial-deconvolution-adata",  default=None,
        help="Path to a pre-computed <sample_id>.deconvoluted.h5ad from a prior Stage 2 run.\n"
             "Required when Stage 3 is in --run-stages without Stage 2.")

    multistage = parser.add_argument_group("Multi-stage (Stage 1 or 2, and Stage 2 or 3)")
    multistage.add_argument("--sample-id",                default=None,
        help="Spatial sample identifier. 'ASAP_PMBDS_' is auto-corrected to 'ASAP_PMDBS_'.\n"
             "Required when Stage 2 and/or 3 is in --run-stages.")
    multistage.add_argument("--accelerator",              default=_ACCELERATOR,
        help=f"PyTorch Lightning accelerator: 'mps' (Apple Silicon), 'gpu' (CUDA), 'cpu'.\n"
             f"Used by Stage 1 and/or Stage 2 training. Default: {_ACCELERATOR}")
    multistage.add_argument("--posterior-num-samples",    type=int, default=_REF_N_SAMPLES,
        help=f"Posterior samples drawn during export_posterior. Used by Stage 1 and/or Stage 2.\n"
             f"Reduce to lower peak RAM (e.g. 200–300) at the cost of noisier estimates.\n"
             f"Default: {_REF_N_SAMPLES}")
    multistage.add_argument("--posterior-batch-size",     type=int, default=_REF_BATCH_SIZE,
        help=f"Batch size for posterior sampling. Used by Stage 1 and/or Stage 2.\n"
             f"Reduce to lower peak RAM (e.g. 500–1000). Default: {_REF_BATCH_SIZE}")
    multistage.add_argument("--keep-checkpoints",         action="store_true", default=False,
        help="Retain model checkpoint directories and intermediate posterior h5ad files\n"
             "after each stage. Used by Stage 1 and/or Stage 2. Useful for debugging.\n"
             "Default: remove on success.")

    return parser.parse_args()


if __name__ == "__main__":
    main()
