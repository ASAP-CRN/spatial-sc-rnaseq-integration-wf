#!/usr/bin/env python3
"""
Merges cached H5ADs from four scRNA-seq sources (Kamath, Menon, Siletti, ASAP),
relevant for Parkinson's disease, in particular for identifying dopaminergic neurons (DA).
trains scVI/scANVI, uses Kamath labels for cell type annotation and computes UMAP.
Results are are stored under {output_dir}/{suffix}/, including an integrated_{suffix}.h5ad.
See reference_building/da_neurons_identification_analysis/README.md for more details.

This is a one-time analysis, not part of the spatial deconvolution WDL pipeline.

Relationship to other scripts:
- Run integrate_kamath_menon_siletti_asap.py after downsample_kamath_menon_siletti_asap.py
  with the same values for the following parameters:
  --kamath-non-da-per-type (or --kamath-all-random / --kamath-cell-type-min, both paired with
  --downsample-kamath), --downsample-menon,
  --siletti-scale and --downsample-asap-cohort so integrate_kamath_menon_siletti_asap.py
  finds the right cached files.
  These downsample_* and integrate_* scripts were split to:
  a) control the type of machine needed CPU (cheaper) vs. GPU (expensive)
  b) cache intermediate files to restart the run if needed
- Run integrate_menon_siletti_asap.py instead of integrate_kamath_menon_siletti_asap.py
  for a run without Kamath's data and using Siletti (Allen) labels as ground truth.

Integration strategy
--------------------
1. Load each source's cached AnnData, either preQC or postQC
2. Select genes for scVI training, either highly variable genes (HVGs), HVGs with or without 
   DA markers, or a custom gene list provided by the user.
3. Train scVI for a batch-corrected latent space.
4. Train scANVI semi-supervised: Kamath labels are used as ground truth for cell type annotation.
   Each other source's native labels are preserved in the H5AD for post-hoc agreement comparison.
5. Compute UMAP on the scANVI latent.


Usage example
-------------
A typical integrate_kamath_menon_siletti_asap.py run, with ~100k cells per dataset (total ~400k cells) could
run in a MacBook Pro with 36GB RAM and an Apple M3 Pro chip in 7-9 hours.

```
python3 integrate_kamath_menon_siletti_asap.py \
    --output-dir /path/to/reference_building/da_neurons_identification_analysis/output/ \
    --kamath-non-da-per-type 13000 \
    --downsample-menon 100000 \
    --siletti-scale 100000 \
    --downsample-asap-cohort 100000 \
    --gene-selection hvg_u_da_markers
```

Outputs
-------
All under {output_dir}/{suffix}/. For the example above, suffix will be:
`asap96k_kamath97k_menon99k_siletti100k_hvg_u_da_markers`
See reference_building/da_neurons_identification_analysis/README.md for more details.


"""

import os
import time
import argparse
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scvi
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
from pathlib import Path
from datetime import timedelta

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "utils"))
from _reference_building_utils import (
    _REFERENCE_DATA_DIR,
    _SILETTI_SCALE_PRESETS,
    _SILETTI_DEFAULT_SCALE,
    _ASAP_H5AD_PATH,
    _KAMATH_DA_MARKER_SYMBOLS,
    _DEFAULT_KAMATH_NON_DA_PER_TYPE,
    _DEFAULT_DOWNSAMPLE_MENON,
    _DEFAULT_DOWNSAMPLE_ASAP_COHORT,
    _qc_cache_paths,
    _sn_cache_path,
    _find_menon_h5ad,
    log_run_command,
)

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message=".*copying X to .raw.*")


def _resolve_marker_ensembl_ids(
    feature_name: pd.Series,
    marker_symbols: tuple,
    source_label: str,
) -> list[str]:
    """
    Resolve gene symbols to Ensembl IDs via a source's var['feature_name'].

    Kamath/ASAP's own var has no gene-symbol column post-load, and
    ad.concat(merge='same') doesn't reliably preserve one onto the merged
    object either, so this looks symbols up against one source (Menon or
    Siletti) BEFORE the merge, not against merged.var itself.

    Parameters
    ----------
    feature_name : pd.Series
        A source's `adata.var['feature_name']`, indexed by Ensembl ID.
    marker_symbols : tuple
        Gene symbols to resolve.
    source_label : str
        Label for the warning message (e.g. 'Menon').

    Returns
    -------
    list[str]
        Ensembl IDs found. Unmatched symbols are warned about, not raised on.
    """
    # feature_name is categorical and can have duplicate symbol values
    # (distinct Ensembl IDs sharing one symbol) -- both break .reindex().
    symbols_str = feature_name.astype(str)
    symbol_to_ensembl = pd.Series(symbols_str.index, index=symbols_str.values)
    symbol_to_ensembl = symbol_to_ensembl[~symbol_to_ensembl.index.duplicated(keep="first")]
    found = symbol_to_ensembl.reindex(marker_symbols)
    missing = [s for s, e in zip(marker_symbols, found) if pd.isna(e)]
    if missing:
        print(f"  Warning: {len(missing)}/{len(marker_symbols)} Kamath DA marker symbols "
              f"not found in {source_label}'s var['feature_name']: {missing}")
    return found.dropna().tolist()


def preprocess_merged(adata: ad.AnnData, n_hvgs: int) -> ad.AnnData:
    """
    Normalize, log-transform, and select HVGs (raw counts kept in layers['counts']).

    Matches sc-rnaseq-wf's cohort_analysis `process` task: Pearson residuals
    flavor over raw counts (not seurat_v3). Uses batch_key='source' (per
    dataset) rather than sc-rnaseq-wf's 'batch_id' (per sample), since there's
    no shared per-sample key across the four studies.

    Parameters
    ----------
    adata : ad.AnnData
        Merged AnnData with raw counts in layers['counts'].
    n_hvgs : int
        Number of HVGs to select.

    Returns
    -------
    ad.AnnData
        Same object, var['highly_variable'] populated, X holds log-norm counts.
    """
    print(f"  Normalizing and log-transforming")
    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    print(f"  Selecting {n_hvgs} HVGs (flavor='pearson_residuals', batch_key='source')")
    sc.experimental.pp.highly_variable_genes(
        adata,
        flavor="pearson_residuals",
        n_top_genes=n_hvgs,
        layer="counts",
        batch_key="source",
        check_values=True,
        subset=False,  # keep all genes; HVG flag used downstream for gene selection
    )
    n_hvg = adata.var["highly_variable"].sum()
    print(f"  HVGs selected: {n_hvg:,}")
    return adata


# ---------------------------------------------------------------------------
# Integration
# ---------------------------------------------------------------------------

def run_scvi(
    adata: ad.AnnData,
    n_latent: int,
    n_layers: int,
    max_epochs: int,
    batch_size: int,
    use_gpu: bool,
) -> tuple[ad.AnnData, scvi.model.SCVI]:
    """
    Train scVI on raw counts for a batch-corrected latent embedding.

    Trains on the gene panel `adata` already has -- gene selection
    happens earlier in main(), not here. Stores obsm['X_scVI'].
    """
    accelerator = "gpu" if use_gpu else "cpu"
    print(f"  Setting up scVI (n_latent={n_latent}, n_layers={n_layers}, "
          f"batch_key='batch', accelerator={accelerator})")

    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key="batch",
    )
    model = scvi.model.SCVI(
        adata,
        n_latent=n_latent,
        n_layers=n_layers,
        gene_likelihood="nb",
        dispersion="gene-batch",
    )
    print(model)

    model.train(
        max_epochs=max_epochs,
        batch_size=batch_size,
        accelerator=accelerator,
        early_stopping=True,
        early_stopping_patience=15,
        plan_kwargs={"lr": 1e-3},
    )

    adata.obsm["X_scVI"] = model.get_latent_representation()
    print(f"  scVI latent stored in obsm['X_scVI'] — shape: {adata.obsm['X_scVI'].shape}")
    return adata, model

def run_scanvi(
    adata: ad.AnnData,
    scvi_model: scvi.model.SCVI,
    labels_key: str,
    unlabeled_category: str,
    max_epochs: int,
    batch_size: int,
    use_gpu: bool,
) -> tuple[ad.AnnData, scvi.model.SCANVI]:
    """
    Train scANVI in semi-supervised mode starting from the scVI model.

    Kamath cells provide cell type labels (DA subtypes for dopaminergic neurons,
    broad type for others). Menon cells are unlabeled. scANVI simultaneously
    refines the latent space and predicts labels for unlabeled cells.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData used for scVI training, with obsm['X_scVI'].
    scvi_model : scvi.model.SCVI
        Trained scVI model to initialize scANVI from.
    labels_key : str
        obs column containing training labels for labeled cells.
    unlabeled_category : str
        Value in labels_key that marks unlabeled cells (Menon cells).
    max_epochs : int
        Maximum scANVI training epochs.
    batch_size : int
        Mini-batch size.
    use_gpu : bool
        Whether to use GPU.

    Returns
    -------
    ad.AnnData
        AnnData with obsm['X_scANVI'] and obs['cell_type_pred'] added.
    scvi.model.SCANVI
        Trained scANVI model.
    """
    accelerator = "gpu" if use_gpu else "cpu"
    n_labeled = (adata.obs[labels_key] != unlabeled_category).sum()
    n_unlabeled = (adata.obs[labels_key] == unlabeled_category).sum()
    print(
        f"  Setting up scANVI (labels_key='{labels_key}', "
        f"unlabeled='{unlabeled_category}')\n"
        f"  Labeled: {n_labeled:,}  Unlabeled: {n_unlabeled:,}"
    )

    model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        adata=adata,
        labels_key=labels_key,
        unlabeled_category=unlabeled_category,
    )
    print(model)

    model.train(
        max_epochs=max_epochs,
        batch_size=batch_size,
        accelerator=accelerator,
        early_stopping=True,
        early_stopping_patience=10,
    )

    adata.obsm["X_scANVI"] = model.get_latent_representation()
    adata.obs["cell_type_pred"] = model.predict()
    print(f"  scANVI latent stored in obsm['X_scANVI']")
    print(f"  Label predictions stored in obs['cell_type_pred']")
    return adata, model


# ---------------------------------------------------------------------------
# Post-integration
# ---------------------------------------------------------------------------

def transfer_da_subtypes(adata: ad.AnnData) -> ad.AnnData:
    """
    Populate obs['da_subtype'] for non-Kamath cells whose cell_type_pred
    matches one of Kamath's DA subtype labels. Non-DA predictions stay NaN.
    """
    kamath_da_labels = set(
        adata.obs.loc[adata.obs["source"] == "kamath", "da_subtype"]
        .dropna()
        .unique()
    )
    if not kamath_da_labels:
        print("  Warning: no Kamath DA subtype labels found — da_subtype transfer skipped")
        return adata

    print(f"  Kamath DA subtype labels: {sorted(kamath_da_labels)}")

    non_kamath_mask = adata.obs["source"] != "kamath"
    pred_is_da = adata.obs["cell_type_pred"].isin(kamath_da_labels)
    transfer_mask = non_kamath_mask & pred_is_da

    adata.obs.loc[transfer_mask, "da_subtype"] = adata.obs.loc[
        transfer_mask, "cell_type_pred"
    ]

    for source in adata.obs.loc[non_kamath_mask, "source"].unique():
        src_mask = adata.obs["source"] == source
        n = (src_mask & pred_is_da).sum()
        print(
            f"  DA subtype transferred to {n:,} {source} cells "
            f"({n / src_mask.sum():.1%} of {source} SN cells)"
        )
    return adata

def save_elbo_plot(
    model_scvi: scvi.model.SCVI,
    model_scanvi: scvi.model.SCANVI,
    output_dir: str,
) -> None:
    """
    Plot per-epoch train/validation ELBO for scVI and scANVI (best epoch
    marked) and save as `output_dir`/elbo_convergence.png.
    """
    purple, amber, grey, bg = "#3D1F61", "#D98E04", "#6B6470", "#FAF9FB"
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), facecolor=bg)

    for ax, (name, model) in zip(axes, [("scVI", model_scvi), ("scANVI", model_scanvi)]):
        train = model.history_["elbo_train"].iloc[:, 0].values
        val = model.history_["elbo_validation"].iloc[:, 0].values
        epochs = range(len(train))

        ax.set_facecolor(bg)
        ax.plot(epochs, train, color=purple, linewidth=2, label="Train ELBO")
        ax.plot(epochs, val, color=amber, linewidth=2, label="Validation ELBO")

        best_epoch = int(val.argmin())
        ax.scatter([best_epoch], [val[best_epoch]], color=amber, zorder=5, s=40,
                   edgecolor="white", linewidth=1)
        ax.annotate(
            f"best epoch {best_epoch}",
            xy=(best_epoch, val[best_epoch]),
            xytext=(best_epoch - len(train) * 0.28, val[best_epoch] + (max(val) - min(val)) * 0.14),
            fontsize=9, color=grey,
            arrowprops=dict(arrowstyle="-", color=grey, linewidth=0.8),
        )

        ax.set_title(f"{name} training ({len(train)} epochs)", fontsize=13, fontweight="bold", pad=10)
        ax.set_xlabel("Epoch")
        ax.set_ylabel("ELBO (lower is better)")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", color=grey, alpha=0.15, linewidth=0.8)
        ax.legend(frameon=False, fontsize=9, loc="upper right")

    fig.suptitle(
        "Train vs. validation ELBO per epoch — best epoch marks the minimum validation ELBO",
        fontsize=11, color=grey, y=1.03,
    )
    fig.tight_layout()
    out_path = os.path.join(output_dir, "elbo_convergence.png")
    fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor=bg)
    plt.close(fig)
    print(f"  ELBO convergence plot → {out_path}")

def save_timing_file(timings: list[tuple[str, float]], output_dir: str) -> None:
    """Save a TSV of named phase durations plus their sum to `output_dir`/timing.tsv."""
    rows = [(label, seconds, str(timedelta(seconds=round(seconds)))) for label, seconds in timings]
    total_seconds = sum(seconds for _, seconds in timings)
    rows.append(("Total", total_seconds, str(timedelta(seconds=round(total_seconds)))))

    df = pd.DataFrame(rows, columns=["phase", "duration_seconds", "duration_hms"])
    out_path = os.path.join(output_dir, "timing.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"  Timing table → {out_path}")
    print(df.to_string(index=False))

def compute_umap(adata: ad.AnnData, use_rep: str, n_neighbors: int) -> ad.AnnData:
    """Compute k-NN graph and UMAP on obsm[use_rep], storing obsm['X_umap']."""
    print(f"  Computing k-NN graph (use_rep='{use_rep}', n_neighbors={n_neighbors})")
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors)
    print("  Computing UMAP")
    sc.tl.umap(adata)
    return adata

def _validate_merged_obs(merged: ad.AnnData) -> None:
    """
    Check all required obs columns are present, and that source-specific
    columns have non-null values for their own source's cells (else the
    post-concat re-attach failed).

    Raises
    ------
    KeyError
        If a required column is missing.
    ValueError
        If a source-specific column has no non-null values for its source.
    """
    required = [
        "source", "batch", "cell_type_scanvi", "da_subtype", "condition_id",
        "menon_cell_type", "menon_author_cell_type",
        "siletti_cxg_cell_type", "siletti_cxg_author_cell_type", "siletti_cxg_region", "siletti_allen_da_vglut2",
        "siletti_allen_neurotransmitter", "siletti_allen_subcluster",
        "pmdbs_sc_rnaseq.C_scANVI", "pmdbs_sc_rnaseq.cell_type",
    ]
    missing = [c for c in required if c not in merged.obs.columns]
    if missing:
        raise KeyError(
            f"Merged AnnData is missing required obs columns: {missing}\n"
            f"Available columns: {merged.obs.columns.tolist()}"
        )

    # siletti_cxg_author_cell_type and asap_cohort's own columns are excluded
    # from this check: both are legitimately always-NaN for their source
    # (missing upstream data), not a re-attach failure.
    source_cols = {
        "menon": ("menon_cell_type", "menon_author_cell_type"),
        "siletti": (
            "siletti_cxg_cell_type", "siletti_cxg_region", "siletti_allen_da_vglut2",
            "siletti_allen_neurotransmitter", "siletti_allen_subcluster",
        ),
    }
    for source, cols in source_cols.items():
        mask = merged.obs["source"] == source
        for col in cols:
            n_valid = mask.sum() - merged.obs.loc[mask, col].isna().sum()
            if n_valid == 0:
                raise ValueError(
                    f"Column '{col}' has no non-null values for {source} cells. "
                    "This column should have been re-attached after ad.concat."
                )

    sources = sorted(merged.obs["source"].unique())
    print(
        f"  obs columns validated: {len(required)} required columns present "
        f"({len(merged.obs.columns)} total) | sources: {sources}"
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def _require_cache(path: str) -> None:
    """Raise a clear error if a required preQC/postQC cache is missing."""
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Required cache not found: {path}\n"
            "Run downsample_kamath_menon_siletti_asap.py first, with the "
            "same --kamath-non-da-per-type (or --kamath-all-random / --kamath-cell-type-min, "
            "both paired with --downsample-kamath) / "
            "--downsample-menon / --siletti-scale / --downsample-asap-cohort values."
        )


def main(args: argparse.Namespace) -> None:
    """Merge the four sources' cached preQC/postQC AnnDatas and train scVI/scANVI."""
    os.makedirs(args.output_dir, exist_ok=True)
    scvi.settings.seed = 42
    t_start = time.monotonic()

    if args.kamath_all_random and args.kamath_cell_type_min > 0:
        parser.error("--kamath-all-random and --kamath-cell-type-min are mutually exclusive.")
    if args.kamath_all_random and args.downsample_kamath is None:
        parser.error("--downsample-kamath is required when --kamath-all-random is set.")
    if args.kamath_cell_type_min > 0 and args.downsample_kamath is None:
        parser.error("--downsample-kamath is required when --kamath-cell-type-min is set.")
    if args.downsample_kamath is not None and not args.kamath_all_random and args.kamath_cell_type_min == 0:
        parser.error("--downsample-kamath requires --kamath-all-random or --kamath-cell-type-min.")

    print("\n=== Locating cached preQC/postQC intermediate files ===")
    kamath_dir = os.path.join(_REFERENCE_DATA_DIR, "kamath")
    menon_dir = os.path.join(_REFERENCE_DATA_DIR, "menon")
    siletti_dir = os.path.join(_REFERENCE_DATA_DIR, "siletti")
    menon_path = _find_menon_h5ad(menon_dir)

    qc_label = "postQC" if args.use_qc_filtered else "preQC"
    print(f"  --use-qc-filtered={args.use_qc_filtered} -> using {qc_label} pools for all four sources")

    if args.kamath_all_random:
        kamath_anchor = os.path.join(kamath_dir, "intermediate", "kamath_merged_all_rand")
        kamath_tag = str(args.downsample_kamath)
    elif args.kamath_cell_type_min > 0:
        kamath_anchor = os.path.join(kamath_dir, "intermediate", f"kamath_merged_cell_type_min{args.kamath_cell_type_min}")
        kamath_tag = str(args.downsample_kamath)
    else:
        kamath_anchor = os.path.join(kamath_dir, "intermediate", "kamath_merged_da_all_non_da")
        kamath_tag = str(args.kamath_non_da_per_type)
    kamath_preqc_cache, kamath_postqc_cache = _qc_cache_paths(kamath_anchor, kamath_tag)
    kamath_cache = kamath_postqc_cache if args.use_qc_filtered else kamath_preqc_cache
    _require_cache(kamath_cache)
    kamath = sc.read_h5ad(kamath_cache)
    print(f"  Kamath ({qc_label}): {kamath.n_obs:,} cells")

    menon_tag = f"rand_{args.downsample_menon}" if args.downsample_menon > 0 else "all"
    menon_preqc_cache, menon_postqc_cache = _qc_cache_paths(_sn_cache_path(menon_path), menon_tag)
    menon_cache = menon_postqc_cache if args.use_qc_filtered else menon_preqc_cache
    _require_cache(menon_cache)
    menon = sc.read_h5ad(menon_cache)
    print(f"  Menon ({qc_label}): {menon.n_obs:,} cells")
    # Only needed for --gene-selection hvg_u_da_markers/list; cheap to always grab.
    menon_feature_name = menon.var["feature_name"].copy() if "feature_name" in menon.var.columns else None

    siletti_neuron_target, siletti_cell_target = _SILETTI_SCALE_PRESETS[args.siletti_scale]
    siletti_anchor = os.path.join(
        siletti_dir, "intermediate",
        f"midbrain_da_neurons{siletti_neuron_target}_cells{siletti_cell_target}",
    )
    siletti_preqc_cache, siletti_postqc_cache = _qc_cache_paths(siletti_anchor, "")
    siletti_cache = siletti_postqc_cache if args.use_qc_filtered else siletti_preqc_cache
    _require_cache(siletti_cache)
    siletti = sc.read_h5ad(siletti_cache)
    print(f"  Siletti ({qc_label}): {siletti.n_obs:,} cells")

    asap_tag = f"rand_{args.downsample_asap_cohort}" if args.downsample_asap_cohort > 0 else "all"
    asap_preqc_cache, asap_postqc_cache = _qc_cache_paths(_ASAP_H5AD_PATH, asap_tag)
    asap_cache = asap_postqc_cache if args.use_qc_filtered else asap_preqc_cache
    _require_cache(asap_cache)
    asap = sc.read_h5ad(asap_cache)
    print(f"  ASAP cohort ({qc_label}): {asap.n_obs:,} cells")

    # --- Merge ---
    print("\n=== Merging datasets ===")

    # ad.concat(join='inner') drops obs columns not shared by all inputs;
    # save source-specific columns here and re-attach after concat.
    _source_specific = {}
    for col in ("menon_cell_type", "menon_author_cell_type"):
        _source_specific[col] = menon.obs[[col]].copy()
    for col in (
        "siletti_cxg_cell_type", "siletti_cxg_author_cell_type", "siletti_cxg_region",
        "siletti_allen_da_vglut2", "siletti_allen_neurotransmitter", "siletti_allen_subcluster",
    ):
        _source_specific[col] = siletti.obs[[col]].copy()
    for col in ("pmdbs_sc_rnaseq.C_scANVI", "pmdbs_sc_rnaseq.cell_type"):
        _source_specific[col] = asap.obs[[col]].copy()

    sources_to_merge = [kamath, menon, siletti, asap]
    gene_counts = {
        "Kamath": kamath.n_vars,
        "Menon": menon.n_vars,
        "Siletti": siletti.n_vars,
        "ASAP cohort": asap.n_vars,
    }

    merged = ad.concat(
        sources_to_merge,
        join="inner",           # intersection of gene sets
        merge="same",           # keep shared var metadata
    )

    # Boolean columns default to False (not NaN) -- NaN would force object
    # dtype (mixed bool/float), which the H5AD writer can't serialize.
    for col, df in _source_specific.items():
        merged.obs[col] = False if pd.api.types.is_bool_dtype(df[col]) else np.nan
        valid_idx = df.index.intersection(merged.obs_names)
        merged.obs.loc[valid_idx, col] = df.loc[valid_idx, col].values

    gene_summary = "  |  ".join(f"{k}: {v:,}" for k, v in gene_counts.items())
    print(f"  {gene_summary}  |  Shared (inner join): {merged.n_vars:,}")
    print(f"  Merged total: {merged.n_obs:,} cells")
    print(f"  Source breakdown:\n{merged.obs['source'].value_counts().to_string()}")

    # --- Per-run output subfolder ---
    # e.g. output_dir/asap50k_kamath52k_menon50k_siletti51k_<gene_selection>/.
    # Computed here (post-merge, real counts known) so models/TSVs/plot/
    # timing all land alongside the final H5AD.
    _source_short = {
        "asap_cohort": "asap", "kamath": "kamath", "menon": "menon", "siletti": "siletti",
    }
    _counts = merged.obs["source"].value_counts()
    # "list" alone isn't descriptive -- use the gene-list file's stem instead.
    _gene_selection_tag = (
        f"list_{Path(args.infile_gene_selection).stem}"
        if args.gene_selection == "list" else args.gene_selection
    )
    _suffix = "_".join(
        f"{_source_short.get(src, src)}{round(_counts[src] / 1000)}k"
        for src in sorted(_counts.index)
        if src in _counts
    ) + f"_{_gene_selection_tag}"
    run_dir = os.path.join(args.output_dir, _suffix)
    os.makedirs(run_dir, exist_ok=True)
    print(f"  Run output directory: {run_dir}")

    # Rename is safe while the log is still open (POSIX keeps the same inode).
    if getattr(args, "_log_file_path", None) and os.path.exists(args._log_file_path):
        new_log_path = os.path.join(run_dir, os.path.basename(args._log_file_path))
        os.rename(args._log_file_path, new_log_path)
        args._log_file_path = new_log_path
        print(f"  Moved run log → {new_log_path}")

    # --- Validate merged obs ---
    print("\n=== Validating merged obs columns ===")
    _validate_merged_obs(merged)

    # --- Preprocess ---
    print("\n=== Preprocessing ===")
    merged = preprocess_merged(merged, args.n_hvgs)

    # --- Gene selection: which genes scVI/scANVI actually train on ---
    # merged_train feeds scVI/scANVI; merged is what gets saved. Same object
    # for every mode except "list", where training uses exactly
    # --infile-gene-selection's genes but the saved H5AD keeps the wider
    # HVG ∪ DA-marker panel regardless (so markers stay plottable).
    print(f"\n=== Gene selection: {args.gene_selection} ===")
    if args.gene_selection == "gene_intersection":
        print(f"  No restriction — training on all {merged.n_vars:,} shared (inner-join) genes")
        merged_train = merged
    else:
        if args.gene_selection == "hvg_only":
            train_mask = merged.var["highly_variable"].to_numpy()
        elif args.gene_selection == "hvg_u_da_markers":
            if menon_feature_name is None:
                raise KeyError(
                    "--gene-selection hvg_u_da_markers requires Menon's "
                    "var['feature_name'] column to resolve marker gene symbols "
                    "to Ensembl IDs, but it was not found on the loaded Menon cache."
                )
            marker_ensembl_ids = _resolve_marker_ensembl_ids(
                menon_feature_name, _KAMATH_DA_MARKER_SYMBOLS, "Menon"
            )
            in_markers = merged.var_names.isin(marker_ensembl_ids)
            train_mask = merged.var["highly_variable"].to_numpy() | in_markers
            n_markers_added = int((in_markers & ~merged.var["highly_variable"].to_numpy()).sum())
            print(f"  {len(marker_ensembl_ids)}/{len(_KAMATH_DA_MARKER_SYMBOLS)} Kamath DA markers "
                  f"resolved to Ensembl IDs in the merged panel; {n_markers_added} added beyond the HVG set")
        else:  # "list"
            with open(args.infile_gene_selection) as f:
                list_gene_ids = [line.strip() for line in f if line.strip()]
            train_mask = merged.var_names.isin(list_gene_ids)
            n_found = int(train_mask.sum())
            n_missing = len(list_gene_ids) - n_found
            print(f"  --infile-gene-selection {args.infile_gene_selection}: "
                  f"{len(list_gene_ids):,} gene IDs listed, {n_found:,} found in the merged panel"
                  + (f", {n_missing:,} not found (skipped)" if n_missing else ""))
        n_before = merged.n_vars
        merged_train = merged[:, train_mask].copy()
        print(f"  Training gene panel: {n_before:,} → {merged_train.n_vars:,} genes")

    if args.gene_selection == "list":
        if menon_feature_name is None:
            raise KeyError(
                "--gene-selection list still needs Menon's var['feature_name'] "
                "to resolve the DA marker genes for the SAVED (HVG ∪ markers) "
                "panel, even though training uses --infile-gene-selection instead."
            )
        marker_ensembl_ids = _resolve_marker_ensembl_ids(
            menon_feature_name, _KAMATH_DA_MARKER_SYMBOLS, "Menon"
        )
        save_mask = merged.var["highly_variable"].to_numpy() | merged.var_names.isin(marker_ensembl_ids)
        n_before = merged.n_vars
        merged = merged[:, save_mask].copy()
        print(f"  Save panel (HVG ∪ DA markers, independent of the training list): "
              f"{n_before:,} → {merged.n_vars:,} genes")
    elif args.gene_selection != "gene_intersection":
        merged = merged_train

    # --- scVI ---
    scvi_path = os.path.join(run_dir, "scvi_model")
    scanvi_path = os.path.join(run_dir, "scanvi_model")

    if args.use_cached_models:
        for path in (scvi_path, scanvi_path):
            if not os.path.exists(path):
                raise FileNotFoundError(
                    f"--use-cached-models requested but model directory not found: {path}"
                )
        print("\n=== Loading cached scVI model ===")
        model_scvi = scvi.model.SCVI.load(scvi_path, adata=merged_train)
        merged_train.obsm["X_scVI"] = model_scvi.get_latent_representation()
        print(f"  scVI latent stored in obsm['X_scVI'] — shape: {merged_train.obsm['X_scVI'].shape}")
        t_after_scvi = time.monotonic()

        print("\n=== Loading cached scANVI model ===")
        model_scanvi = scvi.model.SCANVI.load(scanvi_path, adata=merged_train)
        merged_train.obsm["X_scANVI"] = model_scanvi.get_latent_representation()
        merged_train.obs["cell_type_pred"] = model_scanvi.predict()
        print(f"  scANVI latent stored in obsm['X_scANVI']")
        print(f"  Label predictions stored in obs['cell_type_pred']")
        t_after_scanvi = time.monotonic()
    else:
        print("\n=== scVI integration ===")
        merged_train, model_scvi = run_scvi(
            adata=merged_train,
            n_latent=args.n_latent,
            n_layers=args.n_layers,
            max_epochs=args.max_epochs_scvi,
            batch_size=args.batch_size,
            use_gpu=args.use_gpu,
        )
        print("  Saving scVI model")
        model_scvi.save(scvi_path, overwrite=True)
        t_after_scvi = time.monotonic()

        print("\n=== scANVI label transfer ===")
        merged_train, model_scanvi = run_scanvi(
            adata=merged_train,
            scvi_model=model_scvi,
            labels_key="cell_type_scanvi",
            unlabeled_category="Unknown",
            max_epochs=args.max_epochs_scanvi,
            batch_size=args.batch_size,
            use_gpu=args.use_gpu,
        )
        print("  Saving scANVI model")
        model_scanvi.save(scanvi_path, overwrite=True)
        t_after_scanvi = time.monotonic()

    # Copy training results onto the (possibly wider) save panel. No-op
    # when merged_train IS merged (every mode except "list").
    if merged_train is not merged:
        merged.obsm["X_scVI"] = merged_train.obsm["X_scVI"]
        merged.obsm["X_scANVI"] = merged_train.obsm["X_scANVI"]
        merged.obs["cell_type_pred"] = merged_train.obs["cell_type_pred"]

    print("\n=== Saving ELBO convergence plot ===")
    save_elbo_plot(model_scvi, model_scanvi, run_dir)

    # --- DA subtype transfer ---
    print("\n=== Transferring DA subtypes to Menon cells ===")
    merged = transfer_da_subtypes(merged)

    # --- UMAP ---
    if args.compute_umap:
        print("\n=== Computing UMAP (on scANVI latent) ===")
        merged = compute_umap(merged, use_rep="X_scANVI", n_neighbors=args.n_neighbors)

    # --- Summary ---
    print("\n=== Integration summary ===")
    print("Predicted cell type distribution:")
    print(merged.obs["cell_type_pred"].value_counts().to_string())
    da_counts = merged.obs["da_subtype"].dropna().value_counts()
    if not da_counts.empty:
        print("\nDA subtype distribution (Kamath original + Menon/ASAP transferred):")
        print(da_counts.to_string())

    # --- Save ---
    out_path = os.path.join(run_dir, f"integrated_{_suffix}.h5ad")
    print(f"\nSaving integrated AnnData → {out_path}")
    # Drop heavy HVG stats columns to reduce file size.
    for col in ("means", "variances", "dispersions", "variances_norm", "dispersions_norm"):
        if col in merged.var.columns:
            merged.var.drop(columns=[col], inplace=True)
    merged.write_h5ad(out_path)

    # --- Timing ---
    t_end = time.monotonic()
    print("\n=== Timing ===")
    save_timing_file(
        [
            ("Data loading + preprocessing -> scVI training", t_after_scvi - t_start),
            ("scVI saved -> scANVI training", t_after_scanvi - t_after_scvi),
            ("scANVI saved -> label transfer, UMAP, save", t_end - t_after_scanvi),
        ],
        run_dir,
    )

    print("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--use-qc-filtered",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use the QC-filtered (postQC) cache for each of the four sources\n"
             "(built by downsample_kamath_menon_siletti_asap.py). Pass\n"
             "--no-use-qc-filtered to use the preQC (downsampled but\n"
             "unfiltered) cache instead, for a controlled A/B comparison of\n"
             "QC's actual effect on integration — both were downsampled\n"
             "identically (same seed=42 draw) before QC ever ran, so they\n"
             "differ ONLY in whether QC was applied.\n"
             "[Default: True]",
    )
    parser.add_argument(
        "--use-cached-models",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Load scVI and scANVI models from --output-dir instead of training.\n"
             "Requires a prior run to have saved scvi_model/ and scanvi_model/ there.\n"
             "[Default: False]",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for output files (integrated.h5ad, models, agreement TSVs).",
    )
    parser.add_argument(
        "--kamath-non-da-per-type",
        type=int,
        default=_DEFAULT_KAMATH_NON_DA_PER_TYPE,
        help="Must match the value passed to downsample_kamath_menon_siletti_asap.py\n"
             "-- used only to reconstruct the matching Kamath cache filename,\n"
             "not to re-downsample. Ignored when --kamath-all-random is set.\n"
             f"[Default: {_DEFAULT_KAMATH_NON_DA_PER_TYPE}]",
    )
    parser.add_argument(
        "--kamath-all-random",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Load the all-random Kamath pool (kamath_merged_all_rand_N) produced\n"
             "by downsample_kamath_menon_siletti_asap.py --kamath-all-random.\n"
             "Requires --downsample-kamath. Mutually exclusive with --kamath-cell-type-min.\n"
             "[Default: False]",
    )
    parser.add_argument(
        "--kamath-cell-type-min",
        type=int,
        default=0,
        help="Load the per-cell-type-min Kamath pool (kamath_merged_cell_type_minN_M)\n"
             "produced by downsample_kamath_menon_siletti_asap.py --kamath-cell-type-min.\n"
             "Must match the value used there. Requires --downsample-kamath.\n"
             "Mutually exclusive with --kamath-all-random. 0 = disabled.\n"
             "[Default: 0]",
    )
    parser.add_argument(
        "--downsample-kamath",
        type=int,
        default=None,
        help="Must match --downsample-kamath passed to downsample_kamath_menon_siletti_asap.py\n"
             "-- used only to reconstruct the matching Kamath cache filename,\n"
             "not to re-downsample. Required when --kamath-all-random or\n"
             "--kamath-cell-type-min is set.\n"
             "[Default: None]",
    )
    parser.add_argument(
        "--downsample-menon",
        type=int,
        default=_DEFAULT_DOWNSAMPLE_MENON,
        help="Must match the value passed to downsample_kamath_menon_siletti_asap.py\n"
             "-- used only to reconstruct the matching Menon cache filename,\n"
             "not to re-downsample.\n"
             f"[Default: {_DEFAULT_DOWNSAMPLE_MENON}]",
    )
    parser.add_argument(
        "--siletti-scale",
        type=int,
        choices=sorted(_SILETTI_SCALE_PRESETS),
        default=_SILETTI_DEFAULT_SCALE,
        help="Must match the value passed to downsample_kamath_menon_siletti_asap.py\n"
             "-- used only to reconstruct the matching Siletti cache filename,\n"
             "not to re-downsample.\n"
             f"[Default: {_SILETTI_DEFAULT_SCALE}]",
    )
    parser.add_argument(
        "--downsample-asap-cohort",
        type=int,
        default=_DEFAULT_DOWNSAMPLE_ASAP_COHORT,
        help="Must match the value passed to downsample_kamath_menon_siletti_asap.py\n"
             "-- used only to reconstruct the matching ASAP cache filename,\n"
             "not to re-downsample.\n"
             f"[Default: {_DEFAULT_DOWNSAMPLE_ASAP_COHORT}]",
    )
    parser.add_argument(
        "--n-hvgs",
        type=int,
        default=3000,
        help="Number of highly variable genes for scVI training. Matches\n"
             "sc-rnaseq-wf's main.wdl n_top_genes default (its own\n"
             "cohort_analysis.wdl standalone default/docs say 8000, but the\n"
             "real production entrypoint main.wdl -- which is what actually\n"
             "produced asap-cohort.final.h5ad/mmc.h5ad -- defaults to 3000).\n"
             "[Default: 3000]",
    )
    parser.add_argument(
        "--gene-selection",
        choices=["gene_intersection", "hvg_only", "hvg_u_da_markers", "list"],
        default="gene_intersection",
        help="Which genes scVI/scANVI actually train on (var['highly_variable']\n"
             "is always computed regardless, for reference):\n"
             "  gene_intersection - all genes shared by the 4-way inner join\n"
             "    (no HVG restriction applied). Matches prior runs' behavior.\n"
             "  hvg_only          - subset to just the --n-hvgs HVGs, matching\n"
             "    sc-rnaseq-wf's own cohort_analysis process task, which\n"
             "    physically discards non-HVG genes the same way.\n"
             "  hvg_u_da_markers  - union of the HVGs and the 17 Kamath et al.\n"
             "    (Nat Neurosci 2022) DA marker genes (_KAMATH_DA_MARKER_SYMBOLS),\n"
             "    so canonical DA identity/subtype markers are never dropped\n"
             "    from training even if HVG selection alone would miss them.\n"
             "  list              - subset to exactly the Ensembl gene IDs listed\n"
             "    in --infile-gene-selection (required for this mode). For\n"
             "    controlled A/B comparisons against a specific gene panel\n"
             "    (e.g. a prior run's HVG set with markers added/removed).\n"
             "Appended to the run's output subfolder name.\n"
             "[Default: gene_intersection]",
    )
    parser.add_argument(
        "--infile-gene-selection",
        default=None,
        help="Path to a text file with one Ensembl gene ID per row. Required\n"
             "when --gene-selection list is used; ignored otherwise.",
    )
    parser.add_argument(
        "--n-latent",
        type=int,
        default=30,
        help="Dimensionality of the scVI / scANVI latent space.\n"
             "[Default: 30]",
    )
    parser.add_argument(
        "--n-layers",
        type=int,
        default=2,
        help="Number of hidden layers in the scVI encoder/decoder.\n"
             "[Default: 2]",
    )
    parser.add_argument(
        "--max-epochs-scvi",
        type=int,
        default=400,
        help="Maximum training epochs for scVI. Early stopping is enabled\n"
             "(patience=15), so training usually ends before this limit.\n"
             "[Default: 400]",
    )
    parser.add_argument(
        "--max-epochs-scanvi",
        type=int,
        default=200,
        help="Maximum training epochs for scANVI (patience=10).\n"
             "[Default: 200]",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=1024,
        help="Mini-batch size for scVI and scANVI. Reduce to 512 if GPU OOM.\n"
             "[Default: 1024]",
    )
    parser.add_argument(
        "--n-neighbors",
        type=int,
        default=30,
        help="Number of neighbours for the k-NN graph used in UMAP computation.\n"
             "[Default: 30]",
    )
    parser.add_argument(
        "--compute-umap",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Compute UMAP on the scANVI latent representation.\n"
             "Use --no-compute-umap to skip (saves time on very large datasets).\n"
             "[Default: True]",
    )
    parser.add_argument(
        "--use-gpu",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use GPU acceleration for scVI and scANVI.\n"
             "Use --no-use-gpu to run on CPU only (very slow for >1M cells).\n"
             "[Default: True]",
    )
    args = parser.parse_args()

    if args.gene_selection == "list" and not args.infile_gene_selection:
        parser.error("--gene-selection list requires --infile-gene-selection")
    if args.gene_selection != "list" and args.infile_gene_selection:
        print(f"[INFO] --infile-gene-selection is ignored (--gene-selection is "
              f"'{args.gene_selection}', not 'list')")
    if args.infile_gene_selection and not os.path.exists(args.infile_gene_selection):
        parser.error(f"--infile-gene-selection file not found: {args.infile_gene_selection}")

    args.output_dir = os.path.abspath(args.output_dir)

    repo_root = next(p for p in Path(__file__).resolve().parents if (p / ".git").exists())
    if Path.cwd() != repo_root:
        print(f"[INFO] Changing working directory to repo root: {repo_root}")
        os.chdir(repo_root)

    log_run_command(args)
    main(args)
