#!/usr/bin/env python3
"""
Kamath-free variant of 04_integrate_kamath_menon_siletti_asap.py. Uses same logic to
integrate Menon + Siletti + ASAP only, via scVI + scANVI, but without Kamath data.

Since Kamath is the dataset providing the ground truth labels for scANVI, in
04_integrate_kamath_menon_siletti_asap.py, here instead we use the Allen cell type annotations 
of Siletti's dataset (siletti_allen_cluster). This scripts derives those annotations
from siletti_allen_subcluster via the Allen taxonomy's subcluster->cluster mapping.

Notes:
- Cluster Splat_395 is exactly and only the DA VGLUT2 population.
  All 805 DA VGLUT2 cells map to it.
- KAMATH_DA_MARKER_SYMBOLS still applies for the --gene-selection hvg_u_da_markers/list
  parameter, to use the exact same list of HVGs + DA markers as in the Kamath-inclusive run.

Outputs land in {output_dir}/{suffix}_no_kamath/.

See reference_building/da_neurons_identification_analysis/README.md and 04_integrate_kamath_menon_siletti_asap.py help for more details.

"""

import os
import time
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scvi
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "utils"))
from _reference_building_utils import (
    _REFERENCE_DATA_DIR,
    _SILETTI_SCALE_PRESETS,
    _SILETTI_DEFAULT_SCALE,
    _ASAP_H5AD_PATH,
    _DEFAULT_DOWNSAMPLE_MENON,
    _DEFAULT_DOWNSAMPLE_ASAP_COHORT,
    _qc_cache_paths,
    _sn_cache_path,
    _find_menon_h5ad,
    log_run_command,
)
from integrate_kamath_menon_siletti_asap import (
    preprocess_merged,
    run_scvi,
    run_scanvi,
    save_elbo_plot,
    save_timing_file,
    compute_umap,
    _validate_merged_obs,
    _require_cache,
)

_repo_root = next(p for p in Path(__file__).resolve().parents if (p / ".git").exists())
sys.path.insert(0, str(_repo_root))
from utils.h5ad_utils import resolve_marker_ensembl_ids
from utils.genomics_constants import KAMATH_DA_MARKER_SYMBOLS

# Allen taxonomy files.
# See reference_building/da_neurons_identification_analysis/reference_data/allen-brain-cell-atlas/README.md for details
_ALLEN_MEMBERSHIP_PATH = os.path.join(
    _REFERENCE_DATA_DIR, "allen-brain-cell-atlas", "metadata",
    "WHB-taxonomy", "20240330", "cluster_to_cluster_annotation_membership.csv.gz",
)


def _load_siletti_subcluster_to_cluster() -> pd.Series:
    """
    Map Allen subcluster -> parent cluster name (verified stable many-to-one).

    Returns
    -------
    pd.Series
        Indexed by subcluster name, values are the parent cluster name.
    """
    membership = pd.read_csv(_ALLEN_MEMBERSHIP_PATH, keep_default_na=False)
    sub = membership.loc[
        membership["cluster_annotation_term_set_name"] == "subcluster",
        ["cluster_alias", "cluster_annotation_term_name"],
    ].rename(columns={"cluster_annotation_term_name": "subcluster_name"})
    clu = membership.loc[
        membership["cluster_annotation_term_set_name"] == "cluster",
        ["cluster_alias", "cluster_annotation_term_name"],
    ].rename(columns={"cluster_annotation_term_name": "cluster_name"})
    merged = sub.merge(clu, on="cluster_alias").drop_duplicates("subcluster_name")
    return merged.set_index("subcluster_name")["cluster_name"]


def main(args: argparse.Namespace) -> None:
    """Merge Menon/Siletti/ASAP's cached preQC/postQC AnnDatas and train scVI/scANVI."""
    os.makedirs(args.output_dir, exist_ok=True)
    scvi.settings.seed = 42
    t_start = time.monotonic()

    print("\n=== Locating cached preQC/postQC intermediate files ===")
    menon_dir = os.path.join(_REFERENCE_DATA_DIR, "menon")
    siletti_dir = os.path.join(_REFERENCE_DATA_DIR, "siletti")
    menon_path = _find_menon_h5ad(menon_dir)

    qc_label = "postQC" if args.use_qc_filtered else "preQC"
    print(f"  --use-qc-filtered={args.use_qc_filtered} -> using {qc_label} pools for all three sources")

    menon_tag = f"rand_{args.downsample_menon}" if args.downsample_menon > 0 else "all"
    menon_preqc_cache, menon_postqc_cache = _qc_cache_paths(_sn_cache_path(menon_path), menon_tag)
    menon_cache = menon_postqc_cache if args.use_qc_filtered else menon_preqc_cache
    _require_cache(menon_cache)
    menon = sc.read_h5ad(menon_cache)
    print(f"  Menon ({qc_label}): {menon.n_obs:,} cells")
    # Only needed for --gene-selection hvg_u_da_markers/list -- cheap to always grab.
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

    sources_to_merge = [menon, siletti, asap]
    gene_counts = {"Menon": menon.n_vars, "Siletti": siletti.n_vars, "ASAP cohort": asap.n_vars}

    merged = ad.concat(
        sources_to_merge,
        join="inner",           # intersection of gene sets
        merge="same",           # keep shared var metadata
    )

    for col, df in _source_specific.items():
        merged.obs[col] = False if pd.api.types.is_bool_dtype(df[col]) else np.nan
        valid_idx = df.index.intersection(merged.obs_names)
        merged.obs.loc[valid_idx, col] = df.loc[valid_idx, col].values

    gene_summary = "  |  ".join(f"{k}: {v:,}" for k, v in gene_counts.items())
    print(f"  {gene_summary}  |  Shared (inner join): {merged.n_vars:,}")
    print(f"  Merged total: {merged.n_obs:,} cells")
    print(f"  Source breakdown:\n{merged.obs['source'].value_counts().to_string()}")

    # --- scANVI training labels: Siletti's native Allen cluster label, not Kamath ---
    # siletti_allen_cluster isn't cached upstream (only siletti_allen_subcluster
    # is) -- derived here via the same subcluster->cluster mapping the notebook
    # uses. 'Splat_395' is verified (empirically) to be exactly the DA VGLUT2
    # population -- see module docstring.
    print("\n=== Building scANVI training labels from Siletti's native Allen cluster ===")
    subcluster_to_cluster = _load_siletti_subcluster_to_cluster()
    merged.obs["cell_type_scanvi"] = merged.obs["cell_type_scanvi"].astype(str)
    siletti_mask = merged.obs["source"] == "siletti"
    siletti_cluster = merged.obs.loc[siletti_mask, "siletti_allen_subcluster"].map(subcluster_to_cluster)
    n_unmapped = siletti_cluster.isna().sum()
    if n_unmapped:
        print(f"  Warning: {n_unmapped:,}/{siletti_mask.sum():,} Siletti cells had no "
              f"subcluster->cluster mapping (unexpected -- check Allen taxonomy version)")
    merged.obs.loc[siletti_mask, "cell_type_scanvi"] = siletti_cluster
    n_classes = merged.obs.loc[siletti_mask, "cell_type_scanvi"].nunique()
    n_da = (merged.obs.loc[siletti_mask, "cell_type_scanvi"] == "Splat_395").sum()
    print(f"  {n_classes} distinct Siletti clusters used as training labels "
          f"('Splat_395' = DA VGLUT2, {n_da:,} cells)")
    print(merged.obs.loc[siletti_mask, "cell_type_scanvi"].value_counts().head(15).to_string())

    # --- Per-run output subfolder ---
    _source_short = {"asap_cohort": "asap", "menon": "menon", "siletti": "siletti"}
    _counts = merged.obs["source"].value_counts()
    _gene_selection_tag = (
        f"list_{Path(args.infile_gene_selection).stem}"
        if args.gene_selection == "list" else args.gene_selection
    )
    _suffix = "_".join(
        f"{_source_short.get(src, src)}{round(_counts[src] / 1000)}k"
        for src in sorted(_counts.index)
        if src in _counts
    ) + f"_{_gene_selection_tag}_no_kamath"
    run_dir = os.path.join(args.output_dir, _suffix)
    os.makedirs(run_dir, exist_ok=True)
    print(f"  Run output directory: {run_dir}")

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

    # --- Gene selection: which genes scVI/scANVI actually train on (see
    # 04_integrate_kamath_menon_siletti_asap.py for the full rationale) ---
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
            marker_ensembl_ids = resolve_marker_ensembl_ids(
                menon_feature_name, KAMATH_DA_MARKER_SYMBOLS, "Menon"
            )
            in_markers = merged.var_names.isin(marker_ensembl_ids)
            train_mask = merged.var["highly_variable"].to_numpy() | in_markers
            n_markers_added = int((in_markers & ~merged.var["highly_variable"].to_numpy()).sum())
            print(f"  {len(marker_ensembl_ids)}/{len(KAMATH_DA_MARKER_SYMBOLS)} DA markers "
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
        marker_ensembl_ids = resolve_marker_ensembl_ids(
            menon_feature_name, KAMATH_DA_MARKER_SYMBOLS, "Menon"
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

        print("\n=== scANVI label transfer (DA vs non-DA, from Siletti) ===")
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

    if merged_train is not merged:
        merged.obsm["X_scVI"] = merged_train.obsm["X_scVI"]
        merged.obsm["X_scANVI"] = merged_train.obsm["X_scANVI"]
        merged.obs["cell_type_pred"] = merged_train.obs["cell_type_pred"]

    print("\n=== Saving ELBO convergence plot ===")
    save_elbo_plot(model_scvi, model_scanvi, run_dir)

    # No transfer_da_subtypes() -- Kamath-specific (see module docstring).
    # Use the notebook's label_agreement_widget() for post-hoc validation.

    # --- UMAP ---
    if args.compute_umap:
        print("\n=== Computing UMAP (on scANVI latent) ===")
        merged = compute_umap(merged, use_rep="X_scANVI", n_neighbors=args.n_neighbors)

    # --- Summary ---
    print("\n=== Integration summary ===")
    print("Predicted cell type distribution (DA vs non-DA, transferred from Siletti):")
    print(merged.obs["cell_type_pred"].value_counts().to_string())
    print("\nBy source:")
    print(merged.obs.groupby("source")["cell_type_pred"].value_counts().to_string())

    # --- Save ---
    out_path = os.path.join(run_dir, f"integrated_{_suffix}.h5ad")
    print(f"\nSaving integrated AnnData → {out_path}")
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
            ("scANVI saved -> UMAP, save", t_end - t_after_scanvi),
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
        help="Use the QC-filtered (postQC) cache for each of the three sources\n"
             "(built by 03_downsample_kamath_menon_siletti_asap.py). Pass\n"
             "--no-use-qc-filtered to use the preQC (downsampled but\n"
             "unfiltered) cache instead.\n"
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
        help="Directory for output files (integrated.h5ad, models, ELBO plot).",
    )
    parser.add_argument(
        "--downsample-menon",
        type=int,
        default=_DEFAULT_DOWNSAMPLE_MENON,
        help="Must match the value passed to 03_downsample_kamath_menon_siletti_asap.py\n"
             "-- used only to reconstruct the matching Menon cache filename,\n"
             "not to re-downsample.\n"
             f"[Default: {_DEFAULT_DOWNSAMPLE_MENON}]",
    )
    parser.add_argument(
        "--siletti-scale",
        type=int,
        choices=sorted(_SILETTI_SCALE_PRESETS),
        default=_SILETTI_DEFAULT_SCALE,
        help="Must match the value passed to 03_downsample_kamath_menon_siletti_asap.py\n"
             "-- used only to reconstruct the matching Siletti cache filename,\n"
             "not to re-downsample.\n"
             f"[Default: {_SILETTI_DEFAULT_SCALE}]",
    )
    parser.add_argument(
        "--downsample-asap-cohort",
        type=int,
        default=_DEFAULT_DOWNSAMPLE_ASAP_COHORT,
        help="Must match the value passed to 03_downsample_kamath_menon_siletti_asap.py\n"
             "-- used only to reconstruct the matching ASAP cache filename,\n"
             "not to re-downsample.\n"
             f"[Default: {_DEFAULT_DOWNSAMPLE_ASAP_COHORT}]",
    )
    parser.add_argument(
        "--n-hvgs",
        type=int,
        default=3000,
        help="Number of highly variable genes for scVI training.\n"
             "[Default: 3000]",
    )
    parser.add_argument(
        "--gene-selection",
        choices=["gene_intersection", "hvg_only", "hvg_u_da_markers", "list"],
        default="gene_intersection",
        help="Which genes scVI/scANVI actually train on -- same semantics as\n"
             "04_integrate_kamath_menon_siletti_asap.py's --gene-selection.\n"
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
             "[Default: True]",
    )
    parser.add_argument(
        "--use-gpu",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use GPU acceleration for scVI and scANVI.\n"
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
