#!/usr/bin/env python3
"""
Build a multi-source 'reference' H5AD from an 'integrated' H5AD (e.g., including Siletti, Kamath, ASAP cohort),
which has only HVG's, and their corresponding all-gene H5ADs.

The resulting 'reference' will have all-genes intersecting across sources and cells (Ensembl IDs).

Usage (dry-run, default)
------------------------
python3 scripts/prepare_deconvolution_reference_external.py \\
    --integrated-h5ad /path/to/integrated.h5ad \\
    --source-preqc-h5ads \\
        asap_cohort=/path/to/asap_cohort_preQC.h5ad \\
        kamath=/path/to/kamath_preQC.h5ad \\
        siletti=/path/to/siletti_preQC.h5ad \\
    --output-h5ad /path/to/deconv_reference.h5ad
"""

import argparse
import os
import sys
import numpy as np
import anndata as ad
import h5py
import pandas as pd
from pathlib import Path
from scipy.sparse import vstack

repo_root = next(p for p in Path(__file__).resolve().parents if (p / ".git").exists())
sys.path.insert(0, str(repo_root))

from utils.h5ad_utils import extract_raw_counts, resolve_marker_ensembl_ids
from utils.logging_extra import log_run_command
from utils.genomics_constants import KAMATH_DA_MARKER_SYMBOLS


# obs columns to rename in the output reference H5AD (old → new)
_RENAME_OBS_COLS = {
    "batch": "batch_id",
}

_BATCH_COL = "batch_id"  # obs column used for per-batch cell-count filtering


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _obs_index_key(h5ad_path: str) -> str:
    """Return the obs index column name stored in the H5AD file attrs."""
    with h5py.File(h5ad_path, "r") as f:
        return f["obs"].attrs.get("_index", "index")


def _read_barcodes(h5ad_path: str) -> list:
    """Read obs index values without loading the full H5AD."""
    key = _obs_index_key(h5ad_path)
    with h5py.File(h5ad_path, "r") as f:
        raw = f["obs"][key][:]
    return [v.decode() if isinstance(v, bytes) else v for v in raw]


def _read_var_names(h5ad_path: str) -> list:
    """Read var index values without loading the full H5AD."""
    with h5py.File(h5ad_path, "r") as f:
        key = f["var"].attrs.get("_index", "index")
        raw = f["var"][key][:]
    return [v.decode() if isinstance(v, bytes) else v for v in raw]


def _sample_per_source(
    integrated: ad.AnnData,
    source_col: str,
    subset_sources: list,
    n_cells: int,
    seed: int,
) -> dict:
    """
    Sample `n_cells` barcodes per source from the integrated obs table.

    Uses contiguous-block sampling within each source group to minimise
    random HDF5 seeks when the caller later reads preQC rows.

    Parameters
    ----------
    integrated : ad.AnnData
        Backed integrated H5AD; only .obs is accessed.
    source_col : str
        Column in .obs that identifies the source dataset.
    subset_sources : list
        Sources to include (must be values in obs[source_col]).
    n_cells : int
        Target cells per source.
    seed : int
        Random seed for block-start selection.

    Returns
    -------
    dict
        source → list of barcode strings
    """
    rng = np.random.default_rng(seed)
    result = {}
    obs = integrated.obs

    missing = [s for s in subset_sources if s not in obs[source_col].values]
    if missing:
        raise ValueError(
            f"Source(s) not found in obs['{source_col}']: {missing}. "
            f"Available: {sorted(obs[source_col].unique())}"
        )

    for src in subset_sources:
        mask = obs[source_col] == src
        barcodes = obs.index[mask].tolist()
        n = len(barcodes)
        if n_cells is None or n_cells >= n:
            print(f"  [{src}] {n:,} cells — taking all")
            result[src] = barcodes
        else:
            start = int(rng.integers(0, n - n_cells + 1))
            result[src] = barcodes[start:start + n_cells]
            print(f"  [{src}] {n:,} total → sampled {n_cells:,} (block offset {start:,})")

    return result


def _load_source_counts(
    preqc_h5ad: str,
    target_barcodes: list,
    source_name: str,
) -> tuple:
    """
    Load raw counts for `target_barcodes` from a preQC H5AD.

    Parameters
    ----------
    preqc_h5ad : str
        Path to the preQC H5AD with all-gene raw counts.
    target_barcodes : list
        Barcodes to extract (must all be present in the preQC obs index).
    source_name : str
        Display label for messages.

    Returns
    -------
    tuple
        (csr_matrix int32, list of var_names, list of matched barcodes,
        var['feature_name'] pandas.Series or None if not present)
        The matrix rows match `matched_barcodes` in the same order.
    """
    all_barcodes = _read_barcodes(preqc_h5ad)
    bc_to_idx = {bc: i for i, bc in enumerate(all_barcodes)}

    target_set = set(target_barcodes)
    missing = target_set - bc_to_idx.keys()
    if missing:
        examples = sorted(missing)[:5]
        raise ValueError(
            f"[{source_name}] {len(missing):,} barcode(s) from integrated not found "
            f"in preQC obs (e.g. {examples}). Check --integrated-h5ad and "
            "--source-preqc-h5ads point to the same cohort."
        )

    sorted_idx = sorted(bc_to_idx[bc] for bc in target_barcodes)
    matched_barcodes = [all_barcodes[i] for i in sorted_idx]

    print(f"  [{source_name}] loading {len(sorted_idx):,} cells from {preqc_h5ad}")
    preqc = ad.read_h5ad(preqc_h5ad, backed="r")
    subset = preqc[sorted_idx].to_memory()
    var_names = list(subset.var_names)
    feature_name = subset.var["feature_name"] if "feature_name" in subset.var.columns else None
    counts = extract_raw_counts(subset, source_name=source_name)

    return counts, var_names, matched_barcodes, feature_name


def _resolve_da_marker_ensembl_ids(
    source_feature_names: dict,
    subset_sources: list,
) -> list:
    """
    Resolve Kamath's 17 DA marker gene symbols to Ensembl IDs via the first
    source that has a var['feature_name'] column.

    Delegates to `utils.h5ad_utils.resolve_marker_ensembl_ids`; only the "try
    each source in order" part is new here, since callers elsewhere resolve
    against a single, fixed source rather than an arbitrary source list.

    Parameters
    ----------
    source_feature_names : dict
        source_name → var['feature_name'] pandas.Series (indexed by Ensembl ID), or None.
    subset_sources : list
        Source names to try, in order.

    Returns
    -------
    list
        Ensembl IDs resolved for `KAMATH_DA_MARKER_SYMBOLS`. Unresolved symbols are
        warned about, not raised on.
    """
    for src in subset_sources:
        feature_name = source_feature_names.get(src)
        if feature_name is None:
            continue
        return resolve_marker_ensembl_ids(feature_name, KAMATH_DA_MARKER_SYMBOLS, src)

    raise ValueError(
        "--filter-genes cell2location_u_da requires at least one source with a "
        f"var['feature_name'] column to resolve Kamath DA marker symbols to Ensembl IDs. "
        f"None of {subset_sources} had one."
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def prepare_reference(
    integrated_h5ad: str,
    preqc_map: dict,
    subset_sources: list,
    n_cells: int,
    source_col: str,
    obs_cols: list,
    rename_obs_cols: dict,
    min_cells_per_batch: int,
    cell_type_key: str,
    output_h5ad: str,
    seed: int,
    promote: bool,
    filter_genes_mode: str,
) -> None:
    """
    Build a multi-source deconvolution reference H5AD.

    Parameters
    ----------
    integrated_h5ad : str
        Path to the scVI/scANVI integrated H5AD (backed; only obs accessed).
    preqc_map : dict
        source_name → path to preQC H5AD with all-gene raw counts.
    subset_sources : list
        Source names to include (must be keys in `preqc_map` and values in
        obs[source_col] of the integrated H5AD).
    n_cells : int
        Target cells per source.
    source_col : str
        obs column in the integrated H5AD that holds source identifiers.
    obs_cols : list
        obs columns to carry into the output H5AD.
    rename_obs_cols : dict
        old_name → new_name renames applied to obs before writing.
    min_cells_per_batch : int
        Drop cells whose `_BATCH_COL` group has fewer than this many cells.
        0 disables filtering.
    cell_type_key : str
        obs column with scVI/scANVI cell-type labels (must be in obs_cols).
    output_h5ad : str
        Destination path.
    seed : int
        Random seed for sampling.
    promote : bool
        If False, dry-run only.
    filter_genes_mode : str
        Gene set to write beyond the raw cross-source intersection:
        "intersection" : no further filtering.
        "cell2location_default" : apply cell2location.utils.filtering.filter_genes defaults
         (cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12).
        "cell2location_u_da" : cell2location_default, unioned with Kamath's 17 DA marker genes.

    Returns
    -------
    None
    """
    print(f"[INFO] Loading integrated obs from: {integrated_h5ad}")
    integrated = ad.read_h5ad(integrated_h5ad, backed="r")
    print(f"  {integrated.n_obs:,} cells total, sources: {sorted(integrated.obs[source_col].unique())}")

    missing_cols = [c for c in obs_cols if c not in integrated.obs.columns]
    if missing_cols:
        raise ValueError(
            f"obs columns not found in integrated H5AD: {missing_cols}. "
            f"Available: {list(integrated.obs.columns)}"
        )

    n_cells_label = "all" if n_cells is None else f"{n_cells:,}"
    print(f"[INFO] Sampling {n_cells_label} cells per source from: {subset_sources}")
    source_barcodes = _sample_per_source(integrated, source_col, subset_sources, n_cells, seed)

    # Collect per-source metadata (labels) from integrated obs
    source_obs = {}
    for src, barcodes in source_barcodes.items():
        source_obs[src] = integrated.obs.loc[barcodes, obs_cols].copy()

    print("[INFO] Loading preQC counts per source ...")
    source_counts = {}
    source_vars   = {}
    source_barcodes_matched = {}
    source_feature_names = {}
    for src in subset_sources:
        counts, var_names, matched, feature_name = _load_source_counts(
            preqc_map[src], source_barcodes[src], src
        )
        source_counts[src] = counts
        source_vars[src]   = var_names
        source_barcodes_matched[src] = matched
        source_feature_names[src] = feature_name
        print(f"  [{src}] {counts.shape[0]:,} cells × {counts.shape[1]:,} genes")

    print("[INFO] Computing gene intersection ...")
    gene_sets = [set(v) for v in source_vars.values()]
    shared_genes = sorted(gene_sets[0].intersection(*gene_sets[1:]))
    print(f"  {len(shared_genes):,} genes in common across {len(subset_sources)} sources")

    if not promote:
        n_total = sum(len(bc) for bc in source_barcodes_matched.values())
        print(f"[DRY-RUN] Would write: {output_h5ad}")
        print(f"          Shape: {n_total:,} cells × {len(shared_genes):,} genes (raw intersection)")
        if filter_genes_mode != "intersection":
            print(f"          --filter-genes {filter_genes_mode} will reduce the gene count further; "
                  "re-run with --promote to see the final number.")
        print("          Re-run with --promote to write the file.")
        return

    print("[INFO] Subsetting each source to shared gene set ...")
    gene_index = {g: i for i, g in enumerate(shared_genes)}
    all_matrices = []
    all_obs_dfs  = []
    all_barcodes = []

    for src in subset_sources:
        src_gene_idx = np.array([gene_index[g] for g in shared_genes
                                  if g in gene_index], dtype=np.intp)
        var_to_pos = {g: i for i, g in enumerate(source_vars[src])}
        col_idx = np.array([var_to_pos[g] for g in shared_genes], dtype=np.intp)
        mat = source_counts[src][:, col_idx]
        all_matrices.append(mat)
        all_obs_dfs.append(source_obs[src].loc[source_barcodes_matched[src]])
        all_barcodes.extend(source_barcodes_matched[src])

    X_combined = vstack(all_matrices, format="csr").astype(np.int32)
    obs_combined = pd.concat(all_obs_dfs, axis=0)
    if rename_obs_cols:
        obs_combined = obs_combined.rename(columns=rename_obs_cols)
        print(f"[INFO] Renamed obs columns: {rename_obs_cols}")

    if min_cells_per_batch > 0:
        if _BATCH_COL not in obs_combined.columns:
            raise ValueError(
                f"--minimum-cells-per-batch requires '{_BATCH_COL}' in obs, "
                f"but it was not found. Available columns: {list(obs_combined.columns)}"
            )
        batch_counts = obs_combined[_BATCH_COL].value_counts()
        keep_batches = batch_counts[batch_counts >= min_cells_per_batch].index
        n_before = len(obs_combined)
        keep_mask = obs_combined[_BATCH_COL].isin(keep_batches)
        obs_combined = obs_combined[keep_mask]
        X_combined = X_combined[keep_mask.values]
        all_barcodes = [bc for bc, k in zip(all_barcodes, keep_mask) if k]
        n_dropped = n_before - len(obs_combined)
        n_batches_dropped = len(batch_counts) - len(keep_batches)
        print(
            f"[INFO] Dropped {n_dropped:,} cells from {n_batches_dropped} batches "
            f"with < {min_cells_per_batch} cells ('{_BATCH_COL}'). "
            f"{len(obs_combined):,} cells remaining."
        )

    for col in obs_combined.select_dtypes("category").columns:
        obs_combined[col] = obs_combined[col].cat.remove_unused_categories()

    var_df = pd.DataFrame(index=shared_genes)
    var_df.index.name = "gene_id"

    if filter_genes_mode != "intersection":
        print(f"[INFO] Applying --filter-genes {filter_genes_mode} ...")
        import matplotlib
        matplotlib.use("Agg")
        from cell2location.utils.filtering import filter_genes as c2l_filter_genes

        filter_adata = ad.AnnData(X=X_combined.copy(), var=var_df.copy())
        selected = set(c2l_filter_genes(
            filter_adata, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12
        ))
        print(f"  cell2location filter_genes: {len(shared_genes):,} → {len(selected):,} genes")

        if filter_genes_mode == "cell2location_u_da":
            da_marker_ensembl_ids = _resolve_da_marker_ensembl_ids(source_feature_names, subset_sources)
            added = [g for g in da_marker_ensembl_ids if g in gene_index and g not in selected]
            selected.update(added)
            print(f"  + {len(added)} Kamath DA marker gene(s) added back via union "
                  f"(present in the gene intersection but dropped by the cell2location filter)")

        shared_genes = [g for g in shared_genes if g in selected]
        col_idx = np.array([gene_index[g] for g in shared_genes], dtype=np.intp)
        X_combined = X_combined[:, col_idx]
        var_df = pd.DataFrame(index=shared_genes)
        var_df.index.name = "gene_id"
        print(f"  Final gene set: {len(shared_genes):,} genes")

    n_total = X_combined.shape[0]
    print(f"[INFO] Writing {n_total:,} cells × {len(shared_genes):,} genes to {output_h5ad} ...")
    out_dir = os.path.dirname(output_h5ad)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    out = ad.AnnData(X=X_combined, obs=obs_combined, var=var_df)
    out.obs_names = all_barcodes
    out.write_h5ad(output_h5ad, compression="gzip")

    print(f"[INFO] Combined obs dataframe:")
    print(obs_combined)

    label_counts = obs_combined[cell_type_key].value_counts()
    print(f"[INFO] Cell type distribution ({cell_type_key}):")
    for ct, n in label_counts.items():
        print(f"  {ct}: {n:,}")

    print(f"[DONE] {output_h5ad}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--integrated-h5ad",
        required=True,
        help=(
            "scVI/scANVI integrated H5AD. Used only for obs\n"
            "(cell labels, source column); raw counts come from --source-preqc-h5ads.\n"
        ),
    )
    parser.add_argument(
        "--source-preqc-h5ads",
        nargs="+",
        required=True,
        metavar="SOURCE=PATH",
        help=(
            "Per-source preQC H5AD paths with all-gene raw counts.\n"
            "Format: SOURCE_NAME=/path/to/preQC.h5ad\n"
            "Example:\n"
            "  asap_cohort=/path/asap-cohort.merged_cleaned_unfiltered_preQC.h5ad\n"
            "  kamath=/path/kamath_preQC.h5ad\n"
            "  siletti=/path/siletti_preQC.h5ad"
        ),
    )
    parser.add_argument(
        "--n-cells",
        type=int,
        default=None,
        help="Target cells per source. Default: all cells for each source.",
    )
    parser.add_argument(
        "--source-col",
        default="source",
        help=(
            "obs column in --integrated-h5ad identifying the source dataset.\n"
            "Default: source"
        ),
    )
    parser.add_argument(
        "--cell-type-key",
        default="cell_type_pred",
        help=(
            "obs column with scVI/scANVI cell-type labels.\n"
            "Must be present in --obs-cols. Default: cell_type_pred"
        ),
    )
    parser.add_argument(
        "--obs-cols",
        default="source,batch,cell_type_pred,pmdbs_sc_rnaseq.C_scANVI",
        help=(
            "Comma-separated obs columns to carry into the output H5AD.\n"
            "Default: source,batch,cell_type_pred,pmdbs_sc_rnaseq.C_scANVI"
        ),
    )
    parser.add_argument(
        "--minimum-cells-per-batch",
        type=int,
        default=0,
        metavar="N",
        help=(
            "Drop cells whose batch group has fewer than N cells.\n"
            "Applied to --batch-col after renaming. Default: 0 (disabled)."
        ),
    )
    parser.add_argument(
        "--filter-genes",
        choices=["intersection", "cell2location_default", "cell2location_u_da"],
        default="intersection",
        help=(
            "Gene set to write, beyond the raw cross-source intersection:\n"
            "  intersection           : no further filtering (default; current behaviour).\n"
            "  cell2location_default : apply cell2location.utils.filtering.filter_genes\n"
            "                           with the tutorial's defaults (cell_count_cutoff=5,\n"
            "                           cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12).\n"
            "  cell2location_u_da    : cell2location_default, unioned with Kamath's 17 DA\n"
            "                           marker genes."
        ),
    )
    parser.add_argument(
        "--output-h5ad",
        required=True,
        help="Destination path for the deconvolution reference H5AD.",
    )
    parser.add_argument(
        "--promote",
        action="store_true",
        default=False,
        help="Write the output file. Default is dry-run.",
    )

    args = parser.parse_args()

    preqc_map = {}
    for kv in args.source_preqc_h5ads:
        if "=" not in kv:
            parser.error(
                f"--source-preqc-h5ads entries must be SOURCE=PATH, got: {kv!r}"
            )
        src, path = kv.split("=", 1)
        preqc_map[src.strip()] = os.path.abspath(path.strip())

    subset_sources = list(preqc_map.keys())
    obs_cols = [c.strip() for c in args.obs_cols.split(",") if c.strip()]

    args.integrated_h5ad = os.path.abspath(args.integrated_h5ad)
    args.output_h5ad     = os.path.abspath(args.output_h5ad)
    args.output_dir      = os.path.dirname(args.output_h5ad)
    log_run_command(Path(args.output_dir) / f"{Path(sys.argv[0]).stem}.log")

    prepare_reference(
        integrated_h5ad=args.integrated_h5ad,
        preqc_map=preqc_map,
        subset_sources=subset_sources,
        n_cells=args.n_cells,
        source_col=args.source_col,
        obs_cols=obs_cols,
        rename_obs_cols=_RENAME_OBS_COLS,
        min_cells_per_batch=args.minimum_cells_per_batch,
        cell_type_key=args.cell_type_key,
        output_h5ad=args.output_h5ad,
        seed=42,
        promote=args.promote,
        filter_genes_mode=args.filter_genes,
    )
