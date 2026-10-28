#!/usr/bin/env python3
"""Build a multi-source deconvolution reference H5AD from an integrated (HVG-only) obs table
(H5AD or tabular) plus each source's all-gene preQC H5AD, intersecting genes (Ensembl IDs) and
cells across sources. The resulting 'reference' will have all cells in the integrated obs table
and all-genes intersecting across sources."""

import argparse
import json
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
from utils.genomics_constants import KAMATH_DA_MARKER_SYMBOLS, KAMATH_DA_SUBTYPE_LABELS


# obs columns always renamed in the output reference H5AD (old → new); merged
# with, and overridable by, the JSON's optional rename_obs_cols field.
_DEFAULT_RENAME_OBS_COLS = {
    "batch": "batch_id",
}

_BATCH_COL = "batch_id"  # obs column used for per-batch cell-count filtering

# sc-rnaseq-wf's spike-in convention: source is only set for spike-in cells
# (e.g. "kamath", "siletti"); native ASAP cohort cells are left NaN. Treated
# as this literal source name throughout.
_NAN_SOURCE_LABEL = "asap_cohort"


def _resolve_repo_path(path: str) -> str:
    """
    Resolve `path` against the repo root if relative, else return it unchanged.

    Lets JSON configs use repo-relative paths (portable across machines/clones)
    while still accepting absolute paths for anything outside the repo.

    Parameters
    ----------
    path : str
        Path from a JSON config; relative (to the repo root) or absolute.

    Returns
    -------
    str
        Absolute path.
    """
    p = Path(path)
    return str(p if p.is_absolute() else (repo_root / p).resolve())


def _stringify_mixed_type_cols(obs: pd.DataFrame, cols: list) -> None:
    """
    Coerce mixed-type object columns in `cols` to strings dtypes.
    NaN values are left untouched (not coerced to "nan").

    Parameters
    ----------
    obs : pandas.DataFrame
        Obs table, modified in place.
    cols : list
        Column names to check/fix (duplicates ignored).
    """
    for col in dict.fromkeys(cols):
        if obs[col].dtype != object:
            continue
        types = obs[col].map(lambda v: type(v).__name__ if not pd.isna(v) else None).dropna()
        distinct = sorted(types.unique())
        if len(distinct) > 1:
            print(f"[INFO] obs['{col}'] has mixed types ({distinct}); "
                  "coercing non-null values to str.")
            obs[col] = obs[col].map(lambda v: v if pd.isna(v) else str(v))


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


def _detect_file_type(path: str) -> dict:
    """
    Infer file kind, delimiter, and compression from `path`'s filename.

    Parameters
    ----------
    path : str
        Path ending in .h5ad, .csv[.gz|.bz2], or .tsv[.gz|.bz2].

    Returns
    -------
    dict
        kind : str
            "h5ad" or "tabular".
        sep : str or None
            "," for .csv, "\\t" for .tsv; None when kind is "h5ad".
        compression : str or None
            "gzip", "bz2", or None; always None when kind is "h5ad".
    """
    suffixes = [s.lower() for s in Path(path).suffixes]

    compression = None
    if suffixes and suffixes[-1] in (".gz", ".bz2"):
        compression = {".gz": "gzip", ".bz2": "bz2"}[suffixes[-1]]
        suffixes = suffixes[:-1]

    ext = suffixes[-1] if suffixes else ""
    if ext == ".h5ad":
        if compression is not None:
            raise ValueError(f"Unexpected compression suffix on an .h5ad file: {path}")
        return {"kind": "h5ad", "sep": None, "compression": None}
    if ext == ".csv":
        return {"kind": "tabular", "sep": ",", "compression": compression}
    if ext == ".tsv":
        return {"kind": "tabular", "sep": "\t", "compression": compression}

    raise ValueError(
        f"Unsupported file extension for --integrated-file: {path}. "
        "Expected .h5ad, .csv[.gz|.bz2], or .tsv[.gz|.bz2]."
    )


def _load_integrated_obs(path: str) -> pd.DataFrame:
    """
    Load the obs/metadata table from `path`, auto-detecting H5AD vs tabular input.

    For an H5AD, only .obs is read (backed, raw counts are never loaded).
    For a tabular file, the first column is used as the barcode index.

    Parameters
    ----------
    path : str
        Path to a .h5ad, .csv[.gz|.bz2], or .tsv[.gz|.bz2] file.

    Returns
    -------
    pandas.DataFrame
        Obs table indexed by barcode.
    """
    file_type = _detect_file_type(path)
    if file_type["kind"] == "h5ad":
        return ad.read_h5ad(path, backed="r").obs
    return pd.read_csv(path, sep=file_type["sep"], compression=file_type["compression"], index_col=0)


def _sample_per_source(
    integrated_obs: pd.DataFrame,
    source_col: str,
    subset_sources: list,
    n_cells_map: dict,
    seed: int,
) -> dict:
    """
    Sample barcodes per source from the integrated obs table, using a
    per-source target cell count.

    Uses contiguous-block sampling within each source group to minimise
    random HDF5 seeks when the caller later reads preQC rows.

    Parameters
    ----------
    integrated_obs : pandas.DataFrame
        Obs table (barcode-indexed) loaded via `_load_integrated_obs`.
    source_col : str
        Column identifying the source dataset.
    subset_sources : list
        Sources to include (must be values in obs[source_col]).
    n_cells_map : dict
        source → target cell count. A source absent from this dict, or
        mapped to None, takes all of that source's cells.
    seed : int
        Random seed for block-start selection.

    Returns
    -------
    dict
        source → list of barcode strings
    """
    rng = np.random.default_rng(seed)
    result = {}
    obs = integrated_obs

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
        n_cells = n_cells_map.get(src)
        if n_cells is None or n_cells >= n:
            print(f"  [{src}] {n:,} cells — taking all")
            result[src] = barcodes
        else:
            start = int(rng.integers(0, n - n_cells + 1))
            result[src] = barcodes[start:start + n_cells]
            print(f"  [{src}] {n:,} total → sampled {n_cells:,} (block offset {start:,})")

    return result


def _top_up_rare_cell_types(
    integrated_obs: pd.DataFrame,
    source_barcodes: dict,
    source_col: str,
    subset_sources: list,
    cell_type_key: str,
    min_cells_per_type: int,
    seed: int,
) -> dict:
    """
    Top up cell types under-represented in `source_barcodes`, drawing
    extra cells from integrated_obs to top up min_cells_per_type.

    Parameters
    ----------
    integrated_obs : pandas.DataFrame
        Full obs table (barcode-indexed).
    source_barcodes : dict
        source → list of already-selected barcodes (not mutated).
    source_col : str
        Column identifying the source dataset.
    subset_sources : list
        Sources eligible as top-up candidates (must have a preQC H5AD).
    cell_type_key : str
        obs column with cell-type labels to balance.
    min_cells_per_type : int
        Minimum cells desired per distinct cell_type_key value.
    seed : int
        Random seed for selecting which extra cells to add.

    Returns
    -------
    dict
        source → list of barcodes, with any needed top-ups merged in.
    """
    rng = np.random.default_rng(seed)
    result = {src: list(bcs) for src, bcs in source_barcodes.items()}
    selected = set(bc for bcs in result.values() for bc in bcs)

    eligible = integrated_obs[integrated_obs[source_col].isin(subset_sources)]
    current_counts = eligible.loc[eligible.index.isin(selected), cell_type_key].value_counts()

    for ct, group in eligible.groupby(cell_type_key, observed=True):
        have = current_counts.get(ct, 0)
        if have >= min_cells_per_type:
            continue
        candidates = group.index[~group.index.isin(selected)]
        n_take = min(min_cells_per_type - have, len(candidates))
        if n_take <= 0:
            continue
        chosen = rng.choice(candidates.to_numpy(), size=n_take, replace=False)
        for bc in chosen:
            result.setdefault(eligible.at[bc, source_col], []).append(bc)
        selected.update(chosen)
        reached = have + n_take
        status = "reached" if reached >= min_cells_per_type else f"capped, only {have + len(candidates)} available"
        print(f"[INFO] Topped up '{ct}' from {have:,} to {reached:,} cells ({status})")

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
            f"in preQC obs (e.g. {examples}). Check --integrated-file and "
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


def _collapse_kamath_da_subtypes(labels: pd.Series) -> pd.Series:
    """
    Collapse Kamath et al. (2022) DA subtype labels (e.g. SOX6_AGTR1, CALB1_GEM)
    back into a single "Dopaminergic" category. All other values are unchanged.

    Parameters
    ----------
    labels : pandas.Series
        Cell-type labels (e.g. obs[cell_type_key]).

    Returns
    -------
    pandas.Series
        Copy of `labels`, object dtype, with Kamath DA subtype values replaced
        by "Dopaminergic".
    """
    labels = labels.astype(object)
    is_da_subtype = labels.isin(KAMATH_DA_SUBTYPE_LABELS)
    return labels.where(~is_da_subtype, "Dopaminergic")


def prepare_reference(
    integrated_file: str,
    preqc_map: dict,
    subset_sources: list,
    n_cells_map: dict,
    source_col: str,
    obs_cols: list,
    rename_obs_cols: dict,
    min_cells_per_batch: int,
    min_cells_per_type: int,
    cell_type_key: str,
    collapse_kamath_da_subtypes: bool,
    output_h5ad: str,
    seed: int,
    promote: bool,
    filter_genes_mode: str,
) -> None:
    """
    Build a multi-source deconvolution reference H5AD.

    Parameters
    ----------
    integrated_file : str
        scVI/scANVI integrated obs table: .h5ad (backed; only obs accessed),
        or tabular (.csv/.tsv, optionally .gz/.bz2) with barcodes in the first column.
    preqc_map : dict
        source_name → path to preQC H5AD with all-gene raw counts.
    subset_sources : list
        Source names to include (must be keys in `preqc_map` and values in obs[source_col]).
    n_cells_map : dict
        source → target cell count. A source absent from this dict, or mapped to
        None, takes all of that source's cells.
    source_col : str
        obs column in the integrated H5AD that holds source identifiers.
    obs_cols : list
        obs columns to carry into the output H5AD.
    rename_obs_cols : dict
        old_name → new_name renames applied to obs before writing.
    min_cells_per_batch : int
        Drop cells whose `_BATCH_COL` group has fewer than this many cells. 0 disables filtering.
    min_cells_per_type : int
        Top up any cell_type_key value under this many cells, drawing extras from
        integrated_obs not already selected (capped at however many actually exist).
        0 disables this.
    cell_type_key : str
        obs column with cell-type labels (must be in obs_cols).
    collapse_kamath_da_subtypes : bool
        If True, obs[<cell_type_key>_final] collapses Kamath DA subtype labels
        in `cell_type_key` (e.g. SOX6_AGTR1, CALB1_GEM) into "Dopaminergic".
        If False, obs[<cell_type_key>_final] is a verbatim copy of `cell_type_key`.
    output_h5ad : str
        Destination path.
    seed : int
        Random seed for sampling.
    promote : bool
        If False, dry-run only.
    filter_genes_mode : str
        Gene set to write beyond the raw cross-source intersection; see --filter-genes
        argparse help for the "intersection" / "cell2location_default" / "cell2location_u_da" semantics.

    Returns
    -------
    None
    """
    print(f"[INFO] Loading integrated obs from: {integrated_file}")
    integrated_obs = _load_integrated_obs(integrated_file)
    n_nan_source = integrated_obs[source_col].isna().sum()
    if n_nan_source:
        integrated_obs[source_col] = integrated_obs[source_col].fillna(_NAN_SOURCE_LABEL)
        print(f"[INFO] Filled {n_nan_source:,} NaN '{source_col}' value(s) as '{_NAN_SOURCE_LABEL}' "
              "(native cohort cells are unlabelled; only spike-ins carry a source).")
    print(f"  {len(integrated_obs):,} cells total, sources: {sorted(integrated_obs[source_col].unique())}")

    missing_cols = [c for c in obs_cols if c not in integrated_obs.columns]
    if missing_cols:
        raise ValueError(
            f"obs columns not found in integrated file: {missing_cols}. "
            f"Available: {list(integrated_obs.columns)}"
        )

    _stringify_mixed_type_cols(integrated_obs, obs_cols + [source_col])

    n_cells_label = ", ".join(
        f"{src}={n_cells_map[src]:,}" if n_cells_map.get(src) is not None else f"{src}=all"
        for src in subset_sources
    )
    print(f"[INFO] Sampling cells per source: {n_cells_label}")
    source_barcodes = _sample_per_source(integrated_obs, source_col, subset_sources, n_cells_map, seed)

    if min_cells_per_type > 0:
        source_barcodes = _top_up_rare_cell_types(
            integrated_obs, source_barcodes, source_col, subset_sources,
            cell_type_key, min_cells_per_type, seed,
        )

    source_obs = {}
    for src, barcodes in source_barcodes.items():
        source_obs[src] = integrated_obs.loc[barcodes, obs_cols].copy()

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

    final_col = f"{cell_type_key}_final"
    if collapse_kamath_da_subtypes:
        obs_combined[final_col] = _collapse_kamath_da_subtypes(obs_combined[cell_type_key])
        print(f"[INFO] Collapsed Kamath DA subtype labels in '{cell_type_key}' into "
              f"'Dopaminergic' for obs['{final_col}']")
    else:
        obs_combined[final_col] = obs_combined[cell_type_key]

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
    out.write_h5ad(output_h5ad, compression="lzf")

    print(f"[INFO] Combined obs dataframe:")
    print(obs_combined)

    label_counts = obs_combined[cell_type_key].value_counts()
    print(f"[INFO] Cell type distribution ({cell_type_key}):")
    for ct, n in label_counts.items():
        print(f"  {ct}: {n:,}")

    final_label_counts = obs_combined[final_col].value_counts()
    print(f"[INFO] Cell type distribution ({final_col}):")
    for ct, n in final_label_counts.items():
        print(f"  {ct}: {n:,}")

    print(f"[DONE] {output_h5ad}")

# _DEFAULT_OBS_COLS = ["source", "batch", "cell_type_pred", "pmdbs_sc_rnaseq.C_scANVI"] # used in releases < v5.0.0
_DEFAULT_OBS_COLS = ["source", "batch", "ground_truth_cell_type", "cell_type", "C_scANVI"]
_FILTER_GENES_CHOICES = ("intersection", "cell2location_default", "cell2location_u_da")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input-json",
        required=True,
        help=(
            "Path to a JSON config file, e.g.:\n"
            "\n"
            "{\n"
            '  "integrated_file": "/path/to/asap-cohort.final_metadata.csv.gz",\n'
            '  "source_preqc_h5ads": {\n'
            '    "asap_cohort": "/path/to/asap-cohort.merged_cleaned_unfiltered_preQC.h5ad",\n'
            '    "kamath": "/path/to/kamath_preQC.h5ad",\n'
            '    "siletti": "/path/to/siletti_preQC.h5ad"\n'
            "  },\n"
            '  "n_cells": {"kamath": 20000, "siletti": 50000, "asap_cohort": 100000},\n'
            '  "output_h5ad": "/path/to/deconvolution_reference.h5ad"\n'
            "}\n"
            "\n"
            "All paths may be relative to the repo root or absolute.\n"
            "\n"
            "Fields:\n"
            "  integrated_file (required)\n"
            "      scVI/scANVI integrated obs table: .h5ad (obs only; raw counts come\n"
            "      from source_preqc_h5ads), or tabular (.csv/.tsv, optionally .gz/.bz2)\n"
            "      with barcodes in the first column. Repo-relative or absolute path.\n"
            "\n"
            "  source_preqc_h5ads (required)\n"
            "      source -> preQC H5AD path (repo-relative or absolute) with all-gene\n"
            "      raw counts. Defines which sources are included; every other\n"
            "      source-keyed field (e.g. n_cells) must use these same source names.\n"
            "\n"
            "  n_cells (optional)\n"
            "      source -> target cell count: either i) an integer to downsample that many\n"
            "      cells from the integrated_file, or ii) \"all_integrated\" to take every\n"
            "      cell present for that source in integrated_file. Cells must also exist\n"
            "      in the source's preQC H5AD. A source omitted here, also takes all cells\n"
            "      for that source.\n"
            "  source_col (optional)\n"
            "      obs column in integrated_file identifying the source dataset.\n"
            "      NaN values are treated as \"asap_cohort\" because sc-rnaseq-wf\n"
            "      only indicates the source for spike-in cells.\n"
            "      Default: \"source\"\n"
            "\n"
            "  cell_type_key (optional)\n"
            "      obs column with scVI/scANVI cell-type labels. Must be present in\n"
            "      obs_cols. Default: \"C_scANVI\"\n"
            "\n"
            "  collapse_kamath_da_subtypes (optional)\n"
            "      Collapse Kamath et al. (2022) DA subtype labels (e.g. SOX6_AGTR1,\n"
            "      CALB1_GEM) in cell_type_key back into \"Dopaminergic\" when writing\n"
            "      obs[<cell_type_key>_final]. False writes obs[<cell_type_key>_final]\n"
            "      as a verbatim copy of cell_type_key instead. Default: true\n"
            "\n"
            "  obs_cols (optional)\n"
            "      obs columns to carry into the output H5AD.\n"
            f"      Default: {_DEFAULT_OBS_COLS}\n"
            "\n"
            "  rename_obs_cols (optional)\n"
            "      old_name -> new_name renames applied to obs_cols in the output H5AD.\n"
            f"      Merged with, and able to override, the always-applied {_DEFAULT_RENAME_OBS_COLS}.\n"
            "      Example: {\"cell_type\": \"cell_type_MMC\"}\n"
            "\n"
            "  minimum_cells_per_batch (optional)\n"
            "      Drop cells whose batch group has fewer than N cells. Applied to\n"
            "      batch_col after renaming. Default: 0 (disabled)\n"
            "\n"
            "  minimum_cells_per_type (optional)\n"
            "      Opposite of minimum_cells_per_batch: tops up (not drops) any\n"
            "      cell_type_key value sampled below N cells, drawing extras from\n"
            "      integrated_file not already selected. A type with fewer than N\n"
            "      cells available in total is capped at however many actually exist\n"
            "      (no error). Default: 0 (disabled)\n"
            "\n"
            "  filter_genes (optional)\n"
            "      Gene set to write, beyond the raw cross-source intersection:\n"
            "        intersection          : no further filtering (default).\n"
            "        cell2location_default : cell2location.utils.filtering.filter_genes\n"
            "                                with the tutorial's defaults (cell_count_cutoff=5,\n"
            "                                cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12).\n"
            "        cell2location_u_da    : cell2location_default, unioned with Kamath's\n"
            "                                17 DA marker genes.\n"
            "\n"
            "  output_h5ad (required)\n"
            "      Destination path (repo-relative or absolute) for the deconvolution\n"
            "      reference H5AD.\n"
        ),
    )
    parser.add_argument(
        "--promote",
        action="store_true",
        default=False,
        help="Write the output file. Default is dry-run.",
    )

    args = parser.parse_args()

    try:
        with open(args.input_json, "r") as f:
            config = json.load(f)
    except FileNotFoundError:
        parser.error(f"{args.input_json} JSON file does not exist.")

    try:
        integrated_file = config["integrated_file"]
        preqc_map_raw = config["source_preqc_h5ads"]
        output_h5ad = config["output_h5ad"]
    except KeyError as e:
        parser.error(f"Missing required field {e} in {args.input_json}.")

    preqc_map = {src: _resolve_repo_path(path) for src, path in preqc_map_raw.items()}
    subset_sources = list(preqc_map.keys())

    n_cells_map = {}
    for src, n in config.get("n_cells", {}).items():
        if src not in preqc_map:
            parser.error(
                f"n_cells source {src!r} not found in source_preqc_h5ads: {sorted(preqc_map)}"
            )
        if n == "all_integrated":
            n_cells_map[src] = None
        elif isinstance(n, int) and not isinstance(n, bool):
            n_cells_map[src] = n
        else:
            parser.error(
                f'n_cells value for {src!r} must be an integer or "all_integrated", got: {n!r}'
            )

    filter_genes_mode = config.get("filter_genes", "intersection")
    if filter_genes_mode not in _FILTER_GENES_CHOICES:
        parser.error(f"filter_genes must be one of {_FILTER_GENES_CHOICES}, got: {filter_genes_mode!r}")

    integrated_file = _resolve_repo_path(integrated_file)
    output_h5ad     = _resolve_repo_path(output_h5ad)
    output_dir      = os.path.dirname(output_h5ad)
    log_run_command(Path(output_dir) / f"{Path(output_h5ad).stem}.log")

    prepare_reference(
        integrated_file=integrated_file,
        preqc_map=preqc_map,
        subset_sources=subset_sources,
        n_cells_map=n_cells_map,
        source_col=config.get("source_col", "source"),
        obs_cols=config.get("obs_cols", _DEFAULT_OBS_COLS),
        rename_obs_cols={**_DEFAULT_RENAME_OBS_COLS, **config.get("rename_obs_cols", {})},
        min_cells_per_batch=config.get("minimum_cells_per_batch", 0),
        min_cells_per_type=config.get("minimum_cells_per_type", 0),
        cell_type_key=config.get("cell_type_key", "C_scANVI"),
        collapse_kamath_da_subtypes=config.get("collapse_kamath_da_subtypes", True),
        output_h5ad=output_h5ad,
        seed=42,
        promote=args.promote,
        filter_genes_mode=filter_genes_mode,
    )
