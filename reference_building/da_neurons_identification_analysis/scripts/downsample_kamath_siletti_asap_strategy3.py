#!/usr/bin/env python3
"""
Build a scRNA-seq deconvolution reference (Strategy 3) using ASAP cells labeled by the sc-rnaseq-wf v5.0.0,
and the same Kamath/Siletti cells that were used in the Strategy 1 reference (which used the sc-rnaseq-wf v3.1.0 ASAP cells)

Source   Counts source (raw UMIs)                                 Label source                       Label handling
-------  -----------------------------------------------------    ---------------------------------  --------------------------------
Kamath   kamath/intermediate/kamath_merged_cell_type_min200_      cell_type_scanvi (same file)        Mapped via KAMATH_TO_SILETTI_ALLEN
         100000_preQC.h5ad (33,295 genes; 99,329/99,329
         barcodes match the Strategy 1 reference's Kamath cells)
Siletti  siletti/intermediate/midbrain_da_neurons15000_           supercluster_term (same file)       Used as-is (native Siletti/Allen
         cells47643_preQC.h5ad (58,232 genes; all 99,553                                            annotation, not a prediction)
         Strategy 1-reference Siletti barcodes present)
ASAP     v5.0.0's asap-cohort.merged_cleaned_unfiltered.h5ad      C_scANVI (v5.0.0's                  Used as-is (sc-rnaseq-wf's own
         (38,606 genes, full pre-QC-filter 5.0.0 cohort)          final_metadata.csv.gz)              scANVI prediction)

ASAP cell selection (--n-cells-asap, default 100,000):
  1. Of the Strategy 1 reference's 94,833 ASAP cells, keep whichever re-identify
     1:1 in v5.0.0's final (QC-filtered) metadata (barcode matching logic
     shared with map_ASAP_barcodes_v3_vs_v5_releases.py).
  2. Top up with any native ASAP cells (source is NaN in the 5.0.0 metadata,
     i.e. not a Kamath/Siletti spike-in) whose C_scANVI is one of the 10
     Kamath DA subtypes, not already selected in step 1.
  3. Random-fill with native ASAP cells to reach --n-cells-asap.

Outputs
-------
<outdir>/<stem>.h5ad
<outdir>/<stem>.log
"""

import sys
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

repo_root = next(p for p in Path(__file__).resolve().parents if (p / ".git").exists())
sys.path.insert(0, str(repo_root))

from utils.logging_extra import log_run_command
from utils.cell_type_maps import KAMATH_TO_SILETTI_ALLEN
from utils.genomics_constants import KAMATH_DA_MARKER_SYMBOLS, KAMATH_DA_SUBTYPE_LABELS
from utils.h5ad_utils import apply_cell2location_gene_filter, resolve_marker_ensembl_ids

sys.path.insert(0, str(repo_root / "reference_building/benchmark/reference_data/ma/scripts"))
from map_ASAP_barcodes_v3_vs_v5_releases import load_barcode_index, parse_v3_barcodes, parse_v5_barcodes, build_mapping

_STEM = "deconvolution_reference_strategy3"
_N_CELLS_ASAP_DEFAULT = 100_000

_CELL_TYPE_COLUMN = "harmonized_cell_type_siletti_allen"

_FILTER_GENES_CHOICES = ("intersection", "cell2location_default", "cell2location_u_da")
_FILTER_GENES_DEFAULT = "intersection"


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


def select_asap_barcodes(
    v3_asap_barcodes: pd.Index,
    v5_final_metadata_csv: str,
    n_cells_asap: int,
    min_cells_per_type: int,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """
    Select v5.0.0 ASAP barcodes: re-identified Strategy 1 (sc-rnaseq-wf's v3.1.0) cells, DA top-up,
    rare-cell-type top-up, random fill.

    Parameters
    ----------
    v3_asap_barcodes : pandas.Index
        Strategy 1 reference's ASAP barcodes (obs_names where source == "asap_cohort").
    v5_final_metadata_csv : str
        Path to v5.0.0's final (QC-filtered) metadata CSV.
    n_cells_asap : int
        Target total ASAP cell count.
    min_cells_per_type : int
        Top up any C_scANVI value sampled below this many cells.
    rng : numpy.random.Generator
        Random generator for the fill step.

    Returns
    -------
    pandas.DataFrame
        Indexed by v5.0.0 barcode: "C_scANVI", "selection_reason" (one of
        "mapped_from_v3.1.0", "da_topup", "type_topup", "random_fill").
    """
    print(f"  {len(v3_asap_barcodes):,} Strategy 1 (sc-rnaseq-wf's v3.1.0) ASAP barcode(s)")
    v3_df = parse_v3_barcodes(v3_asap_barcodes)
    print(f"  {len(v3_df):,} Strategy 1 (sc-rnaseq-wf's v3.1.0) ASAP barcode(s) parsed")

    print(f"[INFO] Loading v5.0.0 ASAP metadata from {v5_final_metadata_csv}")
    v5_cols = ["source", "is_spike_in", "C_scANVI", "team", "dataset", "batch"]
    header = pd.read_csv(v5_final_metadata_csv, nrows=0)
    idx_col = header.columns[0]
    v5_meta = pd.read_csv(v5_final_metadata_csv, usecols=[idx_col] + v5_cols, index_col=idx_col, low_memory=False)
    native_mask = v5_meta["source"].isna()
    v5_native = v5_meta.loc[native_mask]
    print(f"  {len(v5_meta):,} total rows; {len(v5_native):,} native ASAP cell(s) (source is NaN)")

    # batch_id: <team>-<batch> like `sc-rnaseq-wf to avoid collisions using "batch" alone (e.g. "Batch1") across datasets.
    v5_native = v5_native.copy()
    v5_native["batch_id"] = (
        v5_native["team"].astype(str) + "-" + v5_native["dataset"].astype(str) + "-" + v5_native["batch"].astype(str)
    )

    print("[INFO] Mapping Strategy 1 (sc-rnaseq-wf's v3.1.0) ASAP barcodes to v5.0.0 ...")
    v5_df = parse_v5_barcodes(pd.Index(v5_native.index))
    mapping = build_mapping(v3_df, v5_df)
    mapped_v5_barcodes = set(mapping.loc[mapping["status"] == "1:1_vs_v5.0.0_final", "v5.0.0_barcode"])
    print(f"  {len(mapped_v5_barcodes):,} v3.1.0 ASAP barcode(s) re-identified 1:1 in v5.0.0")

    is_da = v5_native["C_scANVI"].isin(KAMATH_DA_SUBTYPE_LABELS)
    da_barcodes = set(v5_native.index[is_da]) - mapped_v5_barcodes
    print(f"  {len(da_barcodes):,} additional native ASAP DA cell(s) available for top-up")

    reason_by_barcode = {bc: "mapped_from_v3.1.0" for bc in mapped_v5_barcodes}
    reason_by_barcode.update({bc: "da_topup" for bc in da_barcodes})

    if min_cells_per_type > 0:
        current_counts = v5_native.loc[list(reason_by_barcode), "C_scANVI"].value_counts()
        for ct, group in v5_native.groupby("C_scANVI", observed=True):
            have = current_counts.get(ct, 0)
            if have >= min_cells_per_type:
                continue
            candidates = group.index[~group.index.isin(reason_by_barcode)]
            n_take = min(min_cells_per_type - have, len(candidates))
            if n_take > 0:
                take = rng.choice(candidates, size=n_take, replace=False)
                reason_by_barcode.update({bc: "type_topup" for bc in take})
            reached = have + n_take
            status = "reached" if reached >= min_cells_per_type else f"capped, only {have + len(candidates)} available"
            print(f"[INFO] C_scANVI={ct!r}: had {have:,}, topped up {n_take:,} ({status})")

    n_fill = n_cells_asap - len(reason_by_barcode)
    if n_fill > 0:
        remaining_pool = np.array(sorted(set(v5_native.index) - set(reason_by_barcode)))
        if n_fill >= len(remaining_pool):
            fill = remaining_pool
            print(f"  Warning: only {len(remaining_pool):,} native ASAP cell(s) left in the pool, "
                  f"fewer than the {n_fill:,} needed for random fill")
        else:
            fill = rng.choice(remaining_pool, size=n_fill, replace=False)
        reason_by_barcode.update({bc: "random_fill" for bc in fill})
        print(f"  Random-filled {len(fill):,} more native ASAP cell(s) to reach {n_cells_asap:,}")
    elif n_fill < 0:
        print(f"  Note: mapped + DA top-up already {len(reason_by_barcode):,} cells, "
              f"exceeding target {n_cells_asap:,} — keeping all (no trimming)")

    selected = pd.DataFrame({
        "C_scANVI": v5_native.loc[list(reason_by_barcode), "C_scANVI"],
        "batch_id": v5_native.loc[list(reason_by_barcode), "batch_id"],
        "selection_reason": pd.Series(reason_by_barcode),
    })
    print(f"[INFO] {len(selected):,} ASAP barcode(s) selected total "
          f"({(selected['selection_reason'] == 'mapped_from_v3.1.0').sum():,} mapped_from_v3.1.0, "
          f"{(selected['selection_reason'] == 'da_topup').sum():,} da_topup, "
          f"{(selected['selection_reason'] == 'type_topup').sum():,} type_topup, "
          f"{(selected['selection_reason'] == 'random_fill').sum():,} random_fill)")
    return selected


def build_kamath_block(v3_kamath_barcodes: pd.Index, kamath_preqc_h5ad: str) -> ad.AnnData:
    """
    Load raw-UMI Kamath counts for the Strategy 1 (sc-rnaseq-wf's v3.1.0) reference's Kamath cells.

    Parameters
    ----------
    v3_kamath_barcodes : pandas.Index
        Strategy 1 reference's Kamath barcodes (obs_names where source == "kamath").
    kamath_preqc_h5ad : str
        Path to Kamath's preQC intermediate H5AD (full gene panel, raw counts).

    Returns
    -------
    ad.AnnData
        Kamath cells, X = raw UMI counts, obs[_CELL_TYPE_COLUMN] harmonized
        via KAMATH_TO_SILETTI_ALLEN, obs["source"] = "kamath",
        obs["batch_id"] = "kamath-" + obs["batch"].

    Raises
    ------
    KeyError
        A Kamath barcode is absent from `kamath_preqc_h5ad`, or its
        cell_type_scanvi value has no entry in KAMATH_TO_SILETTI_ALLEN.
    """
    print(f"  {len(v3_kamath_barcodes):,} Strategy 1 (sc-rnaseq-wf's v3.1.0) Kamath barcode(s)")

    print(f"[INFO] Loading raw Kamath counts from {kamath_preqc_h5ad}")
    kamath = ad.read_h5ad(kamath_preqc_h5ad, backed="r")
    missing = set(v3_kamath_barcodes) - set(kamath.obs_names)
    if missing:
        raise KeyError(f"{len(missing):,} v3.1.0 Kamath barcode(s) not found in {kamath_preqc_h5ad}, "
                        f"e.g. {sorted(missing)[:3]}")
    kamath = kamath[v3_kamath_barcodes].to_memory()
    kamath.X = kamath.layers["counts"]

    unmapped = sorted(set(kamath.obs["cell_type_scanvi"]) - set(KAMATH_TO_SILETTI_ALLEN))
    if unmapped:
        raise KeyError(f"cell_type_scanvi value(s) not in KAMATH_TO_SILETTI_ALLEN: {unmapped}")
    kamath.obs[_CELL_TYPE_COLUMN] = kamath.obs["cell_type_scanvi"].map(KAMATH_TO_SILETTI_ALLEN)
    kamath.obs["source"] = "kamath"
    kamath.obs["batch_id"] = "kamath-" + kamath.obs["batch"].astype(str)
    print(f"  {kamath.n_obs:,} cells x {kamath.n_vars:,} genes")
    return kamath


def build_siletti_block(v3_siletti_barcodes: pd.Index, siletti_preqc_h5ad: str) -> ad.AnnData:
    """
    Load raw-UMI Siletti counts for the Strategy 1 (sc-rnaseq-wf's v3.1.0) reference's Siletti cells.

    Parameters
    ----------
    v3_siletti_barcodes : pandas.Index
        Strategy 1 reference's Siletti barcodes (obs_names where source == "siletti").
    siletti_preqc_h5ad : str
        Path to Siletti's preQC intermediate H5AD (full gene panel, raw counts,
        native supercluster_term annotation).

    Returns
    -------
    ad.AnnData
        Siletti cells, X = raw UMI counts, obs[_CELL_TYPE_COLUMN] = native
        supercluster_term, obs["source"] = "siletti",
        obs["batch_id"] = "siletti-" + obs["batch"].

    Raises
    ------
    KeyError
        A v3.1.0 Siletti barcode is absent from `siletti_preqc_h5ad`.
    """
    print(f"  {len(v3_siletti_barcodes):,} Strategy 1 (sc-rnaseq-wf's v3.1.0) Siletti barcode(s)")

    print(f"[INFO] Loading raw Siletti counts from {siletti_preqc_h5ad}")
    siletti = ad.read_h5ad(siletti_preqc_h5ad, backed="r")
    missing = set(v3_siletti_barcodes) - set(siletti.obs_names)
    if missing:
        raise KeyError(f"{len(missing):,} v3.1.0 Siletti barcode(s) not found in {siletti_preqc_h5ad}, "
                        f"e.g. {sorted(missing)[:3]}")
    siletti = siletti[v3_siletti_barcodes].to_memory()
    siletti.X = siletti.layers["counts"]

    siletti.obs[_CELL_TYPE_COLUMN] = siletti.obs["supercluster_term"].astype(str)
    siletti.obs["source"] = "siletti"
    siletti.obs["batch_id"] = "siletti-" + siletti.obs["batch"].astype(str)
    print(f"  {siletti.n_obs:,} cells x {siletti.n_vars:,} genes")
    return siletti


def build_asap_block(selected: pd.DataFrame, asap_v5_unfiltered_h5ad: str) -> ad.AnnData:
    """
    Load raw-UMI v5.0.0 ASAP counts for the selected barcodes.

    Parameters
    ----------
    selected : pandas.DataFrame
        From `select_asap_barcodes` (indexed by v5.0.0 barcode).
    asap_v5_unfiltered_h5ad : str
        Path to v5.0.0's pre-QC-filter merged_cleaned_unfiltered H5AD.

    Returns
    -------
    ad.AnnData
        Selected ASAP cells, X = raw UMI counts, obs[_CELL_TYPE_COLUMN] =
        C_scANVI, obs["source"] = "asap_cohort", obs["selection_reason"] and
        obs["batch_id"] ("{team}-{dataset}-{batch}") carried through from
        `selected`.

    Raises
    ------
    KeyError
        A selected barcode is absent from `asap_v5_unfiltered_h5ad`.
    """
    print(f"[INFO] Loading raw ASAP counts from {asap_v5_unfiltered_h5ad}")
    asap = ad.read_h5ad(asap_v5_unfiltered_h5ad, backed="r")
    missing = set(selected.index) - set(asap.obs_names)
    if missing:
        raise KeyError(f"{len(missing):,} selected ASAP barcode(s) not found in {asap_v5_unfiltered_h5ad}, "
                        f"e.g. {sorted(missing)[:3]}")
    asap = asap[selected.index].to_memory()
    counts_layer = "counts" if "counts" in asap.layers else None
    if counts_layer is not None:
        asap.X = asap.layers[counts_layer]

    # var_names here are mostly gene symbols Ensembl ID fallback.
    if "gene_id" in asap.var.columns:
        asap.var_names = asap.var["gene_id"].values
    else:
        print("  Warning: 'gene_id' not in ASAP var columns — var_names left as-is. "
              "The inner join with Kamath/Siletti will likely produce far fewer shared genes.")

    asap.obs[_CELL_TYPE_COLUMN] = selected.loc[asap.obs_names, "C_scANVI"].values
    asap.obs["selection_reason"] = selected.loc[asap.obs_names, "selection_reason"].values
    asap.obs["batch_id"] = selected.loc[asap.obs_names, "batch_id"].values
    asap.obs["source"] = "asap_cohort"
    print(f"  {asap.n_obs:,} cells x {asap.n_vars:,} genes")
    return asap


def main(
    v3_reference_h5ad: str,
    kamath_preqc_h5ad: str,
    siletti_preqc_h5ad: str,
    asap_v5_unfiltered_h5ad: str,
    asap_v5_final_metadata_csv: str,
    n_cells_asap: int,
    min_cells_per_type: int,
    min_cells_per_batch: int,
    outdir: str,
    filter_genes_mode: str,
    promote: bool,
) -> None:
    """
    Build and write the ASAP-5.0.0 counterpart reference H5AD.

    Parameters
    ----------
    v3_reference_h5ad : str
        Path to the Strategy 1 (sc-rnaseq-wf's v3.1.0) deconvolution reference H5AD.
    kamath_preqc_h5ad : str
        Path to Kamath's preQC intermediate H5AD.
    siletti_preqc_h5ad : str
        Path to Siletti's preQC intermediate H5AD.
    asap_v5_unfiltered_h5ad : str
        Path to v5.0.0's pre-QC-filter merged_cleaned_unfiltered H5AD.
    asap_v5_final_metadata_csv : str
        Path to v5.0.0's final (QC-filtered) metadata CSV.
    n_cells_asap : int
        Target total ASAP cell count.
    min_cells_per_type : int
        Top up any ASAP C_scANVI value sampled below this many cells.
    min_cells_per_batch : int
        Drop cells whose batch_id group has fewer than this many cells.
    outdir : str
        Destination directory. Output is written as `{_STEM}.h5ad`.
    filter_genes_mode : str
        Gene set to write beyond the raw cross-source intersection.
    promote : bool
        If False, dry-run only (no file written).

    Returns
    -------
    None
    """
    rng = np.random.default_rng(42)

    print(f"[INFO] Loading Strategy 1 (sc-rnaseq-wf's v3.1.0) reference barcodes from {v3_reference_h5ad}")
    v3_adata = ad.read_h5ad(v3_reference_h5ad, backed="r")
    v3_kamath_barcodes = v3_adata.obs_names[v3_adata.obs["source"] == "kamath"]
    v3_siletti_barcodes = v3_adata.obs_names[v3_adata.obs["source"] == "siletti"]
    v3_asap_barcodes = v3_adata.obs_names[v3_adata.obs["source"] == "asap_cohort"]

    kamath = build_kamath_block(v3_kamath_barcodes, kamath_preqc_h5ad)
    siletti = build_siletti_block(v3_siletti_barcodes, siletti_preqc_h5ad)
    siletti_feature_name = siletti.var["feature_name"] if "feature_name" in siletti.var.columns else None
    selected_asap = select_asap_barcodes(
        v3_asap_barcodes, asap_v5_final_metadata_csv, n_cells_asap, min_cells_per_type, rng
    )
    asap = build_asap_block(selected_asap, asap_v5_unfiltered_h5ad)

    print("[INFO] Concatenating Kamath + Siletti + ASAP (inner join on genes) ...")
    combined = ad.concat(
        {"kamath": kamath, "siletti": siletti, "asap_cohort": asap},
        join="inner", label="_concat_source", index_unique=None,
    )
    combined.obs = combined.obs.drop(columns=["_concat_source"])
    print(f"[INFO] Combined: {combined.n_obs:,} cells x {combined.n_vars:,} genes")
    print(combined.obs["source"].value_counts())
    print(combined.obs[_CELL_TYPE_COLUMN].value_counts())

    if min_cells_per_batch > 0:
        batch_counts = combined.obs["batch_id"].value_counts()
        keep_batches = batch_counts[batch_counts >= min_cells_per_batch].index
        n_before = combined.n_obs
        keep_mask = combined.obs["batch_id"].isin(keep_batches)
        combined = combined[keep_mask.values].copy()
        n_dropped = n_before - combined.n_obs
        n_batches_dropped = len(batch_counts) - len(keep_batches)
        print(
            f"[INFO] Dropped {n_dropped:,} cells from {n_batches_dropped} batch(es) "
            f"with < {min_cells_per_batch} cells ('batch_id'). "
            f"{combined.n_obs:,} cells remaining."
        )

    if filter_genes_mode != "intersection":
        print(f"[INFO] Applying --filter-genes {filter_genes_mode} ...")
        da_marker_ensembl_ids = None
        if filter_genes_mode == "cell2location_u_da":
            if siletti_feature_name is None:
                raise ValueError(
                    "--filter-genes cell2location_u_da requires Siletti's "
                    "var['feature_name'] to resolve Kamath DA marker symbols to Ensembl IDs, "
                    "but it was not found."
                )
            da_marker_ensembl_ids = resolve_marker_ensembl_ids(
                siletti_feature_name, KAMATH_DA_MARKER_SYMBOLS, "siletti"
            )
        combined = apply_cell2location_gene_filter(combined, filter_genes_mode, da_marker_ensembl_ids)
        print(f"[INFO] After gene filtering: {combined.n_obs:,} cells x {combined.n_vars:,} genes")

    out_path = Path(outdir) / f"{_STEM}.h5ad"
    if not promote:
        print(f"[DRY-RUN] Would write: {out_path}")
        print(f"          Shape: {combined.n_obs:,} cells x {combined.n_vars:,} genes")
        print("          Re-run with --promote to write the file.")
        return

    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)
    combined.write_h5ad(out_path, compression="gzip")
    print(f"[DONE] Written {out_path}")


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
            '  "v3_reference_h5ad": "/path/to/deconvolution_reference_300k_c2lUda_Kamath200Rand.h5ad",\n'
            '  "kamath_preqc_h5ad": "/path/to/kamath_merged_cell_type_min200_100000_preQC.h5ad",\n'
            '  "siletti_preqc_h5ad": "/path/to/midbrain_da_neurons15000_cells47643_preQC.h5ad",\n'
            '  "asap_v5_unfiltered_h5ad": "/path/to/asap-cohort.merged_cleaned_unfiltered.h5ad",\n'
            '  "asap_v5_final_metadata_csv": "/path/to/asap-cohort.final_metadata.csv.gz",\n'
            '  "outdir": "/path/to/output_dir",\n'
            f'  "n_cells_asap": {_N_CELLS_ASAP_DEFAULT},\n'
            f'  "filter_genes": "{_FILTER_GENES_DEFAULT}"\n'
            "}\n"
            "\n"
            "Paths may be relative to the repo root or absolute.\n"
            "\n"
            "Fields:\n"
            "  v3_reference_h5ad (required)\n"
            "      Path to the Strategy 1 (sc-rnaseq-wf's v3.1.0) deconvolution reference H5AD.\n"
            "\n"
            "  kamath_preqc_h5ad (required)\n"
            "      Path to Kamath's preQC intermediate H5AD.\n"
            "\n"
            "  siletti_preqc_h5ad (required)\n"
            "      Path to Siletti's preQC intermediate H5AD.\n"
            "\n"
            "  asap_v5_unfiltered_h5ad (required)\n"
            "      Path to v5.0.0's pre-QC-filter merged_cleaned_unfiltered H5AD.\n"
            "\n"
            "  asap_v5_final_metadata_csv (required)\n"
            "      Path to v5.0.0's final (QC-filtered) metadata CSV.\n"
            "\n"
            "  outdir (required)\n"
            f"      Destination directory. Output is written as {_STEM}.h5ad.\n"
            "\n"
            "  n_cells_asap (optional)\n"
            "      Target total ASAP cell count.\n"
            f"      Default: {_N_CELLS_ASAP_DEFAULT}\n"
            "\n"
            "  minimum_cells_per_type (optional)\n"
            "      Top up any ASAP C_scANVI cell type sampled below this many cells.\n"
            "      Default: 0 (disabled)\n"
            "\n"
            "  minimum_cells_per_batch (optional)\n"
            "      Drop cells whose batch_id group has fewer than this many cells\n"
            "      Default: 0 (disabled)\n"
            "\n"
            "  filter_genes (optional)\n"
            "      Gene set to write, beyond the raw cross-source intersection:\n"
            "        intersection          : no further filtering (default).\n"
            "        cell2location_default : cell2location.utils.filtering.filter_genes\n"
            "                                with the tutorial's defaults (cell_count_cutoff=5,\n"
            "                                cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12).\n"
            "        cell2location_u_da    : cell2location_default, unioned with Kamath's\n"
            "                                17 DA marker genes.\n"
            f"      Default: {_FILTER_GENES_DEFAULT}\n"
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
        v3_reference_h5ad = config["v3_reference_h5ad"]
        kamath_preqc_h5ad = config["kamath_preqc_h5ad"]
        siletti_preqc_h5ad = config["siletti_preqc_h5ad"]
        asap_v5_unfiltered_h5ad = config["asap_v5_unfiltered_h5ad"]
        asap_v5_final_metadata_csv = config["asap_v5_final_metadata_csv"]
        outdir = config["outdir"]
    except KeyError as e:
        parser.error(f"Missing required field {e} in {args.input_json}.")

    filter_genes_mode = config.get("filter_genes", _FILTER_GENES_DEFAULT)
    if filter_genes_mode not in _FILTER_GENES_CHOICES:
        parser.error(f"filter_genes must be one of {_FILTER_GENES_CHOICES}, got: {filter_genes_mode!r}")

    v3_reference_h5ad = _resolve_repo_path(v3_reference_h5ad)
    kamath_preqc_h5ad = _resolve_repo_path(kamath_preqc_h5ad)
    siletti_preqc_h5ad = _resolve_repo_path(siletti_preqc_h5ad)
    asap_v5_unfiltered_h5ad = _resolve_repo_path(asap_v5_unfiltered_h5ad)
    asap_v5_final_metadata_csv = _resolve_repo_path(asap_v5_final_metadata_csv)
    n_cells_asap = config.get("n_cells_asap", _N_CELLS_ASAP_DEFAULT)
    min_cells_per_type = config.get("minimum_cells_per_type", 0)
    min_cells_per_batch = config.get("minimum_cells_per_batch", 0)
    outdir = _resolve_repo_path(outdir)

    log_run_command(Path(outdir) / f"{_STEM}.log")
    main(
        v3_reference_h5ad,
        kamath_preqc_h5ad,
        siletti_preqc_h5ad,
        asap_v5_unfiltered_h5ad,
        asap_v5_final_metadata_csv,
        n_cells_asap,
        min_cells_per_type,
        min_cells_per_batch,
        outdir,
        filter_genes_mode,
        args.promote,
    )
