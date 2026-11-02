#!/usr/bin/env python3
"""
Build a scRNA-seq deconvolution reference (Strategy 3) using ASAP cells labeled by the sc-rnaseq-wf v5.0.0,
and the same Kamath/Siletti cells that were used in the Strategy 1 reference (which used the sc-rnaseq-wf v3.1.0 ASAP cells)

Source   Counts source (raw UMIs)                               Label source                       Label handling
-------  -----------------------------------------------------  ---------------------------------  --------------------------------
Kamath   kamath/intermediate/kamath_merged_cell_type_min200_    cell_type_scanvi (same file)        Mapped via KAMATH_TO_SILETTI_ALLEN
         100000_preQC.h5ad (33,295 genes; 99,329/99,329
         barcodes match the v3.1.0 reference's Kamath cells)
Siletti  siletti/intermediate/midbrain_da_neurons15000_         supercluster_term (same file)       Used as-is (native Siletti/Allen
         cells47643_preQC.h5ad (58,232 genes; all 99,553                                            annotation, not a prediction)
         v3.1.0-reference Siletti barcodes present)
ASAP     v5.0.0's asap-cohort.merged_cleaned_unfiltered.h5ad    C_scANVI (v5.0.0's                  Used as-is (sc-rnaseq-wf's own
         (38,606 genes, full pre-QC-filter 5.0.0 cohort)        final_metadata.csv.gz)              scANVI prediction)

ASAP cell selection (--n-cells-asap, default 100,000):
  1. Of the v3.1.0 reference's 94,833 ASAP cells, keep whichever re-identify
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
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

repo_root = next(p for p in Path(__file__).resolve().parents if (p / ".git").exists())
sys.path.insert(0, str(repo_root))

from utils.logging_extra import log_run_command
from utils.cell_type_maps import KAMATH_TO_SILETTI_ALLEN
from utils.genomics_constants import KAMATH_DA_SUBTYPE_LABELS

sys.path.insert(0, str(repo_root / "reference_building/benchmark/reference_data/ma/scripts"))
from map_ASAP_barcodes_v3_vs_v5_releases import load_barcode_index, parse_v3_barcodes, parse_v5_barcodes, build_mapping

_V3_REFERENCE_H5AD = "reference_building/reference_h5ad/deconvolution_reference_300k_c2lUda_Kamath200Rand.h5ad"
_KAMATH_PREQC_H5AD = (
    "reference_building/da_neurons_identification_analysis/reference_data/kamath/intermediate/"
    "kamath_merged_cell_type_min200_100000_preQC.h5ad"
)
_SILETTI_PREQC_H5AD = (
    "reference_building/da_neurons_identification_analysis/reference_data/siletti/intermediate/"
    "midbrain_da_neurons15000_cells47643_preQC.h5ad"
)
_ASAP_V5_UNFILTERED_H5AD = (
    "reference_building/da_neurons_identification_analysis/reference_data/asap_cohort/"
    "pmdbs_sc_rnaseq/cohort_analysis/5.0.0/2026-10-19T14-59-48Z/asap-cohort.merged_cleaned_unfiltered.h5ad"
)
_ASAP_V5_FINAL_METADATA_CSV = (
    "reference_building/da_neurons_identification_analysis/reference_data/asap_cohort/"
    "pmdbs_sc_rnaseq/cohort_analysis/5.0.0/2026-10-19T14-59-48Z/asap-cohort.final_metadata.csv.gz"
)
_OUTDIR = "reference_building/reference_h5ad"
_STEM = "deconvolution_reference_strategy3"
_N_CELLS_ASAP_DEFAULT = 100_000

_CELL_TYPE_COLUMN = "harmonized_cell_type_siletti_allen"


def select_asap_barcodes(
    v3_reference_h5ad: str,
    v5_final_metadata_csv: str,
    n_cells_asap: int,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """
    Select v5.0.0 ASAP barcodes: re-identified v3.1.0 cells, DA top-up, random fill.

    Parameters
    ----------
    v3_reference_h5ad : str
        Path to the v3.1.0-era deconvolution reference H5AD.
    v5_final_metadata_csv : str
        Path to v5.0.0's final (QC-filtered) metadata CSV.
    n_cells_asap : int
        Target total ASAP cell count.
    rng : numpy.random.Generator
        Random generator for the fill step.

    Returns
    -------
    pandas.DataFrame
        Indexed by v5.0.0 barcode: "C_scANVI", "selection_reason" (one of
        "mapped_from_v3.1.0", "da_topup", "random_fill").
    """
    print(f"[INFO] Loading v3.1.0 reference ASAP barcodes from {v3_reference_h5ad}")
    v3_adata = ad.read_h5ad(v3_reference_h5ad, backed="r")
    v3_asap_barcodes = v3_adata.obs_names[v3_adata.obs["source"] == "asap_cohort"]
    v3_df = parse_v3_barcodes(v3_asap_barcodes)
    print(f"  {len(v3_df):,} v3.1.0 ASAP barcode(s) parsed")

    print(f"[INFO] Loading v5.0.0 ASAP metadata from {v5_final_metadata_csv}")
    v5_cols = ["source", "is_spike_in", "C_scANVI"]
    header = pd.read_csv(v5_final_metadata_csv, nrows=0)
    idx_col = header.columns[0]
    v5_meta = pd.read_csv(v5_final_metadata_csv, usecols=[idx_col] + v5_cols, index_col=idx_col, low_memory=False)
    native_mask = v5_meta["source"].isna()
    v5_native = v5_meta.loc[native_mask]
    print(f"  {len(v5_meta):,} total rows; {len(v5_native):,} native ASAP cell(s) (source is NaN)")

    print("[INFO] Mapping v3.1.0 ASAP barcodes to v5.0.0 ...")
    v5_df = parse_v5_barcodes(pd.Index(v5_native.index))
    mapping = build_mapping(v3_df, v5_df)
    mapped_v5_barcodes = set(mapping.loc[mapping["status"] == "1:1_vs_v5.0.0_final", "v5.0.0_barcode"])
    print(f"  {len(mapped_v5_barcodes):,} v3.1.0 ASAP barcode(s) re-identified 1:1 in v5.0.0")

    is_da = v5_native["C_scANVI"].isin(KAMATH_DA_SUBTYPE_LABELS)
    da_barcodes = set(v5_native.index[is_da]) - mapped_v5_barcodes
    print(f"  {len(da_barcodes):,} additional native ASAP DA cell(s) available for top-up")

    reason_by_barcode = {bc: "mapped_from_v3.1.0" for bc in mapped_v5_barcodes}
    reason_by_barcode.update({bc: "da_topup" for bc in da_barcodes})

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
        "selection_reason": pd.Series(reason_by_barcode),
    })
    print(f"[INFO] {len(selected):,} ASAP barcode(s) selected total "
          f"({(selected['selection_reason'] == 'mapped_from_v3.1.0').sum():,} mapped_from_v3.1.0, "
          f"{(selected['selection_reason'] == 'da_topup').sum():,} da_topup, "
          f"{(selected['selection_reason'] == 'random_fill').sum():,} random_fill)")
    return selected


def build_kamath_block(v3_reference_h5ad: str, kamath_preqc_h5ad: str) -> ad.AnnData:
    """
    Load raw-UMI Kamath counts for the v3.1.0 reference's Kamath cells.

    Parameters
    ----------
    v3_reference_h5ad : str
        Path to the v3.1.0-era deconvolution reference H5AD.
    kamath_preqc_h5ad : str
        Path to Kamath's preQC intermediate H5AD (full gene panel, raw counts).

    Returns
    -------
    ad.AnnData
        Kamath cells, X = raw UMI counts, obs[_CELL_TYPE_COLUMN] harmonized
        via KAMATH_TO_SILETTI_ALLEN, obs["source"] = "kamath".

    Raises
    ------
    KeyError
        A v3.1.0 Kamath barcode is absent from `kamath_preqc_h5ad`, or its
        cell_type_scanvi value has no entry in KAMATH_TO_SILETTI_ALLEN.
    """
    print(f"[INFO] Loading v3.1.0 reference Kamath barcodes from {v3_reference_h5ad}")
    v3_adata = ad.read_h5ad(v3_reference_h5ad, backed="r")
    v3_kamath_barcodes = v3_adata.obs_names[v3_adata.obs["source"] == "kamath"]
    print(f"  {len(v3_kamath_barcodes):,} v3.1.0 Kamath barcode(s)")

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
    print(f"  {kamath.n_obs:,} cells x {kamath.n_vars:,} genes")
    return kamath


def build_siletti_block(v3_reference_h5ad: str, siletti_preqc_h5ad: str) -> ad.AnnData:
    """
    Load raw-UMI Siletti counts for the v3.1.0 reference's Siletti cells.

    Parameters
    ----------
    v3_reference_h5ad : str
        Path to the v3.1.0-era deconvolution reference H5AD.
    siletti_preqc_h5ad : str
        Path to Siletti's preQC intermediate H5AD (full gene panel, raw counts,
        native supercluster_term annotation).

    Returns
    -------
    ad.AnnData
        Siletti cells, X = raw UMI counts, obs[_CELL_TYPE_COLUMN] = native
        supercluster_term, obs["source"] = "siletti".

    Raises
    ------
    KeyError
        A v3.1.0 Siletti barcode is absent from `siletti_preqc_h5ad`.
    """
    print(f"[INFO] Loading v3.1.0 reference Siletti barcodes from {v3_reference_h5ad}")
    v3_adata = ad.read_h5ad(v3_reference_h5ad, backed="r")
    v3_siletti_barcodes = v3_adata.obs_names[v3_adata.obs["source"] == "siletti"]
    print(f"  {len(v3_siletti_barcodes):,} v3.1.0 Siletti barcode(s)")

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
        C_scANVI, obs["source"] = "asap_cohort", obs["selection_reason"]
        carried through from `selected`.

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

    # var_names here are mostly gene symbols (Ensembl ID only as a fallback
    # when no symbol exists) -- re-index to Ensembl IDs so the inner join
    # against Kamath/Siletti's Ensembl-ID var_names actually matches, same as
    # downsample_kamath_menon_siletti_asap.py's load_asap_cohort does.
    if "gene_id" in asap.var.columns:
        asap.var_names = asap.var["gene_id"].values
    else:
        print("  Warning: 'gene_id' not in ASAP var columns — var_names left as-is. "
              "The inner join with Kamath/Siletti will likely produce far fewer shared genes.")

    asap.obs[_CELL_TYPE_COLUMN] = selected.loc[asap.obs_names, "C_scANVI"].values
    asap.obs["selection_reason"] = selected.loc[asap.obs_names, "selection_reason"].values
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
    outdir: str,
) -> None:
    """
    Build and write the ASAP-5.0.0 counterpart reference H5AD.

    Parameters
    ----------
    v3_reference_h5ad : str
        Path to the v3.1.0-era deconvolution reference H5AD.
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
    outdir : str
        Destination directory. Output is written as `{_STEM}.h5ad`.

    Returns
    -------
    None
    """
    rng = np.random.default_rng(42)

    kamath = build_kamath_block(v3_reference_h5ad, kamath_preqc_h5ad)
    siletti = build_siletti_block(v3_reference_h5ad, siletti_preqc_h5ad)
    selected_asap = select_asap_barcodes(v3_reference_h5ad, asap_v5_final_metadata_csv, n_cells_asap, rng)
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

    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)
    out_path = outdir_path / f"{_STEM}.h5ad"
    combined.write_h5ad(out_path, compression="gzip")
    print(f"[DONE] Written {out_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--v3-reference-h5ad",
        default=_V3_REFERENCE_H5AD,
        help=f"Path to the v3.1.0-era deconvolution reference H5AD.\n[Default: {_V3_REFERENCE_H5AD}]",
    )
    parser.add_argument(
        "--kamath-preqc-h5ad",
        default=_KAMATH_PREQC_H5AD,
        help=f"Path to Kamath's preQC intermediate H5AD.\n[Default: {_KAMATH_PREQC_H5AD}]",
    )
    parser.add_argument(
        "--siletti-preqc-h5ad",
        default=_SILETTI_PREQC_H5AD,
        help=f"Path to Siletti's preQC intermediate H5AD.\n[Default: {_SILETTI_PREQC_H5AD}]",
    )
    parser.add_argument(
        "--asap-v5-unfiltered-h5ad",
        default=_ASAP_V5_UNFILTERED_H5AD,
        help=(
            "Path to v5.0.0's pre-QC-filter merged_cleaned_unfiltered H5AD.\n"
            f"[Default: {_ASAP_V5_UNFILTERED_H5AD}]"
        ),
    )
    parser.add_argument(
        "--asap-v5-final-metadata-csv",
        default=_ASAP_V5_FINAL_METADATA_CSV,
        help=f"Path to v5.0.0's final (QC-filtered) metadata CSV.\n[Default: {_ASAP_V5_FINAL_METADATA_CSV}]",
    )
    parser.add_argument(
        "--n-cells-asap",
        type=int,
        default=_N_CELLS_ASAP_DEFAULT,
        help=f"Target total ASAP cell count.\n[Default: {_N_CELLS_ASAP_DEFAULT}]",
    )
    parser.add_argument(
        "--outdir",
        default=_OUTDIR,
        help=f"Destination directory. Output is written as {_STEM}.h5ad.\n[Default: {_OUTDIR}]",
    )
    args = parser.parse_args()

    log_run_command(Path(args.outdir) / f"{_STEM}.log")
    main(
        args.v3_reference_h5ad,
        args.kamath_preqc_h5ad,
        args.siletti_preqc_h5ad,
        args.asap_v5_unfiltered_h5ad,
        args.asap_v5_final_metadata_csv,
        args.n_cells_asap,
        args.outdir,
    )
