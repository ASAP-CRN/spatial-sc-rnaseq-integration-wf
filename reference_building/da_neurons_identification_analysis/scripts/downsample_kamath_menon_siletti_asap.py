#!/usr/bin/env python3
"""
Build and cache each source's preQC and postQC intermediate AnnData for each
of the following sources:
- Kamath 2022
- Menon 2026
- Siletti 2023
- ASAP v3.1.0/2025-12-01T22-36-30Z

Phase 1 of a two-script pipeline (CPU-only: scrublet, filtering,
downsampling). Run this script first, then integrate_kamath_menon_siletti_asap.py
(GPU-bound scVI, scANVI) with the SAME --kamath-non-da-per-type /
--kamath-all-random / --downsample-kamath / --downsample-menon / --siletti-scale /
--downsample-asap-cohort values, so it locates the caches this script wrote.

- Kamath et al. (Nat Neurosci 2022): 22,048 DA neurons (10 subtypes) + non-DA
  cell types from SNpc.
- Menon et al. (bioRxiv 2026): ~4.5M nuclei across 5 regions; filtered here to
  substantia nigra only.
- Siletti et al. (Science 2023): whole-brain atlas; only the 3 midbrain
  dissections used. Its own 'cell_type' has no dopaminergic label (see
  load_siletti — DA VGLUT2 comes from an Allen Brain Cell Atlas cross-ref).
- ASAP: PMDBS cohort, sourced from sc-rnaseq-wf's unfiltered output.

Every source is downsampled BEFORE QC (so scrublet only touches the target
scale), producing subsampled pair (preQC, postQC) that have the same cells
except for cells filtered out by QC. See each load_*() docstring for source-specific 
sampling logic. ASAP is the one exception where preQC isn't simply "unfiltered"
— see load_asap_cohort.

QC (mt%, doublet score, total counts, genes detected and UMI counts — see
_apply_qc_filters) matches sc-rnaseq-wf's cohort_analysis `filter` task
defaults unless otherwise specified.

Inputs
------
1) H5AD's for each source.
   Get Kamath, Menon, Siletti H5AD files by running download_kamath_menon_siletti.py.
   ASAP H5AD files need to be manually downloaded.
2) Requires Siletti re-annotation manually-downloaded Allen Brain Cell Atlas files at
   reference_building/da_neurons_identification_analysis/reference_data/allen-brain-cell-atlas/metadata/ (see that directory's README).

Outputs
-------
reference_building/da_neurons_identification_analysis/reference_data/
├── downsample_kamath_menon_siletti_asap.log
├── {kamath,menon,siletti}/
│   └── intermediate/
│       └── *_{preQC,postQC}.h5ad
├── asap_cohort/
│   └── .../intermediate/
│       └── *_{preQC,postQC}.h5ad
└── <dataset>/
    └── .../intermediate/
        └── *_qc_summary.tsv
"""

import os
import re
import glob
import tempfile
import argparse
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "utils"))
from _reference_building_utils import (
    _REFERENCE_DATA_DIR,
    _SILETTI_SCALE_PRESETS,
    _SILETTI_DEFAULT_SCALE,
    _ASAP_H5AD_PATH,
    _ASAP_SAMPLE_METADATA_PATH,
    _ASAP_FINAL_H5AD_PATH,
    _DEFAULT_KAMATH_NON_DA_PER_TYPE,
    _DEFAULT_DOWNSAMPLE_MENON,
    _DEFAULT_DOWNSAMPLE_ASAP_COHORT,
    _qc_cache_paths,
    _qc_summary_path,
    _sn_cache_path,
    _find_menon_h5ad,
    extract_raw_counts,
    log_run_command,
)

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message=".*copying X to .raw.*")

# The 3 midbrain dissections carrying ~81% of Siletti's Allen-annotated "DA
# VGLUT2" cells (SN-RN 58%, SN 19%, RN 4%). label -> dissection file stem.
_SILETTI_MIDBRAIN_REGIONS = {
    "SN-RN": "dissection_midbrain_m_substantia_nigra_red_nucleus_and_nearby_nuclei_sn_rn",
    "SN": "dissection_midbrain_m_substantia_nigra_sn",
    "RN": "dissection_midbrain_rn_red_nucleus_rn",
}

# Siletti's whole-atlas file, used only to identify which barcodes are neurons.
_SILETTI_ALL_NEURONS_FILE = "all_neurons"

# Allen Brain Cell Atlas annotation of Siletti's dataset include DA cells
# Siletti's own obs deposited at CELLxGENE has no such label). Not auto-downloaded — see metadata/README.md.
_ALLEN_METADATA_DIR = os.path.join(_REFERENCE_DATA_DIR, "allen-brain-cell-atlas", "metadata")
_ALLEN_CELL_METADATA_PATH = os.path.join(_ALLEN_METADATA_DIR, "WHB-10Xv3", "20241115", "cell_metadata.csv.gz")
_ALLEN_MEMBERSHIP_PATH = os.path.join(
    _ALLEN_METADATA_DIR, "WHB-taxonomy", "20240330", "cluster_to_cluster_annotation_membership.csv.gz"
)

# CELLxGENE 'disease' values not in this map (e.g. Kamath's Lewy body
# dementia) are harmonized to 'Other', not dropped.
_CELLXGENE_DISEASE_TO_CONDITION_ID = {
    "normal": "Control",
    "Parkinson disease": "PD",
}

_SN_TISSUE_TERMS = {
    "substantia nigra",
    "substantia nigra pars compacta",
    "UBERON:0002038",
    "UBERON:0002450",
}

# ASAP sc-rnaseq-wf's cohort_analysis QC filter defaults, applied here too so all
# four sources go through comparable QC.
_QC_PCT_COUNTS_MT_MAX = 10
_QC_DOUBLET_SCORE_MAX = 0.2
# Raised from sc-rnaseq-wf's ASAP-calibrated 100,000: real Kamath DA neurons
# reach up to 317,511 UMIs and aren't disproportionately doublets.
_QC_TOTAL_COUNTS_LIMITS = (100, 150_000)
# Raised from sc-rnaseq-wf's ASAP-calibrated 10,000: real Kamath DA neurons
# reach up to 13,368 genes and aren't disproportionately doublets.
_QC_N_GENES_LIMITS = (100, 20_000)

# The 13 mitochondrial protein-coding genes (GRCh38 Ensembl IDs) — used
# uniformly across all sources since Kamath has no gene-symbol var column.
_MT_GENE_ENSEMBL_IDS = {
    "ENSG00000198888",  # MT-ND1
    "ENSG00000198763",  # MT-ND2
    "ENSG00000198804",  # MT-CO1
    "ENSG00000198712",  # MT-CO2
    "ENSG00000228253",  # MT-ATP8
    "ENSG00000198899",  # MT-ATP6
    "ENSG00000198938",  # MT-CO3
    "ENSG00000198840",  # MT-ND3
    "ENSG00000212907",  # MT-ND4L
    "ENSG00000198886",  # MT-ND4
    "ENSG00000198786",  # MT-ND5
    "ENSG00000198695",  # MT-ND6
    "ENSG00000198727",  # MT-CYB
}

# ASAP's QC loss on the unfiltered file is far higher than the other three
# (~52% vs <3.5%) — the raw draw is oversampled by this factor to compensate
# (2.1x calibrated against ~47.6% observed postQC survival), then topped up
# with a controlled quota of doublet-only-flagged cells (see load_asap_cohort).
_ASAP_OVERSAMPLE_FACTOR = 2.1
_ASAP_DOUBLET_REINTRODUCE_FRAC = 0.035


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Accumulates rows from _apply_qc_filters since last flush; each load_*()
# flushes to its own qc_summary.tsv rather than one combined file.
_QC_SUMMARY_ROWS: list[dict] = []


def _write_qc_summary(path: str) -> None:
    """Write accumulated QC filter stats to a TSV, creating its parent directory if needed, then clear the buffer."""
    if not _QC_SUMMARY_ROWS:
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    pd.DataFrame(_QC_SUMMARY_ROWS).to_csv(path, sep="\t", index=False)
    print(f"  Wrote QC summary ({len(_QC_SUMMARY_ROWS)} rows) → {path}")
    _QC_SUMMARY_ROWS.clear()


def _log_qc_pair(
    source: str,
    preqc: ad.AnnData,
    postqc: ad.AnnData,
    preqc_note: str = "",
    postqc_note: str = "",
) -> None:
    """
    Print a single combined preQC / postQC summary line.

    Parameters
    ----------
    source : str
        Human-readable source name (e.g. "Kamath", "Menon").
    preqc : ad.AnnData
        preQC AnnData — may be backed or in-memory.
    postqc : ad.AnnData
        postQC AnnData — may be backed or in-memory.
    preqc_note : str
        Optional parenthetical after the preQC gene count.
    postqc_note : str
        Optional parenthetical after the postQC gene count.
    """
    pre_note  = f" ({preqc_note})"  if preqc_note  else ""
    post_note = f" ({postqc_note})" if postqc_note else ""
    print(
        f"  preQC {preqc.n_obs:,} cells x {preqc.n_vars:,} genes{pre_note}"
        f" / postQC {postqc.n_obs:,} cells x {postqc.n_vars:,} genes{post_note}"
    )


def _log_cached_qc(source: str, tag: str, preqc: ad.AnnData, postqc: ad.AnnData) -> None:
    """
    Print the standard cache-hit block for a source's preQC/postQC pair.

    Parameters
    ----------
    source : str
        Human-readable source name (e.g. "Kamath", "Menon").
    tag : str
        Parameter tag shown in the header (e.g. "rand_100000", "scale_100000").
    preqc : ad.AnnData
        preQC AnnData — may be backed or in-memory.
    postqc : ad.AnnData
        postQC AnnData — may be backed or in-memory.
    """
    print(f"  Loading cached {source} preQC/postQC ({tag})")
    _log_qc_pair(source, preqc, postqc)
    print(preqc.obs.head())
    print(postqc.obs.head())


def _log_input_file(path: str) -> None:
    """Print a Found/Missing line for a raw input H5AD."""
    if os.path.exists(path):
        size_mb = os.path.getsize(path) / 1e6
        print(f"  Found ({size_mb:.0f} MB): {path}")
    else:
        print(f"  Missing: {path}")


def _apply_qc_filters(
    adata: ad.AnnData,
    label: str,
    pct_counts_mt_max: float = _QC_PCT_COUNTS_MT_MAX,
    doublet_score_max: float = _QC_DOUBLET_SCORE_MAX,
    total_counts_limits: tuple = _QC_TOTAL_COUNTS_LIMITS,
    n_genes_limits: tuple = _QC_N_GENES_LIMITS,
    doublet_batch_key: str | None = None,
) -> ad.AnnData:
    """
    Apply sc-rnaseq-wf's cohort_analysis QC filters (mt%, doublet, total
    counts, genes detected) to `adata`. Requires layers['counts'] to be set.

    Doublet scoring runs per `doublet_batch_key` (donor_id/sample_id) via a
    manual per-batch scrublet loop, not scanpy's own `batch_key=` — that
    crashes on real data (it filters near-empty cells per-batch internally,
    then does a strict reindex against the original obs_names → KeyError).
    Matches the physical reality that doublets only form within one
    droplet-generation run, and is far lower peak memory than one joint run.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData with layers['counts'] set.
    label : str
        Label for log lines (e.g. filename or region name).
    doublet_batch_key : str or None
        obs column for the finest sample/donor/library unit. None runs
        scrublet jointly (only safe for small inputs).

    Returns
    -------
    ad.AnnData
        Filtered AnnData (independent copy).
    """
    counts = adata.layers["counts"]
    total_counts = np.asarray(counts.sum(axis=1)).ravel()
    n_genes = np.asarray((counts > 0).sum(axis=1)).ravel()

    mt_present = adata.var_names.isin(_MT_GENE_ENSEMBL_IDS)
    if mt_present.any():
        mt_counts = np.asarray(counts[:, mt_present].sum(axis=1)).ravel()
        pct_counts_mt = np.where(total_counts > 0, mt_counts / total_counts * 100, 0.0)
    else:
        print(
            f"  Warning: none of the {len(_MT_GENE_ENSEMBL_IDS)} MT gene IDs found "
            f"in {label} var_names — pct_counts_mt filter skipped (treated as 0%)."
        )
        pct_counts_mt = np.zeros(adata.n_obs)

    # Restrict to a reduced HVG panel before scrublet — the recommended way to
    # run it, and avoids OOM on large files (full-panel run OOM-killed on
    # Kamath's largest file).
    n_top_genes = min(2000, adata.n_vars)
    try:
        hvg = sc.pp.highly_variable_genes(
            adata, layer="counts", flavor="seurat_v3", n_top_genes=n_top_genes, inplace=False,
        )["highly_variable"].values
    except ValueError:
        # seurat_v3's LOESS fit can fail on very small (post-downsample) cell
        # counts; fall back to the full panel since OOM isn't a concern there.
        print(f"  Warning: HVG selection failed (too few cells for a stable seurat_v3 fit) "
              f"in {label} — running scrublet on the full {adata.n_vars:,}-gene panel instead.")
        hvg = np.ones(adata.n_vars, dtype=bool)
    adata_hvg = adata[:, hvg].copy()
    adata_hvg.X = adata_hvg.layers["counts"]
    if doublet_batch_key is not None and doublet_batch_key not in adata_hvg.obs.columns:
        print(f"  Warning: doublet_batch_key={doublet_batch_key!r} not found in {label} obs "
              f"— running scrublet jointly instead of per-batch.")
        doublet_batch_key = None

    if doublet_batch_key is not None:
        # Batches <5 cells skipped (too small for doublet simulation).
        # scrublet's own internal min_genes filter can still shrink a batch
        # further before PCA runs, so wrap in try/except rather than
        # predicting the failure exactly; any batch it can't handle is left
        # NaN, which fails the doublet QC filter below.
        doublet_score = pd.Series(np.nan, index=adata_hvg.obs_names)
        n_skipped_small = 0
        n_failed = 0
        for batch in np.unique(adata_hvg.obs[doublet_batch_key]):
            batch_mask = (adata_hvg.obs[doublet_batch_key] == batch).values
            batch_size = int(batch_mask.sum())
            if batch_size < 5:
                n_skipped_small += 1
                continue
            batch_n_prin_comps = min(30, max(1, batch_size - 1))
            adata_batch = adata_hvg[batch_mask].copy()
            try:
                sc.pp.scrublet(adata_batch, n_prin_comps=batch_n_prin_comps, verbose=False)
                doublet_score.loc[adata_batch.obs_names] = adata_batch.obs["doublet_score"].values
            except ValueError:
                n_failed += 1
        n_dropped = int(doublet_score.isna().sum())
        if n_dropped:
            print(f"  Note: {n_dropped} {label} cells had no doublet_score "
                  f"({n_skipped_small} batches too small to attempt, {n_failed} batches "
                  f"failed scrublet internally, plus scrublet's own min_genes dropouts) "
                  f"— treated as failing the doublet QC filter.")
        adata.obs["doublet_score"] = doublet_score.values
    else:
        sc.pp.scrublet(adata_hvg, verbose=False)
        adata.obs["doublet_score"] = adata_hvg.obs["doublet_score"].values
    doublet_score = adata.obs["doublet_score"].values

    pass_mt = pct_counts_mt <= pct_counts_mt_max
    pass_doublet = doublet_score <= doublet_score_max
    pass_total_counts = (total_counts >= total_counts_limits[0]) & (total_counts <= total_counts_limits[1])
    pass_n_genes = (n_genes >= n_genes_limits[0]) & (n_genes <= n_genes_limits[1])
    keep = pass_mt & pass_doublet & pass_total_counts & pass_n_genes

    print(
        f"  QC filter ({label}): {adata.n_obs:,} -> {int(keep.sum()):,} cells "
        f"(fail mt%>{pct_counts_mt_max}: {int((~pass_mt).sum()):,}, "
        f"fail doublet>{doublet_score_max}: {int((~pass_doublet).sum()):,}, "
        f"fail counts outside {total_counts_limits}: {int((~pass_total_counts).sum()):,}, "
        f"fail genes outside {n_genes_limits}: {int((~pass_n_genes).sum()):,}; "
        f"failure categories can overlap)"
    )
    _QC_SUMMARY_ROWS.append({
        "label": label,
        "n_before": adata.n_obs,
        "n_after": int(keep.sum()),
        "n_removed": int((~keep).sum()),
        "fail_pct_counts_mt": int((~pass_mt).sum()),
        "fail_doublet_score": int((~pass_doublet).sum()),
        "fail_total_counts": int((~pass_total_counts).sum()),
        "fail_n_genes": int((~pass_n_genes).sum()),
        "pct_counts_mt_max": pct_counts_mt_max,
        "doublet_score_max": doublet_score_max,
        "total_counts_limits": str(total_counts_limits),
        "n_genes_limits": str(n_genes_limits),
    })
    return adata[keep].copy()

def subsample_adata(adata: ad.AnnData, n: int, label: str) -> ad.AnnData:
    """
    Randomly subsample to `n` cells (fixed seed). n=0 or n>=n_obs: no-op.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData to subsample.
    n : int
        Target cell count.
    label : str
        Label for progress messages.

    Returns
    -------
    ad.AnnData
        Subsampled copy, or the original if n is out of range.

    Raises
    ------
    ValueError
        If n is negative.
    """
    if n < 0:
        raise ValueError(f"{label}: requested {n:,} < 0 — must be >= 0")
    if n == 0:
        print(f"  {label}: {adata.n_obs:,} cells — no downsampling applied")
        return adata
    if n >= adata.n_obs:
        print(f"  {label}: requested {n:,} >= actual {adata.n_obs:,} — no downsampling applied")
        return adata
    rng = np.random.default_rng(42)
    idx = rng.choice(adata.n_obs, size=n, replace=False)
    adata = adata[np.sort(idx)].copy()
    print(f"  {label}: downsampled to {adata.n_obs:,} cells")
    return adata

def _write_h5ad_cache(adata: ad.AnnData, path: str) -> None:
    """Write `adata` to `path`, creating its parent directory if needed."""
    print(adata.obs.head())
    os.makedirs(os.path.dirname(path), exist_ok=True)
    adata.write_h5ad(path)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _harmonize_cellxgene_disease(disease: pd.Series) -> pd.Series:
    """Map a CELLxGENE 'disease' column to condition_id ("PD"/"Control"/"Other")."""
    # astype(str) first: .map() on a categorical Series would otherwise keep
    # only mapped values as categories, breaking a later fillna("Other").
    return disease.astype(str).map(_CELLXGENE_DISEASE_TO_CONDITION_ID).fillna("Other")


def _annotate_kamath_obs(adata: ad.AnnData, basename: str, da_subtype_column: str) -> None:
    """
    Set standard obs columns on a Kamath AnnData in place.

    Parameters
    ----------
    adata : ad.AnnData
        Kamath AnnData (in memory, not backed).
    basename : str
        Filename, used in error messages.
    da_subtype_column : str
        obs column with fine-grained DA subtype labels.
    """
    adata.layers["counts"] = extract_raw_counts(adata)
    adata.obs["source"] = "kamath"
    if "donor_id" not in adata.obs.columns:
        raise KeyError(f"'donor_id' column not found in {basename}.")
    adata.obs["batch"] = "kamath_" + adata.obs["donor_id"].astype(str)
    if da_subtype_column in adata.obs.columns:
        adata.obs["cell_type_scanvi"] = adata.obs[da_subtype_column].astype(str)
    else:
        adata.obs["cell_type_scanvi"] = adata.obs["cell_type"].astype(str)
    is_da = adata.obs["cell_type"].str.lower().str.contains("dopaminergic")
    if da_subtype_column in adata.obs.columns and is_da.any():
        adata.obs["da_subtype"] = adata.obs[da_subtype_column].where(is_da, other=np.nan)
    else:
        adata.obs["da_subtype"] = np.nan
    if "disease" not in adata.obs.columns:
        raise KeyError(f"'disease' column not found in {basename}.")
    adata.obs["condition_id"] = _harmonize_cellxgene_disease(adata.obs["disease"])


def load_kamath(
    kamath_h5ad_dir: str,
    exclude_files: list[str],
    da_subtype_column: str,
    n_non_da_per_type: int = 0,
    all_random: bool = False,
    n_total: int = 0,
    qc_pct_counts_mt_max: float = _QC_PCT_COUNTS_MT_MAX,
    qc_doublet_score_max: float = _QC_DOUBLET_SCORE_MAX,
    qc_total_counts_limits: tuple = _QC_TOTAL_COUNTS_LIMITS,
    qc_n_genes_limits: tuple = _QC_N_GENES_LIMITS,
) -> tuple[ad.AnnData, ad.AnnData]:
    """
    Load and merge Kamath cell-type H5AD files into (preQC, postQC).

    Default mode: each file is downsampled BEFORE QC (DA neurons kept in full;
    other cell types capped to `n_non_da_per_type`), then QC'd and concatenated.

    When `all_random=True`: all non-excluded files are loaded in full in backed mode,
    concatenated, randomly subsampled to `n_total` cells (BEFORE QC), then QC'd once 
    on the small pool. `n_non_da_per_type` is ignored.

    Both Nurr1-sorted files must be in `exclude_files` — barcode overlap
    confirmed their cells already appear in the 7 cell-type files (see
    kamath_overlap.py); including them would duplicate ~184k cells.

    Parameters
    ----------
    kamath_h5ad_dir : str
        Directory containing the Kamath CELLxGENE H5ADs.
    exclude_files : list[str]
        Filenames (basenames) to skip — must include both Nurr files.
    da_subtype_column : str
        obs column with fine-grained DA subtype labels.
    n_non_da_per_type : int
        Max cells per non-DA file (0 = no subsampling). DA neurons file
        always kept in full. Ignored when `all_random=True`.
    all_random : bool
        If True, pool all non-excluded files and draw `n_total` cells at random
        (no per-file cap, no DA-first priority). `n_total` must be > 0.
    n_total : int
        Total cells to draw when `all_random=True`. Must be > 0 in that mode.

    Returns
    -------
    tuple[ad.AnnData, ad.AnnData]
        (preQC, postQC) concatenated Kamath AnnDatas.
    """
    paths = sorted(glob.glob(os.path.join(kamath_h5ad_dir, "*.h5ad")))
    if not paths:
        raise FileNotFoundError(
            f"No H5AD files found in {kamath_h5ad_dir}. "
            "Run download_kamath_menon_siletti.py first."
        )

    exclude_basenames = {os.path.basename(f) for f in exclude_files}

    def _valid_paths() -> list[tuple[str, str]]:
        return [
            (p, os.path.basename(p)) for p in paths
            if os.path.basename(p) not in exclude_basenames
            and not os.path.basename(p).startswith("kamath_merged")
        ]

    if all_random:
        # ── Step 1: backed obs scan — collect (path, basename, n_obs) ──────────
        valid = []
        for path, basename in _valid_paths():
            backed = sc.read_h5ad(path, backed="r")
            n = backed.n_obs
            backed.file.close()
            valid.append((path, basename, n))
            print(f"  Backed scan: {basename} ({n:,} cells)")

        if not valid:
            raise RuntimeError("No Kamath files found after exclusions.")

        # ── Step 2: sample n_total flat indices across the concatenated pool ───
        total_available = sum(n for _, _, n in valid)
        n_draw = min(n_total, total_available)
        print(f"  Global pool: {total_available:,} cells → sampling {n_draw:,}")
        rng = np.random.default_rng(0)
        flat_idx = np.sort(rng.choice(total_available, size=n_draw, replace=False))

        # Map flat indices back to per-file integer row positions.
        indices_by_path: dict[str, tuple[str, np.ndarray]] = {}
        offset = 0
        for path, basename, n in valid:
            in_file = flat_idx[(flat_idx >= offset) & (flat_idx < offset + n)] - offset
            if len(in_file) > 0:
                indices_by_path[path] = (basename, in_file)
            offset += n

        # ── Step 3: load only sampled rows from each file ─────────────────────
        subsets = []
        for path, (basename, file_idx) in indices_by_path.items():
            print(f"  Loading {len(file_idx):,} sampled cells from: {basename}")
            backed = sc.read_h5ad(path, backed="r")
            chunk = backed[np.sort(file_idx)].to_memory()
            backed.file.close()
            _annotate_kamath_obs(chunk, basename, da_subtype_column)
            print(f"    cell types: {sorted(chunk.obs['cell_type_scanvi'].unique())}")
            subsets.append(chunk)
            del chunk

        # ── Step 4: concat the small subsets, then QC once ───────────────────
        kamath_preqc = ad.concat(subsets, join="inner")
        del subsets
        kamath_preqc.obs["original_barcode"] = kamath_preqc.obs_names.astype(str)
        kamath_preqc.obs_names = "kamath_" + kamath_preqc.obs_names.astype(str)
        kamath_postqc = _apply_qc_filters(
            kamath_preqc, "Kamath all random",
            pct_counts_mt_max=qc_pct_counts_mt_max, doublet_score_max=qc_doublet_score_max,
            total_counts_limits=qc_total_counts_limits, n_genes_limits=qc_n_genes_limits,
            doublet_batch_key="donor_id",
        )
        _log_qc_pair("Kamath", kamath_preqc, kamath_postqc)
        tag = str(n_total)

    else:
        # ── Default mode: per-file subsample + QC, temp files keep peak RAM low ─
        # "da_neurons" is a substring of "non_da_neurons" — match exact filename.
        tmp_dir = tempfile.mkdtemp(prefix="kamath_qc_")
        tmp_paths_preqc: list[str] = []
        tmp_paths_postqc: list[str] = []

        for path, basename in _valid_paths():
            print(f"  Loading: {basename}")
            adata = sc.read_h5ad(path)
            _annotate_kamath_obs(adata, basename, da_subtype_column)
            print(
                f"    {adata.n_obs:,} cells | "
                f"cell types: {sorted(adata.obs['cell_type_scanvi'].unique())}"
            )

            preqc_path = os.path.join(tmp_dir, f"preQC_{basename}")
            adata.write_h5ad(preqc_path)
            tmp_paths_preqc.append(preqc_path)

            is_da_file = basename == "human_da_neurons_10x_scrna_seq.h5ad"
            if is_da_file:
                print(f"    DA neurons file — keeping all {adata.n_obs:,} cells")
            elif n_non_da_per_type > 0:
                adata = subsample_adata(adata, n_non_da_per_type, f"Kamath non-DA ({basename})")

            adata = _apply_qc_filters(
                adata, basename,
                pct_counts_mt_max=qc_pct_counts_mt_max, doublet_score_max=qc_doublet_score_max,
                total_counts_limits=qc_total_counts_limits, n_genes_limits=qc_n_genes_limits,
                doublet_batch_key="donor_id",
            )
            postqc_path = os.path.join(tmp_dir, f"postQC_{basename}")
            adata.write_h5ad(postqc_path)
            tmp_paths_postqc.append(postqc_path)
            del adata

        if not tmp_paths_preqc:
            raise RuntimeError("No Kamath files were loaded. Check directory and exclude list.")

        print(f"  Concatenating {len(tmp_paths_preqc)} Kamath files (preQC)")
        kamath_preqc = ad.concat((sc.read_h5ad(p) for p in tmp_paths_preqc), join="inner")
        kamath_preqc.obs["original_barcode"] = kamath_preqc.obs_names.astype(str)
        kamath_preqc.obs_names = "kamath_" + kamath_preqc.obs_names.astype(str)
        print(f"  Concatenating {len(tmp_paths_postqc)} Kamath files (postQC)")
        kamath_postqc = ad.concat((sc.read_h5ad(p) for p in tmp_paths_postqc), join="inner")
        kamath_postqc.obs["original_barcode"] = kamath_postqc.obs_names.astype(str)
        kamath_postqc.obs_names = "kamath_" + kamath_postqc.obs_names.astype(str)
        _log_qc_pair("Kamath", kamath_preqc, kamath_postqc)

        for p in tmp_paths_preqc + tmp_paths_postqc:
            os.remove(p)
        os.rmdir(tmp_dir)
        tag = str(n_non_da_per_type)

    anchor_name = "kamath_merged_all_rand" if all_random else "kamath_merged_da_all_non_da"
    kamath_anchor = os.path.join(kamath_h5ad_dir, "intermediate", anchor_name)
    _write_qc_summary(_qc_summary_path(kamath_anchor, tag))

    return kamath_preqc, kamath_postqc

def load_menon_sn(
    menon_h5ad: str,
    n_cells: int = 0,
    use_cached: bool = False,
    qc_pct_counts_mt_max: float = _QC_PCT_COUNTS_MT_MAX,
    qc_doublet_score_max: float = _QC_DOUBLET_SCORE_MAX,
    qc_total_counts_limits: tuple = _QC_TOTAL_COUNTS_LIMITS,
    qc_n_genes_limits: tuple = _QC_N_GENES_LIMITS,
) -> tuple[ad.AnnData, ad.AnnData]:
    """
    Load Menon's full H5AD, filter to SN cells, downsample, then QC.

    Native labels are preserved under renamed columns ('menon_cell_type',
    'menon_author_cell_type'); cell_type_scanvi is 'Unknown' so scANVI label
    transfer runs from Kamath to Menon.

    Opened in backed mode so only obs loads into RAM. The SN filter and
    the `n_cells` downsample are both applied to the barcode index while
    still backed, so `.to_memory()` loads only the target rows — peak RAM
    is `n_cells` rather than all SN cells. When `use_cached=True` and both
    caches exist, loads them directly, skipping the 4.5M-cell read.

    Parameters
    ----------
    menon_h5ad : str
        Path to the Menon CELLxGENE H5AD (4.5M cells).
    n_cells : int
        Cells to downsample the SN pool to, before QC. 0 = no downsampling.
    use_cached : bool
        If True and both caches exist, load them directly.

    Returns
    -------
    tuple[ad.AnnData, ad.AnnData]
        (preQC, postQC) Menon AnnDatas.

    Raises
    ------
    KeyError
        Required obs column absent.
    ValueError
        No SN cells found, or non-integer counts.
    """
    tag = f"rand_{n_cells}" if n_cells > 0 else "all"
    preqc_cache, postqc_cache = _qc_cache_paths(_sn_cache_path(menon_h5ad), tag)
    if use_cached and os.path.exists(preqc_cache) and os.path.exists(postqc_cache):
        menon_preqc = sc.read_h5ad(preqc_cache)
        menon_postqc = sc.read_h5ad(postqc_cache)
        _log_cached_qc("Menon", tag, menon_preqc, menon_postqc)
        return menon_preqc, menon_postqc

    print(f"  Loading Menon H5AD in backed mode (obs only, matrix stays on disk): {menon_h5ad}")
    menon = sc.read_h5ad(menon_h5ad, backed='r')
    print(f"  Opened (backed): {menon.n_obs:,} cells")

    tissue_col = None
    for col in ("tissue", "tissue_ontology_term_id"):
        if col in menon.obs.columns:
            tissue_col = col
            break
    if tissue_col is None:
        raise KeyError("No 'tissue' or 'tissue_ontology_term_id' column found in Menon obs.")

    sn_mask = menon.obs[tissue_col].isin(_SN_TISSUE_TERMS)
    if not sn_mask.any():
        available = sorted(menon.obs[tissue_col].unique())
        raise ValueError(
            f"No substantia nigra cells found using column '{tissue_col}'. "
            f"Available values: {available}"
        )
    # Sample barcodes while still backed (X stays on disk).
    # Combining the SN filter and the downsample here means .to_memory() loads
    # only the target rows — peak RAM is n_cells rather than all SN cells.
    sn_idx = np.where(sn_mask)[0]
    print(f"  SN barcodes (backed): {len(sn_idx):,}")
    if n_cells > 0 and n_cells < len(sn_idx):
        sn_idx = np.sort(np.random.default_rng(0).choice(sn_idx, size=n_cells, replace=False))
        print(f"  Sampled {len(sn_idx):,} barcodes before to_memory()")
    menon = menon[sn_idx].to_memory()
    print(f"  Loaded into memory: {menon.n_obs:,} cells")

    menon.layers["counts"] = extract_raw_counts(menon)

    # batch encodes donor_id + region for per-library resolution (region is
    # always 'SN' post-filter, but future-proofs the key if more are added).
    menon.obs["source"] = "menon"
    for col in ("donor_id", "region"):
        if col not in menon.obs.columns:
            raise KeyError(f"'{col}' column not found in Menon obs.")
    menon.obs["batch"] = (
        "menon_" + menon.obs["donor_id"].astype(str)
        + "_" + menon.obs["region"].astype(str)
    )

    for col, dest in (("cell_type", "menon_cell_type"),
                      ("author_cell_type", "menon_author_cell_type")):
        if col not in menon.obs.columns:
            raise KeyError(
                f"Required column '{col}' not found in Menon obs. "
                f"Available columns: {menon.obs.columns.tolist()}"
            )
        menon.obs[dest] = menon.obs[col].astype(str)
        print(f"  Preserved {dest}: {menon.obs[dest].nunique()} unique labels")
    print(f"    Top menon_cell_type labels: "
          f"{menon.obs['menon_cell_type'].value_counts().head(10).to_dict()}")

    if "disease" not in menon.obs.columns:
        raise KeyError("'disease' column not found in Menon obs.")
    menon.obs["condition_id"] = _harmonize_cellxgene_disease(menon.obs["disease"])
    print(f"  condition_id: {menon.obs['condition_id'].value_counts().to_dict()}")

    menon.obs["cell_type_scanvi"] = "Unknown"
    menon.obs["da_subtype"] = np.nan

    menon.obs["original_barcode"] = menon.obs_names.astype(str)
    menon.obs_names = "menon_" + menon.obs_names.astype(str)
    print(f"  Writing Menon preQC cache → {preqc_cache}")
    _write_h5ad_cache(menon, preqc_cache)

    if "donor_id" not in menon.obs.columns:
        raise KeyError("'donor_id' column not found in Menon obs.")
    menon_postqc = _apply_qc_filters(
        menon, "Menon SN",
        pct_counts_mt_max=qc_pct_counts_mt_max, doublet_score_max=qc_doublet_score_max,
        total_counts_limits=qc_total_counts_limits, n_genes_limits=qc_n_genes_limits,
        doublet_batch_key="donor_id",
    )
    _log_qc_pair("Menon", menon, menon_postqc)
    print(f"  Writing Menon postQC cache → {postqc_cache}")
    _write_h5ad_cache(menon_postqc, postqc_cache)

    _write_qc_summary(_qc_summary_path(_sn_cache_path(menon_h5ad), tag))

    return menon, menon_postqc

def _random_subset(pool: set, n: int, rng: np.random.Generator) -> set:
    """Random subset of `pool` of size `n` (all of `pool` if n >= len(pool))."""
    pool = sorted(pool)
    if n >= len(pool):
        return set(pool)
    idx = rng.choice(len(pool), size=n, replace=False)
    return {pool[i] for i in idx}

def _load_allen_annotations(
    cell_metadata_path: str = _ALLEN_CELL_METADATA_PATH,
    membership_path: str = _ALLEN_MEMBERSHIP_PATH,
) -> pd.DataFrame:
    """
    Per-cell Allen Brain Cell Atlas annotations, indexed by barcode (format
    confirmed identical to Siletti's obs_names).

    Returns
    -------
    pd.DataFrame
        Indexed by barcode: neurotransmitter, subcluster (Allen's terms).

    Raises
    ------
    FileNotFoundError
        Either Allen file missing.
    """
    for path in (cell_metadata_path, membership_path):
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"{path} not found. See reference_building/da_neurons_identification_analysis/"
                "reference_data/allen-brain-cell-atlas/metadata/README.md to fetch the Allen Brain Cell Atlas files."
            )
    # keep_default_na=False: "None" is a real neurotransmitter category here,
    # not a missing value — pandas would otherwise silently NaN it.
    membership = pd.read_csv(membership_path, keep_default_na=False)
    neurotransmitter_by_alias = membership.loc[
        membership["cluster_annotation_term_set_name"] == "neurotransmitter",
        ["cluster_alias", "cluster_annotation_term_name"],
    ].set_index("cluster_alias")["cluster_annotation_term_name"]
    subcluster_by_alias = membership.loc[
        membership["cluster_annotation_term_set_name"] == "subcluster",
        ["cluster_alias", "cluster_annotation_term_name"],
    ].set_index("cluster_alias")["cluster_annotation_term_name"]

    cell_meta = pd.read_csv(cell_metadata_path, usecols=["cell_label", "cluster_alias"])
    annotations = pd.DataFrame(index=pd.Index(cell_meta["cell_label"].values, name="barcode"))
    annotations["neurotransmitter"] = cell_meta["cluster_alias"].map(neurotransmitter_by_alias).values
    annotations["subcluster"] = cell_meta["cluster_alias"].map(subcluster_by_alias).values
    print(f"  Loaded Allen annotations for {len(annotations):,} barcodes")
    return annotations

def _load_all_neuron_barcodes(siletti_dir: str) -> set:
    """Barcodes present in Siletti's whole-atlas 'all neurons' file."""
    path = os.path.join(siletti_dir, f"{_SILETTI_ALL_NEURONS_FILE}.h5ad")
    if not os.path.exists(path):
        raise FileNotFoundError(f"{path} not found. Run download_kamath_menon_siletti.py first.")
    a = sc.read_h5ad(path, backed="r")
    return set(a.obs_names)

def load_siletti(
    siletti_dir: str,
    allen_annotations: pd.DataFrame,
    neuron_barcodes: set,
    neuron_target: int = _SILETTI_SCALE_PRESETS[_SILETTI_DEFAULT_SCALE][0],
    cell_target: int = _SILETTI_SCALE_PRESETS[_SILETTI_DEFAULT_SCALE][1],
    qc_pct_counts_mt_max: float = _QC_PCT_COUNTS_MT_MAX,
    qc_doublet_score_max: float = _QC_DOUBLET_SCORE_MAX,
    qc_total_counts_limits: tuple = _QC_TOTAL_COUNTS_LIMITS,
    qc_n_genes_limits: tuple = _QC_N_GENES_LIMITS,
) -> tuple[ad.AnnData, ad.AnnData]:
    """
    Load Siletti's 3 midbrain dissections with DA-aware per-region sampling.

    Per region: every DA VGLUT2 cell is kept, topped up with random neurons
    to `neuron_target`, then random cells to `cell_target` (both capped by
    what the region has). This selection IS the preQC pool; QC then runs
    per region (a few forced-included DA cells can be removed — logged).

    Parameters
    ----------
    siletti_dir : str
        Directory containing the Siletti dissection H5ADs.
    allen_annotations : pd.DataFrame
        Per-barcode Allen annotations (see _load_allen_annotations).
    neuron_barcodes : set
        Barcodes considered neurons (see _load_all_neuron_barcodes).

    Returns
    -------
    tuple[ad.AnnData, ad.AnnData]
        (preQC, postQC) concatenated Siletti AnnDatas.

    Raises
    ------
    FileNotFoundError
        A region's dissection H5AD is missing.
    KeyError
        Required obs column absent.
    """
    da_barcodes = set(allen_annotations.index[allen_annotations["neurotransmitter"] == "DA VGLUT2"])

    rng = np.random.default_rng(42)
    tmp_dir = tempfile.mkdtemp(prefix="siletti_qc_")
    tmp_paths_preqc = []
    tmp_paths_postqc = []
    for region_label, file_stem in _SILETTI_MIDBRAIN_REGIONS.items():
        h5ad_path = os.path.join(siletti_dir, f"{file_stem}.h5ad")
        if not os.path.exists(h5ad_path):
            raise FileNotFoundError(f"{h5ad_path} not found. Run download_kamath_menon_siletti.py first.")

        backed = sc.read_h5ad(h5ad_path, backed="r")
        obs_names = set(backed.obs_names)
        region_neurons = obs_names & neuron_barcodes

        must_include = obs_names & da_barcodes
        neuron_pool = region_neurons - must_include
        extra_neurons = _random_subset(neuron_pool, max(0, neuron_target - len(must_include)), rng)

        selected = must_include | extra_neurons
        non_neuron_pool = obs_names - region_neurons
        extra_cells = _random_subset(non_neuron_pool, max(0, cell_target - len(selected)), rng)
        selected |= extra_cells

        print(
            f"  {region_label}: {len(must_include)} DA (forced) + {len(extra_neurons)} extra neurons "
            f"+ {len(extra_cells)} extra non-neuronal = {len(selected):,} cells "
            f"(region has {len(obs_names):,} total, {len(region_neurons):,} neurons)"
        )

        chunk = backed[sorted(selected)].to_memory()
        chunk.layers["counts"] = extract_raw_counts(chunk)
        chunk.obs["siletti_cxg_region"] = region_label
        chunk_allen = allen_annotations.reindex(chunk.obs_names)
        chunk.obs["siletti_allen_neurotransmitter"] = chunk_allen["neurotransmitter"].values
        chunk.obs["siletti_allen_subcluster"] = chunk_allen["subcluster"].values
        chunk.obs["siletti_allen_da_vglut2"] = chunk.obs_names.isin(da_barcodes)

        preqc_path = os.path.join(tmp_dir, f"preQC_{region_label}.h5ad")
        chunk.write_h5ad(preqc_path)
        tmp_paths_preqc.append(preqc_path)

        if "sample_id" not in chunk.obs.columns:
            raise KeyError(f"'sample_id' column not found in Siletti {region_label} obs.")
        n_da_before = int(chunk.obs["siletti_allen_da_vglut2"].sum())
        chunk = _apply_qc_filters(
            chunk, f"Siletti {region_label}",
            pct_counts_mt_max=qc_pct_counts_mt_max, doublet_score_max=qc_doublet_score_max,
            total_counts_limits=qc_total_counts_limits, n_genes_limits=qc_n_genes_limits,
            doublet_batch_key="sample_id",
        )
        n_da_after = int(chunk.obs["siletti_allen_da_vglut2"].sum())
        if n_da_after != n_da_before:
            print(
                f"    Note: {n_da_before - n_da_after} of {n_da_before} forced-included "
                f"DA VGLUT2 cells in {region_label} removed by QC filters."
            )
        postqc_path = os.path.join(tmp_dir, f"postQC_{region_label}.h5ad")
        chunk.write_h5ad(postqc_path)
        tmp_paths_postqc.append(postqc_path)
        del chunk

    print(f"  Concatenating {len(tmp_paths_preqc)} midbrain regions (preQC)")
    siletti_preqc = ad.concat((sc.read_h5ad(p) for p in tmp_paths_preqc), join="inner", merge="same")
    print(f"  Concatenating {len(tmp_paths_postqc)} midbrain regions (postQC)")
    siletti_postqc = ad.concat((sc.read_h5ad(p) for p in tmp_paths_postqc), join="inner", merge="same")
    for p in tmp_paths_preqc + tmp_paths_postqc:
        os.remove(p)
    os.rmdir(tmp_dir)

    def _finalize_siletti_obs(siletti: ad.AnnData) -> None:
        """Set source/batch/native-label/condition_id/scANVI obs columns in place."""
        siletti.obs["source"] = "siletti"
        if "sample_id" not in siletti.obs.columns:
            raise KeyError(
                f"'sample_id' column not found in Siletti obs. "
                f"Available columns: {siletti.obs.columns.tolist()}"
            )
        siletti.obs["batch"] = "siletti_" + siletti.obs["sample_id"].astype(str)

        if "cell_type" not in siletti.obs.columns:
            raise KeyError(
                f"'cell_type' column not found in Siletti obs. "
                f"Available columns: {siletti.obs.columns.tolist()}"
            )
        siletti.obs["siletti_cxg_cell_type"] = siletti.obs["cell_type"].astype(str)

        # Siletti has no 'author_cell_type'; closest analog is the 'supercluster_term'.
        if "supercluster_term" in siletti.obs.columns:
            siletti.obs["siletti_cxg_author_cell_type"] = siletti.obs["supercluster_term"].astype(str)
        else:
            siletti.obs["siletti_cxg_author_cell_type"] = np.nan

        if "disease" not in siletti.obs.columns:
            raise KeyError("'disease' column not found in Siletti obs.")
        siletti.obs["condition_id"] = _harmonize_cellxgene_disease(siletti.obs["disease"])

        siletti.obs["cell_type_scanvi"] = "Unknown"
        siletti.obs["da_subtype"] = np.nan

        siletti.obs["original_barcode"] = siletti.obs_names.astype(str)
        siletti.obs_names = "siletti_" + siletti.obs_names.astype(str)

    _finalize_siletti_obs(siletti_preqc)
    _finalize_siletti_obs(siletti_postqc)
    n_da_pre  = int(siletti_preqc.obs["siletti_allen_da_vglut2"].sum())
    n_da_post = int(siletti_postqc.obs["siletti_allen_da_vglut2"].sum())
    _log_qc_pair("Siletti", siletti_preqc, siletti_postqc,
                 f"{n_da_pre} DA VGLUT2", f"{n_da_post} DA VGLUT2")

    siletti_anchor = os.path.join(
        siletti_dir, "intermediate", f"midbrain_da_neurons{neuron_target}_cells{cell_target}",
    )
    _write_qc_summary(_qc_summary_path(siletti_anchor, ""))

    return siletti_preqc, siletti_postqc

def _normalize_asap_sample_id(sample: str) -> str:
    """Normalize an ASAP H5AD 'sample' value to SAMPLE.csv's ASAP_sample_id (fixes a 'PMBDS' typo, strips replicate suffix)."""
    sample = sample.replace("PMBDS", "PMDBS")
    return re.sub(r"_(?:rep)?\d+$", "", sample, flags=re.IGNORECASE)

def _load_asap_condition_id(sample_metadata_path: str, samples: pd.Series) -> pd.Series:
    """
    Harmonize ASAP condition_id (PD/Control/Prodromal) from SAMPLE.csv.

    A few subjects have inconsistent condition_id across their own samples
    (different ASAP teams label the same subject Control vs Prodromal
    differently) — every sample of an affected subject is harmonized to
    'Prodromal', the more specific label.

    Parameters
    ----------
    sample_metadata_path : str
        Path to the ASAP cohort's SAMPLE.csv.
    samples : pd.Series
        Raw obs['sample'] values from the ASAP H5AD.

    Returns
    -------
    pd.Series
        condition_id aligned to `samples`.

    Raises
    ------
    KeyError
        A normalized sample ID has no match in SAMPLE.csv.
    """
    meta = pd.read_csv(sample_metadata_path, keep_default_na=False)
    condition_by_sample_id = meta.set_index("ASAP_sample_id")["condition_id"].to_dict()

    subject_condition_counts = meta.groupby("ASAP_subject_id")["condition_id"].nunique()
    conflicted_subjects = set(subject_condition_counts[subject_condition_counts > 1].index)
    if conflicted_subjects:
        print(
            f"  {len(conflicted_subjects)} ASAP subjects have inconsistent "
            f"condition_id across samples — harmonized to 'Prodromal': "
            f"{sorted(conflicted_subjects)}"
        )
        conflicted_sample_ids = meta.loc[
            meta["ASAP_subject_id"].isin(conflicted_subjects), "ASAP_sample_id"
        ]
        for sample_id in conflicted_sample_ids:
            condition_by_sample_id[sample_id] = "Prodromal"

    normalized = samples.map(_normalize_asap_sample_id)
    unmatched = sorted(set(normalized) - set(condition_by_sample_id))
    if unmatched:
        raise KeyError(
            f"{len(unmatched)} normalized ASAP sample IDs not found in "
            f"SAMPLE.csv: {unmatched[:10]}"
        )
    return normalized.map(condition_by_sample_id)

def load_asap_cohort(
    asap_h5ad: str,
    asap_sample_metadata_csv: str,
    n_cells: int = 0,
    use_cached: bool = False,
    qc_pct_counts_mt_max: float = _QC_PCT_COUNTS_MT_MAX,
    qc_doublet_score_max: float = _QC_DOUBLET_SCORE_MAX,
    qc_total_counts_limits: tuple = _QC_TOTAL_COUNTS_LIMITS,
    qc_n_genes_limits: tuple = _QC_N_GENES_LIMITS,
    oversample_factor: float = _ASAP_OVERSAMPLE_FACTOR,
    doublet_reintroduce_frac: float = _ASAP_DOUBLET_REINTRODUCE_FRAC,
) -> tuple[ad.AnnData, ad.AnnData]:
    """
    Load the ASAP cohort's unfiltered H5AD, downsample (backed mode), QC.

    `asap_h5ad` must be the genuinely unfiltered merged_cleaned_unfiltered.h5ad
    (36,601 genes incl. all 17 Kamath DA markers) — not sc-rnaseq-wf's own
    HVG-subset final.h5ad (3,000 genes, 0/17 markers). Downsampling happens
    BEFORE QC so scrublet only touches the target-scale pool.

    ASAP's QC loss here is far higher than the other three sources (~52% vs
    <3.5%, driven by genuinely empty/low-complexity droplets). To land close
    to `n_cells` after QC:
      1. Raw draw is `n_cells * oversample_factor` (default 2.1x) — an
         intermediate pool, never itself cached.
      2. QC runs on that pool. The strict survivor set is trimmed to
         `n_cells - n_reintroduce` if it exceeds that quota — this becomes `postQC`.
      3. A quota (`doublet_reintroduce_frac`, default 3.5%) of cells that
         are clean except for a failed doublet flag (valid, non-NaN
         scrublet score) is drawn from the same oversampled pool and
         concatenated onto postQC, forming `preQC`. So preQC >= postQC for
         ASAP too, same relationship as the other three sources — though
         neither is the full 2.1x working pool itself, which is discarded.

    var_names are re-indexed to Ensembl IDs (from var['gene_id']) so the
    4-way inner join can match genes. This unfiltered file has no
    'C_scANVI'/'cell_type' (added later in sc-rnaseq-wf's own pipeline) —
    pmdbs_sc_rnaseq.C_scANVI/cell_type are NaN for all cells.

    Parameters
    ----------
    asap_h5ad : str
        Path to the unfiltered ASAP cohort H5AD.
    asap_sample_metadata_csv : str
        Cohort's SAMPLE.csv, joined on obs['sample'] for condition_id.
    n_cells : int
        Target cell count after QC. 0 = no downsampling/oversampling.
    use_cached : bool
        If True and both caches exist, load them directly.
    oversample_factor : float
        Multiplier on `n_cells` for the raw draw. See _ASAP_OVERSAMPLE_FACTOR.
    doublet_reintroduce_frac : float
        Fraction of `n_cells` reintroduced from doublet-only failures.

    Returns
    -------
    tuple[ad.AnnData, ad.AnnData]
        (preQC, postQC) ASAP cohort AnnDatas.

    Raises
    ------
    FileNotFoundError
        Input file missing.
    KeyError
        Required obs column absent.
    """
    for path in (asap_h5ad, asap_sample_metadata_csv):
        if not os.path.exists(path):
            raise FileNotFoundError(f"ASAP cohort file not found: {path}")

    tag = f"rand_{n_cells}" if n_cells > 0 else "all"
    preqc_cache, postqc_cache = _qc_cache_paths(asap_h5ad, tag)
    if use_cached and os.path.exists(preqc_cache) and os.path.exists(postqc_cache):
        asap_preqc = sc.read_h5ad(preqc_cache)
        asap_postqc = sc.read_h5ad(postqc_cache)
        _log_cached_qc("ASAP", tag, asap_preqc, asap_postqc)
        return asap_preqc, asap_postqc

    print(f"  Opening H5AD in backed mode: {asap_h5ad}")
    asap = sc.read_h5ad(asap_h5ad, backed='r')
    print(f"  Opened (backed): {asap.n_obs:,} cells, {asap.n_vars:,} genes")

    if n_cells > 0 and n_cells < asap.n_obs:
        draw_n = min(round(n_cells * oversample_factor), asap.n_obs)
        rng = np.random.default_rng(42)
        chosen_idx = np.sort(rng.choice(asap.n_obs, size=draw_n, replace=False))
        print(f"  Subsampled {draw_n:,} barcodes ({oversample_factor}x the {n_cells:,} "
              f"target, before QC — compensates for ASAP's high QC loss rate)")
    else:
        chosen_idx = np.arange(asap.n_obs)
        print(f"  Using all {asap.n_obs:,} barcodes")

    asap = asap[chosen_idx].to_memory()
    print(f"  Loaded into memory: {asap.n_obs:,} cells, {asap.n_vars:,} genes")

    if "gene_id" in asap.var.columns:
        asap.var_names = asap.var["gene_id"].values
        print(f"  Re-indexed var_names to Ensembl IDs (from var['gene_id'])")
    else:
        print(
            "  Warning: 'gene_id' not in var columns — var_names left as gene symbols. "
            "The inner join with Kamath/Menon/Siletti will likely produce 0 shared genes."
        )

    asap.layers["counts"] = extract_raw_counts(asap)

    # batch = 'asap_{team}_{dataset}_{sample}': the 'batch' column alone
    # (B1, B2, ...) is reused inconsistently across ASAP teams, so this is
    # the finest key that's actually globally unique.
    asap.obs["source"] = "asap_cohort"
    if all(c in asap.obs.columns for c in ("team", "dataset", "sample")):
        asap.obs["batch"] = (
            "asap_" + asap.obs["team"].astype(str)
            + "_" + asap.obs["dataset"].astype(str)
            + "_" + asap.obs["sample"].astype(str)
        )
        n_batches = asap.obs["batch"].nunique()
        print(f"  Batch key: asap_{{team}}_{{dataset}}_{{sample}} — {n_batches} unique batches")
    else:
        fallback_col = next(
            (c for c in ("donor_id", "sample", "batch_id") if c in asap.obs.columns),
            None,
        )
        if fallback_col is None:
            raise KeyError(
                "Neither 'team'+'dataset'+'sample' nor any of 'donor_id', 'sample', "
                "'batch_id' found in ASAP cohort obs. "
                f"Available columns: {asap.obs.columns.tolist()}"
            )
        print(
            f"  Warning: 'team', 'dataset', or 'sample' absent — falling back to "
            f"'{fallback_col}' as batch key"
        )
        asap.obs["batch"] = "asap_" + asap.obs[fallback_col].astype(str)

    # 'pmdbs_sc_rnaseq.' prefix distinguishes from Kamath/Menon labels and
    # our own scANVI predictions. This unfiltered file has neither column yet.
    for src_col, dest_col in (
        ("C_scANVI", "pmdbs_sc_rnaseq.C_scANVI"),
        ("cell_type", "pmdbs_sc_rnaseq.cell_type"),
    ):
        if src_col not in asap.obs.columns:
            print(f"  Warning: '{src_col}' not found in ASAP obs — "
                  f"'{dest_col}' will be NaN")
            asap.obs[dest_col] = np.nan
        else:
            asap.obs[dest_col] = asap.obs[src_col].astype(str)
            n_unique = asap.obs[dest_col].nunique()
            print(f"  Preserved {dest_col}: {n_unique} unique labels")

    if "sample" not in asap.obs.columns:
        raise KeyError(
            "'sample' column not found in ASAP obs — required to join "
            "SAMPLE.csv for condition_id."
        )
    asap.obs["condition_id"] = _load_asap_condition_id(asap_sample_metadata_csv, asap.obs["sample"])
    print(f"  condition_id: {asap.obs['condition_id'].value_counts().to_dict()}")

    asap.obs["cell_type_scanvi"] = "Unknown"
    asap.obs["da_subtype"] = np.nan

    # ASAP barcodes already embed sample/dataset and are globally unique —
    # no prefix needed, just record them.
    asap.obs["original_barcode"] = asap.obs_names.astype(str)

    # Read-only join of sc-rnaseq-wf's own QC'd/annotated obs (final.h5ad),
    # by barcode. Only ~47% of cells here match (final.h5ad is a strict QC
    # subset); unmatched cells get NaN. 'wdl_final_' prefix avoids colliding
    # with this script's own obs columns.
    if os.path.exists(_ASAP_FINAL_H5AD_PATH):
        final_obs = sc.read_h5ad(_ASAP_FINAL_H5AD_PATH, backed="r").obs
        in_wdl_final = asap.obs_names.isin(final_obs.index)
        final_obs_aligned = final_obs.reindex(asap.obs_names)
        for col in final_obs_aligned.columns:
            asap.obs[f"wdl_final_{col}"] = final_obs_aligned[col].values
        asap.obs["in_wdl_final"] = in_wdl_final
        print(f"  Joined {final_obs_aligned.shape[1]} columns from asap-cohort.final.h5ad "
              f"(wdl_final_ prefix); {in_wdl_final.sum():,}/{asap.n_obs:,} cells matched")
    else:
        print(f"  Warning: {_ASAP_FINAL_H5AD_PATH} not found — skipping wdl_final_ join")
        asap.obs["in_wdl_final"] = False

    qc_survivors = _apply_qc_filters(
        asap, "ASAP cohort",
        pct_counts_mt_max=qc_pct_counts_mt_max, doublet_score_max=qc_doublet_score_max,
        total_counts_limits=qc_total_counts_limits, n_genes_limits=qc_n_genes_limits,
        doublet_batch_key="batch",
    )
    print(f"  ASAP QC survivors (untrimmed): {qc_survivors.n_obs:,} cells "
          f"(target was ~{n_cells:,} after {oversample_factor}x oversampling)")

    # postQC = strict QC survivors, trimmed to n_cells - n_reintroduce.
    # preQC = postQC + reintroduced doublet-only cells (~n_cells total).
    if n_cells > 0:
        n_reintroduce = round(n_cells * doublet_reintroduce_frac)
        main_quota = n_cells - n_reintroduce

        if qc_survivors.n_obs > main_quota:
            asap_postqc = subsample_adata(qc_survivors, main_quota, "ASAP postQC (QC survivors, trimmed to quota)")
        else:
            asap_postqc = qc_survivors
        counts = asap.layers["counts"]
        total_counts = np.asarray(counts.sum(axis=1)).ravel()
        n_genes_arr = np.asarray((counts > 0).sum(axis=1)).ravel()
        mt_present = asap.var_names.isin(_MT_GENE_ENSEMBL_IDS)
        if mt_present.any():
            mt_counts = np.asarray(counts[:, mt_present].sum(axis=1)).ravel()
            pct_counts_mt = np.where(total_counts > 0, mt_counts / total_counts * 100, 0.0)
        else:
            pct_counts_mt = np.zeros(asap.n_obs)

        doublet_score = asap.obs["doublet_score"].values
        pass_mt = pct_counts_mt <= qc_pct_counts_mt_max
        pass_total_counts = (
            (total_counts >= qc_total_counts_limits[0]) & (total_counts <= qc_total_counts_limits[1])
        )
        pass_n_genes = (n_genes_arr >= qc_n_genes_limits[0]) & (n_genes_arr <= qc_n_genes_limits[1])
        has_valid_doublet_score = ~pd.isna(doublet_score)
        fails_only_doublet = (
            pass_mt & pass_total_counts & pass_n_genes
            & has_valid_doublet_score & (doublet_score > qc_doublet_score_max)
        )
        candidate_idx = np.where(fails_only_doublet)[0]

        rng2 = np.random.default_rng(42)
        if len(candidate_idx) > n_reintroduce:
            reintroduce_idx = rng2.choice(candidate_idx, size=n_reintroduce, replace=False)
        else:
            reintroduce_idx = candidate_idx
            if len(candidate_idx) < n_reintroduce:
                print(f"  Warning: only {len(candidate_idx):,} clean-except-doublet cells "
                      f"available, fewer than the {n_reintroduce:,} requested "
                      f"({doublet_reintroduce_frac:.1%} of {n_cells:,})")

        reintroduced = asap[reintroduce_idx].copy()
        print(f"  Reintroducing {reintroduced.n_obs:,} clean-except-doublet cells "
              f"({doublet_reintroduce_frac:.1%} of {n_cells:,} target, from "
              f"{len(candidate_idx):,} eligible)")

        asap_preqc = ad.concat([asap_postqc, reintroduced], join="inner")
        _log_qc_pair("ASAP", asap_preqc, asap_postqc,
                     f"postQC + {reintroduced.n_obs:,} doublet-only cells reintroduced",
                     "QC survivors")
    else:
        # No downsampling: no oversample/reintroduce quotas, fall back to
        # the plain preQC (full pool) >= postQC (survivors) relationship.
        asap_preqc = asap
        asap_postqc = qc_survivors
        _log_qc_pair("ASAP", asap_preqc, asap_postqc, postqc_note="QC survivors")

    print(f"  Writing ASAP preQC cache → {preqc_cache}")
    _write_h5ad_cache(asap_preqc, preqc_cache)

    print(f"  Writing ASAP postQC cache → {postqc_cache}")
    _write_h5ad_cache(asap_postqc, postqc_cache)

    _write_qc_summary(_qc_summary_path(asap_h5ad, tag))

    return asap_preqc, asap_postqc


# ---------------------------------------------------------------------------
# Preprocessing
# ---------------------------------------------------------------------------


def main(args: argparse.Namespace) -> None:
    """Build and cache each of the four sources' (preQC, postQC) intermediate AnnData."""
    print("\n=== Locating cached input files ===")
    kamath_dir = os.path.join(_REFERENCE_DATA_DIR, "kamath")
    menon_dir = os.path.join(_REFERENCE_DATA_DIR, "menon")
    siletti_dir = os.path.join(_REFERENCE_DATA_DIR, "siletti")

    menon_path = _find_menon_h5ad(menon_dir)

    kamath_input_paths = sorted(
        p for p in glob.glob(os.path.join(kamath_dir, "*.h5ad"))
        if not os.path.basename(p).startswith("kamath_merged")
    )
    for p in kamath_input_paths:
        _log_input_file(p)

    for file_stem in _SILETTI_MIDBRAIN_REGIONS.values():
        _log_input_file(os.path.join(siletti_dir, f"{file_stem}.h5ad"))
    _log_input_file(os.path.join(siletti_dir, f"{_SILETTI_ALL_NEURONS_FILE}.h5ad"))

    _log_input_file(_ASAP_H5AD_PATH)

    exclude = [f.strip() for f in args.exclude_kamath_files.split(",") if f.strip()]

    print("\n=== Loading Kamath datasets ===")
    if args.kamath_all_random and args.downsample_kamath is None:
        parser.error("--downsample-kamath is required when --kamath-all-random is set.")
    if not args.kamath_all_random and args.downsample_kamath is not None:
        parser.error("--downsample-kamath is only valid when --kamath-all-random is set.")

    kamath_anchor_name = "kamath_merged_all_rand" if args.kamath_all_random else "kamath_merged_da_all_non_da"
    kamath_anchor = os.path.join(kamath_dir, "intermediate", kamath_anchor_name)
    kamath_tag = str(args.downsample_kamath) if args.kamath_all_random else str(args.kamath_non_da_per_type)
    kamath_preqc_cache, kamath_postqc_cache = _qc_cache_paths(kamath_anchor, kamath_tag)
    if args.use_cached_intermediate_h5ad and os.path.exists(kamath_preqc_cache) and os.path.exists(kamath_postqc_cache):
        _kp = sc.read_h5ad(kamath_preqc_cache, backed="r")
        _kq = sc.read_h5ad(kamath_postqc_cache, backed="r")
        _log_cached_qc("Kamath", kamath_tag, _kp, _kq)
        _kp.file.close(); _kq.file.close()
    else:
        kamath_preqc, kamath_postqc = load_kamath(
            kamath_dir, exclude, args.da_subtype_column,
            n_non_da_per_type=args.kamath_non_da_per_type,
            all_random=args.kamath_all_random,
            n_total=args.downsample_kamath or 0,
            qc_pct_counts_mt_max=args.qc_pct_counts_mt_max,
            qc_doublet_score_max=args.qc_doublet_score_max,
            qc_total_counts_limits=tuple(args.qc_total_counts_limits),
            qc_n_genes_limits=tuple(args.qc_n_genes_limits),
        )
        print(f"  Writing Kamath preQC/postQC caches → {kamath_preqc_cache} / {kamath_postqc_cache}")
        _write_h5ad_cache(kamath_preqc, kamath_preqc_cache)
        _write_h5ad_cache(kamath_postqc, kamath_postqc_cache)

    print("\n=== Loading Menon dataset (SN only) ===")
    load_menon_sn(
        menon_path, n_cells=args.downsample_menon,
        use_cached=args.use_cached_intermediate_h5ad,
        qc_pct_counts_mt_max=args.qc_pct_counts_mt_max,
        qc_doublet_score_max=args.qc_doublet_score_max,
        qc_total_counts_limits=tuple(args.qc_total_counts_limits),
        qc_n_genes_limits=tuple(args.qc_n_genes_limits),
    )

    print("\n=== Loading Siletti dataset (DA-aware midbrain sampling) ===")
    siletti_neuron_target, siletti_cell_target = _SILETTI_SCALE_PRESETS[args.siletti_scale]
    print(f"  Scale {args.siletti_scale:,}: "
          f"neuron_target={siletti_neuron_target}, cell_target={siletti_cell_target}")
    siletti_anchor = os.path.join(
        siletti_dir, "intermediate",
        f"midbrain_da_neurons{siletti_neuron_target}_cells{siletti_cell_target}",
    )
    siletti_preqc_cache, siletti_postqc_cache = _qc_cache_paths(siletti_anchor, "")
    if args.use_cached_intermediate_h5ad and os.path.exists(siletti_preqc_cache) and os.path.exists(siletti_postqc_cache):
        _sp = sc.read_h5ad(siletti_preqc_cache, backed="r")
        _sq = sc.read_h5ad(siletti_postqc_cache, backed="r")
        _log_cached_qc("Siletti", f"scale_{args.siletti_scale}", _sp, _sq)
        _sp.file.close(); _sq.file.close()
    else:
        allen_annotations = _load_allen_annotations()
        neuron_barcodes = _load_all_neuron_barcodes(siletti_dir)
        siletti_preqc, siletti_postqc = load_siletti(
            siletti_dir, allen_annotations, neuron_barcodes,
            neuron_target=siletti_neuron_target,
            cell_target=siletti_cell_target,
            qc_pct_counts_mt_max=args.qc_pct_counts_mt_max,
            qc_doublet_score_max=args.qc_doublet_score_max,
            qc_total_counts_limits=tuple(args.qc_total_counts_limits),
            qc_n_genes_limits=tuple(args.qc_n_genes_limits),
        )
        print(f"  Writing Siletti preQC/postQC caches → {siletti_preqc_cache} / {siletti_postqc_cache}")
        _write_h5ad_cache(siletti_preqc, siletti_preqc_cache)
        _write_h5ad_cache(siletti_postqc, siletti_postqc_cache)

    print("\n=== Loading ASAP PMDBS scRNA-seq cohort ===")
    load_asap_cohort(
        asap_h5ad=_ASAP_H5AD_PATH,
        asap_sample_metadata_csv=_ASAP_SAMPLE_METADATA_PATH,
        n_cells=args.downsample_asap_cohort,
        use_cached=args.use_cached_intermediate_h5ad,
        qc_pct_counts_mt_max=args.qc_pct_counts_mt_max,
        qc_doublet_score_max=args.qc_doublet_score_max,
        qc_total_counts_limits=tuple(args.qc_total_counts_limits),
        qc_n_genes_limits=tuple(args.qc_n_genes_limits),
    )

    print(
        "\nDone. Intermediate preQC/postQC caches (and each dataset's own "
        "qc_summary.tsv, next to its caches) are ready under each source's "
        "intermediate/ directory. Run "
        "integrate_kamath_menon_siletti_asap.py next, with the SAME "
        "--kamath-non-da-per-type / --kamath-all-random / --downsample-menon / "
        "--siletti-scale / --downsample-asap-cohort values, to locate them and "
        "train the scVI/scANVI reference."
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--use-cached-intermediate-h5ad",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Skip a source entirely if its (preQC, postQC) cache pair\n"
             "already exists. Default False (always regenerate).\n"
             "[Default: False]",
    )
    parser.add_argument(
        "--da-subtype-column",
        default="author_cell_type",
        help="obs column in the Kamath DA Neurons H5AD carrying the 10 DA subtype\n"
             "labels (e.g. SOX6_AGTR1, CALB1_GEM).\n"
             "[Default: author_cell_type]",
    )
    parser.add_argument(
        "--exclude-kamath-files",
        default=(
            "human_nurr_positive_nuclei_10x_scrna_seq.h5ad,"
            "human_nurr_negative_nuclei_10x_scrna_seq.h5ad"
        ),
        help="Comma-separated Kamath H5AD filenames to skip. Both Nurr files\n"
             "excluded by default (see kamath_overlap.py).\n"
             "[Default: human_nurr_positive_nuclei_10x_scrna_seq.h5ad,\n"
             "          human_nurr_negative_nuclei_10x_scrna_seq.h5ad]",
    )
    parser.add_argument(
        "--kamath-non-da-per-type",
        type=int,
        default=_DEFAULT_KAMATH_NON_DA_PER_TYPE,
        help="Cells to retain from each non-DA Kamath cell-type file, downsampled\n"
             "BEFORE that file's own QC pass. DA neurons file always kept in full.\n"
             "Ignored when --kamath-all-random is set.\n"
             f"[Default: {_DEFAULT_KAMATH_NON_DA_PER_TYPE}]",
    )
    parser.add_argument(
        "--kamath-all-random",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Pool all non-excluded Kamath files and draw --downsample-kamath\n"
             "cells at random (no per-file cap, no DA-first priority).\n"
             "Requires --downsample-kamath. Incompatible without it.\n"
             "[Default: False]",
    )
    parser.add_argument(
        "--downsample-kamath",
        type=int,
        default=None,
        help="Total cells to draw when --kamath-all-random is True. The global\n"
             "random draw happens BEFORE QC. Required when --kamath-all-random\n"
             "is set; incompatible when it is not.\n"
             "[Default: None]",
    )
    parser.add_argument(
        "--downsample-menon",
        type=int,
        default=_DEFAULT_DOWNSAMPLE_MENON,
        help="Randomly subsample the Menon SN dataset to N cells, BEFORE QC.\n"
             "0 = no downsampling.\n"
             f"[Default: {_DEFAULT_DOWNSAMPLE_MENON}]",
    )
    parser.add_argument(
        "--siletti-scale",
        type=int,
        choices=sorted(_SILETTI_SCALE_PRESETS),
        default=_SILETTI_DEFAULT_SCALE,
        help="Siletti downsampling scale (approximate grand total cell count),\n"
             "applied before QC (see load_siletti and _SILETTI_SCALE_PRESETS).\n"
             "  50000  -> neuron_target=15000, cell_target=23000  (~50.7k total)\n"
             "  100000 -> neuron_target=15000, cell_target=47643  (100k total)\n"
             f"[Default: {_SILETTI_DEFAULT_SCALE}]",
    )
    parser.add_argument(
        "--qc-pct-counts-mt-max",
        type=float,
        default=_QC_PCT_COUNTS_MT_MAX,
        help="Maximum %% mitochondrial counts per cell, applied to all four sources.\n"
             f"[Default: {_QC_PCT_COUNTS_MT_MAX}]",
    )
    parser.add_argument(
        "--qc-doublet-score-max",
        type=float,
        default=_QC_DOUBLET_SCORE_MAX,
        help="Maximum scrublet doublet score per cell, applied to all four sources.\n"
             f"[Default: {_QC_DOUBLET_SCORE_MAX}]",
    )
    parser.add_argument(
        "--qc-total-counts-limits",
        type=int,
        nargs=2,
        metavar=("MIN", "MAX"),
        default=list(_QC_TOTAL_COUNTS_LIMITS),
        help="Minimum and maximum total UMI counts per cell, applied to all four sources.\n"
             f"[Default: {list(_QC_TOTAL_COUNTS_LIMITS)}]",
    )
    parser.add_argument(
        "--qc-n-genes-limits",
        type=int,
        nargs=2,
        metavar=("MIN", "MAX"),
        default=list(_QC_N_GENES_LIMITS),
        help="Minimum and maximum genes detected per cell, applied to all four sources.\n"
             f"[Default: {list(_QC_N_GENES_LIMITS)}]",
    )
    parser.add_argument(
        "--downsample-asap-cohort",
        type=int,
        default=_DEFAULT_DOWNSAMPLE_ASAP_COHORT,
        help="Randomly subsample the ASAP cohort to N cells, BEFORE QC. 0 = no\n"
             "downsampling (not recommended, QC would run on all 7.55M cells).\n"
             f"[Default: {_DEFAULT_DOWNSAMPLE_ASAP_COHORT}]",
    )
    args = parser.parse_args()

    repo_root = next(p for p in Path(__file__).resolve().parents if (p / ".git").exists())
    if Path.cwd() != repo_root:
        print(f"[INFO] Changing working directory to repo root: {repo_root}")
        os.chdir(repo_root)

    log_run_command(args, log_dir=_REFERENCE_DATA_DIR)
    main(args)
