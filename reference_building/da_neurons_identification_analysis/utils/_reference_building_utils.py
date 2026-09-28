"""
Shared path/config constants and helpers for downsample_kamath_menon_siletti_asap.py,
integrate_kamath_menon_siletti_asap.py, and build_deconv_reference_h5ad.py.
"""

import os
import sys
import glob
import argparse
from pathlib import Path
from datetime import datetime

import numpy as np
from scipy.sparse import csr_matrix, issparse

# Expected local cache directory name
_REFERENCE_DATA_DIR = "reference_building/da_neurons_identification_analysis/reference_data"

# Every DA VGLUT2 Siletti cell is kept, then topped up to neuron_target, then
# cell_target -- RN's small size (~4,714 cells) is why cell_target roughly
# doubles rather than exactly doubles the grand total. Keyed by grand total.
_SILETTI_SCALE_PRESETS = {
    50_000: (15_000, 23_000),
    100_000: (15_000, 47_643),
}
_SILETTI_DEFAULT_SCALE = 50_000

# Current pmdbs_sc_rnaseq cohort_analysis output. Hardcoded, not a CLI arg.
_ASAP_COHORT_DIR = os.path.join(
    _REFERENCE_DATA_DIR, "asap_cohort", "pmdbs_sc_rnaseq", "cohort_analysis",
    "3.1.0", "2025-12-01T22-36-30Z",
)
# Unfiltered, full-gene-panel ASAP file (36,601 genes incl. all 17 Kamath DA
# markers) -- not asap-cohort.final.h5ad (sc-rnaseq-wf's own HVG subset, only
# 3,000 genes, 0/17 markers). Barcode selection is via backed-mode obs access.
_ASAP_H5AD_PATH = os.path.join(_ASAP_COHORT_DIR, "asap-cohort.merged_cleaned_unfiltered.h5ad")
_ASAP_SAMPLE_METADATA_PATH = os.path.join(_ASAP_COHORT_DIR, "metadata", "release", "v5.0.0", "SAMPLE.csv")
# sc-rnaseq-wf's own QC'd/annotated output, joined read-only onto our ASAP
# preQC/postQC obs under a 'wdl_final_' prefix (see load_asap_cohort).
_ASAP_FINAL_H5AD_PATH = os.path.join(_ASAP_COHORT_DIR, "asap-cohort.final.h5ad")

# Kamath et al. (Nat Neurosci 2022) DA marker gene symbols. Resolved to
# Ensembl IDs at runtime via a source's var['feature_name'] (see
# _resolve_marker_ensembl_ids in the integrate script).
_KAMATH_DA_MARKER_SYMBOLS = (
    # Pan-DA identity
    "TH", "SLC6A3", "SLC18A2", "DDC",
    # Branch markers (the two main DA subdivisions)
    "SOX6", "CALB1",
    # Subtype-defining markers (e.g. SOX6_AGTR1, CALB1_GEM)
    "AGTR1", "DDT", "GFRA2", "PART1", "CALCR", "CRYM", "CCDC68", "GEM", "PPP1R17", "RBP4", "TRHR",
)

# CLI defaults shared between both scripts' argparse blocks, so they can't
# drift out of sync -- downsample uses them to downsample, integrate only to
# reconstruct matching cache filenames.
_DEFAULT_KAMATH_NON_DA_PER_TYPE = 5000
_DEFAULT_DOWNSAMPLE_MENON = 0
_DEFAULT_DOWNSAMPLE_ASAP_COHORT = 50000


# ---------------------------------------------------------------------------
# Cache-path construction
# ---------------------------------------------------------------------------

def _intermediate_dir(h5ad_path: str) -> Path:
    """'intermediate/' dir for derived caches of a dataset. Chain-safe (won't nest)."""
    parent = Path(h5ad_path).parent
    return parent if parent.name == "intermediate" else parent / "intermediate"

def _sn_cache_path(menon_h5ad: str) -> str:
    """Conventional path for the SN-filtered Menon cache file."""
    p = Path(menon_h5ad)
    return str(_intermediate_dir(menon_h5ad) / (p.stem + "_sn_only.h5ad"))

def _qc_cache_paths(anchor_h5ad: str, tag: str) -> tuple[str, str]:
    """(preQC, postQC) cache paths for `anchor_h5ad`'s dataset, tagged (e.g. 'rand_100000'; '' if unique)."""
    d = _intermediate_dir(anchor_h5ad)
    stem = Path(anchor_h5ad).stem
    suffix = f"_{tag}" if tag else ""
    return str(d / f"{stem}{suffix}_preQC.h5ad"), str(d / f"{stem}{suffix}_postQC.h5ad")

def _qc_summary_path(anchor_h5ad: str, tag: str) -> str:
    """qc_summary.tsv path, next to that dataset's own preQC/postQC caches."""
    d = _intermediate_dir(anchor_h5ad)
    stem = Path(anchor_h5ad).stem
    suffix = f"_{tag}" if tag else ""
    return str(d / f"{stem}{suffix}_qc_summary.tsv")

def _find_menon_h5ad(menon_dir: str) -> str:
    """Locate the already-downloaded Menon H5AD in `menon_dir`."""
    existing = sorted(glob.glob(os.path.join(menon_dir, "*.h5ad")))
    if not existing:
        raise FileNotFoundError(
            f"No H5AD files found in {menon_dir}. "
            "Run download_kamath_menon_siletti.py first."
        )
    if len(existing) > 1:
        print(f"  Warning: {len(existing)} H5ADs in {menon_dir}; using first: {existing[0]}")
    else:
        size_mb = os.path.getsize(existing[0]) / 1e6
        print(f"  Found ({size_mb:.0f} MB): {existing[0]}")
    return existing[0]


# ---------------------------------------------------------------------------
# Raw count extraction and validation
# ---------------------------------------------------------------------------

def extract_raw_counts(adata, source_name: str = "") -> csr_matrix:
    """
    Extract raw integer count matrix, checking layers['counts'], adata.raw.X,
    then adata.X in order. Validates the full matrix is integer-valued
    (whole-number floats like 1.0 pass; 1.5 raises) and casts to int32.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData to extract counts from.
    source_name : str, optional
        Label for error/info messages (e.g. file basename).

    Returns
    -------
    csr_matrix
        Sparse int32 count matrix (cells × genes).

    Raises
    ------
    ValueError
        No raw count layer found, or matrix contains non-integer values.
    """
    if "counts" in adata.layers:
        mat = adata.layers["counts"]
    elif adata.raw is not None:
        mat = adata.raw.X
    else:
        _pfx = f"{source_name}: " if source_name else ""
        print(
            f"  [INFO] {_pfx}no layers['counts'] or adata.raw — "
            "using adata.X (integer-valued check happens below)"
        )
        mat = adata.X

    mat = csr_matrix(mat) if not issparse(mat) else mat.tocsr()

    if not np.issubdtype(mat.dtype, np.integer) and mat.data.size > 0:
        non_int_mask = mat.data % 1 != 0
        if non_int_mask.any():
            examples = mat.data[non_int_mask][:5].tolist()
            _pfx = f"{source_name}: " if source_name else ""
            raise ValueError(
                f"{_pfx}count matrix contains non-integer values "
                f"(e.g. {examples}). Raw counts must be integer-valued. "
                "Use a raw H5AD download or check --counts-layer."
            )

    if mat.dtype != np.int32:
        mat = mat.astype(np.int32)

    return mat


# ---------------------------------------------------------------------------
# Run logging
# ---------------------------------------------------------------------------

class _Tee:
    """Mirror writes to both the original stream and a log file (delegates unknown attrs)."""

    def __init__(self, original, log_file) -> None:
        self._original = original
        self._log = log_file

    def write(self, data: str) -> int:
        self._original.write(data)
        self._log.write(data)
        return len(data)

    def flush(self) -> None:
        self._original.flush()
        self._log.flush()

    def __getattr__(self, name):
        return getattr(self._original, name)

def log_run_command(args: argparse.Namespace, log_dir: str = None) -> None:
    """
    Log the invocation and tee stdout/stderr to `{log_dir or args.output_dir}/{script_stem}.log`
    for the process lifetime. Uses sys.argv[0] (not __file__) for script_stem --
    __file__ here would resolve to this shared module, not the calling script.

    Parameters
    ----------
    args : argparse.Namespace
        Must have an `output_dir` attribute unless `log_dir` is given.
    log_dir : str, optional
        Directory for the log file. Defaults to `args.output_dir` -- pass
        this explicitly for scripts whose log shouldn't live under
        --output-dir (e.g. downsample_kamath_menon_siletti_asap.py's log
        lives under reference_building/da_neurons_identification_analysis/reference_data/, next to the caches
        it actually produces).

    Returns
    -------
    None
    """
    log_file = Path(log_dir if log_dir is not None else args.output_dir) / f"{Path(sys.argv[0]).stem}.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    f = open(log_file, "w", encoding="utf-8")  # noqa: SIM115 — intentionally kept open
    f.write(f"[{timestamp}] {' '.join(sys.argv)}\n\n")
    sys.stdout = _Tee(sys.stdout, f)
    sys.stderr = _Tee(sys.stderr, f)
    args._log_file_path = str(log_file)
