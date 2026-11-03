#!/usr/bin/env python3
"""
Format Ma et al. (2025, Nat Commun) Visium Substantia Nigra spot data into
per-sample h5ad files for this repo's Cell2Location pipeline
(docker/cell2location_tools/scripts/cell2location_deconv.py --spatial-h5ad).

Each `*_result-*.zip` is one sample's full Space Ranger `outs/` dir.

Uses `squidpy.read.visium()` (not scanpy's) since its output already
matches the structure `cell2location_deconv.py`'s `_load_spatial()`
expects. Every gene is kept unfiltered -- that script intersects against
whatever reference it's given at Stage 2, so no filtering happens here.
"""

import os
import re
import sys
import warnings
import argparse
import tempfile
import zipfile
from pathlib import Path
from datetime import datetime

import openpyxl
import numpy as np
import pandas as pd
import squidpy as sq
import anndata as ad
import scipy.sparse as sp

_DATA_DIR = "reference_building/benchmark/reference_data/ma/data"
_ORIGINAL_DIR = f"{_DATA_DIR}/original"
_SUPP_XLSX = f"{_ORIGINAL_DIR}/41467_2025_62478_MOESM3_ESM.xlsx"
_OUTPUT_DIR = f"{_DATA_DIR}/reformatted/spatial_h5ad"

_SN_SHEET_A = "A"
_ZIP_GLOB = "*_result-*.zip"
_ZIP_TB_RE = re.compile(r"_T(\d{4})_result")


# ---------------------------------------------------------------------------
# Logging (mirrors 03_build_benchmark.py / 04_plot_sn_dopa_gaba_meg3_pd_vs_control.py)
# ---------------------------------------------------------------------------

class _Tee:
    """Mirror writes to both the original stream and a log file (used for stdout/stderr)."""

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


def log_run_command(args: argparse.Namespace) -> None:
    """Log the invocation to {_DATA_DIR}/{script_stem}.log and tee all subsequent output there."""
    log_file = Path(_DATA_DIR) / f"{Path(__file__).stem}.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    f = open(log_file, "w", encoding="utf-8")  # noqa: SIM115 — intentionally kept open
    f.write(f"[{timestamp}] {' '.join(sys.argv)}\n\n")
    sys.stdout = _Tee(sys.stdout, f)
    sys.stderr = _Tee(sys.stderr, f)


# ---------------------------------------------------------------------------
# Tab A metadata
# ---------------------------------------------------------------------------

def load_tb_metadata(xlsx_path: str) -> pd.DataFrame:
    """
    Load per-sample Condition/Neuropath_Diagnosis from tab A, keyed by TB #.

    Parameters
    ----------
    xlsx_path : str
        Path to 41467_2025_62478_MOESM3_ESM.xlsx.

    Returns
    -------
    pandas.DataFrame
        One row per TB #, columns `TB #`, `Condition`, `Neuropath_Diagnosis`.
    """
    wb = openpyxl.load_workbook(xlsx_path, read_only=True, data_only=True)
    rows = list(wb[_SN_SHEET_A].iter_rows(values_only=True))
    header = rows[0]
    meta = pd.DataFrame(rows[1:], columns=header)
    meta = meta[["TB #", "Condition", "Neuropath_Diagnosis"]].dropna(subset=["TB #"])
    return meta.drop_duplicates(subset=["TB #"])


# ---------------------------------------------------------------------------
# Visium zip -> AnnData
# ---------------------------------------------------------------------------

def find_visium_zips(original_dir: str) -> dict:
    """
    Find each Space Ranger result zip and parse its TB # from the filename.

    Parameters
    ----------
    original_dir : str
        Directory holding the `*_result-*.zip` files.

    Returns
    -------
    dict
        Maps TB # (e.g. "T-4424") to the zip's `Path`.
    """
    zips = {}
    for zip_path in sorted(Path(original_dir).glob(_ZIP_GLOB)):
        m = _ZIP_TB_RE.search(zip_path.name)
        if not m:
            print(f"  WARNING: could not parse a TB # out of {zip_path.name}, skipping")
            continue
        zips[f"T-{m.group(1)}"] = zip_path
    return zips


def extract_outs_dir(zip_path: Path, extract_root: str) -> str:
    """
    Extract a Space Ranger zip's `outs/` subtree into `extract_root`.

    Parameters
    ----------
    zip_path : pathlib.Path
        Path to the `*_result-*.zip` file.
    extract_root : str
        Directory to extract into (a fresh temporary directory per sample).

    Returns
    -------
    str
        Path to the extracted `.../outs` directory.
    """
    with zipfile.ZipFile(zip_path) as zf:
        members = [m for m in zf.namelist() if "/outs/" in m or m.endswith("/outs")]
        zf.extractall(extract_root, members=members)

    outs_dirs = list(Path(extract_root).glob("*/outs"))
    if len(outs_dirs) != 1:
        raise ValueError(f"Expected exactly one outs/ dir under {extract_root}, found {outs_dirs}")
    return str(outs_dirs[0])


def build_sample_adata(outs_dir: str, tb_id: str, meta_row: pd.Series | None) -> ad.AnnData:
    """
    Load a Space Ranger `outs/` directory and annotate it for this sample.

    Parameters
    ----------
    outs_dir : str
        Path to the extracted `outs/` directory (filtered_feature_bc_matrix.h5
        + spatial/).
    tb_id : str
        This sample's TB # (e.g. "T-4424").
    meta_row : pandas.Series or None
        The sample's row from `load_tb_metadata`, or None if no tab A match.

    Returns
    -------
    anndata.AnnData
        Spots x genes, raw counts in `.X`, Ensembl IDs in `.var['gene_ids']`,
        `.obs['sample_id'/'Condition'/'Neuropath_Diagnosis']` added.
    """
    with warnings.catch_warnings():
        # squidpy warns about duplicate var_names from its scanpy-based h5 loader --
        # var_names_make_unique() below already fixes it, so this is just noise.
        warnings.filterwarnings("ignore", message="Variable names are not unique")
        adata = sq.read.visium(outs_dir)
    adata.var_names_make_unique()

    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    sample_vals = adata.X[: min(500, adata.n_obs)].toarray()
    if not np.all(sample_vals % 1 == 0):
        raise ValueError(f"{tb_id}: .X contains non-integer values -- not raw counts?")

    adata.obs["sample_id"] = tb_id
    adata.obs["Condition"] = meta_row["Condition"] if meta_row is not None else pd.NA
    adata.obs["Neuropath_Diagnosis"] = meta_row["Neuropath_Diagnosis"] if meta_row is not None else pd.NA
    return adata


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(args: argparse.Namespace) -> None:
    """Convert every sample's Space Ranger zip into a Cell2Location-ready spatial h5ad."""
    all_zips = find_visium_zips(args.original_dir)
    print(f"Found {len(all_zips)} Visium sample zip(s) in {args.original_dir}: {sorted(all_zips)}")

    tb_meta = load_tb_metadata(args.xlsx_path)
    tb_meta_by_id = {row["TB #"]: row for _, row in tb_meta.iterrows()}

    os.makedirs(args.output_dir, exist_ok=True)

    for tb_id, zip_path in sorted(all_zips.items()):
        meta_row = tb_meta_by_id.get(tb_id)
        if meta_row is None:
            print(f"  WARNING: {tb_id} has no tab A metadata match (Condition/Neuropath_Diagnosis will be NA)")

        with tempfile.TemporaryDirectory(prefix=f"visium_{tb_id}_") as tmp_dir:
            outs_dir = extract_outs_dir(zip_path, tmp_dir)
            adata = build_sample_adata(outs_dir, tb_id, meta_row)

        condition = meta_row["Condition"] if meta_row is not None else "NA"
        print(f"  {tb_id} ({condition}): {adata.n_obs:,} spots x {adata.n_vars:,} genes")

        out_path = os.path.join(args.output_dir, f"{tb_id}.spatial.h5ad")
        adata.write_h5ad(out_path, compression="lzf")
        print(f"    Written {out_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--original-dir",
        default=None,
        help="Directory holding the Space Ranger *_result-*.zip files. A relative path is\n"
             "resolved against the directory this command is run from (not the repo root).\n"
             f"[Default: {_ORIGINAL_DIR}, relative to the repo root]",
    )
    parser.add_argument(
        "--xlsx-path",
        default=None,
        help="Ma et al. 2025 supplementary workbook (for Condition/Neuropath_Diagnosis).\n"
             "A relative path is resolved against the directory this command is run from.\n"
             f"[Default: {_SUPP_XLSX}, relative to the repo root]",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory for the per-sample spatial h5ad files. A relative path is\n"
             "resolved against the directory this command is run from.\n"
             f"[Default: {_OUTPUT_DIR}, relative to the repo root]",
    )
    args = parser.parse_args()

    invocation_cwd = Path.cwd()
    for attr in ("original_dir", "xlsx_path", "output_dir"):
        value = getattr(args, attr)
        if value is not None:
            setattr(args, attr, str((invocation_cwd / value).resolve()))

    repo_root = next(p for p in Path(__file__).resolve().parents if (p / ".git").exists())
    if Path.cwd() != repo_root:
        print(f"[INFO] Changing working directory to repo root: {repo_root}")
        os.chdir(repo_root)

    if args.original_dir is None:
        args.original_dir = _ORIGINAL_DIR
    if args.xlsx_path is None:
        args.xlsx_path = _SUPP_XLSX
    if args.output_dir is None:
        args.output_dir = _OUTPUT_DIR

    log_run_command(args)
    main(args)
