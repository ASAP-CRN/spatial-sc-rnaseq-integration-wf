#!/usr/bin/env python3
"""
Precompute cell2location's per-spot cell-type fraction correlation (Pearson,
"co-occurrence") tables for Strategy 1, one table per sample plus one pooled
across all samples -- so the notebook can just load a CSV instead of
recomputing `.corr()` live per sample/strategy selection.

This is the Strategy-1 counterpart to 08_cell_type_co_occurrence_reference_RefStrategy2.py
and RefStrategy3.py (point --cell2location-dir at the matching cell2location_StrategyN dir).
"""

import os
import sys
import argparse
from pathlib import Path
from datetime import datetime

import pandas as pd
import anndata as ad

_DATA_DIR = "reference_building/benchmark/reference_data/ma/data"
_CELL2LOCATION_DIR = f"{_DATA_DIR}/reformatted/cell2location_Strategy1"
_OUTPUT_DIR = _CELL2LOCATION_DIR
_H5AD_SUFFIX = ".deconvoluted_clustered.h5ad"
_COOCCURRENCE_OUT = "cell_type_cooccurrence.csv.gz"

# Sample value standing in for the pooled-across-all-samples row; the
# notebook's "Sample" dropdown must use this exact string for its "All
# samples" option so it matches rows written here.
ALL_SAMPLES_LABEL = "All samples (pooled)"


# ---------------------------------------------------------------------------
# Logging
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
# Loading
# ---------------------------------------------------------------------------

def _strip_ct_prefix(name: str) -> str:
    """Strip cell2location's `q05_cell_abundance_w_sf` obsm/obs column prefix."""
    for prefix in ("frac_q05cell_abundance_w_sf_", "q05cell_abundance_w_sf_"):
        if name.startswith(prefix):
            return name[len(prefix):]
    return name


def load_sample_fractions(h5ad_path: Path) -> pd.DataFrame:
    """
    Load a deconvoluted sample's per-spot fine-grained cell-type fractions.

    Parameters
    ----------
    h5ad_path : pathlib.Path
        Path to a <sample_id>.deconvoluted_clustered.h5ad from cell2location_deconv.py.

    Returns
    -------
    pandas.DataFrame
        Spots x fine cell_type_pred labels; each row sums to 1.
    """
    adata = ad.read_h5ad(h5ad_path)
    cell_types = list(adata.obsm["q05_cell_abundance_w_sf"].columns)
    frac_df = adata.obsm["cell_type_fractions"][[f"frac_{ct}" for ct in cell_types]].copy()
    frac_df.columns = [_strip_ct_prefix(ct) for ct in cell_types]
    frac_df.index = adata.obs_names
    return frac_df


# ---------------------------------------------------------------------------
# Co-occurrence
# ---------------------------------------------------------------------------

def _corr_long(frac_df: pd.DataFrame, sample_label: str) -> pd.DataFrame:
    """Pearson correlation matrix for one sample's (or the pooled) fractions, in long/tidy form."""
    corr = frac_df.corr()
    long = corr.rename_axis("cell_type_a").reset_index().melt(
        id_vars="cell_type_a", var_name="cell_type_b", value_name="pcc"
    )
    long.insert(0, "sample", sample_label)
    return long


def build_cooccurrence_tables(cell2location_dir: str, sample_ids=None) -> pd.DataFrame:
    """
    Per-sample and pooled-across-samples cell-type fraction correlation tables.

    Parameters
    ----------
    cell2location_dir : str
        Directory of <sample_id>.deconvoluted_clustered.h5ad files.
    sample_ids : list of str or None
        Restrict to these sample IDs (e.g. "T-0081"); None processes all found.

    Returns
    -------
    pandas.DataFrame
        Long/tidy table: sample, cell_type_a, cell_type_b, pcc. `sample` is
        either a real sample ID or `ALL_SAMPLES_LABEL` for the pooled rows
        (correlation computed on all samples' spots concatenated together).
    """
    h5ad_paths = sorted(Path(cell2location_dir).glob(f"*{_H5AD_SUFFIX}"))
    if sample_ids:
        h5ad_paths = [p for p in h5ad_paths if p.name.replace(_H5AD_SUFFIX, "") in sample_ids]

    tables, frac_dfs = [], []
    for h5ad_path in h5ad_paths:
        sample_name = h5ad_path.name.replace(_H5AD_SUFFIX, "")
        print(f"  Loading {sample_name}")
        frac_df = load_sample_fractions(h5ad_path)
        tables.append(_corr_long(frac_df, sample_name))
        frac_dfs.append(frac_df)

    pooled_frac = pd.concat(frac_dfs)
    tables.append(_corr_long(pooled_frac, ALL_SAMPLES_LABEL))
    return pd.concat(tables, ignore_index=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(args: argparse.Namespace) -> None:
    """Build and write the per-sample + pooled cell-type co-occurrence tables."""
    print(f"Loading cell2location outputs from {args.cell2location_dir}")
    cooccurrence = build_cooccurrence_tables(args.cell2location_dir, args.samples)

    n_samples = cooccurrence.loc[cooccurrence["sample"] != ALL_SAMPLES_LABEL, "sample"].nunique()
    print(f"\nBuilt co-occurrence tables for {n_samples} sample(s), plus pooled ('{ALL_SAMPLES_LABEL}')")

    os.makedirs(args.output_dir, exist_ok=True)
    out_path = os.path.join(args.output_dir, _COOCCURRENCE_OUT)
    cooccurrence.to_csv(out_path, index=False, compression="gzip")
    print(f"\nWritten {out_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--cell2location-dir",
        default=None,
        help="Directory of <sample_id>.deconvoluted_clustered.h5ad files -- point this at\n"
             "cell2location_Strategy1. A relative path is\n"
             "resolved against the directory this command is run from.\n"
             f"[Default: {_CELL2LOCATION_DIR}, relative to the repo root]",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory for the co-occurrence table. A relative\n"
             "path is resolved against the directory this command is run from.\n"
             f"[Default: {_OUTPUT_DIR}, relative to the repo root]",
    )
    parser.add_argument(
        "--samples",
        nargs="+",
        default=None,
        help="Restrict to these sample IDs (e.g. T-0081 T-2465). Default: all samples found.",
    )
    args = parser.parse_args()

    invocation_cwd = Path.cwd()
    for attr in ("cell2location_dir", "output_dir"):
        value = getattr(args, attr)
        if value is not None:
            setattr(args, attr, str((invocation_cwd / value).resolve()))

    repo_root = next(p for p in Path(__file__).resolve().parents if (p / ".git").exists())
    if Path.cwd() != repo_root:
        print(f"[INFO] Changing working directory to repo root: {repo_root}")
        os.chdir(repo_root)

    if args.cell2location_dir is None:
        args.cell2location_dir = _CELL2LOCATION_DIR
    if args.output_dir is None:
        args.output_dir = _OUTPUT_DIR

    log_run_command(args)
    main(args)
