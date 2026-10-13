"""Downsample an h5ad to n_cells via contiguous per-group blocks (minimises HDF5 seeks).
Output: {input_dir}/{input_stem}_rand_{n}.h5ad."""

import sys
import argparse
import numpy as np
import anndata as ad
from pathlib import Path

repo_root = next(p for p in Path(__file__).resolve().parents if (p / ".git").exists())
sys.path.insert(0, str(repo_root))
from utils.logging_extra import log_run_command

_PREFERRED_KEYWORDS = ("batch", "sample", "donor", "id")
_MAX_CATEGORIES = 500


def _detect_grouping_col(obs):
    """
    Detect a suitable categorical column for grouping cells.
    Prefers columns with dtype "category".

    Parameters
    ----------
    obs : pandas.DataFrame
        The .obs table from an AnnData object.

    Returns
    -------
    str or None
        Column name, or None if no suitable column is found.
    """
    candidates = [
        c for c in obs.columns
        if hasattr(obs[c], "cat") and 2 <= obs[c].nunique() <= _MAX_CATEGORIES
    ]
    if not candidates:
        # fall back to object columns
        candidates = [
            c for c in obs.columns
            if obs[c].dtype == object and 2 <= obs[c].nunique() <= _MAX_CATEGORIES
        ]
    if not candidates:
        return None

    def _score(col):
        return any(kw in col.lower() for kw in _PREFERRED_KEYWORDS)

    return sorted(candidates, key=_score, reverse=True)[0]


def _sample_idx(obs, group_col, n_target, seed):
    """
    Sample cell indices as contiguous blocks within each group.

    Parameters
    ----------
    obs : pandas.DataFrame
        Obs table with a "_i" integer position column added.
    group_col : str or None
        Column to group by; if None, draws a random sample with no grouping.
    n_target : int
        Total number of cells to select.
    seed : int
        Random seed.

    Returns
    -------
    numpy.ndarray
        Sorted 1-D array of integer row positions, length <= n_target.
    """
    rng = np.random.default_rng(seed)

    if group_col is None:
        idx = rng.choice(len(obs), size=n_target, replace=False)
        return np.sort(idx)

    parts = []
    for _, grp in obs.groupby(group_col, observed=True):
        n = max(1, round(len(grp) / len(obs) * n_target))
        pos = grp["_i"].values
        start = rng.integers(0, max(1, len(pos) - n))
        parts.append(pos[start : start + n])

    return np.sort(np.concatenate(parts))[:n_target]


def downsample(infile, n_cells, seed, promote):
    """
    Downsample an h5ad file to n_cells using contiguous per-group blocks.

    Parameters
    ----------
    infile : Path
        Input h5ad path.
    n_cells : int
        Target number of cells.
    seed : int
        Random seed.
    promote : bool
        Write output file when True; dry-run when False.

    Returns
    -------
    None
    """
    outfile = infile.parent / f"{infile.stem}_rand_{n_cells}.h5ad"

    print(f"[in]  {infile}")
    print(f"[out] {outfile}")

    adata = ad.read_h5ad(infile, backed="r")
    print(f"[info] shape: {adata.shape}")

    obs = adata.obs.copy()
    obs["_i"] = np.arange(adata.n_obs)

    group_col = _detect_grouping_col(obs)
    if group_col:
        n_groups = obs[group_col].nunique()
        print(f"[info] grouping by '{group_col}' ({n_groups} groups)")
    else:
        print("[warn] no suitable categorical column found; using ungrouped random sample")

    idx = _sample_idx(obs, group_col, n_cells, seed)
    print(f"[info] selected {len(idx):,} cells")

    if not promote:
        print(f"[dry-run] would write {outfile}  (pass --promote to write)")
        adata.file.close()
        return

    print("[write] writing (gzip)...")
    adata[idx].write_h5ad(outfile, compression="gzip")
    adata.file.close()
    print(f"[done] {outfile}")


def _parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-i", "--infile-h5ad",
        required=True,
        help="Input h5ad file to downsample.",
    )
    parser.add_argument(
        "-n", "--number-of-cells",
        type=int,
        default=100_000,
        help="Target number of cells (default: 100000).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed (default: 42).",
    )
    parser.add_argument(
        "--promote",
        action="store_true",
        default=False,
        help="Write the output file. Without this flag, runs in dry-run mode.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    log_run_command(Path(args.infile_h5ad).parent / f"{Path(sys.argv[0]).stem}.log")
    downsample(
        infile=Path(args.infile_h5ad),
        n_cells=args.number_of_cells,
        seed=args.seed,
        promote=args.promote,
    )
