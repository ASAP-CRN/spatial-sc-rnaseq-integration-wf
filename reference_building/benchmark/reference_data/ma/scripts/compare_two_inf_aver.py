#!/usr/bin/env python3
"""
Compare two Cell2location inf_aver.csv.gz reference signature files.

Reports row/column header overlap between the two files (rows = genes,
columns = cell types), then compares values at shared row/column positions:
counts exactly identical cells (after parsing, so differing text precision
in the source CSVs doesn't matter), and -- if any cell differs -- writes a
scatter plot of all shared-cell values (file 1 vs file 2, concatenated across
every row/column pair) annotated with the Pearson correlation coefficient
(PCC).

Outputs
-------
<outdir>/<label-1>__<label-2>.compare_two_inf_aver_summary.txt
<outdir>/<label-1>__<label-2>.compare_two_inf_aver_scatter.png   -- only written if any shared cell differs
<outdir>/compare_two_inf_aver.log
"""

import sys
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

repo_root = next(p for p in Path(__file__).resolve().parents if (p / ".git").exists())
sys.path.insert(0, str(repo_root))

from utils.logging_extra import log_run_command

_SCATTER_PNG_SUFFIX = "compare_two_inf_aver_scatter.png"
_SUMMARY_TXT_SUFFIX = "compare_two_inf_aver_summary.txt"
_MAX_PRINTED_LABELS = 20


def _load_inf_aver(infile_csv: str) -> pd.DataFrame:
    """
    Load an inf_aver.csv.gz reference signature file.

    Parameters
    ----------
    infile_csv : str
        Path to a Cell2location inf_aver.csv.gz file (rows = genes, indexed by
        the first column; columns = cell types).

    Returns
    -------
    pandas.DataFrame
        Rows x columns of float values, indexed by gene ID.
    """
    return pd.read_csv(infile_csv, index_col=0)


def _format_label_list(labels: list) -> str:
    """
    Format a list of row/column labels for display, truncating if long.

    Parameters
    ----------
    labels : list
        Labels to format.

    Returns
    -------
    str
        Comma-separated labels, truncated with a "... and N more" suffix
        beyond `_MAX_PRINTED_LABELS`.
    """
    if len(labels) <= _MAX_PRINTED_LABELS:
        return ", ".join(labels)
    shown = ", ".join(labels[:_MAX_PRINTED_LABELS])
    return f"{shown}, ... and {len(labels) - _MAX_PRINTED_LABELS} more"


def _summarize_headers(df_1: pd.DataFrame, df_2: pd.DataFrame) -> dict:
    """
    Summarize row/column header overlap between two inf_aver DataFrames.

    Parameters
    ----------
    df_1 : pandas.DataFrame
        First inf_aver table.
    df_2 : pandas.DataFrame
        Second inf_aver table.

    Returns
    -------
    dict
        n_rows_1 : int
            Row count in `df_1`.
        n_rows_2 : int
            Row count in `df_2`.
        rows_only_1 : list
            Row labels present only in `df_1`.
        rows_only_2 : list
            Row labels present only in `df_2`.
        shared_rows : list
            Row labels present in both, sorted.
        n_cols_1 : int
            Column count in `df_1`.
        n_cols_2 : int
            Column count in `df_2`.
        cols_only_1 : list
            Column labels present only in `df_1`.
        cols_only_2 : list
            Column labels present only in `df_2`.
        shared_cols : list
            Column labels present in both, sorted.
    """
    rows_1, rows_2 = set(df_1.index), set(df_2.index)
    cols_1, cols_2 = set(df_1.columns), set(df_2.columns)
    return {
        "n_rows_1": len(rows_1),
        "n_rows_2": len(rows_2),
        "rows_only_1": sorted(rows_1 - rows_2),
        "rows_only_2": sorted(rows_2 - rows_1),
        "shared_rows": sorted(rows_1 & rows_2),
        "n_cols_1": len(cols_1),
        "n_cols_2": len(cols_2),
        "cols_only_1": sorted(cols_1 - cols_2),
        "cols_only_2": sorted(cols_2 - cols_1),
        "shared_cols": sorted(cols_1 & cols_2),
    }


def _compare_shared_cells(df_1: pd.DataFrame, df_2: pd.DataFrame, shared_rows: list, shared_cols: list) -> dict:
    """
    Compare values at shared row/column positions between two inf_aver DataFrames.

    Parameters
    ----------
    df_1 : pandas.DataFrame
        First inf_aver table.
    df_2 : pandas.DataFrame
        Second inf_aver table.
    shared_rows : list
        Row labels present in both tables.
    shared_cols : list
        Column labels present in both tables.

    Returns
    -------
    dict
        n_cells : int
            Total shared-row x shared-column cells compared.
        n_identical : int
            Cells with an exactly identical (post-parsing) float value in both tables.
        values_1 : numpy.ndarray
            Flattened `df_1` values at every shared cell, in matched order.
        values_2 : numpy.ndarray
            Flattened `df_2` values at every shared cell, in the same order as `values_1`.
    """
    aligned_1 = df_1.loc[shared_rows, shared_cols]
    aligned_2 = df_2.loc[shared_rows, shared_cols]

    values_1 = aligned_1.to_numpy().ravel()
    values_2 = aligned_2.to_numpy().ravel()

    return {
        "n_cells": values_1.size,
        "n_identical": int((values_1 == values_2).sum()),
        "values_1": values_1,
        "values_2": values_2,
    }


def _plot_scatter(values_1: np.ndarray, values_2: np.ndarray, label_1: str, label_2: str, outfile_png: str) -> float:
    """
    Scatter-plot shared-cell values from both tables and annotate with PCC.

    Parameters
    ----------
    values_1 : numpy.ndarray
        Flattened values from the first table.
    values_2 : numpy.ndarray
        Flattened values from the second table, in the same cell order as `values_1`.
    label_1 : str
        Axis label for `values_1`.
    label_2 : str
        Axis label for `values_2`.
    outfile_png : str
        Path to write the scatter plot PNG to.

    Returns
    -------
    float
        Pearson correlation coefficient between `values_1` and `values_2`.
    """
    pcc, _ = pearsonr(values_1, values_2)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(values_1, values_2, s=4, alpha=0.3, linewidths=0)
    lims = [min(values_1.min(), values_2.min()), max(values_1.max(), values_2.max())]
    ax.plot(lims, lims, color="grey", linestyle="--", linewidth=1)
    ax.set_xlabel(label_1)
    ax.set_ylabel(label_2)
    ax.set_title(f"inf_aver shared-cell values (PCC = {pcc:.4f})")
    fig.tight_layout()
    fig.savefig(outfile_png, dpi=150)
    plt.close(fig)

    return pcc


def main(infile_csv_1: str, infile_csv_2: str, label_1: str, label_2: str, outdir: str) -> None:
    """
    Compare two Cell2location inf_aver.csv.gz files and report/plot the result.

    Parameters
    ----------
    infile_csv_1 : str
        Path to the first inf_aver.csv.gz file.
    infile_csv_2 : str
        Path to the second inf_aver.csv.gz file.
    label_1 : str
        Display label for `infile_csv_1` (plot axis / summary text).
    label_2 : str
        Display label for `infile_csv_2`.
    outdir : str
        Directory to write the summary text and scatter plot PNG to.

    Returns
    -------
    None
    """
    df_1 = _load_inf_aver(infile_csv_1)
    df_2 = _load_inf_aver(infile_csv_2)

    header_summary = _summarize_headers(df_1, df_2)
    cell_summary = _compare_shared_cells(
        df_1, df_2, header_summary["shared_rows"], header_summary["shared_cols"]
    )

    lines = [
        f"File 1 ({label_1}): {infile_csv_1}",
        f"File 2 ({label_2}): {infile_csv_2}",
        "",
        f"Rows    -- file 1: {header_summary['n_rows_1']:,}  file 2: {header_summary['n_rows_2']:,}  "
        f"shared: {len(header_summary['shared_rows']):,}  only-in-1: {len(header_summary['rows_only_1']):,}  "
        f"only-in-2: {len(header_summary['rows_only_2']):,}",
        f"Columns -- file 1: {header_summary['n_cols_1']:,}  file 2: {header_summary['n_cols_2']:,}  "
        f"shared: {len(header_summary['shared_cols']):,}  only-in-1: {len(header_summary['cols_only_1']):,}  "
        f"only-in-2: {len(header_summary['cols_only_2']):,}",
    ]
    if header_summary["rows_only_1"]:
        lines.append(f"  Rows only in file 1: {_format_label_list(header_summary['rows_only_1'])}")
    if header_summary["rows_only_2"]:
        lines.append(f"  Rows only in file 2: {_format_label_list(header_summary['rows_only_2'])}")
    if header_summary["cols_only_1"]:
        lines.append(f"  Columns only in file 1: {_format_label_list(header_summary['cols_only_1'])}")
    if header_summary["cols_only_2"]:
        lines.append(f"  Columns only in file 2: {_format_label_list(header_summary['cols_only_2'])}")

    lines.append("")
    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)
    stem = f"{label_1}__{label_2}"

    if cell_summary["n_cells"] == 0:
        lines.append("No shared rows/columns -- nothing to compare.")
    else:
        pct_identical = cell_summary["n_identical"] / cell_summary["n_cells"]
        lines.append(
            f"Shared-cell comparison (shared rows x shared columns): {cell_summary['n_cells']:,} cells total, "
            f"{cell_summary['n_identical']:,} identical ({pct_identical:.2%})"
        )
        if cell_summary["n_identical"] < cell_summary["n_cells"]:
            outfile_png = outdir_path / f"{stem}.{_SCATTER_PNG_SUFFIX}"
            pcc = _plot_scatter(cell_summary["values_1"], cell_summary["values_2"], label_1, label_2, str(outfile_png))
            lines.append(f"Not all shared cells are identical -- scatter plot written to {outfile_png} (PCC = {pcc:.4f})")
        else:
            lines.append("All shared cells are identical -- no scatter plot needed.")

    summary_text = "\n".join(lines)
    print(summary_text)

    outfile_summary = outdir_path / f"{stem}.{_SUMMARY_TXT_SUFFIX}"
    outfile_summary.write_text(summary_text + "\n")
    print(f"\nSummary written to {outfile_summary}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--infile-csv-1",
        required=True,
        help="Path to the first inf_aver.csv.gz file.",
    )
    parser.add_argument(
        "--infile-csv-2",
        required=True,
        help="Path to the second inf_aver.csv.gz file.",
    )
    parser.add_argument(
        "--label-1",
        default=None,
        help="Display label for --infile-csv-1 (plot axis / summary text).\n"
             "[Default: derived from --infile-csv-1's filename]",
    )
    parser.add_argument(
        "--label-2",
        default=None,
        help="Display label for --infile-csv-2 (plot axis / summary text).\n"
             "[Default: derived from --infile-csv-2's filename]",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Directory to write the summary text and scatter plot PNG to.",
    )
    args = parser.parse_args()

    label_1 = args.label_1 or Path(args.infile_csv_1).name
    label_2 = args.label_2 or Path(args.infile_csv_2).name

    log_run_command(Path(args.outdir) / "compare_two_inf_aver.log")
    main(
        infile_csv_1=args.infile_csv_1,
        infile_csv_2=args.infile_csv_2,
        label_1=label_1,
        label_2=label_2,
        outdir=args.outdir,
    )
