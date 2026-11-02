#!/usr/bin/env python3
"""
Build the per-spot Ma et al. SN Visium benchmark table: joins per-spot
BayesSpace cluster assignments (from 01_fetch_al_dalahmah_zenodo_data.py /
02_export_spot_metadata.R) onto the paper's per-cluster RCTD proportions
(Supplementary Data 1, sheet D).

CAVEAT: RCTD produces per-spot weights, but only this per-BayesSpace-cluster
average was deposited publicly -- every spot sharing a cluster gets an
identical proportion value below. Treat as cluster-resolution ground truth
(44 clusters across 10 samples), not true per-spot.

Note: in the workflow, we DON'T assign the per-spot winner across the Ma (2025)
13 cell-type columns. We instead use the continuous proportion values for
downstream analysis. The reason for this is because Dopaminergic_Neuron and
Fibrous_Astrocyte proportions are large enough in nearly every cluster to
swamp the rest, so they would always "win" if we took the argmax, which
would lead to flawed results.

Usage:
  python3 03_build_benchmark.py
"""

import os
import csv
import gzip
import argparse
import sys
from pathlib import Path
from datetime import datetime

import openpyxl

_DATA_DIR = "reference_building/benchmark/reference_data/ma/data"
_SPOT_METADATA_CSV = f"{_DATA_DIR}/reformatted/spot_metadata_all_samples.csv.gz"
_SUPP_XLSX = f"{_DATA_DIR}/original/41467_2025_62478_MOESM3_ESM.xlsx"
_SUPP_SHEET = "D"
_OUTPUT_CSV = f"{_DATA_DIR}/reformatted/per_spot_benchmark_al_dalahmah_sn.csv.gz"


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


def extract_bayes_cluster_table(xlsx_path: str, sheet_name: str) -> list[dict]:
    """
    Pull the "deconvolution of BayesSpace Clusters" table from Supplementary
    Data 1, sheet D.

    Parameters
    ----------
    xlsx_path : str
        Path to 41467_2025_62478_MOESM3_ESM.xlsx.
    sheet_name : str
        Worksheet holding the table (sheet "D").

    Returns
    -------
    list[dict]
        One dict per BayesSpace cluster, keyed by the table's own columns
        (Sample, Bayes_Sample, Number_of_Spots, the 13 cell-type
        proportions, Avg_TH_Expression).
    """
    wb = openpyxl.load_workbook(xlsx_path, read_only=True, data_only=True)
    rows = list(wb[sheet_name].iter_rows(values_only=True))

    header_row_idx = next(
        i for i, row in enumerate(rows) if row[0] == "Sample" and row[1] == "Bayes_Sample"
    )
    header = [h for h in rows[header_row_idx] if h is not None]

    cluster_rows = []
    for row in rows[header_row_idx + 1:]:
        if row[0] is None or row[1] is None or "_" not in str(row[1]):
            break
        cluster_rows.append(dict(zip(header, row[: len(header)])))
    return cluster_rows


def open_for_write(path: str):
    """Open `path` for text writing, gzip-compressed if it ends in `.gz`."""
    if path.endswith(".gz"):
        return gzip.open(path, "wt", newline="")
    return open(path, "w", newline="")


def open_for_read(path: str):
    """Open `path` for text reading, transparently gzip-decompressed if it ends in `.gz`."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, newline="")


def join_spots_to_clusters(spot_metadata: list[dict], cluster_rows: list[dict]) -> list[dict]:
    """Left-join each spot onto its BayesSpace cluster's RCTD proportions via `bayes_sample_key`."""
    decon_by_key = {r["Bayes_Sample"]: r for r in cluster_rows}
    decon_cols = [c for c in cluster_rows[0] if c not in ("Sample", "Bayes_Sample")]

    joined = []
    unmatched = 0
    for spot in spot_metadata:
        cluster = decon_by_key.get(spot["bayes_sample_key"])
        row = dict(spot)
        if cluster is None:
            unmatched += 1
            row.update({c: "" for c in decon_cols})
        else:
            row.update({c: cluster[c] for c in decon_cols})
        joined.append(row)

    if unmatched:
        print(f"  WARNING: {unmatched} / {len(spot_metadata)} spots had no matching cluster row")
    return joined


def main(args: argparse.Namespace) -> None:
    """Join per-spot metadata to the per-cluster RCTD deconvolution table and write the benchmark CSV."""
    with open_for_read(args.spot_metadata_csv) as fh:
        spot_metadata = list(csv.DictReader(fh))
    print(f"Loaded {len(spot_metadata)} spots from {args.spot_metadata_csv}")

    cluster_rows = extract_bayes_cluster_table(args.supplementary_xlsx, _SUPP_SHEET)
    print(f"Loaded {len(cluster_rows)} BayesSpace cluster rows from {args.supplementary_xlsx}")

    joined = join_spots_to_clusters(spot_metadata, cluster_rows)

    region_da = {}
    for row in joined:
        value = row.get("Dopaminergic_Neuron", "")
        if value in ("", "NA"):
            continue
        region_da.setdefault(row["region"], []).append(float(value))
    print("\nSanity check -- mean Dopaminergic_Neuron proportion by region:")
    for region, values in region_da.items():
        print(f"  {region}: n_spots={len(values)}, mean={sum(values) / len(values):.4f}")

    os.makedirs(os.path.dirname(_OUTPUT_CSV), exist_ok=True)
    with open_for_write(_OUTPUT_CSV) as fh:
        writer = csv.DictWriter(fh, fieldnames=list(joined[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(joined)
    print(f"\nWritten {_OUTPUT_CSV}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--spot-metadata-csv",
        default=_SPOT_METADATA_CSV,
        help=f"Per-spot metadata CSV from 02_export_spot_metadata.R.\n[Default: {_SPOT_METADATA_CSV}]",
    )
    parser.add_argument(
        "--supplementary-xlsx",
        default=_SUPP_XLSX,
        help=f"Supplementary Data 1 workbook fetched by 01_fetch_al_dalahmah_zenodo_data.py.\n"
             f"[Default: {_SUPP_XLSX}]",
    )
    args = parser.parse_args()

    repo_root = next(p for p in Path(__file__).resolve().parents if (p / ".git").exists())
    if Path.cwd() != repo_root:
        print(f"[INFO] Changing working directory to repo root: {repo_root}")
        os.chdir(repo_root)

    log_run_command(args)
    main(args)
