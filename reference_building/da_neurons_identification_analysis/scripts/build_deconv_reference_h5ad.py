#!/usr/bin/env python3
"""
Build a deconvolution-ready reference H5AD: .obs from a QC-filtered H5AD
(good cells, labels) + .X/.var from an all-gene H5AD (raw counts). HVG-only
references are inadequate for deconvolution since each method (Cell2location,
RCTD) selects its own marker genes from the full transcriptome.

Usage (dry-run, default)
------------------------
python3 reference_building/da_neurons_identification_analysis/scripts/build_deconv_reference_h5ad.py \\
    --final-h5ad  /path/to/asap-cohort.final.h5ad \\
    --full-h5ad   /path/to/asap-cohort.merged_cleaned_unfiltered.h5ad \\
    --output-h5ad /path/to/asap-cohort.deconv_reference.h5ad

Usage (write output)
--------------------
python3 reference_building/da_neurons_identification_analysis/scripts/build_deconv_reference_h5ad.py \\
    --final-h5ad  /path/to/asap-cohort.final.h5ad \\
    --full-h5ad   /path/to/asap-cohort.merged_cleaned_unfiltered.h5ad \\
    --output-h5ad /path/to/asap-cohort.deconv_reference.h5ad \\
    --promote
"""

import argparse
import os
import numpy as np
import anndata as ad
import h5py
import sys
from pathlib import Path
from scipy.sparse import csr_matrix

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "utils"))
from _reference_building_utils import extract_raw_counts, log_run_command

_CHUNK_SIZE = 10_000  # cells per backed-read slice


def _write_x_incremental(
    output_h5ad: str,
    full: ad.AnnData,
    target_indices: list,
    full_name: str,
) -> None:
    """
    Stream CSR data from `full` into an existing H5AD's X group in chunks
    (no chunk accumulation in memory). Expects an empty CSR shell already
    written by anndata; deletes and replaces it with real data.
    """
    n_cells = len(target_indices)
    n_genes = full.n_vars
    n_chunks = (n_cells + _CHUNK_SIZE - 1) // _CHUNK_SIZE

    with h5py.File(output_h5ad, "r+") as f:
        del f["X/data"]
        del f["X/indices"]
        del f["X/indptr"]

        ds_data    = f["X"].create_dataset("data",    shape=(0,), maxshape=(None,), dtype=np.int32, chunks=(1_000_000,), compression="gzip")
        ds_indices = f["X"].create_dataset("indices", shape=(0,), maxshape=(None,), dtype=np.int32, chunks=(1_000_000,), compression="gzip")

        indptr = [np.int64(0)]

        for chunk_num, start in enumerate(range(0, n_cells, _CHUNK_SIZE), 1):
            idx_chunk = target_indices[start:start + _CHUNK_SIZE]
            mat = extract_raw_counts(full[idx_chunk], source_name=full_name)

            nnz = mat.nnz
            old_size = ds_data.shape[0]
            if nnz > 0:
                ds_data.resize(old_size + nnz, axis=0)
                ds_data[old_size:]    = mat.data
                ds_indices.resize(old_size + nnz, axis=0)
                ds_indices[old_size:] = mat.indices

            # int64 to avoid overflow on large files
            offset = indptr[-1]
            indptr.extend((mat.indptr[1:].astype(np.int64) + offset).tolist())
            print(f"  chunk {chunk_num}/{n_chunks} — {min(start + _CHUNK_SIZE, n_cells):,}/{n_cells:,} cells written")

        f["X"].create_dataset("indptr", data=np.array(indptr, dtype=np.int64), compression="gzip")
        f["X"].attrs["shape"] = np.array([n_cells, n_genes], dtype=np.int64)


def build_reference(
    final_h5ad: str,
    full_h5ad: str,
    output_h5ad: str,
    promote: bool,
) -> None:
    """
    Combine QC-filtered obs with all-gene counts into a deconvolution reference.

    Parameters
    ----------
    final_h5ad : str
        QC-filtered H5AD; source of .obs and the authoritative barcode list.
    full_h5ad : str
        All-gene H5AD with raw counts; source of .X/.var. May include
        poor-QC cells, dropped via final_h5ad's barcode list.
    output_h5ad : str
        Destination path.
    promote : bool
        If False, dry-run only.

    Returns
    -------
    None
    """
    print(f"[INFO] Loading QC-filtered cell list from: {final_h5ad}")
    final = ad.read_h5ad(final_h5ad, backed="r")
    n_final = final.n_obs
    print(f"  {n_final:,} QC-filtered cells, {final.n_vars:,} genes (HVG subset)")

    print(f"[INFO] Loading full-gene matrix from: {full_h5ad}")
    full = ad.read_h5ad(full_h5ad, backed="r")
    print(f"  {full.n_obs:,} total cells, {full.n_vars:,} genes")

    final_barcodes = set(final.obs_names)
    full_barcodes  = set(full.obs_names)
    missing = final_barcodes - full_barcodes
    if missing:
        examples = sorted(missing)[:5]
        raise ValueError(
            f"{len(missing):,} cell barcode(s) from {os.path.basename(final_h5ad)} "
            f"are absent in {os.path.basename(full_h5ad)} "
            f"(e.g. {examples}). The two files must share the same cell universe."
        )
    print(f"  {n_final:,} / {n_final:,} cells matched")

    # Sorted indices for sequential HDF5 access (avoids random seeks)
    full_obs_list  = list(full.obs_names)
    full_obs_index = {bc: i for i, bc in enumerate(full_obs_list)}
    target_indices  = sorted(full_obs_index[bc] for bc in final_barcodes)
    target_barcodes = [full_obs_list[i] for i in target_indices]

    print(f"[INFO] Output shape: {n_final:,} cells × {full.n_vars:,} genes")

    if not promote:
        print(f"[DRY-RUN] Would write: {output_h5ad}")
        print("          Re-run with --promote to write the file.")
        return

    out_dir = os.path.dirname(output_h5ad)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Shell AnnData with an empty X, writes obs/var; X streamed in after
    print(f"[INFO] Writing obs/var shell to {output_h5ad} ...")
    shell = ad.AnnData(
        X=csr_matrix((n_final, full.n_vars), dtype=np.int32),
        obs=final.obs.loc[target_barcodes],
        var=full.var.copy(),
    )
    shell.write_h5ad(output_h5ad, compression="gzip")
    del shell

    full_name = os.path.basename(full_h5ad)
    n_chunks  = (n_final + _CHUNK_SIZE - 1) // _CHUNK_SIZE
    print(f"[INFO] Streaming X in {n_chunks} chunks into {output_h5ad} ...")
    _write_x_incremental(output_h5ad, full, target_indices, full_name)

    print(f"[DONE] {output_h5ad}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--final-h5ad",
        required=True,
        help=(
            "QC-filtered H5AD (good cells, HVG-only gene set).\n"
            "Provides .obs metadata and the authoritative barcode list.\n"
            "Example: /path/to/asap-cohort.final.h5ad"
        ),
    )
    parser.add_argument(
        "--full-h5ad",
        required=True,
        help=(
            "All-gene H5AD with raw integer counts (may include poor-QC cells).\n"
            "Provides .X and .var; poor-QC cells are filtered via --final-h5ad barcodes.\n"
            "Example: /path/to/asap-cohort.merged_cleaned_unfiltered.h5ad"
        ),
    )
    parser.add_argument(
        "--output-h5ad",
        required=True,
        help=(
            "Output path for the combined deconvolution reference H5AD\n"
            "(QC-filtered cells × all expressed genes, raw counts).\n"
            "Example: /path/to/asap-cohort.deconv_reference.h5ad"
        ),
    )
    parser.add_argument(
        "--promote",
        action="store_true",
        default=False,
        help="Write the output file. Default is dry-run (prints what would happen).",
    )

    args = parser.parse_args()

    args.final_h5ad  = os.path.abspath(args.final_h5ad)
    args.full_h5ad   = os.path.abspath(args.full_h5ad)
    args.output_h5ad = os.path.abspath(args.output_h5ad)
    args.output_dir  = os.path.dirname(args.output_h5ad)  # required by log_run_command
    log_run_command(args)

    repo_root_cwd = next(
        p for p in Path(__file__).resolve().parents if (p / ".git").exists()
    )
    if Path.cwd() != repo_root_cwd:
        print(f"[INFO] Changing working directory to repo root: {repo_root_cwd}")
        os.chdir(repo_root_cwd)

    build_reference(
        final_h5ad=args.final_h5ad,
        full_h5ad=args.full_h5ad,
        output_h5ad=args.output_h5ad,
        promote=args.promote,
    )
