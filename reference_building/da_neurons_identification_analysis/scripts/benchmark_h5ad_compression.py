#!/usr/bin/env python3
"""Benchmark write time, read time, and file size for H5AD compression filters (gzip, lzf)
and the Zarr store format.

A test run with this ~300k cells H5AD:
gs://asap-wdl-ci/workflow_execution/spatial-sc-rnaseq-integration-wf/reference_building/reference_h5ad/deconvolution_reference_300k_c2lUda_Kamath200Rand.h5ad
and these commands:
`python3 /path/to/benchmark_h5ad_compression.py --input-h5ad /path/to/deconvolution_reference_300k_c2lUda_Kamath200Rand.h5ad --outdir /path/to/outdir/`
Produced these results in a MacBook Pro with 36GB RAM:
┌─────────────┬───────────┬──────────┬─────────┐
│ compression │ write (s) │ read (s) │  size   │
├─────────────┼───────────┼──────────┼─────────┤
│ gzip_level0 │ 3.3       │ 3.2      │ 6.24 GB │
├─────────────┼───────────┼──────────┼─────────┤
│ gzip_level4 │ 87.6      │ 15.0     │ 2.02 GB │
├─────────────┼───────────┼──────────┼─────────┤
│ gzip_level9 │ 2162.2    │ 13.9     │ 1.86 GB │
├─────────────┼───────────┼──────────┼─────────┤
│ lzf         │ 9.7       │ 5.7      │ 4.09 GB │
├─────────────┼───────────┼──────────┼─────────┤
│ zarr        │ 3.4       │ 2.2      │ 1.53 GB │
└─────────────┴───────────┴──────────┴─────────┘

"""

import argparse
import os
import shutil
import sys
import time
from pathlib import Path

import anndata as ad

repo_root = next(p for p in Path(__file__).resolve().parents if (p / ".git").exists())
sys.path.insert(0, str(repo_root))
from utils.logging_extra import log_run_command


def _human_size(num_bytes: float) -> str:
    """Format a byte count as a human-readable string (e.g. "1.23 GB")."""
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if num_bytes < 1024:
            return f"{num_bytes:.2f} {unit}"
        num_bytes /= 1024
    return f"{num_bytes:.2f} PB"


def _path_size(path: str) -> int:
    """Return total size in bytes; sums recursively if `path` is a directory (e.g. a Zarr store)."""
    if os.path.isfile(path):
        return os.path.getsize(path)
    return sum(
        os.path.getsize(os.path.join(dirpath, f))
        for dirpath, _, filenames in os.walk(path)
        for f in filenames
    )


def _remove_path(path: str) -> None:
    """Remove `path`, whether it's a file or a directory (e.g. a Zarr store)."""
    if os.path.isdir(path):
        shutil.rmtree(path)
    else:
        os.remove(path)


def benchmark_compression(
    input_h5ad: str,
    outdir: str,
    compressions: list,
    gzip_levels: list,
    n_read_repeats: int,
    keep_outputs: bool,
) -> list:
    """
    Re-write `input_h5ad` under each compression filter and measure write time,
    read time, and resulting file size.

    Parameters
    ----------
    input_h5ad : str
        Source H5AD to re-write under each compression filter.
    outdir : str
        Directory to write benchmark output files.
    compressions : list of str or None
        Filters to test (e.g. ["gzip", "lzf", "zarr"]). Each "gzip" entry is expanded
        into one run per level in `gzip_levels`. "zarr" writes/reads a Zarr store (via
        `write_zarr`/`read_zarr`) instead of an H5AD; its "file size" is the Zarr store
        directory's total size.
    gzip_levels : list of int
        compression_opts levels (0-9). h5py's gzip default is 4).
    n_read_repeats : int
        Read passes per compression, averaged to reduce OS-cache noise.
    keep_outputs : bool
        If False, delete each written file/store after benchmarking it.

    Returns
    -------
    list of dict
        One result dict per run: compression, write_seconds, read_seconds, file_size_bytes.
    """
    print(f"[INFO] Loading source: {input_h5ad}")
    t0 = time.perf_counter()
    adata = ad.read_h5ad(input_h5ad)
    print(f"  Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes in {time.perf_counter() - t0:.1f}s")

    os.makedirs(outdir, exist_ok=True)

    # Expand each "gzip" entry into one (label, compression, compression_opts) spec per level;
    # other entries ("lzf", "zarr", None) get a single spec with compression_opts=None.
    specs = []
    for compression in compressions:
        is_gzip = isinstance(compression, str) and compression.lower() == "gzip"
        if is_gzip:
            for level in gzip_levels:
                specs.append((f"gzip_level{level}", "gzip", level))
        else:
            label = "zarr" if (isinstance(compression, str) and compression.lower() == "zarr") else (compression or "none")
            specs.append((label, compression, None))

    results = []

    for label, compression, compression_opts in specs:
        is_zarr = compression == "zarr"
        ext = ".zarr" if is_zarr else ".h5ad"
        out_path = os.path.join(outdir, f"benchmark_{label}{ext}")

        print(f"\n[INFO] {label}")
        t0 = time.perf_counter()
        if is_zarr:
            adata.write_zarr(out_path)
        elif compression_opts is not None:
            adata.write_h5ad(out_path, compression=compression, compression_opts=compression_opts)
        else:
            adata.write_h5ad(out_path, compression=compression)
        write_seconds = time.perf_counter() - t0
        file_size_bytes = _path_size(out_path)
        print(f"  Write: {write_seconds:.1f}s, size: {_human_size(file_size_bytes)}")

        read_times = []
        for _ in range(n_read_repeats):
            t0 = time.perf_counter()
            ad.read_zarr(out_path) if is_zarr else ad.read_h5ad(out_path)
            read_times.append(time.perf_counter() - t0)
        read_seconds = sum(read_times) / len(read_times)
        formatted_times = ", ".join(f"{t:.1f}" for t in read_times)
        print(f"  Read ({n_read_repeats}x avg): {read_seconds:.1f}s  (individual: {formatted_times})")

        results.append({
            "compression": label,
            "write_seconds": write_seconds,
            "read_seconds": read_seconds,
            "file_size_bytes": file_size_bytes,
        })

        if not keep_outputs:
            _remove_path(out_path)
            print(f"  [cleanup] removed {out_path}")

    print("\n[SUMMARY]")
    print(f"{'compression':<12} {'write (s)':<12} {'read (s)':<12} {'size':<12}")
    for r in results:
        print(f"{r['compression']:<12} {r['write_seconds']:<12.1f} {r['read_seconds']:<12.1f} "
              f"{_human_size(r['file_size_bytes']):<12}")

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--input-h5ad",
        default="reference_building/reference_h5ad/deconvolution_reference_300k_c2lUda_Kamath200Rand.h5ad",
        help="Source H5AD to benchmark.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Directory to write benchmark output files (deleted afterward unless --keep-outputs).",
    )
    parser.add_argument(
        "--compressions",
        nargs="+",
        default=["gzip", "lzf", "zarr"],
        help="Compression filters to test, space-separated. Use 'none' for no compression,\n"
             "'zarr' to benchmark the Zarr store format instead (write_zarr/read_zarr).\n"
             "Default: gzip lzf zarr",
    )
    parser.add_argument(
        "--gzip-levels",
        type=int,
        nargs="+",
        default=[0, 4, 9],
        help="compression_opts levels (0-9) to benchmark for each 'gzip' entry in\n"
             "--compressions. h5py's gzip default (no compression_opts) is 4.\n"
             "Default: 0 4 9",
    )
    parser.add_argument(
        "--n-read-repeats",
        type=int,
        default=3,
        help="Read passes per compression, averaged to reduce OS-cache noise. Default: 3",
    )
    parser.add_argument(
        "--keep-outputs",
        action="store_true",
        default=False,
        help="Keep the written benchmark files instead of deleting them after measurement.",
    )
    args = parser.parse_args()

    compressions = [None if c.lower() == "none" else c for c in args.compressions]

    log_run_command(Path(args.outdir) / f"{Path(sys.argv[0]).stem}.log")

    benchmark_compression(
        input_h5ad=args.input_h5ad,
        outdir=args.outdir,
        compressions=compressions,
        gzip_levels=args.gzip_levels,
        n_read_repeats=args.n_read_repeats,
        keep_outputs=args.keep_outputs,
    )
