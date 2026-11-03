#!/usr/bin/env python3
"""
Fetch the two Zenodo artifacts (Ma et al. 2025 processed Seurat
objects + Supplementary Data 1) needed for downstream benchmark.
Avoids downloading either record's full multi-GB archive.
"""

import os
import sys
import zipfile
import argparse
import urllib.request
from pathlib import Path
from datetime import datetime

_DATA_DIR = "reference_building/benchmark/reference_data/ma/data"
_OUTPUT_DIR = f"{_DATA_DIR}/original"

_ST_OBJ_ZIP_URL = "https://zenodo.org/api/records/19454950/files/data.zip/content"
_ST_OBJ_ENTRY = "data/processed_st_obj.rds"

_SUPP_ZIP_URL = "https://zenodo.org/api/records/19457732/files/SupplementaryData.zip/content"
_SUPP_ENTRY = "Supplementary Data/41467_2025_62478_MOESM3_ESM.xlsx"


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


class _HTTPRangeFile:
    """
    Minimal seekable, read-only file-like object backed by HTTP range requests --
    lets `zipfile.ZipFile` extract entries from a remote zip without downloading it.
    """

    def __init__(self, url: str) -> None:
        self.url = url
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req) as resp:
            self.size = int(resp.headers["Content-Length"])
        self.pos = 0

    def seek(self, offset: int, whence: int = 0) -> int:
        if whence == 0:
            self.pos = offset
        elif whence == 1:
            self.pos += offset
        elif whence == 2:
            self.pos = self.size + offset
        return self.pos

    def tell(self) -> int:
        return self.pos

    def read(self, n: int = -1) -> bytes:
        end = (self.size - 1) if n < 0 else min(self.pos + n, self.size) - 1
        if end < self.pos:
            return b""
        req = urllib.request.Request(self.url, headers={"Range": f"bytes={self.pos}-{end}"})
        with urllib.request.urlopen(req) as resp:
            data = resp.read()
        self.pos += len(data)
        return data

    def seekable(self) -> bool:
        return True


def fetch_zip_entry(zip_url: str, entry_name: str, dest_path: str, use_cached: bool) -> None:
    """Extract a single named entry from a remote zip via range requests, writing it to dest_path."""
    if use_cached and os.path.exists(dest_path):
        size_mb = os.path.getsize(dest_path) / 1e6
        print(f"  Cached ({size_mb:.0f} MB): {dest_path}")
        return

    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    remote = _HTTPRangeFile(zip_url)
    with zipfile.ZipFile(remote) as zf:
        info = zf.getinfo(entry_name)
        print(
            f"  Extracting {entry_name} "
            f"({info.file_size / 1e6:.0f} MB uncompressed) → {dest_path}"
        )
        with zf.open(info) as src, open(dest_path, "wb") as dst:
            while True:
                chunk = src.read(8 * 1024 * 1024)
                if not chunk:
                    break
                dst.write(chunk)
    print(f"  Done: {os.path.getsize(dest_path) / 1e6:.0f} MB")


def main(args: argparse.Namespace) -> None:
    """Fetch Ma et al.'s processed Seurat object list and Supplementary Data 1 workbook."""
    st_obj_path = os.path.join(args.output_dir, "processed_st_obj.rds")
    supp_path = os.path.join(args.output_dir, "41467_2025_62478_MOESM3_ESM.xlsx")

    print(f"\n=== Fetching {_ST_OBJ_ENTRY} (Zenodo 19454950) ===")
    fetch_zip_entry(_ST_OBJ_ZIP_URL, _ST_OBJ_ENTRY, st_obj_path, args.use_cached)

    print(f"\n=== Fetching {_SUPP_ENTRY} (Zenodo 19457732) ===")
    fetch_zip_entry(_SUPP_ZIP_URL, _SUPP_ENTRY, supp_path, args.use_cached)

    print(
        f"\nDone.\n  Seurat object list: {st_obj_path}\n  Supplementary Data 1: {supp_path}\n"
        "Next: run 02_export_spot_metadata.R, then 03_build_benchmark.py."
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--output-dir",
        default=_OUTPUT_DIR,
        help="Directory to write the two fetched files into.\n"
             f"[Default: {_OUTPUT_DIR}]",
    )
    parser.add_argument(
        "--use-cached",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Skip re-downloading a file that already exists at its destination.\n"
             "Use --no-use-cached to force re-fetching both files.\n"
             "[Default: True]",
    )
    args = parser.parse_args()

    repo_root = next(p for p in Path(__file__).resolve().parents if (p / ".git").exists())
    if Path.cwd() != repo_root:
        print(f"[INFO] Changing working directory to repo root: {repo_root}")
        os.chdir(repo_root)

    log_run_command(args)
    main(args)
