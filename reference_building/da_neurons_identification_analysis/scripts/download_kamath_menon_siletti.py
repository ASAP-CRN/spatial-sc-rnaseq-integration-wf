#!/usr/bin/env python3
"""
Download and cache the raw Kamath 2022, Menon 2026, and Siletti 2023 H5ADs
from CELLxGENE, for later use by integrate_kamath_menon_siletti_asap.py.

Stores files under reference_building/da_neurons_identification_analysis/reference_data/ -- no filtering,
subsampling, or merging. ASAP PMDBS cohort is not downloaded here; the user
must provide its H5ADs and metadata CSV directly to the integration script.

Kamath: 9 files (7 cell-type files + 2 Nurr1-sorted fractions), from
https://cellxgene.cziscience.com/collections/b0f0b447-ac37-45b0-b1bf-5c0b7d871120

Menon: 1 file, selected via --menon-dataset-name (default 'Complete RNA'), from
https://cellxgene.cziscience.com/collections/16876983-d454-43db-a8ec-77fafca9bf38

Siletti: only the 4 files actually used downstream, from
https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443
- 3 "Dissection: <region>" files covering the midbrain (SN-RN, SN, RN) --
  see _SILETTI_NEEDED_FILE_STEMS here and downsample_kamath_menon_siletti_asap.py's
  _SILETTI_MIDBRAIN_REGIONS (keep both in sync if the region selection ever changes).
- "All neurons", the whole-atlas file used only to test barcode membership
  (neuron vs non-neuron) in load_siletti() -- never pooled into the modeled
  dataset, so this doesn't incur the "triple-counting" risk of downloading
  multiple full-coverage partitions (e.g. the 31 "Supercluster:" files, which
  together with the dissections and "All neurons"/"All non-neuronal cells"
  cover the identical ~3.37M-cell pool from different angles).

--use-cached-initial-h5ad (default True) skips re-downloading a source
whose cache directory already has file(s).
"""

import os
import re
import sys
import glob
import argparse
import requests
from pathlib import Path
from datetime import datetime

_CELLXGENE_API = "https://api.cellxgene.cziscience.com/curation/v1"

# Local cache directory for downloaded H5ADs
_REFERENCE_DATA_DIR = "reference_building/da_neurons_identification_analysis/reference_data"

# CELLxGENE collection IDs
_KAMATH_COLLECTION_ID  = "b0f0b447-ac37-45b0-b1bf-5c0b7d871120"
_MENON_COLLECTION_ID   = "16876983-d454-43db-a8ec-77fafca9bf38"
_SILETTI_COLLECTION_ID = "283d65eb-dd53-496d-adb7-7570c7caa443"

# Exactly the Siletti files load_siletti()/_load_all_neuron_barcodes() in
# downsample_kamath_menon_siletti_asap.py actually read -- see module
# docstring for why only these. Normalized-title stems (via _normalize_label).
_SILETTI_NEEDED_FILE_STEMS = {
    "dissection_midbrain_m_substantia_nigra_red_nucleus_and_nearby_nuclei_sn_rn",
    "dissection_midbrain_m_substantia_nigra_sn",
    "dissection_midbrain_rn_red_nucleus_rn",
    "all_neurons",
}


# ---------------------------------------------------------------------------
# Helpers
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
    """Log the invocation to {_REFERENCE_DATA_DIR}/{script_stem}.log and tee all subsequent output there."""
    log_file = Path(_REFERENCE_DATA_DIR) / f"{Path(__file__).stem}.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    f = open(log_file, "w", encoding="utf-8")  # noqa: SIM115 — intentionally kept open
    f.write(f"[{timestamp}] {' '.join(sys.argv)}\n\n")
    sys.stdout = _Tee(sys.stdout, f)
    sys.stderr = _Tee(sys.stderr, f)


def _normalize_label(label: str) -> str:
    """Normalize a CELLxGENE dataset title into a safe filename stem."""
    return re.sub(r"[^a-z0-9]+", "_", label.lower()).strip("_")


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

def _fetch_collection_h5ad_assets(collection_id: str) -> list[dict]:
    """Query the CELLxGENE curation API for a collection's H5AD assets (dataset_id, title, url, filesize)."""
    api_url = f"{_CELLXGENE_API}/collections/{collection_id}"
    resp = requests.get(api_url, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    assets = []
    for dataset in data.get("datasets", []):
        dataset_id = dataset.get("dataset_id", dataset.get("id", "unknown"))
        title = dataset.get("title", dataset_id)
        # API has used both 'assets' and 'dataset_assets' across versions.
        asset_list = dataset.get("assets", dataset.get("dataset_assets", []))
        for asset in asset_list:
            if asset.get("filetype", "").upper() == "H5AD":
                url = asset.get("url", asset.get("download_url", ""))
                assets.append({
                    "dataset_id": dataset_id,
                    "title": title,
                    "url": url,
                    "filesize": asset.get("filesize", 0),
                })
    return assets


def _download_file(url: str, dest_path: str, label: str, use_cached: bool) -> None:
    """Download a file with streaming progress; skip if already cached and use_cached=True."""
    if use_cached and os.path.exists(dest_path):
        size_mb = os.path.getsize(dest_path) / 1e6
        print(f"  Cached ({size_mb:.0f} MB): {dest_path}")
        return

    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    print(f"  Downloading {label} → {dest_path}")
    resp = requests.get(url, stream=True, timeout=60)
    resp.raise_for_status()

    total = int(resp.headers.get("content-length", 0))
    downloaded = 0
    chunk_size = 8 * 1024 * 1024  # 8 MB

    with open(dest_path, "wb") as fh:
        for chunk in resp.iter_content(chunk_size=chunk_size):
            if chunk:
                fh.write(chunk)
                downloaded += len(chunk)
                if total:
                    pct = downloaded / total * 100
                    print(
                        f"\r    {downloaded / 1e6:.0f} / {total / 1e6:.0f} MB"
                        f" ({pct:.0f}%)",
                        end="",
                        flush=True,
                    )
    if total:
        print()
    print(f"  Done: {downloaded / 1e6:.0f} MB")


def download_kamath_h5ads(use_cached: bool) -> str:
    """Download all Kamath H5ADs from CELLxGENE; files named by dataset_id."""
    kamath_dir = os.path.join(_REFERENCE_DATA_DIR, "kamath")
    os.makedirs(kamath_dir, exist_ok=True)

    existing = sorted(glob.glob(os.path.join(kamath_dir, "*.h5ad")))
    if use_cached and existing:
        print(f"  Cached: {len(existing)} files in {kamath_dir} — skipping API call")
        return kamath_dir

    print(f"  Fetching Kamath collection metadata (id: {_KAMATH_COLLECTION_ID})")
    assets = _fetch_collection_h5ad_assets(_KAMATH_COLLECTION_ID)
    if not assets:
        raise RuntimeError(
            f"No H5AD assets found in Kamath collection {_KAMATH_COLLECTION_ID}. "
            "Check the collection ID or the CELLxGENE API."
        )
    print(f"  Found {len(assets)} H5AD assets in Kamath collection")

    for asset in assets:
        filename = f"{_normalize_label(asset['title'])}.h5ad"
        dest = os.path.join(kamath_dir, filename)
        size_mb = asset["filesize"] / 1e6 if asset["filesize"] else 0
        label = f"{asset['title']} ({size_mb:.0f} MB)" if size_mb else asset["title"]
        _download_file(asset["url"], dest, label, use_cached)

    return kamath_dir


def download_menon_h5ad(
    use_cached: bool,
    dataset_name: str,
) -> str:
    """Download the Menon H5AD whose title contains `dataset_name` (case-insensitive)."""
    menon_dir = os.path.join(_REFERENCE_DATA_DIR, "menon")
    os.makedirs(menon_dir, exist_ok=True)

    existing = sorted(glob.glob(os.path.join(menon_dir, "*.h5ad")))
    if use_cached and existing:
        if len(existing) > 1:
            print(f"  Warning: {len(existing)} H5ADs in {menon_dir}; using first: {existing[0]}")
        else:
            size_mb = os.path.getsize(existing[0]) / 1e6
            print(f"  Cached ({size_mb:.0f} MB): {existing[0]}")
        return existing[0]

    print(f"  Fetching Menon collection metadata (id: {_MENON_COLLECTION_ID})")
    assets = _fetch_collection_h5ad_assets(_MENON_COLLECTION_ID)
    if not assets:
        raise RuntimeError(
            f"No H5AD assets found in Menon collection {_MENON_COLLECTION_ID}. "
            "Check the collection ID or the CELLxGENE API."
        )

    needle = dataset_name.lower()
    matches = [a for a in assets if needle in a["title"].lower()]
    if not matches:
        available = [a["title"] for a in assets]
        raise RuntimeError(
            f"No Menon dataset title contains '{dataset_name}'. "
            f"Available titles: {available}\n"
            "Use --menon-dataset-name to adjust the filter."
        )
    if len(matches) > 1:
        print(
            f"  Warning: {len(matches)} Menon datasets match '{dataset_name}'; "
            f"using first: '{matches[0]['title']}'"
        )
    asset = matches[0]

    filename = f"{_normalize_label(asset['title'])}.h5ad"
    dest = os.path.join(menon_dir, filename)
    size_mb = asset["filesize"] / 1e6 if asset["filesize"] else 0
    label = f"{asset['title']} ({size_mb:.0f} MB)" if size_mb else asset["title"]
    _download_file(asset["url"], dest, label, use_cached)
    return dest


def download_siletti_h5ads(use_cached: bool) -> str:
    """Download only the Siletti H5ADs this pipeline actually uses (see _SILETTI_NEEDED_FILE_STEMS)."""
    siletti_dir = os.path.join(_REFERENCE_DATA_DIR, "siletti")
    os.makedirs(siletti_dir, exist_ok=True)

    expected = {os.path.join(siletti_dir, f"{stem}.h5ad") for stem in _SILETTI_NEEDED_FILE_STEMS}
    existing = set(glob.glob(os.path.join(siletti_dir, "*.h5ad")))
    if use_cached and expected <= existing:
        print(f"  Cached: all {len(expected)} needed Siletti files already in {siletti_dir} — skipping API call")
        return siletti_dir

    print(f"  Fetching Siletti collection metadata (id: {_SILETTI_COLLECTION_ID})")
    assets = _fetch_collection_h5ad_assets(_SILETTI_COLLECTION_ID)
    if not assets:
        raise RuntimeError(
            f"No H5AD assets found in Siletti collection {_SILETTI_COLLECTION_ID}. "
            "Check the collection ID or the CELLxGENE API."
        )

    by_stem = {_normalize_label(a["title"]): a for a in assets}
    missing = _SILETTI_NEEDED_FILE_STEMS - by_stem.keys()
    if missing:
        raise RuntimeError(
            f"{len(missing)} needed Siletti dataset(s) not found in the collection: "
            f"{sorted(missing)}. CELLxGENE title formatting may have changed -- check "
            "_SILETTI_NEEDED_FILE_STEMS against the collection's current dataset titles."
        )
    print(
        f"  Found {len(assets)} H5AD assets in Siletti collection; "
        f"downloading the {len(_SILETTI_NEEDED_FILE_STEMS)} needed"
    )

    for stem in sorted(_SILETTI_NEEDED_FILE_STEMS):
        asset = by_stem[stem]
        dest = os.path.join(siletti_dir, f"{stem}.h5ad")
        size_mb = asset["filesize"] / 1e6 if asset["filesize"] else 0
        label = f"{asset['title']} ({size_mb:.0f} MB)" if size_mb else asset["title"]
        _download_file(asset["url"], dest, label, use_cached)

    return siletti_dir


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(args: argparse.Namespace) -> None:
    """Download and cache Kamath, Menon, and Siletti raw H5ADs."""
    print("\n=== Downloading Kamath H5ADs (9 files) ===")
    kamath_dir = download_kamath_h5ads(args.use_cached_initial_h5ad)
    print(f"  → {kamath_dir}")

    print("\n=== Downloading Menon H5AD ===")
    menon_path = download_menon_h5ad(args.use_cached_initial_h5ad, args.menon_dataset_name)
    print(f"  → {menon_path}")

    print("\n=== Downloading Siletti H5ADs (3 midbrain dissections + all_neurons) ===")
    siletti_dir = download_siletti_h5ads(args.use_cached_initial_h5ad)
    print(f"  → {siletti_dir}")

    print(
        "\nDone. Run integrate_kamath_menon_siletti_asap.py next to build "
        "intermediate/derived files and run the integration."
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--use-cached-initial-h5ad",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Skip re-downloading a source whose cache directory already has\n"
             "file(s). Use --no-use-cached-initial-h5ad to force re-downloading\n"
             "everything.\n"
             "[Default: True]",
    )
    parser.add_argument(
        "--menon-dataset-name",
        default="Complete RNA",
        help="Substring matched (case-insensitive) against dataset titles in the\n"
             "Menon CELLxGENE collection to select which H5AD to download.\n"
             "[Default: Complete RNA]",
    )
    args = parser.parse_args()

    repo_root = next(p for p in Path(__file__).resolve().parents if (p / ".git").exists())
    if Path.cwd() != repo_root:
        print(f"[INFO] Changing working directory to repo root: {repo_root}")
        os.chdir(repo_root)

    log_run_command(args)
    main(args)
