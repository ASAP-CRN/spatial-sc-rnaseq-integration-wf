"""
Shared input-location definitions for spatial-sc-rnaseq-integration-wf notebooks.

Centralizes dataset names, GCS/local path conventions, and filename suffixes
that recur across notebooks (ASAP raw/promoted H5ADs, Ma et al. 2025 benchmark H5ADs).
"""

# ----- ASAP (team-scherzer-pmdbs-spatial-visium-mtg) -----

ASAP_DATASET_NAME = "team-scherzer-pmdbs-spatial-visium-mtg"

# Raw/all-genes H5ADs
ASAP_RAW_H5AD_SUFFIX = ".cleaned_unfiltered.h5ad"

# cell2location deconvolution stage 3 outputs.
DECONVOLUTED_H5AD_SUFFIX = ".deconvoluted_clustered.h5ad"


def asap_raw_gcs_dir(dataset_name=ASAP_DATASET_NAME):
    """
    GCS directory for `dataset_name`'s raw (all-genes) H5ADs.

    Parameters
    ----------
    dataset_name : str
        ASAP dataset name.

    Returns
    -------
    str
    """
    return f"gs://asap-raw-{dataset_name}/workflow_execution/spatial_visium/preprocess/counts_to_adata/1.0.0"


def asap_raw_local_cache_dir(repo_root, dataset_name=ASAP_DATASET_NAME):
    """
    Local cache directory mirroring `asap_raw_gcs_dir`'s layout.

    Parameters
    ----------
    repo_root : pathlib.Path
        Repo root, from `get_repo_root`.
    dataset_name : str
        ASAP dataset name.

    Returns
    -------
    str
    """
    return str(repo_root / "spatial_samples" / f"asap-raw-{dataset_name}"
               / "workflow_execution" / "spatial_visium" / "preprocess" / "counts_to_adata" / "1.0.0")


def asap_release_gcs_dir(bucket_env, release_version, dataset_name=ASAP_DATASET_NAME):
    """
    GCS directory for `dataset_name`'s promoted `deconvolution_analysis` outputs.

    Parameters
    ----------
    bucket_env : str
        Bucket environment -- "dev", "uat", or "curated".
    release_version : str
        ASAP-CRN release version (e.g. "v6.0.0").
    dataset_name : str
        ASAP dataset name.

    Returns
    -------
    str
    """
    return f"gs://asap-{bucket_env}-{dataset_name}/spatial_visium/release/{release_version}/deconvolution_analysis"


def asap_strategy3_local_cache_dir(repo_root, dataset_name=ASAP_DATASET_NAME):
    """
    Local cache directory for `dataset_name`'s Strategy 3 deconvolution outputs.

    Parameters
    ----------
    repo_root : pathlib.Path
        Repo root, from `get_repo_root`.
    dataset_name : str
        ASAP dataset name.

    Returns
    -------
    str
    """
    return str(repo_root / "spatial_samples" / f"asap-raw-{dataset_name}"
               / "workflow_execution" / "spatial_visium" / "cell2location" / "Strategy3")


# ----- Ma et al. 2025 benchmark -----

MA_SPATIAL_H5AD_SUFFIX = ".spatial.h5ad"
MA_STRATEGIES = (1, 2, 3)


def ma_base_dir(repo_root):
    """
    Base directory for the Ma (2025) benchmark's reformatted data.

    Parameters
    ----------
    repo_root : pathlib.Path
        Repo root, from `get_repo_root`.

    Returns
    -------
    pathlib.Path
    """
    return repo_root / "reference_building" / "benchmark" / "reference_data" / "ma" / "data" / "reformatted"


def ma_spatial_h5ad_dir(repo_root):
    """
    Directory of Ma (2025)'s raw (all-genes) spatial H5ADs.

    Parameters
    ----------
    repo_root : pathlib.Path
        Repo root, from `get_repo_root`.

    Returns
    -------
    pathlib.Path
    """
    return ma_base_dir(repo_root) / "spatial_h5ad"


def ma_strategy_dirname(strategy_n):
    """
    Subdirectory name for one reference-build strategy's Ma deconvolution outputs.

    Parameters
    ----------
    strategy_n : int
        Strategy number (1, 2, or 3).

    Returns
    -------
    str
    """
    return f"cell2location_Strategy{strategy_n}"


def ma_benchmark_csv_path(repo_root):
    """
    Path to Ma (2025)'s per-spot RCTD benchmark CSV.

    Parameters
    ----------
    repo_root : pathlib.Path
        Repo root, from `get_repo_root`.

    Returns
    -------
    str
    """
    return str(ma_base_dir(repo_root) / "per_spot_benchmark_al_dalahmah_sn.csv.gz")
