#!/usr/bin/env python3
"""
Submit Cell2location Stages 2+3 (spatial mapping + clustering) for one or more
spatial samples against an already-built Stage 1 reference signature, in
chunks, via the omics Workbench CLI.

All sample IDs are read from --infile-sample-ids-all. Samples already processed
(done) are subtracted from either --infile-sample-ids-done or automatically from
the gs:// bucket looking for and validating the presence of expected meta.csv sidecar files.
The remainder sample ID's are split into chunks of --chunk-size.

_MAX_CHUNK_SIZE caps the chunk size to the observed GPU quota (see that constant's docstring for why).

A workflow-params JSON is always written for every chunk, under
<outdir>/<run_timestamp>/chunk_NNN.json; --submit-chunks selects which of
those chunks are actually submitted via `omics workbench runs submit`.

--promote gates whether that submission actually happens (default is
dry-run: every chunk JSON is written, nothing is submitted).

Outputs
-------
<outdir>/<run_timestamp>/chunk_NNN.json              -- one per chunk
<outdir>/submit_cell2location_2and3_prod.log
"""

import sys
import argparse
import json
import re
import shutil
import subprocess
from pathlib import Path
from datetime import datetime, timezone

repo_root = next(p for p in Path(__file__).resolve().parents if (p / ".git").exists())
sys.path.insert(0, str(repo_root))

from utils.logging_extra import log_run_command

_CHUNK_FILENAME_TEMPLATE = "chunk_{index:03d}.json"
_DEFAULT_OUTDIR = "prod_runs"
_DEFAULT_BATCH_KEY = "batch_id"
_DEFAULT_CELL_TYPE_KEY = "harmonized_cell_type_siletti_allen"
_DEFAULT_SUBMIT_CHUNKS = "all"

# GCP on-demand NVIDIA_T4_GPUS quota observed at 16 in us-central1 (2026-11-05);
# cell2location_deconvolute_visium requests 1 T4 per sample.
_MAX_CHUNK_SIZE = 15

# Matches deconvolution.wdl' PMBDS > PMDBS fix
_PMBDS_PATTERN = re.compile(r"^ASAP_PMBDS_")

# Workflow version for outdir like:
# gs://asap-raw-team-scherzer-pmdbs-spatial-visium-mtg/workflow_execution/spatial_visium/deconvolution_analysis/{build_reference|cluster_visium|deconvolute_visium}/1.0.0/
_DEFAULT_CELL2LOCATION_VERSION = "1.0.0"

# Stage 3 (cluster_visium) expected outfiles.
# To count a sample as done, the corresponding *.meta.tsv files of these
# must exist and have the expected content.
_CLUSTER_VISIUM_SUFFIXES = (
    "deconvoluted_clustered.h5ad",
    "deconvolution_clustering_metadata.csv.gz",
    "timing.tsv",
)
_EXPECTED_WORKFLOW_NAME = "spatial_sc_rnaseq_integration"


def _resolve_repo_path(path: str) -> str:
    """
    Resolve `path` against the repo root if relative, else return it unchanged.

    Parameters
    ----------
    path : str
        Path from a JSON config; relative (to the repo root) or absolute.

    Returns
    -------
    str
        Absolute path.
    """
    p = Path(path)
    return str(p if p.is_absolute() else (repo_root / p).resolve())


def _read_sample_ids(infile_sample_ids: str) -> list:
    """
    Read sample IDs, one per line; blank lines and '#'-comments are skipped.

    Parameters
    ----------
    infile_sample_ids : str
        Path to a plain-text file of sample IDs.

    Returns
    -------
    list
        Sample IDs in file order; duplicates are dropped (first occurrence kept).
    """
    seen = set()
    sample_ids = []
    with open(infile_sample_ids) as f:
        for line in f:
            sample_id = line.strip()
            if not sample_id or sample_id.startswith("#"):
                continue
            if sample_id in seen:
                print(f"[INFO] Duplicate sample_id {sample_id!r} in {infile_sample_ids}; keeping first occurrence.")
                continue
            seen.add(sample_id)
            sample_ids.append(sample_id)
    return sample_ids


def _chunk(items: list, chunk_size: int) -> list:
    """Split `items` into consecutive lists of at most `chunk_size` each."""
    return [items[i:i + chunk_size] for i in range(0, len(items), chunk_size)]


def _normalize_sample_id(sample_id: str) -> str:
    """
    Normalize a sample ID to the corrected 'PMDBS' spelling for comparison.

    Parameters
    ----------
    sample_id : str
        Sample ID, possibly using the original (uncorrected) 'PMBDS' spelling.

    Returns
    -------
    str
        `sample_id` with any leading 'ASAP_PMBDS_' replaced by 'ASAP_PMDBS_';
        unchanged if it doesn't start with 'ASAP_PMBDS_'.
    """
    return _PMBDS_PATTERN.sub("ASAP_PMDBS_", sample_id)


def _resolve_remaining_sample_ids(all_ids: list, done_ids: set) -> list:
    """
    Subtract already-done sample IDs from the full list, preserving order.

    Comparison is done on the PMBDS/PMDBS-normalized form.

    Parameters
    ----------
    all_ids : list
        Every sample ID to eventually run, in order.
    done_ids : set
        Sample IDs already completed (any spelling); excluded from the result.

    Returns
    -------
    list
        `all_ids`, in original spelling and order, minus anything matching
        (after normalization) an entry in `done_ids`.
    """
    done_ids_normalized = {_normalize_sample_id(sample_id) for sample_id in done_ids}
    remaining_ids = [
        sample_id for sample_id in all_ids
        if _normalize_sample_id(sample_id) not in done_ids_normalized
    ]
    print(f"[INFO] {len(all_ids):,} total sample ID(s); {len(done_ids):,} already done; "
          f"{len(remaining_ids):,} remaining")
    return remaining_ids


def _cluster_visium_gcs_dir(config: dict) -> str:
    """
    Build the GCS directory for cluster_visium (Stage 3) outputs, matching
    deconvolution.wdl's own raw_data_path construction.

    Parameters
    ----------
    config : dict
        Parsed --input-json config.

    Returns
    -------
    str
        GCS directory path, no trailing slash.
    """
    cell2location_version = config.get("cell2location_version", _DEFAULT_CELL2LOCATION_VERSION)
    bucket = config["cohort_raw_data_bucket_spatial_sc_rnaseq"].rstrip("/")
    return f"{bucket}/workflow_execution/spatial_visium/deconvolution_analysis/cluster_visium/{cell2location_version}"


def _list_bucket_objects(gcs_dir: str) -> set:
    """
    List object basenames directly under a GCS directory via `gcloud storage ls`.

    Parameters
    ----------
    gcs_dir : str
        GCS directory path (no trailing slash required).

    Returns
    -------
    set
        Basenames (not full paths) of objects found; empty if the directory
        doesn't exist or has no objects.
    """
    result = subprocess.run(
        ["gcloud", "storage", "ls", gcs_dir.rstrip("/") + "/*"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        return set()
    return {line.strip().rsplit("/", 1)[-1] for line in result.stdout.splitlines() if line.strip()}


def _validate_meta_tsv_lines(gcs_paths: list, expected_filenames: list) -> list:
    """
    Fetch multiple .meta.tsv sidecars in one `gcloud storage cat` call and
    check each has well-formed upload_outputs content.

    Parameters
    ----------
    gcs_paths : list
        Full GCS paths to .meta.tsv files, in order.
    expected_filenames : list
        Same length as `gcs_paths`; each file's own first TSV field must
        match the corresponding entry.

    Returns
    -------
    list
        Parallel list of bool, one per `gcs_paths` entry; True if that file
        has exactly 6 tab-separated fields, its filename field matches, and
        its workflow-name field matches `_EXPECTED_WORKFLOW_NAME`.
    """
    if not gcs_paths:
        return []

    result = subprocess.run(["gcloud", "storage", "cat"] + gcs_paths, capture_output=True, text=True)
    lines = result.stdout.splitlines()
    if len(lines) != len(gcs_paths):
        print(f"  [WARNING] Bulk .meta.tsv fetch returned {len(lines)} line(s) for {len(gcs_paths)} "
              "file(s); falling back to one call per file.")
        lines = []
        for path in gcs_paths:
            single = subprocess.run(["gcloud", "storage", "cat", path], capture_output=True, text=True)
            lines.append(single.stdout.strip() if single.returncode == 0 else "")

    valid = []
    for line, expected_filename in zip(lines, expected_filenames):
        fields = line.strip().split("\t")
        if len(fields) != 6:
            valid.append(False)
            continue
        filename, md5_hash, timestamp, workflow_name, workflow_version, workflow_release = fields
        valid.append(
            filename == expected_filename
            and bool(md5_hash)
            and bool(timestamp)
            and workflow_name == _EXPECTED_WORKFLOW_NAME
            and bool(workflow_release)
        )
    return valid


def _find_done_sample_ids_from_bucket(cluster_visium_gcs_dir: str) -> list:
    """
    Determine which sample IDs are done by checking 
    cluster_visium's .meta.tsv sidecars directly in the raw bucket.

    A sample counts as done only if all 3 expected .meta.tsv files
    (`_CLUSTER_VISIUM_SUFFIXES`) exist and each has well-formed content
    (see `_validate_meta_tsv_lines`).

    Parameters
    ----------
    cluster_visium_gcs_dir : str
        GCS directory for Stage 3 outputs (see `_cluster_visium_gcs_dir`).

    Returns
    -------
    list
        Sample IDs (corrected PMDBS spelling, matching cluster_visium's own
        output filenames) confirmed done.
    """
    print(f"[INFO] Listing {cluster_visium_gcs_dir} ...")
    objects = _list_bucket_objects(cluster_visium_gcs_dir)
    print(f"  {len(objects):,} object(s) found")

    candidate_sample_ids = set()
    for obj in objects:
        for suffix in _CLUSTER_VISIUM_SUFFIXES:
            meta_suffix = f".{suffix}.meta.tsv"
            if obj.endswith(meta_suffix):
                candidate_sample_ids.add(obj[:-len(meta_suffix)])

    # Only validate candidates where all 3 .meta.tsv files are present at all.
    complete_candidates = [
        sample_id for sample_id in sorted(candidate_sample_ids)
        if all(f"{sample_id}.{suffix}.meta.tsv" in objects for suffix in _CLUSTER_VISIUM_SUFFIXES)
    ]

    gcs_paths = []
    expected_filenames = []
    for sample_id in complete_candidates:
        for suffix in _CLUSTER_VISIUM_SUFFIXES:
            filename = f"{sample_id}.{suffix}"
            gcs_paths.append(f"{cluster_visium_gcs_dir}/{filename}.meta.tsv")
            expected_filenames.append(filename)

    print(f"[INFO] Validating {len(gcs_paths):,} .meta.tsv file(s) for {len(complete_candidates):,} "
          "candidate(s) in one batch ...")
    validity = _validate_meta_tsv_lines(gcs_paths, expected_filenames)

    n_per_sample = len(_CLUSTER_VISIUM_SUFFIXES)
    done_sample_ids = []
    for i, sample_id in enumerate(complete_candidates):
        sample_validity = validity[i * n_per_sample:(i + 1) * n_per_sample]
        if all(sample_validity):
            done_sample_ids.append(sample_id)
        else:
            print(f"  [WARNING] {sample_id}: at least one .meta.tsv failed content validation")

    print(f"[INFO] {len(done_sample_ids):,} sample(s) confirmed done via bucket .meta.tsv validation")
    return done_sample_ids


def _parse_submit_chunks(submit_chunks: str, n_chunks: int) -> set:
    """
    Parse --submit-chunks into a 1-indexed set of chunk numbers to submit.

    Parameters
    ----------
    submit_chunks : str
        "all", or a comma-separated list of 1-indexed chunk numbers (e.g. "1" or "1,2").
    n_chunks : int
        Total number of chunks actually produced, for range validation.

    Returns
    -------
    set
        1-indexed chunk numbers to submit.

    Raises
    ------
    ValueError
        A parsed chunk number is out of range [1, n_chunks].
    """
    if submit_chunks == _DEFAULT_SUBMIT_CHUNKS:
        return set(range(1, n_chunks + 1))

    chunk_numbers = set()
    for token in submit_chunks.split(","):
        n = int(token.strip())
        if not (1 <= n <= n_chunks):
            raise ValueError(f"--submit-chunks value {n} is out of range [1, {n_chunks}]")
        chunk_numbers.add(n)
    return chunk_numbers


def _build_chunk_params(sample_ids: list, config: dict) -> dict:
    """
    Build the spatial_sc_rnaseq_integration workflow-params dict for one chunk.

    Parameters
    ----------
    sample_ids : list
        Sample IDs to include in this chunk.
    config : dict
        Parsed --input-json config; see argparse help for required/optional fields.

    Returns
    -------
    dict
        Ready to json.dump as a workflow-params file.
    """
    samples = [
        {
            "sample_id": sample_id,
            "qc_adata_object": config["qc_adata_object_template"].format(sample_id=sample_id),
        }
        for sample_id in sample_ids
    ]
    params = {
        "spatial_sc_rnaseq_integration.sc_rnaseq_reference_cell2location_model":
            config["sc_rnaseq_reference_cell2location_model"],
        "spatial_sc_rnaseq_integration.spatial_visium_projects": [
            {
                "asap_team_id": config["asap_team_id"],
                "asap_dataset_id": config["asap_dataset_id"],
                "asap_dataset_doi_url": config["asap_dataset_doi_url"],
                "samples": samples,
                "raw_data_bucket": config["raw_data_bucket"],
            }
        ],
        "spatial_sc_rnaseq_integration.batch_key": config.get("batch_key", _DEFAULT_BATCH_KEY),
        "spatial_sc_rnaseq_integration.cell_type_key": config.get("cell_type_key", _DEFAULT_CELL_TYPE_KEY),
        "spatial_sc_rnaseq_integration.run_stages": "2,3",
        "spatial_sc_rnaseq_integration.cohort_raw_data_bucket_spatial_sc_rnaseq":
            config["cohort_raw_data_bucket_spatial_sc_rnaseq"],
        "spatial_sc_rnaseq_integration.container_registry": config["container_registry"],
    }
    if "preemptible" in config:
        params["spatial_sc_rnaseq_integration.preemptible"] = config["preemptible"]
    return params


def submit_chunk(chunk_json: Path, workflow_id: str, workflow_version_name: str, promote: bool) -> None:
    """
    Submit one chunk's workflow-params JSON via the omics Workbench CLI.

    Parameters
    ----------
    chunk_json : pathlib.Path
        Path to the chunk's workflow-params JSON (already written to disk).
    workflow_id : str
        Workbench workflow ID (e.g. from `omics workbench workflows create`).
    workflow_version_name : str
        Workbench workflow registration version (e.g. "v1.0.0").
    promote : bool
        If False, print the command instead of running it.

    Returns
    -------
    None
    """
    cmd = [
        "omics", "workbench", "runs", "submit",
        "--workflow", workflow_id,
        "--version", workflow_version_name,
        "--workflow-params", f"@{chunk_json}",
    ]
    if not promote:
        print(f"[DRY-RUN] Would run: {' '.join(cmd)}")
        return

    if shutil.which("omics") is None:
        raise RuntimeError("'omics' CLI not found on PATH. Activate the venv with the Workbench CLI installed.")

    print(f"[INFO] Submitting {chunk_json} ...")
    subprocess.run(cmd, check=True)


def main(
    infile_sample_ids_all: str,
    infile_sample_ids_done: str,
    validate_bucket_meta_tsv: bool,
    chunk_size: int,
    submit_chunks: str,
    config: dict,
    outdir: str,
    workflow_id: str,
    workflow_version_name: str,
    promote: bool,
) -> None:
    """
    Write one workflow-params JSON per chunk of sample IDs, submitting the
    chunks selected by `submit_chunks`.

    Parameters
    ----------
    infile_sample_ids_all : str
        Path to a plain-text file of every sample ID to eventually run.
    infile_sample_ids_done : str or None
        Subtracted from `infile_sample_ids_all`; mutually exclusive with `validate_bucket_meta_tsv`.
    validate_bucket_meta_tsv : bool
        Derive done samples from cluster_visium's bucket .meta.tsv sidecars instead.
    chunk_size : int
        Maximum number of samples per chunk/submission.
    submit_chunks : str
        "all", or a comma-separated list of 1-indexed chunk numbers to submit.
    config : dict
        Parsed --input-json config.
    outdir : str
        Directory chunk JSONs are written under (a run-timestamp subdir is added).
    workflow_id : str
        Workbench workflow ID to submit against.
    workflow_version_name : str
        Workbench workflow version to submit against.
    promote : bool
        If False, dry-run only (JSONs are still written, nothing is submitted).

    Returns
    -------
    None
    """
    all_ids = _read_sample_ids(infile_sample_ids_all)

    if validate_bucket_meta_tsv:
        if shutil.which("gcloud") is None:
            print("[ERROR] 'gcloud' CLI not found on PATH; required for --validate-bucket-meta-tsv.")
            sys.exit(1)
        done_ids = set(_find_done_sample_ids_from_bucket(_cluster_visium_gcs_dir(config)))
    elif infile_sample_ids_done:
        done_ids = set(_read_sample_ids(infile_sample_ids_done))
    else:
        done_ids = set()

    sample_ids = _resolve_remaining_sample_ids(all_ids, done_ids)

    chunks = _chunk(sample_ids, chunk_size)
    print(f"[INFO] Split into {len(chunks)} chunk(s) of up to {chunk_size} sample(s) each")

    try:
        selected_chunk_numbers = _parse_submit_chunks(submit_chunks, len(chunks))
    except ValueError as e:
        print(f"[ERROR] {e}")
        sys.exit(1)

    run_timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_dir = Path(outdir) / run_timestamp
    run_dir.mkdir(parents=True, exist_ok=True)

    for i, chunk_sample_ids in enumerate(chunks, start=1):
        chunk_json = run_dir / _CHUNK_FILENAME_TEMPLATE.format(index=i)
        params = _build_chunk_params(chunk_sample_ids, config)
        with open(chunk_json, "w") as f:
            json.dump(params, f, indent=2)
        print(f"[INFO] Chunk {i}/{len(chunks)} ({len(chunk_sample_ids)} sample(s)): {chunk_json}")
        print(f"         {chunk_sample_ids}")

        if i not in selected_chunk_numbers:
            print(f"[INFO] Chunk {i} not in --submit-chunks; skipping submission.")
            continue
        submit_chunk(chunk_json, workflow_id, workflow_version_name, promote)

    if not promote:
        print("[DRY-RUN] Re-run with --promote to actually submit.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--infile-sample-ids-all",
        required=True,
        help="Path to a plain-text file of sample IDs, one per line.",
    )
    done_group = parser.add_mutually_exclusive_group()
    done_group.add_argument(
        "--infile-sample-ids-done",
        default=None,
        help="Sample IDs already completed; subtracted from --infile-sample-ids-all.\n"
             "Mutually exclusive with --validate-bucket-meta-tsv.",
    )
    done_group.add_argument(
        "--validate-bucket-meta-tsv",
        action="store_true",
        default=False,
        help="Derive done samples from cluster_visium's bucket .meta.tsv sidecars (all 3\n"
             "must exist with well-formed content). Mutually exclusive with\n"
             "--infile-sample-ids-done.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        required=True,
        help="Maximum number of samples per chunk (one chunk = one JSON = one submitted run).\n"
             f"Must be <= {_MAX_CHUNK_SIZE} (GPU quota headroom; see _MAX_CHUNK_SIZE).",
    )
    parser.add_argument(
        "--submit-chunks",
        default=_DEFAULT_SUBMIT_CHUNKS,
        help="Which chunks to submit (subject to --promote): 'all', or e.g. '1' or '1,2'.\n"
             "Every chunk's JSON is always written regardless.\n"
             f"[Default: {_DEFAULT_SUBMIT_CHUNKS}]",
    )
    parser.add_argument(
        "-i", "--input-json",
        required=True,
        help=(
            "Path to a JSON config file with project/environment settings shared\n"
            "across all chunks, e.g.:\n"
            "\n"
            "{\n"
            '  "asap_team_id": "team-scherzer",\n'
            '  "asap_dataset_id": "PD5D_MTG_Visium10x_ST",\n'
            '  "asap_dataset_doi_url": "https://doi.org/10.5281/zenodo.17242087",\n'
            '  "raw_data_bucket": "gs://asap-raw-team-scherzer-pmdbs-spatial-visium-mtg",\n'
            '  "qc_adata_object_template": "gs://.../counts_to_adata/1.0.0/{sample_id}.cleaned_unfiltered.h5ad",\n'
            '  "sc_rnaseq_reference_cell2location_model": "gs://.../deconvolution_reference_strategy3_inf_aver.csv.gz",\n'
            '  "cohort_raw_data_bucket_spatial_sc_rnaseq": "gs://asap-raw-team-scherzer-pmdbs-spatial-visium-mtg",\n'
            '  "container_registry": "us-central1-docker.pkg.dev/dnastack-asap-parkinsons/workflow-images",\n'
            '  "workflow_id": "865c2967-9af8-499d-9b70-aa10db6c49da",\n'
            f'  "workflow_version_name": "v1.0.0",\n'
            f'  "batch_key": "{_DEFAULT_BATCH_KEY}",\n'
            f'  "cell_type_key": "{_DEFAULT_CELL_TYPE_KEY}",\n'
            f'  "outdir": "{_DEFAULT_OUTDIR}"\n'
            "}\n"
            "\n"
            "qc_adata_object_template must contain a literal '{sample_id}' placeholder.\n"
            "\n"
            "Fields:\n"
            "  asap_team_id (required)\n"
            "  asap_dataset_id (required)\n"
            "  asap_dataset_doi_url (required)\n"
            "  raw_data_bucket (required)\n"
            "      All 4 above describe the single SpatialVisiumProject every chunk shares.\n"
            "\n"
            "  qc_adata_object_template (required)\n"
            "      GCS path template for each sample's post-QC H5AD; must contain '{sample_id}'.\n"
            "\n"
            "  sc_rnaseq_reference_cell2location_model (required)\n"
            "      GCS path to the already-built Stage 1 inf_aver.csv.gz.\n"
            "\n"
            "  cohort_raw_data_bucket_spatial_sc_rnaseq (required)\n"
            "      Raw data bucket root; see workflows/main.wdl parameter_meta.\n"
            "\n"
            "  container_registry (required)\n"
            "\n"
            "  workflow_id (required)\n"
            "      Workbench workflow ID (from `omics workbench workflows create`).\n"
            "\n"
            "  workflow_version_name (required)\n"
            "      Workbench workflow version to submit against. For example, 'v1.0.0'\n"
            "      if used: `omics workbench workflows create ... --version-name v1.0.0`.\n"
            "\n"
            "  batch_key (optional)\n"
            f"      Default: '{_DEFAULT_BATCH_KEY}'\n"
            "\n"
            "  cell_type_key (optional)\n"
            f"      Default: '{_DEFAULT_CELL_TYPE_KEY}'\n"
            "\n"
            "  preemptible (optional)\n"
            "      Overrides the workflow's own preemptible default when set.\n"
            "\n"
            "  outdir (optional)\n"
            "      Directory chunk JSONs are written under (run-timestamp subdir added).\n"
            f"      Default: '{_DEFAULT_OUTDIR}'\n"
            "\n"
            "  cell2location_version (optional)\n"
            "      Must match deconvolution.wdl's cell2location_version; used by\n"
            "      --validate-bucket-meta-tsv to locate cluster_visium's output dir.\n"
            f"      Default: '{_DEFAULT_CELL2LOCATION_VERSION}'\n"
        ),
    )
    parser.add_argument(
        "--promote",
        action="store_true",
        default=False,
        help="Actually submit each chunk via the omics CLI. Default is dry-run "
             "(JSONs are written, nothing is submitted).",
    )
    args = parser.parse_args()

    if args.chunk_size > _MAX_CHUNK_SIZE:
        parser.error(f"--chunk-size must be <= {_MAX_CHUNK_SIZE}, got: {args.chunk_size}")

    try:
        with open(args.input_json, "r") as f:
            config = json.load(f)
    except FileNotFoundError:
        parser.error(f"{args.input_json} JSON file does not exist.")

    required_fields = [
        "asap_team_id", "asap_dataset_id", "asap_dataset_doi_url", "raw_data_bucket",
        "qc_adata_object_template", "sc_rnaseq_reference_cell2location_model",
        "cohort_raw_data_bucket_spatial_sc_rnaseq", "container_registry",
        "workflow_id", "workflow_version_name",
    ]
    missing = [field for field in required_fields if field not in config]
    if missing:
        parser.error(f"Missing required field(s) {missing} in {args.input_json}.")

    if "{sample_id}" not in config["qc_adata_object_template"]:
        parser.error("qc_adata_object_template must contain a literal '{sample_id}' placeholder.")

    outdir = _resolve_repo_path(config.get("outdir", _DEFAULT_OUTDIR))

    log_run_command(Path(outdir) / "submit_cell2location_2and3_prod.log")
    main(
        infile_sample_ids_all=args.infile_sample_ids_all,
        infile_sample_ids_done=args.infile_sample_ids_done,
        validate_bucket_meta_tsv=args.validate_bucket_meta_tsv,
        chunk_size=args.chunk_size,
        submit_chunks=args.submit_chunks,
        config=config,
        outdir=outdir,
        workflow_id=config["workflow_id"],
        workflow_version_name=config["workflow_version_name"],
        promote=args.promote,
    )
