#!/usr/bin/env python3
"""
Submit Cell2location Stage 1 (reference build) via the omics Workbench CLI.

Before submitting, checks whether <ref_stem>_inf_aver.csv.gz already exists
in the raw bucket and prevents overwriting it unless --force is given.

--promote gates whether submission actually happens (default is dry-run).

Outputs
-------
<outdir>/<run_timestamp>/stage1.json   -- the submitted workflow-params JSON
<outdir>/submit_cell2location_1_prod.log
"""

import sys
import argparse
import json
import shutil
import subprocess
from pathlib import Path
from datetime import datetime, timezone

repo_root = next(p for p in Path(__file__).resolve().parents if (p / ".git").exists())
sys.path.insert(0, str(repo_root))

from utils.logging_extra import log_run_command

_DEFAULT_OUTDIR = "prod_runs"
_DEFAULT_CELL2LOCATION_VERSION = "1.0.0"
_WORKFLOW_PARAM_PREFIX = "spatial_sc_rnaseq_integration."


def _resolve_repo_path(path: str) -> str:
    """
    Resolve `path` against the repo root if relative, else return it unchanged.

    Parameters
    ----------
    path : str
        Path from a CLI arg; relative (to the repo root) or absolute.

    Returns
    -------
    str
        Absolute path.
    """
    p = Path(path)
    return str(p if p.is_absolute() else (repo_root / p).resolve())


def _build_reference_gcs_path(config: dict, cell2location_version: str) -> str:
    """
    Build the expected GCS path for Stage 1's inf_aver.csv.gz output, matching
    deconvolution.wdl's own raw_data_path construction.

    Parameters
    ----------
    config : dict
        Parsed --input-json workflow-params (spatial_sc_rnaseq_integration.-prefixed).
    cell2location_version : str
        Must match deconvolution.wdl's own cell2location_version constant.

    Returns
    -------
    str
        GCS path to the expected inf_aver.csv.gz file.
    """
    bucket = config[f"{_WORKFLOW_PARAM_PREFIX}cohort_raw_data_bucket_spatial_sc_rnaseq"].rstrip("/")
    ref_adata = config[f"{_WORKFLOW_PARAM_PREFIX}sc_rnaseq_reference_adata"]
    ref_stem = Path(ref_adata).name.removesuffix(".h5ad")
    return (
        f"{bucket}/workflow_execution/spatial_visium/deconvolution_analysis/"
        f"build_reference/{cell2location_version}/{ref_stem}_inf_aver.csv.gz"
    )


def _gcs_object_exists(gcs_path: str) -> bool:
    """
    Check whether a single GCS object exists via `gcloud storage ls`.

    Parameters
    ----------
    gcs_path : str
        Full GCS object path.

    Returns
    -------
    bool
        True if `gcloud storage ls` finds the object.
    """
    result = subprocess.run(["gcloud", "storage", "ls", gcs_path], capture_output=True, text=True)
    return result.returncode == 0 and gcs_path in result.stdout


def submit_run(params_json: Path, workflow_id: str, workflow_version_name: str, promote: bool) -> None:
    """
    Submit a workflow-params JSON via the omics Workbench CLI.

    Parameters
    ----------
    params_json : pathlib.Path
        Path to the workflow-params JSON (already written to disk).
    workflow_id : str
        Workbench workflow ID (e.g. from `omics workbench workflows create`).
    workflow_version_name : str
        Workbench workflow registration version (e.g. "v1.0.1").
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
        "--workflow-params", f"@{params_json}",
    ]
    if not promote:
        print(f"[DRY-RUN] Would run: {' '.join(cmd)}")
        return

    if shutil.which("omics") is None:
        raise RuntimeError("'omics' CLI not found on PATH. Activate the venv with the Workbench CLI installed.")

    print(f"[INFO] Submitting {params_json} ...")
    subprocess.run(cmd, check=True)


def main(
    config: dict,
    outdir: str,
    workflow_id: str,
    workflow_version_name: str,
    cell2location_version: str,
    force: bool,
    promote: bool,
) -> None:
    """
    Check for an existing reference, then write and submit Stage 1's
    workflow-params JSON.

    Parameters
    ----------
    config : dict
        Parsed --input-json workflow-params (spatial_sc_rnaseq_integration.-prefixed).
    outdir : str
        Directory the submitted JSON is written under (a run-timestamp subdir is added).
    workflow_id : str
        Workbench workflow ID to submit against.
    workflow_version_name : str
        Workbench workflow version to submit against.
    cell2location_version : str
        Must match deconvolution.wdl's cell2location_version; used to locate
        build_reference's expected output path.
    force : bool
        Submit even if the reference already exists in the raw bucket.
    promote : bool
        If False, dry-run only (JSON is still written, nothing is submitted).

    Returns
    -------
    None
    """
    if shutil.which("gcloud") is None:
        print("[ERROR] 'gcloud' CLI not found on PATH; required to check for an existing reference.")
        sys.exit(1)

    expected_gcs_path = _build_reference_gcs_path(config, cell2location_version)
    print(f"[INFO] Checking for an existing reference at {expected_gcs_path} ...")
    if _gcs_object_exists(expected_gcs_path):
        if not force:
            print(f"[ERROR] {expected_gcs_path} already exists. Stage 1 training is expensive -- "
                  "re-run with --force to rebuild it anyway.")
            sys.exit(1)
        print(f"[WARNING] {expected_gcs_path} already exists; --force given, proceeding to rebuild.")
    else:
        print("[INFO] No existing reference found; proceeding.")

    run_timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_dir = Path(outdir) / run_timestamp
    run_dir.mkdir(parents=True, exist_ok=True)
    params_json = run_dir / "stage1.json"
    with open(params_json, "w") as f:
        json.dump(config, f, indent=2)
    print(f"[INFO] Wrote {params_json}")

    submit_run(params_json, workflow_id, workflow_version_name, promote)

    if not promote:
        print("[DRY-RUN] Re-run with --promote to actually submit.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input-json",
        required=True,
        help="Path to a complete spatial_sc_rnaseq_integration workflow-params JSON for\n"
             "Stage 1 (fully-prefixed keys), e.g. prod_runs/submit_cell2location_1_prod_config.json.\n"
             f"Must set {_WORKFLOW_PARAM_PREFIX}run_stages to \"1\".",
    )
    parser.add_argument(
        "--workflow-id",
        required=True,
        help="Workbench workflow ID (from `omics workbench workflows create`).",
    )
    parser.add_argument(
        "--workflow-version-name",
        required=True,
        help="Workbench workflow version to submit against (e.g. 'v1.0.1').",
    )
    parser.add_argument(
        "--outdir",
        default=_DEFAULT_OUTDIR,
        help="Directory the submitted JSON is written under (run-timestamp subdir added).\n"
             f"[Default: {_DEFAULT_OUTDIR}]",
    )
    parser.add_argument(
        "--cell2location-version",
        default=_DEFAULT_CELL2LOCATION_VERSION,
        help="Must match deconvolution.wdl's cell2location_version; used to locate\n"
             f"build_reference's expected output path. [Default: {_DEFAULT_CELL2LOCATION_VERSION}]",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Submit even if the reference already exists in the raw bucket.",
    )
    parser.add_argument(
        "--promote",
        action="store_true",
        default=False,
        help="Actually submit via the omics CLI. Default is dry-run\n"
             "(JSON is written, nothing is submitted).",
    )
    args = parser.parse_args()

    try:
        with open(args.input_json, "r") as f:
            config = json.load(f)
    except FileNotFoundError:
        parser.error(f"{args.input_json} JSON file does not exist.")

    required_fields = [
        f"{_WORKFLOW_PARAM_PREFIX}sc_rnaseq_reference_adata",
        f"{_WORKFLOW_PARAM_PREFIX}cohort_raw_data_bucket_spatial_sc_rnaseq",
        f"{_WORKFLOW_PARAM_PREFIX}container_registry",
    ]
    missing = [field for field in required_fields if field not in config]
    if missing:
        parser.error(f"Missing required field(s) {missing} in {args.input_json}.")

    run_stages = config.get(f"{_WORKFLOW_PARAM_PREFIX}run_stages")
    if run_stages != "1":
        parser.error(
            f"{args.input_json} must set {_WORKFLOW_PARAM_PREFIX}run_stages to \"1\" "
            f"(got: {run_stages!r}) -- this script is Stage 1 only."
        )

    outdir = _resolve_repo_path(args.outdir)
    log_run_command(Path(outdir) / "submit_cell2location_1_prod.log")
    main(
        config=config,
        outdir=outdir,
        workflow_id=args.workflow_id,
        workflow_version_name=args.workflow_version_name,
        cell2location_version=args.cell2location_version,
        force=args.force,
        promote=args.promote,
    )
