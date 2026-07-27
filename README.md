# spatial-sc-rnaseq-integration-wf

Repo for testing and developing a workflow that integrates single-cell/single-nucleus RNA-seq data with spatial transcriptomics data (Nanostring GeoMx and 10x Visium).

This workflow consumes outputs from two upstream sibling workflows:
* [sc-rnaseq-wf](https://github.com/ASAP-CRN/sc-rnaseq-wf)
* [spatial-transcriptomics-wf](https://github.com/ASAP-CRN/spatial-transcriptomics-wf)

Common workflows, tasks, utility scripts, and docker images reused across harmonized ASAP workflows are defined in [the wf-common repository](https://github.com/ASAP-CRN/wf-common).


# Table of contents

- [Workflows](#workflows)
    - [Integration](#integration)
        - [Deconvolution methods](#deconvolution-methods)
- [Inputs](#inputs)
    - [Structs](#structs)
    - [Generating the inputs JSON](#generating-the-inputs-json)
- [Outputs](#outputs)
    - [Output structure](#output-structure)
- [Docker images](#docker-images)
- [wdl-ci](#wdl-ci)


# Workflows

Workflows are defined in [the `workflows` directory](workflows). The Python scripts which process the data at each stage can be found in [the docker/integration_tools/scripts directory](docker/integration_tools/scripts).

Workflow name:
- `spatial_sc_rnaseq_integration`

**Entrypoint**: [workflows/main.wdl](workflows/main.wdl)

**Input template**: [workflows/inputs.json](workflows/inputs.json)

## Integration

The integration workflow receives processed AnnData objects (`.h5ad`) from the upstream sibling workflows (i.e. does not re-run CellRanger, SpaceRanger, or cell type labelling).

The scRNA-seq cohort adata is used as a **reference** for all spatial samples. The spatial side is handled **per sample**: each spatial sample adata carrying raw UMI counts is used before any batch correction, so it can be harmonized independently against the scRNA-seq reference. Samples from all projects are flattened and integrated in parallel.

Upstream inputs consumed per source:

| Source workflow | Level | Output consumed |
| :- | :- | :- |
| sc-rnaseq-wf `cohort_analysis` | cohort (reference) | `final_adata_object`, `final_metadata_csv`, `scanvi_cell_types_parquet` |
| spatial-transcriptomics-wf `spatial_visium preprocess` | per-sample, raw UMI counts | `qc_adata_objects` (one per sample) |
| spatial-transcriptomics-wf `spatial_geomx process_to_adata` | per-sample | `processed_adata_objects` (one per sample) |

The integration sub-workflow ([workflows/integration/integration.wdl](workflows/integration/integration.wdl)) performs a single step per sample:

1. **Integrate** — run the integration method on the per-sample spatial adata against the scRNA-seq cohort reference (`integrate_sample` task).

All samples (across all projects) are scattered flat and integrated in parallel.

### Deconvolution methods

Cell type deconvolution is performed by running two complementary methods per sample, providing both a Python and an R solution:

| Method | Reference | Language | Docker image |
| :- | :- | :- | :- |
| [Cell2location](https://github.com/BayraktarLab/cell2location), Bayesian hierarchical model | Kleshchevnikov *et al.*, *Nat Biotechnol* 2022 | Python — GPU-accelerated | `cell2location_tools` |
| [RCTD](https://github.com/dmcable/spacexr), statistical mixture model | Cable *et al.*, *Nat Biotechnol* 2022 | R — CPU-only | `rctd_tools` |

Both methods use the scRNA-seq cohort as a reference to infer cell type composition per spot. Cell2location estimates absolute cell abundances per spot and provides posterior uncertainty for each prediction. RCTD estimates cell type proportions per spot and is fast and robust across tissue types.

Both methods are designed for spot-based spatial transcriptomics and are applied to Visium samples. GeoMx ROI data may also be processed where ROI size and design make spot-level deconvolution appropriate.

Cell2location is the primary method and requires a GPU. RCTD serves as a robust fallback: because Cell2location can fail on large datasets due to GPU memory constraints, running both methods in parallel ensures deconvolution outputs are always available. Since each method runs as an independent WDL task, GPU resources are provisioned only for Cell2location, keeping RCTD on lower-cost CPU instances.


# Inputs

An input template file can be found at [workflows/inputs.json](workflows/inputs.json).

| Type | Name | Description |
| :- | :- | :- |
| String | cohort_id | Name of the cohort; used to name output files. |
| Array[[ScRnaseqProject](#scrnaseqproject)] | sc_rnaseq_projects | Processed cohort outputs from sc-rnaseq-wf used as the integration reference. Typically one entry representing the cross-team cohort. |
| Array[[SpatialVisiumProject](#spatialvisiumproject)]? | spatial_visium_projects | Processed project-level outputs from spatial-transcriptomics-wf (Visium) to integrate per sample. [[ ]] |
| Array[[SpatialGeomxProject](#spatialgeomxproject)]? | spatial_geomx_projects | Processed project-level outputs from spatial-transcriptomics-wf (GeoMx) to integrate per sample. [[ ]] |
| String? | batch_key | Key in AnnData object for batch information. ['batch_id'] |
| Int? | n_top_genes | Number of highly variable genes to keep. [3000] |
| Int? | n_comps | Number of principal components to compute. [30] |
| Int? | n_neighbors | Size of local neighborhood used for manifold approximation. [15] |
| Array[Float]? | leiden_res | Leiden resolutions controlling the coarseness of clustering. [0.05, 0.1, 0.2, 0.4] |
| String | cohort_raw_data_bucket_spatial_sc_rnaseq | Bucket to upload cohort integration intermediate files to. |
| Array[String] | cohort_staging_data_buckets_spatial_sc_rnaseq | Buckets to upload cohort integration final outputs to. |
| String | container_registry | Container registry where workflow Docker images are hosted. |
| String? | zones | GCP zones where compute will take place. ['us-central1-c us-central1-f'] |

## Structs

### ScRnaseqProject

Wraps the final cohort-level outputs of [sc-rnaseq-wf](https://github.com/ASAP-CRN/sc-rnaseq-wf) `cohort_analysis`. Used as the scRNA-seq reference for all spatial samples.

| Type | Name | Description |
| :- | :- | :- |
| String | asap_team_id | ASAP-generated unique identifier for team; used for naming output files. |
| String | asap_dataset_id | ASAP-generated unique identifier for dataset; used for metadata. |
| String | asap_dataset_doi_url | ASAP-generated Zenodo DOI URL referencing the dataset. |
| File | final_adata_object | `cohort_analysis.final_adata_object` (.h5ad) — Harmony-integrated, cell-type-annotated AnnData. |
| File | final_metadata_csv | `cohort_analysis.final_metadata_csv` — per-cell metadata table. |
| File | scanvi_cell_types_parquet | `cohort_analysis.scanvi_cell_types_parquet` — scANVI cell type predictions. |
| String | raw_data_bucket | Raw data bucket; intermediate output files are stored here. |
| Array[String] | staging_data_buckets | Staging data buckets; final outputs are stored here. |

### SpatialVisiumProject

Groups per-sample Visium adatas from [spatial-transcriptomics-wf](https://github.com/ASAP-CRN/spatial-transcriptomics-wf) `spatial_visium preprocess` for a single team/dataset. Samples are flattened across projects and integrated independently.

| Type | Name | Description |
| :- | :- | :- |
| String | asap_team_id | ASAP-generated unique identifier for team; used for naming output files. |
| String | asap_dataset_id | ASAP-generated unique identifier for dataset; used for metadata. |
| String | asap_dataset_doi_url | ASAP-generated Zenodo DOI URL referencing the dataset. |
| Array[[SpatialVisiumSample](#spatialvisiumsample)] | samples | Per-sample adatas (raw UMI counts, pre-Harmony); each sample is integrated independently. |
| String | raw_data_bucket | Raw data bucket; intermediate output files are stored here. |
| Array[String] | staging_data_buckets | Staging data buckets; final outputs are stored here. |

### SpatialVisiumSample

One Visium sample. Source: `spatial_visium preprocess.qc_adata_objects` from [spatial-transcriptomics-wf](https://github.com/ASAP-CRN/spatial-transcriptomics-wf). The adata carries post-QC raw UMI counts before Harmony batch correction.

| Type | Name | Description |
| :- | :- | :- |
| String | sample_id | Sample identifier; used to name output files. |
| File | qc_adata_object | `spatial_visium preprocess.qc_adata_objects` (.h5ad) — per-sample, post-QC, raw UMI counts. |

### SpatialGeomxProject

Groups per-sample GeoMx adatas from [spatial-transcriptomics-wf](https://github.com/ASAP-CRN/spatial-transcriptomics-wf) `spatial_geomx process_to_adata` for a single team/dataset. Samples are flattened across projects and integrated independently.

| Type | Name | Description |
| :- | :- | :- |
| String | asap_team_id | ASAP-generated unique identifier for team; used for naming output files. |
| String | asap_dataset_id | ASAP-generated unique identifier for dataset; used for metadata. |
| String | asap_dataset_doi_url | ASAP-generated Zenodo DOI URL referencing the dataset. |
| Array[[SpatialGeomxSample](#spatialgeomxsample)] | samples | Per-sample adatas; each sample is integrated independently. |
| String | raw_data_bucket | Raw data bucket; intermediate output files are stored here. |
| Array[String] | staging_data_buckets | Staging data buckets; final outputs are stored here. |

### SpatialGeomxSample

One GeoMx sample. Source: `spatial_geomx process_to_adata.processed_adata_objects` from [spatial-transcriptomics-wf](https://github.com/ASAP-CRN/spatial-transcriptomics-wf). The adata carries Q3-normalized counts (75th-percentile within-sample normalization), before Harmony cross-sample batch correction, which runs later in `spatial_geomx cohort_analysis`.

| Type | Name | Description |
| :- | :- | :- |
| String | sample_id | Sample identifier; used to name output files. |
| File | processed_adata_object | `spatial_geomx process_to_adata.processed_adata_objects` (.h5ad) — per-sample, Q3-normalized, pre-Harmony. |

## Generating the inputs JSON

The inputs JSON may be generated manually. Because this workflow consumes processed outputs from upstream sibling workflows rather than raw FASTQ files, the `sc_rnaseq_projects`, `spatial_visium_projects`, and `spatial_geomx_projects` arrays must reference GCS paths to the `.h5ad` and `.parquet` files already produced by those workflows.

The relevant GCS paths follow the staging bucket structure of the sibling workflows:

```
# sc-rnaseq-wf outputs (cross-team cohort)
asap-dev-{cohort}-pmdbs-sn-rnaseq-{context}/pmdbs_sc_rnaseq/release/${crn_release_version}/cohort_analysis/
  ${cohort_id}.final.h5ad
  ${cohort_id}.final_metadata.csv
  ${cohort_id}.scanvi_cell_types.parquet

# spatial-transcriptomics-wf Visium outputs — per-sample pre-QC adatas (raw UMI counts)
asap-dev-{team-xxyy}-pmdbs-spatial-visium-{context}/spatial_visium/release/${crn_release_version}/preprocess/
  ${sample_id}.qc.h5ad   # one file per sample

# spatial-transcriptomics-wf GeoMx outputs — per-sample processed adatas
asap-dev-{team-xxyy}-pmdbs-spatial-geomx-{context}/spatial_geomx/release/${crn_release_version}/process_to_adata/
  ${sample_id}.processed.h5ad   # one file per sample
```

Example inputs JSON fragment:

```json
{
  "spatial_sc_rnaseq_integration.cohort_id": "asap-cohort",
  "spatial_sc_rnaseq_integration.sc_rnaseq_projects": [
    {
      "asap_team_id": "cohort",
      "asap_dataset_id": "DS_PMDBS_XXXX",
      "asap_dataset_doi_url": "https://doi.org/...",
      "final_adata_object": "gs://asap-dev-cohort-.../cohort_analysis/asap-cohort.final.h5ad",
      "final_metadata_csv": "gs://asap-dev-cohort-.../cohort_analysis/asap-cohort.final_metadata.csv",
      "scanvi_cell_types_parquet": "gs://asap-dev-cohort-.../cohort_analysis/asap-cohort.scanvi_cell_types.parquet",
      "raw_data_bucket": "gs://asap-raw-cohort-...",
      "staging_data_buckets": ["gs://asap-dev-cohort-..."]
    }
  ],
  "spatial_sc_rnaseq_integration.spatial_visium_projects": [
    {
      "asap_team_id": "team-xxyy",
      "asap_dataset_id": "DS_PMDBS_YYYY",
      "asap_dataset_doi_url": "https://doi.org/...",
      "samples": [
        {
          "sample_id": "ASAP_PMDBS_000001_s001_rep1",
          "qc_adata_object": "gs://asap-dev-team-xxyy-.../preprocess/ASAP_PMDBS_000001_s001_rep1.qc.h5ad"
        },
        {
          "sample_id": "ASAP_PMDBS_000002_s001_rep1",
          "qc_adata_object": "gs://asap-dev-team-xxyy-.../preprocess/ASAP_PMDBS_000002_s001_rep1.qc.h5ad"
        }
      ],
      "raw_data_bucket": "gs://asap-raw-team-xxyy-...",
      "staging_data_buckets": ["gs://asap-dev-team-xxyy-..."]
    }
  ],
  "spatial_sc_rnaseq_integration.spatial_geomx_projects": [],
  "spatial_sc_rnaseq_integration.cohort_raw_data_bucket_spatial_sc_rnaseq": "gs://asap-raw-cohort-...",
  "spatial_sc_rnaseq_integration.cohort_staging_data_buckets_spatial_sc_rnaseq": ["gs://asap-dev-cohort-..."],
  "spatial_sc_rnaseq_integration.container_registry": "us-central1-docker.pkg.dev/dnastack-asap-parkinsons/workflow-images"
}
```

# Outputs

## Output structure

- `workflow_run_timestamp`: format: `%Y-%m-%dT%H-%M-%SZ`
- Outputs are written per spatial sample; each file is named with the `sample_id`
- The `MANIFEST.tsv` file in the staging data bucket describes the file name, md5 hash, timestamp, workflow version, workflow name, and workflow release for the run used to generate each file in that directory

### Raw data (intermediate files and final outputs for all runs of the workflow)

The raw data bucket will contain *some* artifacts generated as part of workflow execution. Following successful workflow execution, the artifacts will also be copied into the staging bucket as final outputs.

In the workflow, task outputs are either specified as `String` (final outputs, which will be copied in order to live in raw data buckets and staging buckets) or `File` (intermediate outputs that are periodically cleaned up, which will live in the cromwell-output bucket). This was implemented to reduce storage costs.

```bash
asap-raw-{cohort}-{source}-{modality_flavour}-{context}
└── workflow_execution
    └── spatial_sc_rnaseq_integration
        └── cohort_analysis
            └── ${cohort_analysis_workflow_version}
                └── ${workflow_run_timestamp}
                    └── <per-sample integration outputs>
```

### Staging data (intermediate workflow objects and final workflow outputs for the latest run of the workflow)

Following QC by researchers, the objects in the dev or uat bucket are synced into the curated data buckets, maintaining the same file structure. Curated data buckets are named `asap-curated-{cohort}-{source}-{modality_flavour}-{context}` and `dataset_id` = `{cohort}-{source}-{modality_flavour}-{context}`.

Data may be synced using [the `promote_staging_data` script](#promoting-staging-data).

```bash
asap-dev-{cohort}-{source}-{modality_flavour}-{context}
└── spatial_sc_rnaseq_integration
    └── release
        └── ${crn_release_version}
            ├── cohort_analysis
            │   ├── ${sampleA_id}.integrated.h5ad
            │   ├── ${sampleA_id}.integration_metadata.csv
            │   ├── ${sampleB_id}.integrated.h5ad
            │   ├── ${sampleB_id}.integration_metadata.csv
            │   ├── ...
            │   ├── ${sampleN_id}.integrated.h5ad
            │   ├── ${sampleN_id}.integration_metadata.csv
            │   └── MANIFEST.tsv
            ├── workflow_version  # plain text file
            └── workflow_metadata
                └── ${timestamp}
                    ├── MANIFEST.tsv  # combined
                    └── data_promotion_report.md
```

## Promoting staging data

The [`promote_staging_data` script](https://github.com/ASAP-CRN/wf-common/blob/main/util/promote_staging_data) can be used to promote staging data that has been approved to the curated data bucket.

This script compiles bucket and file information for both the initial (staging) and target (prod) environment. It also runs data integrity tests to ensure staging data can be promoted and generates a Markdown report. It (1) checks that files are not empty and are not less than or equal to 10 bytes (factoring in white space) and (2) checks that files have associated metadata and is present in MANIFEST.tsv.

If data integrity tests pass, this script will upload a combined MANIFEST.tsv and the data promotion Markdown report under a metadata/{timestamp} directory in the staging bucket. Previous manifest files and reports will be kept. Next, it will rsync all files in the staging bucket to the curated bucket's workflow and metadata directories. **Exercise caution when using this script**; files that are not present in the source (staging) bucket will be deleted at the destination (curated) bucket.

If data integrity tests fail, staging data cannot be promoted. The combined `MANIFEST.tsv`, Markdown report, and `promote_staging_data_script.log` will be locally available.

The script defaults to a dry run, printing out the files that would be copied or deleted.

### Options

```
-h  Display this message and exit
-l  List available teams
-w  Workflow name used as a directory in bucket (e.g. 'spatial_sc_rnaseq_integration')
-v  Release version (e.g. v1.0.0)
-p  Promote data. If this option is not selected, data that would be copied or deleted is printed out, but files are not actually changed (dry run)
```

### Usage

```bash
# List available teams
./wf-common/util/promote_staging_data -l -w spatial_sc_rnaseq_integration -v v1.0.0

# Print out the files that would be copied or deleted from the staging bucket to the curated bucket for a specific release version
./wf-common/util/promote_staging_data -w spatial_sc_rnaseq_integration -v v1.0.0

# Promote data for a specific release version
./wf-common/util/promote_staging_data -w spatial_sc_rnaseq_integration -v v1.0.0 -p
```


# Docker images

Docker images are defined in [the `docker` directory](docker). Each image must minimally define a `build.env` file and a `Dockerfile`.

Example directory structure:
```bash
docker
└── integration_tools
    ├── build.env
    ├── Dockerfile
    ├── requirements.txt
    └── scripts
        └── ...
```

## The `build.env` file

Each target image is defined using the `build.env` file, which is used to specify the name and version tag for the corresponding Docker image. It must contain at minimum the following variables:

- `IMAGE_NAME`
- `IMAGE_TAG`

All variables defined in the `build.env` file will be made available as build arguments during Docker image build.

The `DOCKERFILE` variable may be used to specify the path to a Dockerfile if that file is not found alongside the `build.env` file, for example when multiple images use the same base Dockerfile definition.

## Building Docker images

Docker images can be built using the [`build_docker_images`](https://github.com/DNAstack/bioinformatics-scripts/blob/main/scripts/build_docker_images) script.

```bash
# Build a single image (use -x to target linux/amd64)
./build_docker_images -d docker/integration_tools -x

# Build all images in the `docker` directory
./build_docker_images -d docker -x

# Build and push all images in the docker directory, using the `dnastack` container registry
./build_docker_images -d docker -x -c dnastack -p
```

## Tool and library versions

| Image | Major tool versions | Links |
| :- | :- | :- |
| integration_tools | <ul><li>[google-cloud-cli 524.0.0](https://cloud.google.com/sdk/docs/release-notes#52400_2025-05-28)</li><li>[python 3.10.12](https://www.python.org/downloads/release/python-31012/)</li><li>[anndata](https://anndata.readthedocs.io/en/latest/)</li><li>[scanpy](https://scanpy.readthedocs.io/en/stable/)</li></ul> | [Dockerfile](docker/integration_tools) |
| cell2location_tools | <ul><li>[cell2location](https://github.com/BayraktarLab/cell2location)</li><li>[PyTorch](https://pytorch.org/)</li></ul> | [Dockerfile](docker/cell2location_tools) |
| rctd_tools | <ul><li>[spacexr](https://github.com/dmcable/spacexr)</li></ul> | [Dockerfile](docker/rctd_tools) |
| util | <ul><li>[google-cloud-cli 524.0.0](https://cloud.google.com/sdk/docs/release-notes#52400_2025-05-28)</li></ul> | [Dockerfile](https://github.com/ASAP-CRN/wf-common/tree/main/docker/util) |


# wdl-ci

[`wdl-ci`](https://github.com/DNAstack/wdl-ci) provides tools to validate and test workflows and tasks written in [Workflow Description Language (WDL)](https://github.com/openwdl/wdl). In addition to the tests packaged in `wdl-ci`, the [spatial-sc-rnaseq-integration-wdl-ci-custom-test-dir](./spatial-sc-rnaseq-integration-wdl-ci-custom-test-dir) is a directory containing custom WDL-based tests that are used to test workflow tasks. `wdl-ci` in this repository is set up to run on pull request.

In general, `wdl-ci` will use inputs provided in the [wdl-ci.config.json](./wdl-ci.config.json) and compare current outputs and validated outputs based on changed tasks/workflows to ensure outputs are still valid by meeting the criteria in the specified tests. For example, if the `integrate_sample` task was changed, then this task would be submitted and that output would be considered the "current output". When inspecting the integrated AnnData object, there is a test specified in the [wdl-ci.config.json](./wdl-ci.config.json) called `check_hdf5`. The test will compare the "current output" and "validated output" (provided in the [wdl-ci.config.json](./wdl-ci.config.json)) to make sure that the `.h5ad` file is still a valid HDF5 file.
