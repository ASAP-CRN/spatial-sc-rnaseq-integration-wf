# spatial-sc-rnaseq-integration-wf

Repo for testing and developing a workflow that integrates single-cell/single-nucleus RNA-seq data with spatial transcriptomics data. This version (v1.0.0) covers deconvolution of 10x Visium data only.

This workflow consumes outputs from two upstream ASAP workflows:
* [sc-rnaseq-wf](https://github.com/ASAP-CRN/sc-rnaseq-wf)
* [spatial-transcriptomics-wf](https://github.com/ASAP-CRN/spatial-transcriptomics-wf)

Common workflows, tasks, utility scripts, and docker images reused across harmonized ASAP workflows are defined in [the wf-common repository](https://github.com/ASAP-CRN/wf-common).


# Table of contents

- [Workflows](#workflows)
    - [Spatial deconvolution](#spatial-deconvolution)
        - [Deconvolution methods](#deconvolution-methods)
- [Inputs](#inputs)
    - [Structures](#structures)
    - [Generating the inputs JSON](#generating-the-inputs-json)
- [Outputs](#outputs)
    - [Output structure](#output-structure)
- [Docker images](#docker-images)
- [wdl-ci](#wdl-ci)


# Workflows

Workflows are defined in [the `workflows` directory](./workflows/). The Python scripts which process the data at each stage can be found in the `docker/<method>_tools/scripts` directories.

Workflow name:
- `spatial_sc_rnaseq_integration`

**Entrypoint**: [workflows/main.wdl](workflows/main.wdl)

**Input template**: [workflows/inputs.json](workflows/inputs.json)

## Spatial deconvolution

The deconvolution workflow receives processed AnnData objects (`.h5ad`) from the upstream sibling workflows (i.e. does not re-run CellRanger, SpaceRanger, or cell type labelling).

The scRNA-seq cohort adata is used as a **reference** for all spatial samples. The spatial side is handled **per sample**: each spatial sample adata carrying raw UMI counts and all genes--not only HVGs is used before any batch correction, so it can be harmonized independently against the scRNA-seq reference. Samples from all projects are flattened and deconvoluted in parallel.

Upstream inputs consumed per source:

| Source workflow | Level | Output consumed |
| :- | :- | :- |
| sc-rnaseq-wf `cohort_analysis` (downsampled and/or supplemented with other sources) | cohort (reference) | `final_adata_object` |
| spatial-transcriptomics-wf `spatial_visium preprocess` | per-sample, raw UMI counts, all genes--not only HVGs| `qc_adata_objects` (one per sample) |

The deconvolution sub-workflow ([workflows/integration/integration.wdl](workflows/integration/integration.wdl)) runs as 3 independent WDL tasks, gated by `run_stages`:

1. **Stage 1 — NB regression** (skipped if `sc_rnaseq_reference_cell2location_model` is provided): train a Negative Binomial regression model on the scRNA-seq reference to estimate per-cell-type gene expression signatures (`inf_aver`). Runs **once for all spatial samples**.
2. **Stage 2 — Spatial mapping**: run Cell2location to map cell type abundances onto Visium spots. Runs once **per Visium sample**, scattered in parallel.
3. **Stage 3 — Clustering**: PCA → neighbours → UMAP → Leiden clustering on each spot's **deconvolved cell type fractions** (Stage 2's output) — not on gene expression. Runs once **per Visium sample**, scattered in parallel.

All samples (across all projects) are scattered flat and deconvoluted in parallel.

### Deconvolution methods

| Method | Reference | Language | Docker image | Status |
| :- | :- | :- | :- | :- |
| [Cell2location](https://github.com/BayraktarLab/cell2location), Bayesian hierarchical model | [Kleshchevnikov *et al.*, *Nat Biotechnol* 2022](https://pubmed.ncbi.nlm.nih.gov/35027729/) | Python — GPU-accelerated | `cell2location_tools` | Being implemented |
| [RCTD](https://github.com/dmcable/spacexr), statistical mixture model | [Cable *et al.*, *Nat Biotechnol* 2022](https://pubmed.ncbi.nlm.nih.gov/33603203/) | R — CPU-only | `rctd_tools` | To be implemented |

Both methods use the scRNA-seq cohort as a reference to infer cell type composition per spot. Cell2location estimates absolute cell abundances per spot and provides posterior uncertainty for each prediction. RCTD estimates cell type proportions per spot and is fast and robust across tissue types. These methods were reported as top-performers in benchmarks for the task ([Li *et al.*, *Nat Methods* 2022](https://pubmed.ncbi.nlm.nih.gov/35577954/); [Gillespie *et al.*, *Cells* 2025](https://pubmed.ncbi.nlm.nih.gov/40710313/)).

Cell2location is the primary method and requires a GPU. RCTD will serve as a robust fallback when implemented: because Cell2location can fail on large datasets due to GPU memory constraints, running both methods in parallel ensures deconvolution outputs are always available. Since each method runs as an independent WDL task, GPU resources are provisioned only for Cell2location, keeping RCTD on lower-cost CPU instances.


# Inputs

An input template file can be found at [workflows/inputs.json](workflows/inputs.json).

| Type | Name | Description |
| :- | :- | :- |
| File? | sc_rnaseq_reference_adata | scRNA-seq reference AnnData; used to deconvolute all spatial samples. Typically a downsampled version of sc-rnaseq-wf `cohort_analysis`'s `final_adata_object`, either alone or supplemented with other sources, such as Kamath/Siletti. Required when `run_stages` includes Stage 1. |
| Array[[`SpatialVisiumProject`](#spatialvisiumproject)]? | spatial_visium_projects | Processed project-level outputs from spatial-transcriptomics-wf (Visium) to deconvolute per sample. [[ ]] |
| File? | sc_rnaseq_reference_cell2location_model | Pre-computed `<ref_stem>_inf_aver.csv.gz`. Required when `run_stages` excludes Stage 1. |
| String? | batch_key | Column in reference .obs used as batch covariate in the NB regression model. ['batch_id'] |
| String? | cell_type_key | Column in reference .obs to use as the cell type label. ['cell_type_pred'] |
| Int? | n_cells_per_spot | Cell2location N_cells_per_location prior; 8 for brain Visium. [8] |
| Int? | detection_alpha | Cell2location detection rate prior (20–200); higher values assume low technical variability. [100] |
| Int? | max_epochs_ref | NB regression (Stage 1) training epochs. [250] |
| Int? | max_epochs_spatial | Cell2location spatial mapping (Stage 2) training epochs. [30000] |
| Int? | n_comps | PCA components for clustering on deconvolved cell type fractions (Stage 3). [30] |
| Int? | n_neighbors | Neighbourhood size for clustering on deconvolved cell type fractions (Stage 3). [15] |
| Array[Float]? | leiden_res | Leiden resolutions for clustering on deconvolved cell type fractions (Stage 3); one clustering per resolution. [0.2, 0.5, 1.0] |
| String? | run_stages | Comma-separated subset of `{1,2,3}` — which stages to run. `"1"` builds the `inf_aver.csv.gz` reference; `"2"` loads it and deconvolutes spatial samples; `"3"` Leiden-clusters each sample's deconvolved cell type fractions. ['1,2,3'] |
| String | cohort_raw_data_bucket_spatial_sc_rnaseq | Bucket to upload cohort deconvolution intermediate files to. |
| String | container_registry | Container registry where workflow Docker images are hosted. |
| String? | zones | GCP zones where compute will take place. ['us-central1-c us-central1-f'] |

#### Note
For script testing, these inputs can be used:
- sc-RNAseq reference: `gs://asap-wdl-ci/workflow_execution/spatial-sc-rnaseq-integration-wf/reference_building/reference_h5ad/deconvolution_reference_30k.h5ad`
- Spatial: `gs://asap-wdl-ci/workflow_execution/spatial_visium/preprocess/counts_to_adata/1.0.0/ASAP_PMBDS_000003_s007_Rep1.cleaned_unfiltered.h5ad`

## Structures

### `SpatialVisiumProject`

Groups per-sample Visium adatas from [spatial-transcriptomics-wf](https://github.com/ASAP-CRN/spatial-transcriptomics-wf) `spatial_visium preprocess` for a single team/dataset. Samples are flattened across projects and deconvoluted independently.

| Type | Name | Description |
| :- | :- | :- |
| String | asap_team_id | ASAP-generated unique identifier for team; used for naming output files. |
| String | asap_dataset_id | ASAP-generated unique identifier for dataset; used for metadata. |
| String | asap_dataset_doi_url | ASAP-generated Zenodo DOI URL referencing the dataset. |
| Array[[`SpatialVisiumSample`](#spatialvisiumsample)] | samples | Per-sample adatas; each sample is deconvoluted independently. |
| String | raw_data_bucket | Raw data bucket; intermediate output files are stored here. |

### `SpatialVisiumSample`

One Visium sample, with raw UMI counts and all genes--not only HVGs from [spatial-transcriptomics-wf](https://github.com/ASAP-CRN/spatial-transcriptomics-wf) (`*cleaned_unfiltered.h5ad`). Source: `spatial_visium preprocess.qc_adata_objects`. These can be obtained from the dataset buckets: `gs://asap-raw-<dataset_name>/workflow_execution/spatial_visium/preprocess/counts_to_adata/<version>`.

| Type | Name | Description |
| :- | :- | :- |
| String | sample_id | Sample identifier; used to name output files. |
| File | qc_adata_object | `spatial_visium preprocess.qc_adata_objects` (`*cleaned_unfiltered.h5ad`) — per-sample, raw UMI counts, all genes--not only HVGs. |

## Generating the inputs JSON

The inputs JSON may be generated manually. Because this workflow consumes processed outputs from upstream sibling workflows rather than raw FASTQ files, `sc_rnaseq_reference_adata` and the `spatial_visium_projects` array must reference GCS paths to the `.h5ad` files already produced by those workflows.

The relevant GCS paths follow the staging bucket structure of the sibling workflows:

```
# sc-rnaseq-wf outputs (downsampled cohort H5AD) -- `sc_rnaseq_reference_adata` may instead
# point to a custom-built reference (e.g. supplemented with other sources, such as Kamath/Siletti).
asap-dev-{cohort}-pmdbs-sn-rnaseq-{context}/pmdbs_sc_rnaseq/release/${crn_release_version}/cohort_analysis/
  ${cohort_id}.final.h5ad

# spatial-transcriptomics-wf Visium outputs from spatial-transcriptomics-wf (*cleaned_unfiltered.h5ad, with raw UMIs and all genes--not only HVGs).
asap-raw-{dataset_name}/workflow_execution/spatial_visium/preprocess/counts_to_adata/${version}/${sample_id}.cleaned_unfiltered.h5ad   # one file per sample
```

Example inputs JSON fragment:

```json
{
  "spatial_sc_rnaseq_integration.sc_rnaseq_reference_adata": "gs://asap-raw-cohort-.../cohort_analysis/asap-cohort.final.h5ad",
  "spatial_sc_rnaseq_integration.spatial_visium_projects": [
    {
      "asap_team_id": "team-xxyy",
      "asap_dataset_id": "DS_PMDBS_YYYY",
      "asap_dataset_doi_url": "https://doi.org/...",
      "samples": [
        {
          "sample_id": "ASAP_PMDBS_000001_s001_rep1",
          "qc_adata_object": "gs://asap-raw-team-xxyy-.../workflow_execution/spatial_visium/preprocess/counts_to_adata/1.0.0/ASAP_PMDBS_000001_s001_rep1.cleaned_unfiltered.h5ad"
        },
        {
          "sample_id": "ASAP_PMDBS_000002_s001_rep1",
          "qc_adata_object": "gs://asap-raw-team-xxyy-.../workflow_execution/spatial_visium/preprocess/counts_to_adata/1.0.0/ASAP_PMDBS_000002_s001_rep1.cleaned_unfiltered.h5ad"
        }
      ],
      "raw_data_bucket": "gs://asap-raw-team-xxyy-..."
    }
  ],
  "spatial_sc_rnaseq_integration.cohort_raw_data_bucket_spatial_sc_rnaseq": "gs://asap-raw-cohort-...",
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
                    └── <per-sample deconvolution outputs>
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
            │   ├── ${sampleA_id}.deconvoluted.h5ad
            │   ├── ${sampleA_id}.deconvolution_metadata.csv
            │   ├── ${sampleB_id}.deconvoluted.h5ad
            │   ├── ${sampleB_id}.deconvolution_metadata.csv
            │   ├── ...
            │   ├── ${sampleN_id}.deconvoluted.h5ad
            │   ├── ${sampleN_id}.deconvolution_metadata.csv
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
└── cell2location_tools
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
./build_docker_images -d docker/cell2location_tools -x

# Build all images in the `docker` directory
./build_docker_images -d docker -x

# Build and push all images in the docker directory, using the `dnastack` container registry
./build_docker_images -d docker -x -c dnastack -p
```

## Tool and library versions

| Image | Major tool versions | Links |
| :- | :- | :- |
| cell2location_tools | <ul><li>[cell2location](https://github.com/BayraktarLab/cell2location)</li><li>[PyTorch](https://pytorch.org/)</li></ul> | [Dockerfile](docker/cell2location_tools) |
| rctd_tools | <ul><li>[spacexr](https://github.com/dmcable/spacexr)</li></ul> | [Dockerfile](docker/rctd_tools) — *to be implemented* |
| util | <ul><li>[google-cloud-cli 524.0.0](https://cloud.google.com/sdk/docs/release-notes#52400_2025-05-28)</li></ul> | [Dockerfile](https://github.com/ASAP-CRN/wf-common/tree/main/docker/util) |


# wdl-ci

[`wdl-ci`](https://github.com/DNAstack/wdl-ci) provides tools to validate and test workflows and tasks written in [Workflow Description Language (WDL)](https://github.com/openwdl/wdl). In addition to the tests packaged in `wdl-ci`, the [spatial-sc-rnaseq-integration-wdl-ci-custom-test-dir](./spatial-sc-rnaseq-integration-wdl-ci-custom-test-dir) is a directory containing custom WDL-based tests that are used to test workflow tasks. `wdl-ci` in this repository is set up to run on pull request.

In general, `wdl-ci` will use inputs provided in the [wdl-ci.config.json](./wdl-ci.config.json) and compare current outputs and validated outputs based on changed tasks/workflows to ensure outputs are still valid by meeting the criteria in the specified tests. For example, if the `cell2location_deconvolute_visium` task was changed, then this task would be submitted and that output would be considered the "current output". When inspecting the deconvoluted AnnData object, there is a test specified in the [wdl-ci.config.json](./wdl-ci.config.json) called `check_hdf5`. The test will compare the "current output" and "validated output" (provided in the [wdl-ci.config.json](./wdl-ci.config.json)) to make sure that the `.h5ad` file is still a valid HDF5 file.
