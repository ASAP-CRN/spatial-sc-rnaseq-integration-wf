# spatial-sc-rnaseq-integration-wf

Repo for testing and developing a workflow that integrates single-cell/single-nucleus RNA-seq data with spatial transcriptomics data. This version (v1.0.0) covers deconvolution of 10x Visium data using [Cell2Location](https://cell2location.readthedocs.io/en/latest/index.html) only.

This workflow consumes outputs from the following two upstream ASAP workflows and does not re-run CellRanger, SpaceRanger, or cell type labelling.
* [sc-rnaseq-wf](https://github.com/ASAP-CRN/sc-rnaseq-wf)
* [spatial-transcriptomics-wf](https://github.com/ASAP-CRN/spatial-transcriptomics-wf)

Common workflows, tasks, utility scripts, and docker images reused across harmonized ASAP workflows are defined in [the wf-common repository](https://github.com/ASAP-CRN/wf-common).


# Table of contents

- [Workflows](#workflows)
    - [Spatial deconvolution](#spatial-deconvolution)
        - [Cell2Location](#cell2location)
- [Inputs](#inputs)
    - [Structures](#structures)
    - [Input JSON's](#input-jsons)
- [Outputs](#outputs)
    - [Output structure](#output-structure)
    - [Promoting data raw → dev/uat](#promoting-data-raw--devuat)
    - [Promoting data dev/uat → curated](#promoting-data-devuat--curated)
- [Docker images](#docker-images)
- [Tool and library versions](#tool-and-library-versions)
- [wdl-ci test runs](#wdl-ci-test-runs)
- [Registering workflows in Workbench](#registering-workflows-in-workbench)
- [Submitting production runs via omics CLI](#submitting-production-runs-via-omics-cli)


# Workflows

- Workflow name: `spatial_sc_rnaseq_integration`
- Entrypoint: [workflows/main.wdl](workflows/main.wdl)
- Input template: [Input JSON's](#input-jsons)
- Workflow scripts: [workflows/ directory](./workflows/)
- Python scripts: [docker/cell2location_tools/scripts/](docker/cell2location_tools/scripts/)


## Spatial deconvolution

The deconvolution workflow receives processed AnnData objects (`.h5ad`) from the two upstream ASAP workflows:

| Source workflow | Level | Output consumed |
| :- | :- | :- |
| [sc-rnaseq-wf](https://github.com/ASAP-CRN/sc-rnaseq-wf) `cohort_analysis` downsampled + spike-ins from [Kamath et al (2022)](https://www.nature.com/articles/s41593-022-01061-1) and [Siletti et al (2023)](https://www.science.org/doi/10.1126/science.add7046) | cohort H5AD used to build scRNA-seq reference | `final_adata_object` |
| [spatial-transcriptomics-wf](https://github.com/ASAP-CRN/spatial-transcriptomics-wf) `spatial_visium preprocess` | per-sample, raw UMI counts, all genes--not only HVGs| `qc_adata_objects` (one per sample) |

The Cell2Location deconvolution sub-workflow ([workflows/deconvolution/deconvolution.wdl](workflows/deconvolution/deconvolution.wdl)) runs as three independent WDL tasks/stages, gated by `run_stages`:
1. **Stage 1 — Reference build**: train a Negative Binomial regression model on the scRNA-seq reference to estimate per-cell-type gene expression signatures (`inf_aver`). Runs **once for all spatial samples**.
2. **Stage 2 — Spatial mapping**: run Cell2location to map cell type abundances onto Visium spots. Runs once **per Visium sample**, submitted in batches of 15 samples in parallel by [scripts/submit_cell2location_2and3_prod.py](scripts/submit_cell2location_2and3_prod.py) to meet the project on-demand NVIDIA T4 GPU quota. Writes `{sample_id}.deconvoluted.h5ad` / `{sample_id}.deconvolution_metadata.csv.gz`.
3. **Stage 3 — Clustering**: PCA → neighbours → UMAP → Leiden clustering on each spot's **deconvolved cell type fractions** — not on gene expression. Runs once **per Visium sample**, scattered in parallel along stage 2. Writes `{sample_id}.deconvoluted_clustered.h5ad` / `{sample_id}.deconvolution_clustering_metadata.csv.gz`. **Requires Stage 2 in the same submission** — standalone `run_stages` values of `"3"`/`"1,3"` aren't supported.


### Cell2Location

| Method | Reference | Language | Docker image | Hardware requiriments |
| :- | :- | :- | :- | :- |
| [Cell2location](https://github.com/BayraktarLab/cell2location), Bayesian hierarchical model | [Kleshchevnikov *et al.*, *Nat Biotechnol* 2022](https://pubmed.ncbi.nlm.nih.gov/35027729/) | Python | `cell2location_tools` | Stages 1 and 2 require a GPU, stage 3 requires CPU only. This repo runs in a VM with the following specs: CPU = 8, MEM = 52 GB, GPU = 1× nvidia-tesla-t4 |
* See [Registering workflow and submitting production runs via omics CLI](#registering-workflow-and-submitting-production-runs-via-omics-cli) for details.
* An alternative method ([RCTD](https://github.com/dmcable/spacexr)) was used for a [benchmark analysis vs. Ma et al (2025) data](reference_building/benchmark/README.md) but it was not implemented in this repo v1.0.0.


# Inputs

Input templates can be found at [Input JSON's](#input-jsons).

| WDL field | Type | Description | Stage 1 | Stages 2+3 | Data promotion |
| :- | :- | :- | :- | :- | :- |
| run_stages | String? | Comma-separated subset of `{1,2,3}` — which stages to run. `"1"` builds the `inf_aver.csv.gz` reference; `"2"` loads it and deconvolutes spatial samples; `"3"` Leiden-clusters each sample's deconvolved cell type fractions. Stage 3 requires Stage 2 in the same submission — `"3"` and `"1,3"` alone are not currently supported. ['1,2,3'] | ✓ | ✓ | ✓ |
| cohort_raw_data_bucket_spatial_sc_rnaseq | String | Raw data bucket root where outfiles are uploaded. | ✓ | ✓ | ✓ |
| container_registry | String | Container registry where workflow Docker images are hosted. | ✓ | ✓ | ✓ |
| zones | String? | GCP zones where compute takes place. ['us-central1-c us-central1-f'] | ✓ | ✓ | ✓ |
| preemptible | Int? | Number of preemptible-VM attempts before Cromwell falls back to on-demand; `0` forces on-demand. Production runs use `0` because GPU spot availability was unpredictable and kept failing using `3`. [3] | ✓ | ✓ | - |
| batch_key | String? | Column in reference .obs used as batch covariate in the NB regression model. ['batch_id'] | ✓ | - | - |
| cell_type_key | String? | Column in reference .obs to use as the cell type label. ['cell_type_pred'] | ✓ | - | - |
| sc_rnaseq_reference_adata | File? | scRNA-seq reference AnnData; used to deconvolute all spatial samples. Required when `run_stages` includes Stage 1. | ✓ | - | - |
| max_epochs_ref | Int? | NB regression (Stage 1) training epochs. [250] | ✓ | - | - |
| sc_rnaseq_reference_cell2location_model | File? | Pre-computed `<ref_stem>_inf_aver.csv.gz`. Required when `run_stages` excludes Stage 1. | - | ✓ | - |
| n_cells_per_spot | Int? | Cell2location N_cells_per_location prior; 8 for brain Visium. [8] | - | ✓ | - |
| detection_alpha | Int? | Cell2location detection rate prior (20–200); higher values assume low technical variability. [100] | - | ✓ | - |
| max_epochs_spatial | Int? | Cell2location spatial mapping (Stage 2) training epochs. Note early-ending is implemented. [30000] | - | ✓ | - |
| n_comps | Int? | PCA components for clustering on deconvolved cell type fractions (Stage 3). [30] | - | ✓ | - |
| n_neighbors | Int? | Neighbourhood size for clustering on deconvolved cell type fractions (Stage 3). [15] | - | ✓ | - |
| leiden_res | Array[Float]? | Leiden resolutions for clustering on deconvolved cell type fractions (Stage 3); one clustering per resolution. [0.2, 0.5, 1.0] | - | ✓ | - |
| spatial_visium_projects | Array[[`SpatialVisiumProject`](#spatialvisiumproject)]? | Processed project-level outputs from spatial-transcriptomics-wf (Visium) to deconvolute per sample. [[ ]] | - | ✓ | ✓ |
| staging_data_buckets | Array[String]? | Set of buckets to promote completed per-sample deconvolution outputs to (e.g. dev/uat). [[ ]] | - | - | ✓ |
| crn_release_version | String? | Release version (e.g. `v6.0.0`) to stage outputs under; required when `staging_data_buckets` is non-empty. [""] | - | - | ✓ |


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

One Visium sample, with raw UMI counts and all genes--not only HVGs from [spatial-transcriptomics-wf](https://github.com/ASAP-CRN/spatial-transcriptomics-wf) (`*cleaned_unfiltered.h5ad`). Source: `spatial_visium preprocess.qc_adata_objects`. These can be obtained from the dataset buckets: `gs://asap-raw-<dataset_name>/workflow_execution/spatial_visium/preprocess/counts_to_adata/<version>`. Spots are those Space Ranger flagged `in_tissue = 1`. This workflow doesn't apply additional spot-level QC/filtering. The user can obtain the [spatial-transcriptomics-wf](https://github.com/ASAP-CRN/spatial-transcriptomics-wf) own's QC metrics from `gs://asap-raw-<dataset_name>/workflow_execution/spatial_visium/qc/<version>`.

| Type | Name | Description |
| :- | :- | :- |
| String | sample_id | Sample identifier; used to name output files. |
| File | qc_adata_object | `spatial_visium preprocess.qc_adata_objects` (`*cleaned_unfiltered.h5ad`) — per-sample, raw UMI counts, all genes--not only HVGs. |

## Input JSON's

There are three kinds of JSON used in this repo. Those used for the ASAP CRN release v6.0.0 are:

- **Stage 1 — Reference build**: [prod_runs/submit_cell2location_1_prod_config.json](prod_runs/submit_cell2location_1_prod_config.json).
  - Copies of the reference h5ad and inf_aver reference files were stored here:     
  `gs://asap-raw-team-scherzer-pmdbs-spatial-visium-mtg/workflow_execution/spatial_visium/deconvolution_analysis/build_reference/`
- **Stages 2 & 3 — Spatial mapping & clustering**: submit batches of 15 spatial samples: [prod_runs/submit_cell2location_2and3_prod_config.json](prod_runs/submit_cell2location_2and3_prod_config.json).
- **Data promotion**: raw → dev/uat: [prod_runs/promote_raw_to_dev_uat_v6.0.0_all95samples.json](prod_runs/promote_raw_to_dev_uat_v6.0.0_all95samples.json).

# Outputs

## Output structure

- Outputs are written per spatial sample; each file is named with the `sample_id`.
- Each task uploads its own outputs to a CRN release version, like:

### Raw data

```bash
asap-raw-<dataset_name>
└── workflow_execution
    └── spatial_visium
        └── deconvolution_analysis
            ├── build_reference
            │   └── ${cell2location_version}
            │       ├── ${ref_stem}_inf_aver.csv.gz
            │       ├── stage1_ref_training_history.png
            │       └── ${ref_stem}_ref_build.timing.tsv
            ├── deconvolute_visium
            │   └── ${cell2location_version}
            │       ├── ${sample_id}.deconvoluted.h5ad
            │       ├── ${sample_id}.deconvolution_metadata.csv.gz
            │       ├── ${sample_id}.stage2_spatial_training_history.png
            │       └── ${sample_id}.timing.tsv
            └── cluster_visium
                └── ${cell2location_version}
                    ├── ${sample_id}.deconvoluted_clustered.h5ad
                    ├── ${sample_id}.deconvolution_clustering_metadata.csv.gz
                    └── ${sample_id}.timing.tsv
```

### Staging data (dev/uat)

```bash
asap-{dev|uat}-<dataset_name>
└── spatial_visium
    └── release
        └── ${crn_release_version}
            └── deconvolution_analysis
                ├── ${sampleA_id}.deconvoluted_clustered.h5ad
                ├── ${sampleA_id}.deconvolution_clustering_metadata.csv.gz
                ├── ${sampleA_id}.timing.tsv
                ├── ...
                ├── ${sampleN_id}.deconvoluted_clustered.h5ad
                ├── ${sampleN_id}.deconvolution_clustering_metadata.csv.gz
                ├── ${sampleN_id}.timing.tsv
                └── MANIFEST.tsv  # combined, one row per promoted file
```

Only Stage 3's clustered outputs are promoted to dev/uat.

## Promoting data raw → dev/uat

Wired directly in `main.wdl` via `upload_deconvolution_outputs`, gated on `staging_data_buckets` being non-empty. Setting `run_stages=""` promotes already-computed outputs. See [Registering workflow and submitting production runs via omics CLI](#registering-workflow-and-submitting-production-runs-via-omics-cli) for the 95-sample promotion example.

## Promoting data dev/uat → curated

**Not yet exercised for this workflow.** It will be wired for the December release.

# Docker images

Docker images are defined in [the docker/ directory](docker), each with a `build.env`, `Dockerfile`, and `cloudbuild.yaml` (see [docker/cell2location_tools](docker/cell2location_tools)). Images are built and pushed to Artifact Registry via [Google Cloud Build](https://cloud.google.com/build/docs/overview) (`gcloud builds submit --config=docker/cell2location_tools/cloudbuild.yaml`).

Image details as of [PR#4](https://github.com/ASAP-CRN/spatial-sc-rnaseq-integration-wf/pull/4) (workflow v1.0.0 release, Dec 2026)
* `ID: 7a97a551-b395-46f1-ab7c-e84838706ca3`.
* `IMAGE_TAG` remained `1.0.0` across all builds so far.

### Tool and library versions

| Image | Major tool versions | Links |
| :- | :- | :- |
| cell2location_tools | <ul><li>[cell2location](https://github.com/BayraktarLab/cell2location)</li><li>[PyTorch](https://pytorch.org/)</li></ul> | [Dockerfile](docker/cell2location_tools/Dockerfile) |
| util | <ul><li>[google-cloud-cli 524.0.0](https://cloud.google.com/sdk/docs/release-notes#52400_2025-05-28)</li></ul> | [Dockerfile](https://github.com/ASAP-CRN/wf-common/tree/main/docker/util) |


# `wdl-ci` test runs

[`wdl-ci`](https://github.com/DNAstack/wdl-ci) validates and tests this repo's WDL tasks/workflows, using [wdl-ci.config.json](./wdl-ci.config.json). Two ways to run it:
  1. **GitHub Actions** ([.github/workflows/lint-test-workflows.yml](.github/workflows/lint-test-workflows.yml)) — automatic on every PR (`on: pull_request`), skip with a commit message including `[skip ci]`.
  2. **Local CLI** — run five CI subcommands directly against Workbench/Cromwell, no GitHub Actions involved. See [scripts/submit_wdl_ci_local.sh](scripts/submit_wdl_ci_local.sh) (requires `wdl-ci` and a 1Password env-file).

For CI script testing, these inputs can be used:
- 30,000 cells sc-RNAseq reference: `gs://asap-wdl-ci/workflow_execution/spatial_sc_rnaseq_integration/reference_building/reference_h5ad/deconvolution_reference_30k.h5ad`
- One ASAP spatial sample: `gs://asap-wdl-ci/workflow_execution/spatial_visium/preprocess/counts_to_adata/1.0.0/ASAP_PMBDS_000003_s007_Rep1.cleaned_unfiltered.h5ad`


# Registering workflows in Workbench

A new [DNAstack Workbench](https://docs.omics.ai/docs) wokflow registration is made like:
  ```bash
      omics workbench workflows create --entrypoint workflows/main.wdl --version-name v1.0.0
  ```

Workflow bumps are made like:
  ```bash
      omics workbench workflows versions create \
      workflows/main.wdl workflows/structs.wdl workflows/deconvolution/deconvolution.wdl \
      wf-common/wdl/tasks/get_workflow_metadata.wdl wf-common/wdl/tasks/upload_final_outputs.wdl \
      --workflow <workbenck-workflow-id> \
      --name <workbenck-workflow-registration-version> \
      --entrypoint workflows/main.wdl
  ```

This repo workflow details as of [PR#4](https://github.com/ASAP-CRN/spatial-sc-rnaseq-integration-wf/pull/4) (workflow v1.0.0 release, Dec 2026) are:
* `<workbenck-workflow-id>: 865c2967-9af8-499d-9b70-aa10db6c49da`.
* `<workbenck-workflow-registration-version>: v1.0.5`.

# Submitting production runs via omics CLI

After workflows are registered in [DNAstack Workbench](#registering-workflows-in-workbench), production runs are submitted as follows:

- **Stage 1** — builds the `inf_aver.csv.gz` to be re-used for each spatial sample deconvolution in stage 2.
  ```bash
  scripts/submit_cell2location_1_prod.py \
      --input-json prod_runs/submit_cell2location_1_prod_config.json \
      --promote
  ```

- **Stages 2 and 3** — submits chunks of 15 spatial samples to be deconvoluted and clustered in parallel. Repeat this until all samples have been processed.
  ```bash
  scripts/submit_cell2location_2and3_prod.py \
      --infile-sample-ids-all prod_runs/list_of_spatial_samples_to_submit_all.txt \
      --validate-bucket-meta-tsv \
      --chunk-size 15 \
      --submit-chunks 1 \
      --input-json prod_runs/submit_cell2location_2and3_prod_config.json \
      --promote
  ```
 
- **Data promotion** - promotes all 95 spatial samples' outfiles from raw → dev/uat.
  ```bash
  omics workbench runs submit \
      --workflow <workbenck-workflow-id> \
      --version <workbenck-workflow-registration-version> \
      --workflow-params @prod_runs/promote_raw_to_dev_uat_v6.0.0_all95samples.json
  ```
