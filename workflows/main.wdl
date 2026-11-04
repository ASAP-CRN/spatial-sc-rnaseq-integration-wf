version 1.0

# Spatial deconvolution using scRNA-seq data workflow entrypoint.
# Inputs are processed outputs from sc-rnaseq-wf and spatial-transcriptomics-wf;
# this workflow does not handle raw FASTQ or SpaceRanger/CellRanger execution.
#
# The scRNA-seq reference can obtained from a downsampled version of the sc-rnaseq-wf cohort
# alone or supplemented with additional data (e.g. Siletti and/or Kamath cells).
# See `scripts/build_deconv_reference_h5ad.py` for details.
#
# Each spatial sample carries its own per-sample H5AD (raw UMI counts, post-QC) from
# spatial-transcriptomics-wf.
#
# Stage 1: Build the scRNA-seq reference inf_aver file from the H5AD for Cell2location
# Stage 2: Deconvolute spatial samples using the reference inf_aver file with Cell2location
# Stage 3: Use the deconvoluted spatial samples for downstream analysis (e.g., Leiden clustering)
#
# Deconvolution method: Cell2location (Bayesian hierarchical model, GPU-accelerated).

import "structs.wdl"
import "../wf-common/wdl/tasks/get_workflow_metadata.wdl" as GetWorkflowMetadata
import "../wf-common/wdl/tasks/upload_final_outputs.wdl" as UploadFinalOutputs #!UnusedImport -- reserved for planned output-upload wiring, not yet implemented
import "integration/integration.wdl" as Integration

workflow spatial_sc_rnaseq_integration {
	meta {
		name: "scRNA-seq based spatial deconvolution"
		category: "Spatial Transcriptomics"
		description: "Deconvolute per-sample Visium data against an scRNA-seq cohort reference using Cell2location"
		version: "v1.0.0"
		author: "DNAstack Bioinformatics"
		email: "support@dnastack.com"
	}

	input {
		# scRNA-seq reference H5AD file.
		File sc_rnaseq_reference_adata

		# spatial-transcriptomics-wf post-QC spatial H5ADs.
		Array[SpatialVisiumProject] spatial_visium_projects = []

		# Optional pre-computed Cell2location reference signatures (inf_aver.csv.gz).
		# Required when run_stages excludes Stage 1 (e.g. "2,3"); ignored otherwise.
		File? sc_rnaseq_reference_cell2location_model

		# Cell2location deconvolution parameters
		String batch_key = "batch_id"
		String cell_type_key = "cell_type_pred"
		Int n_cells_per_spot = 8
		Int detection_alpha = 100
		Int max_epochs_ref = 250
		Int max_epochs_spatial = 30000

		# Clustering parameters
		Int n_comps = 30
		Int n_neighbors = 15
		Array[Float] leiden_res = [0.2, 0.5, 1.0]

		# Stages to run; see Integration.integration
		String run_stages = "1,2,3"
		Int preemptible = 3

		String cohort_raw_data_bucket_spatial_sc_rnaseq

		String container_registry
		String zones = "us-central1-c us-central1-f"
	}

	call GetWorkflowMetadata.get_workflow_metadata {
		input:
			zones = zones
	}

	# Flatten per-sample data across all projects (WDL 1.0: scatter then flatten)
	scatter (project in spatial_visium_projects) {
		Array[SpatialVisiumSample] project_samples = project.samples
	}
	Array[SpatialVisiumSample] visium_samples = flatten(project_samples)

	call Integration.integration {
		input:
			sc_rnaseq_reference_adata = sc_rnaseq_reference_adata,
			sc_rnaseq_reference_cell2location_model = sc_rnaseq_reference_cell2location_model,
			visium_samples = visium_samples,
			batch_key = batch_key,
			cell_type_key = cell_type_key,
			n_cells_per_spot = n_cells_per_spot,
			detection_alpha = detection_alpha,
			max_epochs_ref = max_epochs_ref,
			max_epochs_spatial = max_epochs_spatial,
			n_comps = n_comps,
			n_neighbors = n_neighbors,
			leiden_res = leiden_res,
			run_stages = run_stages,
			preemptible = preemptible,
			raw_data_path = cohort_raw_data_bucket_spatial_sc_rnaseq,
			billing_project = get_workflow_metadata.billing_project,
			container_registry = container_registry,
			zones = zones
	}

	output {
		# Array[File] — flat over all samples (across all projects). Only meaningful when
		# Stage 2 and/or 3 ran for every sample; see integration.wdl's output{} note.
		Array[File] visium_deconvoluted_adata_objects = integration.visium_deconvoluted_adata_objects
		Array[File] visium_deconvolution_metadata_csvs = integration.visium_deconvolution_metadata_csvs
		Array[Array[File]] visium_inf_aver_csv_gzs = integration.visium_inf_aver_csv_gzs
	}

	parameter_meta {
		# Inputs
		sc_rnaseq_reference_adata: {help: "scRNA-seq reference AnnData (.h5ad) used as the Cell2location deconvolution reference for all spatial samples"}
		spatial_visium_projects: {help: "Processed spatial outputs (from spatial-transcriptomics-wf); samples are flattened across projects and deconvoluted independently. []"}
		sc_rnaseq_reference_cell2location_model: {help: "Pre-computed inf_aver.csv.gz from a previous run; required when run_stages excludes Stage 1 (e.g. '2,3')"}
		batch_key: {help: "Column in reference .obs used as batch covariate in the NB regression model. ['batch_id']"}
		cell_type_key: {help: "Column in reference .obs to use as the cell type label. ['cell_type_pred']"}
		n_cells_per_spot: {help: "Cell2location N_cells_per_location prior; 8 for brain Visium. [8]"}
		detection_alpha: {help: "Cell2location detection rate prior (20–200); higher values assume low technical variability. [100]"}
		max_epochs_ref: {help: "NB regression (Stage 1) training epochs. [250]"}
		max_epochs_spatial: {help: "Cell2location spatial mapping (Stage 2) training epochs. [30000]"}
		n_comps: {help: "Number of PCA components for downstream spot clustering. [30]"}
		n_neighbors: {help: "Neighbourhood size for downstream spot clustering. [15]"}
		leiden_res: {help: "Leiden resolutions controlling the coarseness of clustering. [0.2, 0.5, 1.0]"}
		run_stages: {help: "Comma-separated subset of {1,2,3} — which stages to run: 1 = NB regression (reference signatures), 2 = Cell2location spatial mapping, 3 = downstream Leiden clustering; e.g. '1' builds the reference only, '2,3' loads a pre-computed reference and deconvolutes+clusters, '3' loads a pre-computed deconvolution and clusters only. ['1,2,3']"}
		preemptible: {help: "Number of preemptible-VM attempts before Cromwell falls back to on-demand; 0 forces on-demand from the first attempt; use 0 for long production runs where a preemption-triggered restart-from-scratch would be too costly. [3]"}
		cohort_raw_data_bucket_spatial_sc_rnaseq: {help: "GCS path to raw data bucket; intermediate outputs are stored here"}
		container_registry: {help: "Container registry hosting workflow Docker images"}
		zones: {help: "GCP zones where compute will run. ['us-central1-c us-central1-f']"}
	}
}
