version 1.0

import "../structs.wdl"

# Reference: downsampled H5AD from the sc-rnaseq-wf cohort_analysis (raw UMI counts).
# Spatial input: spatial-transcriptomics-wf Visium per-sample post-QC H5ADs.
# Deconvolution: Cell2location (GPU)

workflow deconvolution {
	meta {
		name: "scRNA-seq based spatial deconvolution"
		category: "Spatial Transcriptomics"
		description: "Deconvolute per-sample Visium data against an scRNA-seq cohort reference using Cell2location"
		version: "v1.0.0"
		author: "DNAstack Bioinformatics"
		email: "support@dnastack.com"
	}

	input {
		# Required when run_stages includes Stage 1
		File? sc_rnaseq_reference_adata

		# Required when run_stages excludes Stage 1
		File? sc_rnaseq_reference_cell2location_model

		# List of per-sample adatas
		Array[SpatialVisiumSample] visium_samples

		# Cell2location deconvolution parameters
		String batch_key
		String cell_type_key
		Int n_cells_per_spot
		Int detection_alpha
		Int max_epochs_ref
		Int max_epochs_spatial

		# Clustering parameters (on deconvolved cell type fractions, not gene expression)
		Int n_comps
		Int n_neighbors
		Array[Float] leiden_res

		# Stages to run; see parameter_meta below
		String run_stages
		Int preemptible

		String raw_data_path_prefix
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	# run_stages (e.g. "2,3"), always single comma-separated digits from {1,2,3}, never multi-digit.
	Boolean include_stage1 = sub(run_stages, "1", "x") != run_stages
	Boolean include_stage2 = sub(run_stages, "2", "x") != run_stages
	Boolean include_stage3 = sub(run_stages, "3", "x") != run_stages

	# This is the cell2location_tools docker image tag (threaded into all 3 tasks below) and the
	# raw-bucket path segment for their outputs -- bump this one value to bump both.
	String cell2location_version = "1.0.0"
	String build_reference_raw_data_path = "~{raw_data_path_prefix}/build_reference/~{cell2location_version}"
	String deconvolute_visium_raw_data_path = "~{raw_data_path_prefix}/deconvolute_visium/~{cell2location_version}"
	String cluster_visium_raw_data_path = "~{raw_data_path_prefix}/cluster_visium/~{cell2location_version}"

	# Stage 1: reference signature, built once to deconvolute all samples -- not per sample.
	if (include_stage1) {
		call cell2location_build_reference {
			input:
				sc_rnaseq_reference_adata = select_first([sc_rnaseq_reference_adata]),
				batch_key = batch_key,
				cell_type_key = cell_type_key,
				max_epochs_ref = max_epochs_ref,
				preemptible = preemptible,
				raw_data_path = build_reference_raw_data_path,
				workflow_info = workflow_info,
				billing_project = billing_project,
				container_registry = container_registry,
				cell2location_version = cell2location_version,
				zones = zones
		}
	}

	# Either the reference signature just built above, or a user-supplied pre-computed inf_aver.csv.gz file.
	File? resolved_reference_model = if include_stage1
		then cell2location_build_reference.inf_aver_csv_gz
		else sc_rnaseq_reference_cell2location_model #!FileCoercion

	# scatter flat to deconvolute all spatial samples
	scatter (sample in visium_samples) {
		if (include_stage2) {
			call cell2location_deconvolute_visium {
				input:
					sample_id = sample.sample_id,
					spatial_adata = sample.qc_adata_object,
					sc_rnaseq_reference_cell2location_model = select_first([resolved_reference_model]),
					n_cells_per_spot = n_cells_per_spot,
					detection_alpha = detection_alpha,
					max_epochs_spatial = max_epochs_spatial,
					preemptible = preemptible,
					raw_data_path = deconvolute_visium_raw_data_path,
					workflow_info = workflow_info,
					billing_project = billing_project,
					container_registry = container_registry,
					cell2location_version = cell2location_version,
					zones = zones
			}
		}

		if (include_stage3) {
			call cell2location_cluster_visium {
				input:
					sample_id = sample.sample_id,
					spatial_deconvolution_adata = select_first([cell2location_deconvolute_visium.deconvoluted_adata_object]), #!FileCoercion
					n_comps = n_comps,
					n_neighbors = n_neighbors,
					leiden_res = leiden_res,
					preemptible = preemptible,
					raw_data_path = cluster_visium_raw_data_path,
					workflow_info = workflow_info,
					billing_project = billing_project,
					container_registry = container_registry,
					cell2location_version = cell2location_version,
					zones = zones
			}
		}
	}

	output {
		# Flat over all samples; null for a sample when run_stages excludes Stage 2.
		Array[String?] visium_deconvoluted_adata_objects = cell2location_deconvolute_visium.deconvoluted_adata_object
		Array[String?] visium_deconvolution_metadata_csvs = cell2location_deconvolute_visium.deconvolution_metadata_csv
		Array[String?] visium_deconvolution_training_history_pngs = cell2location_deconvolute_visium.stage2_spatial_training_history_png
		Array[String?] visium_deconvolution_timing_tsvs = cell2location_deconvolute_visium.timing_tsv
		# Flat over all samples; null for a sample when run_stages excludes Stage 3.
		Array[String?] visium_deconvoluted_clustered_adata_objects = cell2location_cluster_visium.deconvoluted_adata_object
		Array[String?] visium_deconvolution_clustering_metadata_csvs = cell2location_cluster_visium.deconvolution_metadata_csv
		Array[String?] visium_deconvolution_clustering_timing_tsvs = cell2location_cluster_visium.timing_tsv
		# Built once; null when run_stages excludes Stage 1.
		String? sc_rnaseq_reference_inf_aver_csv_gz = cell2location_build_reference.inf_aver_csv_gz
		String? sc_rnaseq_reference_training_history_png = cell2location_build_reference.stage1_ref_training_history_png
		String? sc_rnaseq_reference_build_timing_tsv = cell2location_build_reference.ref_build_timing_tsv
		# Exposed so callers (e.g. main.wdl's promotion wiring) can reconstruct this run's raw
		# output paths without duplicating this constant.
		String cell2location_version_value = cell2location_version
	}

	parameter_meta {
		# Inputs
		sc_rnaseq_reference_adata: {help: "scRNA-seq reference AnnData (.h5ad) used to build the Cell2location inf_aver.csv.gz deconvolution reference; required when run_stages includes Stage 1."}
		sc_rnaseq_reference_cell2location_model: {help: "Pre-computed inf_aver.csv.gz deconvolution reference; required when run_stages includes Stage 2 without Stage 1."}
		visium_samples: {help: "Flat list of Visium per-sample adatas (raw UMI counts, post-QC) to deconvolute."}
		batch_key: {help: "Column in reference .obs used as batch covariate in the NB regression model."}
		cell_type_key: {help: "Column in reference .obs to use as the cell type label."}
		n_cells_per_spot: {help: "Cell2location N_cells_per_location prior; 8 for brain Visium."}
		detection_alpha: {help: "Cell2location detection rate prior (20–200); higher values assume low technical variability."}
		max_epochs_ref: {help: "NB regression (Stage 1) training epochs."}
		max_epochs_spatial: {help: "Cell2location spatial mapping (Stage 2) training epochs."}
		n_comps: {help: "PCA components for clustering on deconvolved cell type fractions (Stage 3)."}
		n_neighbors: {help: "Neighbourhood size for clustering on deconvolved cell type fractions (Stage 3)."}
		leiden_res: {help: "Leiden resolutions for clustering on deconvolved cell type fractions (Stage 3); one clustering per resolution."}
		run_stages: {help: "Comma-separated subset of {1,2,3} — which stages to run: 1 = build Cell2location inf_aver.csv.gz deconvolution reference, 2 = Cell2location spatial mapping (per sample), 3 = Leiden clustering on deconvolved cell type fractions."}
		preemptible: {help: "Number of preemptible-VM attempts before Cromwell falls back to on-demand; 0 forces on-demand from the first attempt."}
		raw_data_path_prefix: {help: "Raw data bucket path prefix; each task's outputs are uploaded to its own versioned subpath below this prefix (`<raw_data_path_prefix>/<task_name>/<cell2location_version>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest alongside each uploaded output."}
		billing_project: {help: "GCP billing project for requester-pays bucket access."}
		container_registry: {help: "Container registry hosting workflow Docker images."}
		zones: {help: "GCP zones where compute will run."}

		# Outputs
		visium_deconvoluted_adata_objects: {help: "GCS paths to Cell2location deconvoluted AnnData objects (.h5ad), one per Visium sample; null for a sample when run_stages excludes Stage 2."}
		visium_deconvolution_metadata_csvs: {help: "GCS paths to per-spot deconvolution metadata CSV cell type fractions; null for a sample when run_stages excludes Stage 2."}
		visium_deconvolution_training_history_pngs: {help: "GCS paths to Cell2location spatial mapping (Stage 2) ELBO training-loss plot, one per Visium sample; null for a sample when run_stages excludes Stage 2."}
		visium_deconvolution_timing_tsvs: {help: "GCS paths to per-step timing breakdown (Stage 2 run), one per Visium sample; null for a sample when run_stages excludes Stage 2."}
		visium_deconvoluted_clustered_adata_objects: {help: "GCS paths to Cell2location deconvoluted AnnData objects (.h5ad) with Leiden cluster columns added, one per Visium sample; null for a sample when run_stages excludes Stage 3."}
		visium_deconvolution_clustering_metadata_csvs: {help: "GCS paths to per-spot deconvolution metadata CSV cell type fractions with Leiden cluster columns added; null for a sample when run_stages excludes Stage 3."}
		visium_deconvolution_clustering_timing_tsvs: {help: "GCS paths to per-step timing breakdown (Stage 3 run), one per Visium sample; null for a sample when run_stages excludes Stage 3."}
		sc_rnaseq_reference_inf_aver_csv_gz: {help: "GCS path to the Cell2location reference signature (<ref_stem>_inf_aver.csv.gz), built once for the whole cohort; null when run_stages excludes Stage 1."}
		sc_rnaseq_reference_training_history_png: {help: "GCS path to the NB regression (Stage 1) ELBO training-loss plot, built once for the whole cohort; null when run_stages excludes Stage 1."}
		sc_rnaseq_reference_build_timing_tsv: {help: "GCS path to the per-step timing breakdown (Stage 1 run), built once for the whole cohort; null when run_stages excludes Stage 1."}
		cell2location_version_value: {help: "Versioned subpath each cell2location task's raw outputs are uploaded under; always populated."}
	}
}

task cell2location_build_reference {
	meta {
		description: "Build the scRNA-seq reference signature (Cell2location inf_aver.csv.gz)"
	}

	input {
		File sc_rnaseq_reference_adata
		String batch_key
		String cell_type_key
		Int max_epochs_ref
		Int preemptible

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String cell2location_version
		String zones
	}

	Int disk_size = ceil(size(sc_rnaseq_reference_adata, "GB") + 50)
	# Matches cell2location_deconv.py's own Path(ref_h5ad).stem naming of out_inf_aver.
	String ref_stem = basename(sc_rnaseq_reference_adata, ".h5ad")

	command <<<
		set -euo pipefail

		/usr/bin/time -v python3 /opt/scripts/cell2location_deconv.py \
			--run-stages 1 \
			--sc-rnaseq-reference-adata ~{sc_rnaseq_reference_adata} \
			--batch-key ~{batch_key} \
			--cell-type-key ~{cell_type_key} \
			--max-epochs-ref ~{max_epochs_ref} \
			--output-dir . \
			--promote

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{ref_stem}_inf_aver.csv.gz" \
			-o "stage1_ref_training_history.png" \
			-o "~{ref_stem}_ref_build.timing.tsv"
	>>>

	output {
		String inf_aver_csv_gz = "~{raw_data_path}/~{ref_stem}_inf_aver.csv.gz"
		String stage1_ref_training_history_png = "~{raw_data_path}/stage1_ref_training_history.png"
		String ref_build_timing_tsv = "~{raw_data_path}/~{ref_stem}_ref_build.timing.tsv"
	}

	runtime {
		docker: "~{container_registry}/cell2location_tools:~{cell2location_version}"
		cpu: 8
		memory: "52 GB"
		disks: "local-disk " + disk_size + " SSD"
		preemptible: preemptible
		bootDiskSizeGb: 40
		gpuType: "nvidia-tesla-t4"
		gpuCount: 1
		zones: zones
	}

	parameter_meta {
		raw_data_path: {help: "Raw data bucket path for Stage 1 outputs; location of raw bucket to upload task outputs to (`<raw_data_path_prefix>/build_reference/<cell2location_version>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest alongside each uploaded output."}
		billing_project: {help: "GCP billing project for requester-pays bucket access."}
		container_registry: {help: "Container registry hosting workflow Docker images."}
		cell2location_version: {help: "cell2location_tools Docker image tag; also the versioned raw-bucket subpath segment."}
		zones: {help: "GCP zones where compute will run."}
	}
}

task cell2location_cluster_visium {
	meta {
		description: "Leiden-cluster one Visium sample's Cell2location deconvoluted cell type fractions"
	}

	input {
		String sample_id
		File spatial_deconvolution_adata

		Int n_comps
		Int n_neighbors
		Array[Float] leiden_res
		Int preemptible

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String cell2location_version
		String zones
	}

	Int disk_size = ceil(size(spatial_deconvolution_adata, "GB") + 20)
	# Fix expected sample_id in outfile name. This is a temporary fix while spatial-transcriptomics-wf is re-run to fix underlying typos.
	String corrected_sample_id = sub(sample_id, "^ASAP_PMBDS_", "ASAP_PMDBS_")

	command <<<
		set -euo pipefail

		/usr/bin/time -v python3 /opt/scripts/cell2location_deconv.py \
			--run-stages 3 \
			--spatial-deconvolution-adata ~{spatial_deconvolution_adata} \
			--sample-id ~{sample_id} \
			--n-comps ~{n_comps} \
			--n-neighbors ~{n_neighbors} \
			--leiden-res ~{sep=" " leiden_res} \
			--output-dir . \
			--promote

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{corrected_sample_id}.deconvoluted_clustered.h5ad" \
			-o "~{corrected_sample_id}.deconvolution_clustering_metadata.csv.gz" \
			-o "~{corrected_sample_id}.timing.tsv"
	>>>

	output {
		String deconvoluted_adata_object = "~{raw_data_path}/~{corrected_sample_id}.deconvoluted_clustered.h5ad"
		String deconvolution_metadata_csv = "~{raw_data_path}/~{corrected_sample_id}.deconvolution_clustering_metadata.csv.gz"
		String timing_tsv = "~{raw_data_path}/~{corrected_sample_id}.timing.tsv"
	}

	runtime {
		docker: "~{container_registry}/cell2location_tools:~{cell2location_version}"
		cpu: 4
		memory: "16 GB"
		disks: "local-disk " + disk_size + " SSD"
		preemptible: preemptible
		zones: zones
	}

	parameter_meta {
		raw_data_path: {help: "Raw data bucket path for Stage 3 outputs; location of raw bucket to upload task outputs to (`<raw_data_path_prefix>/cluster_visium/<cell2location_version>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest alongside each uploaded output."}
		billing_project: {help: "GCP billing project for requester-pays bucket access."}
		container_registry: {help: "Container registry hosting workflow Docker images."}
		cell2location_version: {help: "cell2location_tools Docker image tag; also the versioned raw-bucket subpath segment."}
		zones: {help: "GCP zones where compute will run."}
	}
}

task cell2location_deconvolute_visium {
	meta {
		description: "Deconvolute one Visium sample against a pre-built Cell2location inf_aver.csv.gz deconvolution reference"
	}

	input {
		String sample_id
		File spatial_adata
		File sc_rnaseq_reference_cell2location_model

		Int n_cells_per_spot
		Int detection_alpha
		Int max_epochs_spatial
		Int preemptible

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String cell2location_version
		String zones
	}

	Int disk_size = ceil(size([spatial_adata, sc_rnaseq_reference_cell2location_model], "GB") + 50)
	# Fix expected sample_id in outfile name. This is a temporary fix while spatial-transcriptomics-wf is re-run to fix underlying typos.
	String corrected_sample_id = sub(sample_id, "^ASAP_PMBDS_", "ASAP_PMDBS_")

	command <<<
		set -euo pipefail

		/usr/bin/time -v python3 /opt/scripts/cell2location_deconv.py \
			--run-stages 2 \
			--spatial-h5ad ~{spatial_adata} \
			--sc-rnaseq-reference-cell2location-model ~{sc_rnaseq_reference_cell2location_model} \
			--sample-id ~{sample_id} \
			--n-cells-per-spot ~{n_cells_per_spot} \
			--detection-alpha ~{detection_alpha} \
			--max-epochs-spatial ~{max_epochs_spatial} \
			--output-dir . \
			--promote

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{corrected_sample_id}.deconvoluted.h5ad" \
			-o "~{corrected_sample_id}.deconvolution_metadata.csv.gz" \
			-o "~{corrected_sample_id}.stage2_spatial_training_history.png" \
			-o "~{corrected_sample_id}.timing.tsv"
	>>>

	output {
		String deconvoluted_adata_object = "~{raw_data_path}/~{corrected_sample_id}.deconvoluted.h5ad"
		String deconvolution_metadata_csv = "~{raw_data_path}/~{corrected_sample_id}.deconvolution_metadata.csv.gz"
		String stage2_spatial_training_history_png = "~{raw_data_path}/~{corrected_sample_id}.stage2_spatial_training_history.png"
		String timing_tsv = "~{raw_data_path}/~{corrected_sample_id}.timing.tsv"
	}

	runtime {
		docker: "~{container_registry}/cell2location_tools:~{cell2location_version}"
		cpu: 8
		memory: "52 GB"
		disks: "local-disk " + disk_size + " SSD"
		preemptible: preemptible
		bootDiskSizeGb: 40
		gpuType: "nvidia-tesla-t4"
		gpuCount: 1
		zones: zones
	}

	parameter_meta {
		raw_data_path: {help: "Raw data bucket path for Stage 2 outputs; location of raw bucket to upload task outputs to (`<raw_data_path_prefix>/deconvolute_visium/<cell2location_version>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest alongside each uploaded output."}
		billing_project: {help: "GCP billing project for requester-pays bucket access."}
		container_registry: {help: "Container registry hosting workflow Docker images."}
		cell2location_version: {help: "cell2location_tools Docker image tag; also the versioned raw-bucket subpath segment."}
		zones: {help: "GCP zones where compute will run."}
	}
}
