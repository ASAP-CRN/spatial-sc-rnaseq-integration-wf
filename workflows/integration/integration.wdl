version 1.0

import "../structs.wdl"

# Reference: downsampled H5AD from the sc-rnaseq-wf cohort_analysis (raw UMI counts).
# Spatial input: spatial-transcriptomics-wf Visium per-sample post-QC H5ADs.
# Deconvolution: Cell2location (GPU)

workflow integration {
	meta {
		name: "scRNA-seq based spatial deconvolution"
		category: "Spatial Transcriptomics"
		description: "Deconvolute per-sample Visium data against an scRNA-seq cohort reference using Cell2location"
		version: "v1.0.0"
		author: "DNAstack Bioinformatics"
		email: "support@dnastack.com"
	}

	input {
		File sc_rnaseq_reference_adata

		# Required when run_stages excludes Stage 1; ignored otherwise.
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

		# Clustering parameters
		Int n_comps
		Int n_neighbors
		Array[Float] leiden_res

		# Stages to run; see cell2location_visium
		String run_stages
		Int preemptible

		String raw_data_path
		String billing_project
		String container_registry
		String zones
	}

	# Visium: scatter flat over all samples (across all projects)
	scatter (sample in visium_samples) {
		call cell2location_visium {
			input:
				sample_id = sample.sample_id,
				spatial_adata = sample.qc_adata_object,
				sc_rnaseq_reference_adata = sc_rnaseq_reference_adata,
				sc_rnaseq_reference_cell2location_model = sc_rnaseq_reference_cell2location_model,
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
				raw_data_path = raw_data_path,
				billing_project = billing_project,
				container_registry = container_registry,
				zones = zones
		}
	}

	output {
		# Array[File] — flat over all samples (across all projects). Only meaningful when
		# Stage 2 and/or 3 ran for every sample; see cell2location_visium's output{} note.
		Array[File] visium_deconvoluted_adata_objects = cell2location_visium.deconvoluted_adata_object
		Array[File] visium_deconvolution_metadata_csvs = cell2location_visium.deconvolution_metadata_csv
		Array[Array[File]] visium_inf_aver_csv_gzs = cell2location_visium.inf_aver_csv_gz
	}

	parameter_meta {
		# Inputs
		sc_rnaseq_reference_adata: {help: "scRNA-seq reference AnnData (.h5ad) used as the Cell2location deconvolution reference; required when run_stages includes Stage 1"}
		sc_rnaseq_reference_cell2location_model: {help: "Pre-computed inf_aver.csv.gz from a previous run; required when run_stages includes Stage 2 without Stage 1"}
		visium_samples: {help: "Flat list of Visium per-sample adatas (raw UMI counts, post-QC) to deconvolute independently"}
		batch_key: {help: "Column in reference .obs used as batch covariate in the NB regression model"}
		cell_type_key: {help: "Column in reference .obs to use as the cell type label"}
		n_cells_per_spot: {help: "Cell2location N_cells_per_location prior; 8 for brain Visium"}
		detection_alpha: {help: "Cell2location detection rate prior (20–200); higher values assume low technical variability"}
		max_epochs_ref: {help: "NB regression (Stage 1) training epochs"}
		max_epochs_spatial: {help: "Cell2location spatial mapping (Stage 2) training epochs"}
		n_comps: {help: "Number of PCA components for downstream spot clustering"}
		n_neighbors: {help: "Neighbourhood size for downstream spot clustering"}
		leiden_res: {help: "Leiden resolutions controlling the coarseness of clustering"}
		run_stages: {help: "Comma-separated subset of {1,2,3} — which stages to run: 1 = NB regression (reference signatures), 2 = Cell2location spatial mapping, 3 = downstream Leiden clustering. Default '1,2,3' (full pipeline); e.g. '1' builds the reference only, '2,3' loads a pre-computed reference and deconvolutes+clusters, '3' loads a pre-computed deconvolution and clusters only"}
		preemptible: {help: "Number of preemptible-VM attempts before Cromwell falls back to on-demand; 0 forces on-demand from the first attempt. Default 3; use 0 for long production runs where a preemption-triggered restart-from-scratch would be too costly"}
		raw_data_path: {help: "GCS path to raw data bucket; intermediate outputs are stored here"}
		billing_project: {help: "GCP billing project for requester-pays bucket access"}
		container_registry: {help: "Container registry hosting workflow Docker images"}
		zones: {help: "GCP zones where compute will run"}

		# Outputs
		visium_deconvoluted_adata_objects: {help: "Cell2location deconvoluted AnnData objects (.h5ad), one per Visium sample"}
		visium_deconvolution_metadata_csvs: {help: "Per-spot deconvolution metadata CSV (cell type fractions + Leiden clusters), one per Visium sample"}
		visium_inf_aver_csv_gzs: {help: "Cell2location reference signatures (<ref_stem>_inf_aver.csv.gz), one per sample that ran Stage 1; empty array per sample when run_stages excludes Stage 1"}
	}
}

task cell2location_visium {
	meta {
		description: "Deconvolute one Visium sample against the scRNA-seq cohort reference using Cell2location; run_stages selects which of the 3 stages (reference build, deconvolution, clustering) actually run"
	}

	input {
		String sample_id
		File? spatial_adata
		File? sc_rnaseq_reference_adata
		File? sc_rnaseq_reference_cell2location_model
		File? spatial_deconvolution_adata

		String batch_key
		String cell_type_key
		Int n_cells_per_spot
		Int detection_alpha
		Int max_epochs_ref
		Int max_epochs_spatial
		Int n_comps
		Int n_neighbors
		Array[Float] leiden_res
		String run_stages
		Int preemptible

		String raw_data_path #!UnusedDeclaration -- reserved for planned output-upload wiring, not yet implemented
		String billing_project #!UnusedDeclaration -- reserved for planned output-upload wiring, not yet implemented
		String container_registry
		String zones
	}

	Int disk_size = ceil(size([spatial_adata, sc_rnaseq_reference_adata, sc_rnaseq_reference_cell2location_model, spatial_deconvolution_adata], "GB") + 50)
	# Fix expected sample_id in outfile name.
	String corrected_sample_id = sub(sample_id, "^ASAP_PMBDS_", "ASAP_PMDBS_")

	# run_stages (e.g. "2,3"), always single comma-separated digits from {1,2,3}, never multi-digit.
	Boolean include_stage1 = sub(run_stages, "1", "x") != run_stages
	Boolean include_stage2 = sub(run_stages, "2", "x") != run_stages

	# Gate 1 - scRNA-seq reference: either H5AD or precomputed inf_aver.
	String gate_1_flag = if include_stage1 && defined(sc_rnaseq_reference_adata)
		then "--sc-rnaseq-reference-adata " + select_first([sc_rnaseq_reference_adata])
		else if !include_stage1 && include_stage2 && defined(sc_rnaseq_reference_cell2location_model)
		then "--sc-rnaseq-reference-cell2location-model " + select_first([sc_rnaseq_reference_cell2location_model])
		else ""

	# Gate 2 - Raw spatial H5AD
	String gate_2_flag = if include_stage2 && defined(spatial_adata)
		then "--spatial-h5ad " + select_first([spatial_adata])
		else ""

	# Gate 3 - Deconvoluted spatial H5AD. Forwarded purely on defined() -- cell2location_deconv.py's
	# _validate_run_stages() is the single source of truth for whether that's actually valid given
	# run_stages, not this WDL.
	String gate_3_flag = if defined(spatial_deconvolution_adata)
		then "--spatial-deconvolution-adata " + select_first([spatial_deconvolution_adata])
		else ""

	command <<<
		set -euo pipefail

		/usr/bin/time -v python3 /opt/scripts/cell2location_deconv.py \
			--run-stages ~{run_stages} \
			~{gate_1_flag} \
			~{gate_2_flag} \
			~{gate_3_flag} \
			--sample-id ~{sample_id} \
			--batch-key ~{batch_key} \
			--cell-type-key ~{cell_type_key} \
			--n-cells-per-spot ~{n_cells_per_spot} \
			--detection-alpha ~{detection_alpha} \
			--max-epochs-ref ~{max_epochs_ref} \
			--max-epochs-spatial ~{max_epochs_spatial} \
			--n-comps ~{n_comps} \
			--n-neighbors ~{n_neighbors} \
			--leiden-res ~{sep=" " leiden_res} \
			--output-dir . \
			--promote
	>>>

	output {
		# NOTE: never written when run_stages is "1" only
		File deconvoluted_adata_object = "~{corrected_sample_id}.deconvoluted.h5ad"
		File deconvolution_metadata_csv = "~{corrected_sample_id}.deconvolution_metadata.csv.gz"
		# Non-empty only when run_stages includes Stage 1; empty when a pre-computed model was used.
		Array[File] inf_aver_csv_gz = glob("*_inf_aver.csv.gz")
	}

	runtime {
		docker: "~{container_registry}/cell2location_tools:1.0.0"
		cpu: 8
		memory: "52 GB"
		disks: "local-disk " + disk_size + " SSD"
		preemptible: preemptible
		bootDiskSizeGb: 40
		gpuType: "nvidia-tesla-t4"
		gpuCount: 1
		zones: zones
	}
}
