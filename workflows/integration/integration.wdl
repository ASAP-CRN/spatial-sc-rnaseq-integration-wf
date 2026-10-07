version 1.0

import "../structs.wdl"

# Spatial + scRNA-seq integration sub-workflow.
#
# Deconvolution method: Cell2location (Bayesian hierarchical model, GPU-accelerated).
# RCTD (statistical mixture model, R) is planned for a future version.
#
# For each Visium sample the workflow runs Cell2location against the scRNA-seq cohort
# reference. Samples from all projects are flattened upstream and scattered flat here.
#
# scRNA-seq inputs:  sc-rnaseq-wf cohort_analysis (reference for all spatial samples)
# Spatial inputs:    per-sample post-QC adatas (raw UMI counts) from spatial-transcriptomics-wf
#                    - Visium: spatial_visium preprocess.qc_adata_objects

workflow integration {
	meta {
		description: "Deconvolute per-sample Visium adatas against an scRNA-seq cohort reference using Cell2location"
		version: "v1.0.0"
		author: "DNAstack Bioinformatics Team"
		email: "bioinformatics@dnastack.com"
	}

	input {
		String cohort_id

		# scRNA-seq cohort reference (used for all spatial samples)
		File   sc_rnaseq_reference_adata

		# Optional pre-computed Cell2location reference signatures (inf_aver.csv.gz).
		# When provided, Stage 1 NB regression is skipped. When omitted, Stage 1 is
		# trained from sc_rnaseq_reference_adata (and inf_aver.csv.gz is saved per sample).
		# TODO: add a dedicated build_reference_signatures task so Stage 1 runs once
		#       rather than being re-trained in each scatter when no model is pre-provided.
		File?  sc_rnaseq_reference_cell2location_model

		# Flat list of per-sample adatas (flattened across projects upstream in main.wdl)
		Array[SpatialVisiumSample] visium_samples

		# Cell2location deconvolution parameters
		String batch_key
		String cell_type_key
		Int    n_cells_per_spot
		Int    detection_alpha
		Int    max_epochs_ref
		Int    max_epochs_spatial

		# Downstream clustering parameters
		Int           n_comps
		Int           n_neighbors
		Array[Float]  leiden_res
		Boolean       skip_clustering = false
		Int           preemptible = 3

		String               raw_data_path
		String               billing_project
		String               container_registry
		String               zones
	}

	# Visium: scatter flat over all samples (across all projects)
	scatter (sample in visium_samples) {
		call cell2location_visium {
			input:
				sample_id                               = sample.sample_id,
				spatial_adata                           = sample.qc_adata_object,
				sc_rnaseq_reference_adata               = sc_rnaseq_reference_adata,
				sc_rnaseq_reference_cell2location_model = sc_rnaseq_reference_cell2location_model,
				batch_key                               = batch_key,
				cell_type_key                           = cell_type_key,
				n_cells_per_spot                        = n_cells_per_spot,
				detection_alpha                         = detection_alpha,
				max_epochs_ref                          = max_epochs_ref,
				max_epochs_spatial                      = max_epochs_spatial,
				n_comps                                 = n_comps,
				n_neighbors                             = n_neighbors,
				leiden_res                              = leiden_res,
				skip_clustering                         = skip_clustering,
				preemptible                             = preemptible,
				raw_data_path                           = raw_data_path,
				billing_project                         = billing_project,
				container_registry                      = container_registry,
				zones                                   = zones
		}
	}

	output {
		# Array[File] — flat over all samples (across all projects)
		Array[File]  visium_deconvoluted_adata_objects  = cell2location_visium.deconvoluted_adata_object
		Array[File]  visium_deconvolution_metadata_csvs = cell2location_visium.deconvolution_metadata_csv
		Array[Array[File]] visium_inf_aver_csv_gzs       = cell2location_visium.inf_aver_csv_gz
	}

	parameter_meta {
		# Inputs
		cohort_id:                               {help: "Cohort identifier; used to name output files"}
		sc_rnaseq_reference_adata:               {help: "cohort_analysis.final_adata_object (.h5ad) — Harmony-integrated, cell-type-annotated AnnData used as the deconvolution reference"}
		sc_rnaseq_reference_cell2location_model: {help: "Pre-computed inf_aver.csv.gz from a previous run; skips Stage 1 NB regression when provided"}
		visium_samples:                          {help: "Flat list of Visium per-sample adatas (raw UMI counts, post-QC) to deconvolute independently"}
		batch_key:                               {help: "Column in reference .obs used as batch covariate in the NB regression model"}
		cell_type_key:                           {help: "Column in reference .obs to use as the cell type label"}
		n_cells_per_spot:                        {help: "Cell2location N_cells_per_location prior; 8 for brain Visium"}
		detection_alpha:                         {help: "Cell2location detection rate prior (20–200); higher values assume low technical variability"}
		max_epochs_ref:                          {help: "NB regression (Stage 1) training epochs"}
		max_epochs_spatial:                      {help: "Cell2location spatial mapping (Stage 2) training epochs"}
		n_comps:                                 {help: "Number of PCA components for downstream spot clustering"}
		n_neighbors:                             {help: "Neighbourhood size for downstream spot clustering"}
		leiden_res:                              {help: "Leiden resolutions controlling the coarseness of clustering"}
		skip_clustering:                         {help: "Skip Stage 3 downstream Leiden clustering; Stage 1/2 outputs are still written"}
		preemptible:                             {help: "Number of preemptible-VM attempts before Cromwell falls back to on-demand; 0 forces on-demand from the first attempt. Default 3; use 0 for long production runs where a preemption-triggered restart-from-scratch would be too costly"}
		raw_data_path:                           {help: "GCS path to raw data bucket; intermediate outputs are stored here"}
		billing_project:                         {help: "GCP billing project for requester-pays bucket access"}
		container_registry:                      {help: "Container registry hosting workflow Docker images"}
		zones:                                   {help: "GCP zones where compute will run"}

		# Outputs
		visium_deconvoluted_adata_objects:  {help: "Cell2location deconvoluted AnnData objects (.h5ad), one per Visium sample"}
		visium_deconvolution_metadata_csvs: {help: "Per-spot deconvolution metadata CSV (cell type fractions + Leiden clusters), one per Visium sample"}
		visium_inf_aver_csv_gzs:            {help: "Cell2location reference signatures (<ref_stem>_inf_aver.csv.gz), one per sample that ran Stage 1; empty array per sample when pre-computed model was supplied"}
	}
}

task cell2location_visium {
	meta {
		description: "Deconvolute one Visium sample against the scRNA-seq cohort reference using Cell2location"
	}

	input {
		String sample_id
		File   spatial_adata              # per-sample Visium adata (.h5ad) — raw UMI counts, post-QC

		# Exactly one of the two reference inputs must be provided (validated at runtime):
		#   sc_rnaseq_reference_adata               → trains Stage 1 NB regression, saves inf_aver.csv.gz
		#   sc_rnaseq_reference_cell2location_model → loads pre-computed inf_aver.csv.gz, skips Stage 1
		File?  sc_rnaseq_reference_adata
		File?  sc_rnaseq_reference_cell2location_model

		String batch_key
		String cell_type_key
		Int    n_cells_per_spot
		Int    detection_alpha
		Int    max_epochs_ref
		Int    max_epochs_spatial
		Int    n_comps
		Int    n_neighbors
		Array[Float] leiden_res
		Boolean skip_clustering = false
		Int preemptible = 3

		String               raw_data_path
		String               billing_project
		String               container_registry
		String               zones
	}

	Int disk_size = ceil(size(spatial_adata, "GB") + size(select_first([sc_rnaseq_reference_adata, sc_rnaseq_reference_cell2location_model]), "GB") + 50)
	# Fix expected sample_id in outfile name.
	String corrected_sample_id = sub(sample_id, "^ASAP_PMBDS_", "ASAP_PMDBS_")

	command <<<
		set -euo pipefail

		# utils/ lives alongside the scripts in the container at /opt/scripts/utils/
		# (requires: COPY utils /opt/scripts/utils in the Dockerfile build context)
		export PYTHONPATH=/opt/scripts

		/usr/bin/time -v python3 /opt/scripts/cell2location_deconvolute.py \
			--sample-id ~{sample_id} \
			--spatial-h5ad ~{spatial_adata} \
			~{if defined(sc_rnaseq_reference_cell2location_model)
			  then "--sc-rnaseq-reference-cell2location-model " + select_first([sc_rnaseq_reference_cell2location_model])
			  else "--sc-rnaseq-reference-adata " + select_first([sc_rnaseq_reference_adata])} \
			--batch-key ~{batch_key} \
			--cell-type-key ~{cell_type_key} \
			--n-cells-per-spot ~{n_cells_per_spot} \
			--detection-alpha ~{detection_alpha} \
			--max-epochs-ref ~{max_epochs_ref} \
			--max-epochs-spatial ~{max_epochs_spatial} \
			--n-comps ~{n_comps} \
			--n-neighbors ~{n_neighbors} \
			--leiden-res ~{sep=" " leiden_res} \
			~{if skip_clustering then "--skip-clustering" else ""} \
			--output-dir . \
			--promote
	>>>

	output {
		File  deconvoluted_adata_object  = "~{corrected_sample_id}.deconvoluted.h5ad"
		File  deconvolution_metadata_csv = "~{corrected_sample_id}.deconvolution_metadata.csv.gz"
		# Non-empty only when Stage 1 ran (sc_rnaseq_reference_adata provided); empty in Mode B.
		Array[File] inf_aver_csv_gz      = glob("*_inf_aver.csv.gz")
	}

	runtime {
		docker: "~{container_registry}/cell2location_tools:1.0.0"
		cpu: 8
		memory: "52 GB"
		disks: "local-disk " + disk_size + " SSD"
		bootDiskSizeGb: 40
		gpuType: "nvidia-tesla-t4"
		gpuCount: 1
		nvidiaDriverVersion: "545.23.08" #!UnknownRuntimeKey
		preemptible: preemptible
		zones: zones
	}
}
