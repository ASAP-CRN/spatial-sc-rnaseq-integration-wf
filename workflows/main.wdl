version 1.0

# Spatial deconvolution using scRNA-seq data workflow entrypoint.
# Inputs are processed outputs from sc-rnaseq-wf and spatial-transcriptomics-wf;
# this workflow does not handle raw FASTQ or SpaceRanger/CellRanger execution.
#
# The scRNA-seq reference can obtained from a downsampled version of the sc-rnaseq-wf cohort
# alone or supplemented with additional data (e.g. Siletti and/or Kamath cells).
# See `scripts/build_deconv_reference_h5ad.py` for details.
#
# Each spatial sample carries its own per-sample H5AD (raw UMI counts, all genes--not
# only HVGs, i.e. cleaned_unfiltered) from spatial-transcriptomics-wf.
#
# Stage 1: Build the scRNA-seq reference inf_aver file from the H5AD for Cell2location (once per cohort)
# Stage 2: Deconvolute each spatial sample using the reference inf_aver file with Cell2location (per sample)
# Stage 3: Leiden-cluster each spatial sample on its deconvolved cell type fractions (not on gene expression).
#
# Deconvolution method: Cell2location (Bayesian hierarchical model, GPU-accelerated).

import "structs.wdl"
import "../wf-common/wdl/tasks/get_workflow_metadata.wdl" as GetWorkflowMetadata
import "../wf-common/wdl/tasks/upload_final_outputs.wdl" as UploadFinalOutputs
import "deconvolution/deconvolution.wdl" as Deconvolution

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
		# Required when run_stages includes Stage 1
		File? sc_rnaseq_reference_adata

		# spatial-transcriptomics-wf cleaned_unfiltered spatial H5ADs (raw UMI counts, all genes--not only HVGs)
		Array[SpatialVisiumProject] spatial_visium_projects = [] #!UnverifiedStruct -- empty-array default can't be statically verified against the struct type; known miniwdl limitation

		# Pre-computed inf_aver.csv.gz deconvolution reference
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

		# Stages to run; see Deconvolution.deconvolution
		String run_stages = "1,2,3"
		Int preemptible = 3

		String cohort_raw_data_bucket_spatial_sc_rnaseq

		# Promotion to staging buckets (e.g. DEV/UAT); see parameter_meta below.
		# Skipped entirely when staging_data_buckets is empty.
		Array[String] staging_data_buckets = []
		String crn_release_version = ""

		String container_registry
		String zones = "us-central1-c us-central1-f"
	}

	# Owned by this repo's release process; bump alongside a GitHub release, not at submission time.
	String workflow_name = "spatial_sc_rnaseq_integration"
	String workflow_version = "v1.0.0"
	String workflow_release = "https://github.com/ASAP-CRN/spatial-sc-rnaseq-integration-wf/releases/tag/~{workflow_name}-~{workflow_version}"

	call GetWorkflowMetadata.get_workflow_metadata {
		input:
			zones = zones
	}

	Array[Array[String]] workflow_info = [[get_workflow_metadata.timestamp, workflow_name, workflow_version, workflow_release]]

	# Per the output bucket's team-agreed naming convention, deconvolution outputs are grouped under the
	# existing spatial_visium domain (alongside spatial-transcriptomics-wf's own preprocess outputs),
	# not under this repo's own spatial_sc_rnaseq_integration workflow name.
	String raw_data_path_prefix = "~{cohort_raw_data_bucket_spatial_sc_rnaseq}/workflow_execution/spatial_visium/deconvolution_analysis"

	# Flatten per-sample data across all projects (WDL 1.0: scatter then flatten)
	scatter (project in spatial_visium_projects) {
		Array[SpatialVisiumSample] project_samples = project.samples
	}
	Array[SpatialVisiumSample] visium_samples = flatten(project_samples)

	call validate_staging_combination {
		input:
			run_stages = run_stages,
			staging_data_buckets = staging_data_buckets,
			container_registry = container_registry,
			zones = zones
	}

	call Deconvolution.deconvolution {
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
			# Reading validate_staging_combination.validated_run_stages.
			run_stages = validate_staging_combination.validated_run_stages,
			preemptible = preemptible,
			raw_data_path_prefix = raw_data_path_prefix,
			workflow_info = workflow_info,
			billing_project = get_workflow_metadata.billing_project,
			container_registry = container_registry,
			zones = zones
	}

	# Promote already-completed per-sample deconvolution outputs to one or more staging
	# buckets (e.g. DEV/UAT). Only Stage 3's clustered h5ad is staged, not Stage 2's.
	if (length(staging_data_buckets) > 0) {
		String deconvolution_analysis_release_path = "spatial_visium/release/~{crn_release_version}/deconvolution_analysis"
		String promotion_cluster_visium_path = "~{raw_data_path_prefix}/cluster_visium/~{deconvolution.cell2location_version_value}"

		# corrected_sample_id comes from deconvolution.visium_corrected_sample_ids (single
		# source of truth for the PMBDS->PMDBS fix, see deconvolution.wdl) instead of
		# re-deriving it here.
		scatter (corrected_sample_id in deconvolution.visium_corrected_sample_ids) {
			String deconvoluted_clustered_adata_path = "~{promotion_cluster_visium_path}/~{corrected_sample_id}.deconvoluted_clustered.h5ad"
			String deconvolution_clustering_metadata_csv_path = "~{promotion_cluster_visium_path}/~{corrected_sample_id}.deconvolution_clustering_metadata.csv.gz"
			String deconvolution_clustering_timing_tsv_path = "~{promotion_cluster_visium_path}/~{corrected_sample_id}.timing.tsv"
		}

		call UploadFinalOutputs.upload_final_outputs as upload_deconvolution_outputs {
			input:
				output_file_paths = flatten([
					deconvoluted_clustered_adata_path,
					deconvolution_clustering_metadata_csv_path,
					deconvolution_clustering_timing_tsv_path
				]),
				staging_data_buckets = staging_data_buckets,
				staging_data_path = deconvolution_analysis_release_path,
				billing_project = get_workflow_metadata.billing_project,
				zones = zones
		}
	}

	output {
		# Flat over all samples; null when run_stages excludes Stage 2.
		Array[String?] visium_deconvoluted_adata_objects = deconvolution.visium_deconvoluted_adata_objects
		Array[String?] visium_deconvolution_metadata_csvs = deconvolution.visium_deconvolution_metadata_csvs
		Array[String?] visium_deconvolution_training_history_pngs = deconvolution.visium_deconvolution_training_history_pngs
		Array[String?] visium_deconvolution_timing_tsvs = deconvolution.visium_deconvolution_timing_tsvs
		# Flat over all samples; null when run_stages excludes Stage 3.
		Array[String?] visium_deconvoluted_clustered_adata_objects = deconvolution.visium_deconvoluted_clustered_adata_objects
		Array[String?] visium_deconvolution_clustering_metadata_csvs = deconvolution.visium_deconvolution_clustering_metadata_csvs
		Array[String?] visium_deconvolution_clustering_timing_tsvs = deconvolution.visium_deconvolution_clustering_timing_tsvs
		# Built once; null when run_stages excludes Stage 1.
		String? sc_rnaseq_reference_inf_aver_csv_gz = deconvolution.sc_rnaseq_reference_inf_aver_csv_gz
		String? sc_rnaseq_reference_training_history_png = deconvolution.sc_rnaseq_reference_training_history_png
		String? sc_rnaseq_reference_build_timing_tsv = deconvolution.sc_rnaseq_reference_build_timing_tsv
		# Null when staging_data_buckets is empty (promotion skipped).
		Array[String]? deconvolution_promotion_manifests = upload_deconvolution_outputs.manifests
	}

	parameter_meta {
		# Inputs
		sc_rnaseq_reference_adata: {help: "scRNA-seq reference AnnData (.h5ad) used to build the Cell2location inf_aver.csv.gz deconvolution reference; required when run_stages includes Stage 1."}
		spatial_visium_projects: {help: "Processed spatial outputs (from spatial-transcriptomics-wf); samples are flattened across projects and deconvoluted independently. []"}
		sc_rnaseq_reference_cell2location_model: {help: "Pre-computed inf_aver.csv.gz deconvolution reference; required when run_stages includes Stage 2 without Stage 1."}
		batch_key: {help: "Column in reference .obs used as batch covariate in the NB regression model. ['batch_id']"}
		cell_type_key: {help: "Column in reference .obs to use as the cell type label. ['cell_type_pred']"}
		n_cells_per_spot: {help: "Cell2location N_cells_per_location prior; 8 for brain Visium. [8]"}
		detection_alpha: {help: "Cell2location detection rate prior (20–200); higher values assume low technical variability. [100]"}
		max_epochs_ref: {help: "NB regression (Stage 1) training epochs. [250]"}
		max_epochs_spatial: {help: "Cell2location spatial mapping (Stage 2) training epochs. [30000]"}
		n_comps: {help: "PCA components for clustering on deconvolved cell type fractions (Stage 3). [30]"}
		n_neighbors: {help: "Neighbourhood size for clustering on deconvolved cell type fractions (Stage 3). [15]"}
		leiden_res: {help: "Leiden resolutions for clustering on deconvolved cell type fractions (Stage 3); one clustering per resolution. [0.2, 0.5, 1.0]"}
		run_stages: {help: "Comma-separated subset of {1,2,3} — which stages to run: 1 = build Cell2location inf_aver.csv.gz deconvolution reference, 2 = Cell2location spatial mapping (per sample), 3 = Leiden clustering on deconvolved cell type fractions. Stage 3 requires Stage 2 in the same submission -- standalone \"3\" or \"1,3\" are not supported. ['1,2,3']"}
		preemptible: {help: "Number of preemptible-VM attempts before Cromwell falls back to on-demand; 0 forces on-demand from the first attempt. [3]"}
		cohort_raw_data_bucket_spatial_sc_rnaseq: {help: "GCS path to the raw data bucket root; task outputs are uploaded under `<this>/workflow_execution/spatial_visium/deconvolution_analysis/<task_name>/<cell2location_version>`."}
		staging_data_buckets: {help: "Set of buckets to additionally stage completed per-sample deconvolution outputs to (e.g. DEV and/or UAT). Requires run_stages to be either empty (promote a previous run's already-computed outputs) or include Stage 3 (this submission produces what it promotes) -- any other run_stages value fails validation. Promotion is skipped entirely when empty. []"}
		crn_release_version: {help: "Release version (e.g. 'v5.0.0') to stage outputs under; required when staging_data_buckets is non-empty."}
		container_registry: {help: "Container registry hosting workflow Docker images."}
		zones: {help: "GCP zones where compute will run. ['us-central1-c us-central1-f']"}

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
		deconvolution_promotion_manifests: {help: "GCS paths to the MANIFEST.tsv uploaded to each staging bucket, one per bucket in staging_data_buckets; null when staging_data_buckets is empty."}
	}
}

task validate_staging_combination {
	meta {
		description: "Fail fast if staging_data_buckets is set but run_stages won't produce (or reuse) Stage 3's clustered outputs"
	}

	input {
		String run_stages
		Array[String] staging_data_buckets

		String container_registry
		String zones
	}

	command <<<
		set -euo pipefail

		run_stages="~{run_stages}"
		n_staging_buckets=~{length(staging_data_buckets)}

		if [[ ${n_staging_buckets} -gt 0 && -n "${run_stages}" ]]; then
			case ",${run_stages}," in
				*,3,*)
					;;
				*)
					echo "[ERROR] staging_data_buckets is set but run_stages=\"${run_stages}\" does not include Stage 3 and is not empty." >&2
					echo "[ERROR] Promotion needs Stage 3's clustered outputs -- either set run_stages=\"\" to promote a previous run's outputs, or include Stage 3 in run_stages." >&2
					exit 1
					;;
			esac
		fi
	>>>

	output {
		String validated_run_stages = run_stages
	}

	runtime {
		docker: "~{container_registry}/util:1.2.0"
		cpu: 1
		memory: "1 GB"
		zones: zones
	}

	parameter_meta {
		run_stages: {help: "Same run_stages passed to Deconvolution.deconvolution; validated here against staging_data_buckets."}
		staging_data_buckets: {help: "Same staging_data_buckets passed to upload_final_outputs."}
		container_registry: {help: "Container registry hosting workflow Docker images."}
		zones: {help: "GCP zones where compute will run."}
	}
}
