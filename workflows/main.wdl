version 1.0

# Spatial + scRNA-seq integration workflow entrypoint.
# Inputs are processed outputs from sc-rnaseq-wf and spatial-transcriptomics-wf;
# this workflow does not handle raw FASTQ or SpaceRanger/CellRanger execution.
#
# The scRNA-seq cohort adata is used as a reference for all spatial samples.
# Each spatial sample carries its own per-sample adata (raw UMI counts, post-QC);
# samples from all projects are flattened and deconvoluted in parallel against the reference.
#
# Deconvolution method: Cell2location (Bayesian hierarchical model, GPU-accelerated).
# RCTD support is planned for a future version.

import "structs.wdl"
import "../wf-common/wdl/tasks/get_workflow_metadata.wdl" as GetWorkflowMetadata
import "../wf-common/wdl/tasks/upload_final_outputs.wdl" as UploadFinalOutputs
import "integration/integration.wdl" as Integration

workflow spatial_sc_rnaseq_integration {
	input {
		String cohort_id

		# Processed scRNA-seq cohort outputs (from sc-rnaseq-wf cohort_analysis).
		# Typically one entry representing the cross-team cohort; used as the reference
		# for all spatial samples.
		Array[ScRnaseqProject] sc_rnaseq_projects

		# Processed spatial outputs (from spatial-transcriptomics-wf).
		# Each project carries per-sample adatas (raw UMI counts, post-QC);
		# samples are flattened across projects and deconvoluted independently.
		Array[SpatialVisiumProject] spatial_visium_projects = []

		# Optional pre-computed Cell2location reference signatures (inf_aver.csv.gz).
		# When provided, Stage 1 NB regression is skipped for all samples.
		# When omitted, Stage 1 is trained from sc_rnaseq_projects[0].final_adata_object.
		File? sc_rnaseq_reference_cell2location_model

		# Cell2location deconvolution parameters
		String batch_key        = "batch_id"
		String cell_type_key    = "cell_type_pred"
		Int n_cells_per_spot    = 8
		Int detection_alpha     = 100
		Int max_epochs_ref      = 250
		Int max_epochs_spatial  = 30000

		# Downstream clustering parameters
		Int n_comps             = 30
		Int n_neighbors         = 15
		Array[Float] leiden_res = [0.05, 0.1, 0.2, 0.4]
		Boolean skip_clustering = false

		String cohort_raw_data_bucket_spatial_sc_rnaseq

		String container_registry
		String zones = "us-central1-c us-central1-f"
	}

	String workflow_version = "v1.0.0"
	String workflow_release = "https://github.com/ASAP-CRN/spatial-sc-rnaseq-integration-wf/releases/tag/spatial_sc_rnaseq_integration-~{workflow_version}"

	call GetWorkflowMetadata.get_workflow_metadata {
		input:
			zones = zones
	}

	# scRNA-seq reference: use first project (expected to be the cross-team cohort)
	File sc_rnaseq_reference_adata = sc_rnaseq_projects[0].final_adata_object

	# Flatten per-sample data across all projects (WDL 1.0: scatter then flatten)
	scatter (project in spatial_visium_projects) {
		Array[SpatialVisiumSample] project_samples = project.samples
	}
	Array[SpatialVisiumSample] visium_samples = flatten(project_samples)

	call Integration.integration {
		input:
			cohort_id                                = cohort_id,
			sc_rnaseq_reference_adata                = sc_rnaseq_reference_adata,
			sc_rnaseq_reference_cell2location_model  = sc_rnaseq_reference_cell2location_model,
			visium_samples                           = visium_samples,
			batch_key                                = batch_key,
			cell_type_key                            = cell_type_key,
			n_cells_per_spot                         = n_cells_per_spot,
			detection_alpha                          = detection_alpha,
			max_epochs_ref                           = max_epochs_ref,
			max_epochs_spatial                       = max_epochs_spatial,
			n_comps                                  = n_comps,
			n_neighbors                              = n_neighbors,
			leiden_res                               = leiden_res,
			skip_clustering                          = skip_clustering,
			raw_data_path                            = cohort_raw_data_bucket_spatial_sc_rnaseq,
			billing_project                          = get_workflow_metadata.billing_project,
			container_registry                       = container_registry,
			zones                                    = zones
	}

	output {
		# Array[File] — flat over all samples (across all projects)
		Array[File]  visium_deconvoluted_adata_objects   = integration.visium_deconvoluted_adata_objects
		Array[File]  visium_deconvolution_metadata_csvs  = integration.visium_deconvolution_metadata_csvs
		Array[Array[File]] visium_inf_aver_csv_gzs        = integration.visium_inf_aver_csv_gzs
	}
}
