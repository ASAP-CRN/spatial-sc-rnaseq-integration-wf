version 1.0

# Spatial + scRNA-seq integration workflow entrypoint.
# Inputs are processed outputs from sc-rnaseq-wf and spatial-transcriptomics-wf;
# this workflow does not handle raw FASTQ or SpaceRanger/CellRanger execution.
#
# The scRNA-seq cohort adata is used as a reference for all spatial samples.
# Each spatial sample carries its own per-sample adata (raw UMI counts, pre-Harmony);
# samples from all projects are flattened and integrated in parallel against the reference.

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
		# Each project carries per-sample adatas (raw UMI counts, pre-Harmony);
		# samples are flattened across projects and integrated independently.
		Array[SpatialVisiumProject] spatial_visium_projects = []
		Array[SpatialGeomxProject]  spatial_geomx_projects  = []

		# Integration parameters
		String batch_key            = "batch_id"
		Int n_top_genes             = 3000
		Int n_comps                 = 30
		Int n_neighbors             = 15
		Array[Float] leiden_res     = [0.05, 0.1, 0.2, 0.4]

		String cohort_raw_data_bucket_spatial_sc_rnaseq
		Array[String] cohort_staging_data_buckets_spatial_sc_rnaseq

		String container_registry
		String zones = "us-central1-c us-central1-f"
	}

	String workflow_version = "v1.0.0"
	String workflow_release = "https://github.com/ASAP-CRN/spatial-sc-rnaseq-integration-wf/releases/tag/spatial_sc_rnaseq_integration-~{workflow_version}"

	call GetWorkflowMetadata.get_workflow_metadata {
		input:
			workflow_name    = "spatial_sc_rnaseq_integration",
			workflow_version = workflow_version,
			workflow_release = workflow_release
	}

	# scRNA-seq reference: use first project (expected to be the cross-team cohort)
	File sc_rnaseq_reference_adata        = sc_rnaseq_projects[0].final_adata_object
	File sc_rnaseq_reference_metadata_csv = sc_rnaseq_projects[0].final_metadata_csv
	File sc_rnaseq_cell_types_parquet     = sc_rnaseq_projects[0].scanvi_cell_types_parquet

	# Flatten per-sample data across all projects
	Array[SpatialVisiumSample] visium_samples = flatten(spatial_visium_projects[*].samples)
	Array[SpatialGeomxSample]  geomx_samples  = flatten(spatial_geomx_projects[*].samples)

	call Integration.integration {
		input:
			cohort_id                        = cohort_id,
			sc_rnaseq_reference_adata        = sc_rnaseq_reference_adata,
			sc_rnaseq_reference_metadata_csv = sc_rnaseq_reference_metadata_csv,
			sc_rnaseq_cell_types_parquet     = sc_rnaseq_cell_types_parquet,
			visium_samples                   = visium_samples,
			geomx_samples                    = geomx_samples,
			batch_key                        = batch_key,
			n_top_genes                      = n_top_genes,
			n_comps                          = n_comps,
			n_neighbors                      = n_neighbors,
			leiden_res                       = leiden_res,
			raw_data_path                    = cohort_raw_data_bucket_spatial_sc_rnaseq,
			workflow_info                    = get_workflow_metadata.workflow_info,
			billing_project                  = "dnastack-asap-parkinsons",
			container_registry               = container_registry,
			zones                            = zones
	}

	output {
		# Array[File] — flat over all samples (across all projects)
		Array[File] visium_integrated_adata_objects  = Integration.integration.visium_integrated_adata_objects
		Array[File] visium_integration_metadata_csvs = Integration.integration.visium_integration_metadata_csvs
		Array[File] geomx_integrated_adata_objects   = Integration.integration.geomx_integrated_adata_objects
		Array[File] geomx_integration_metadata_csvs  = Integration.integration.geomx_integration_metadata_csvs
	}
}
