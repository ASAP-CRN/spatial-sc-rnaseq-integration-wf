version 1.0

import "../structs.wdl"

# Spatial + scRNA-seq integration sub-workflow.
#
# For each spatial sample the workflow integrates the per-sample spatial adata
# (raw UMI counts, pre-Harmony) against the scRNA-seq cohort reference.
# Samples from all projects are flattened upstream and scattered flat here.
#
# scRNA-seq inputs:  sc-rnaseq-wf cohort_analysis (reference for all spatial samples)
# Spatial inputs:    per-sample pre-QC adatas (raw UMI counts) from spatial-transcriptomics-wf
#                    - Visium: spatial_visium preprocess.qc_adata_objects
#                    - GeoMx:  spatial_geomx process_to_adata.processed_adata_objects

workflow integration {
	meta {
		description: "Integrate per-sample spatial adatas against an scRNA-seq cohort reference; Visium and GeoMx samples are scattered in parallel"
		version: "v1.0.0"
		author: "DNAstack Bioinformatics Team"
		email: "bioinformatics@dnastack.com"
	}

	input {
		String cohort_id

		# scRNA-seq cohort reference (used for all spatial samples)
		File sc_rnaseq_reference_adata
		File sc_rnaseq_reference_metadata_csv
		File sc_rnaseq_cell_types_parquet

		# Flat lists of per-sample adatas (flattened across projects upstream in main.wdl)
		Array[SpatialVisiumSample] visium_samples
		Array[SpatialGeomxSample]  geomx_samples

		String batch_key
		Int n_top_genes
		Int n_comps
		Int n_neighbors
		Array[Float] leiden_res

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	# Visium: scatter flat over all samples (across all projects)
	scatter (sample in visium_samples) {
		call integrate_sample as integrate_visium {
			input:
				sample_id                        = sample.sample_id,
				spatial_adata                    = sample.qc_adata_object,
				sc_rnaseq_reference_adata        = sc_rnaseq_reference_adata,
				sc_rnaseq_reference_metadata_csv = sc_rnaseq_reference_metadata_csv,
				sc_rnaseq_cell_types_parquet     = sc_rnaseq_cell_types_parquet,
				batch_key                        = batch_key,
				n_top_genes                      = n_top_genes,
				n_comps                          = n_comps,
				n_neighbors                      = n_neighbors,
				leiden_res                       = leiden_res,
				raw_data_path                    = raw_data_path,
				workflow_info                    = workflow_info,
				billing_project                  = billing_project,
				container_registry               = container_registry,
				zones                            = zones
		}
	}

	# GeoMx: same flat scatter pattern
	scatter (sample in geomx_samples) {
		call integrate_sample as integrate_geomx {
			input:
				sample_id                        = sample.sample_id,
				spatial_adata                    = sample.processed_adata_object,
				sc_rnaseq_reference_adata        = sc_rnaseq_reference_adata,
				sc_rnaseq_reference_metadata_csv = sc_rnaseq_reference_metadata_csv,
				sc_rnaseq_cell_types_parquet     = sc_rnaseq_cell_types_parquet,
				batch_key                        = batch_key,
				n_top_genes                      = n_top_genes,
				n_comps                          = n_comps,
				n_neighbors                      = n_neighbors,
				leiden_res                       = leiden_res,
				raw_data_path                    = raw_data_path,
				workflow_info                    = workflow_info,
				billing_project                  = billing_project,
				container_registry               = container_registry,
				zones                            = zones
		}
	}

	output {
		# Array[File] — flat over all samples (across all projects)
		Array[File] visium_integrated_adata_objects  = integrate_visium.integrated_adata_object
		Array[File] visium_integration_metadata_csvs = integrate_visium.integration_metadata_csv
		Array[File] geomx_integrated_adata_objects   = integrate_geomx.integrated_adata_object
		Array[File] geomx_integration_metadata_csvs  = integrate_geomx.integration_metadata_csv
	}

	parameter_meta {
		# Inputs
		cohort_id:                        {help: "Cohort identifier; used to name output files"}
		sc_rnaseq_reference_adata:        {help: "cohort_analysis.final_adata_object (.h5ad) — Harmony-integrated, cell-type-annotated AnnData used as the deconvolution reference"}
		sc_rnaseq_reference_metadata_csv: {help: "cohort_analysis.final_metadata_csv — per-cell metadata table"}
		sc_rnaseq_cell_types_parquet:     {help: "cohort_analysis.scanvi_cell_types_parquet — scANVI cell type predictions"}
		visium_samples:                   {help: "Flat list of Visium per-sample adatas (raw UMI counts, post-QC, pre-Harmony) to integrate independently"}
		geomx_samples:                    {help: "Flat list of GeoMx per-sample adatas (Q3-normalized, pre-Harmony) to integrate independently"}
		batch_key:                        {help: "Key in AnnData object for batch information"}
		n_top_genes:                      {help: "Number of highly variable genes to retain"}
		n_comps:                          {help: "Number of principal components to compute"}
		n_neighbors:                      {help: "Size of local neighbourhood used for manifold approximation"}
		leiden_res:                       {help: "Leiden resolutions controlling the coarseness of clustering"}
		raw_data_path:                    {help: "GCS path to raw data bucket; intermediate outputs are stored here"}
		workflow_info:                    {help: "Workflow name, version, and release metadata from get_workflow_metadata"}
		billing_project:                  {help: "GCP billing project for requester-pays bucket access"}
		container_registry:               {help: "Container registry hosting workflow Docker images"}
		zones:                            {help: "GCP zones where compute will run"}

		# Outputs
		visium_integrated_adata_objects:  {help: "Integrated AnnData objects (.h5ad), one per Visium sample"}
		visium_integration_metadata_csvs: {help: "Integration metadata CSV, one per Visium sample"}
		geomx_integrated_adata_objects:   {help: "Integrated AnnData objects (.h5ad), one per GeoMx sample"}
		geomx_integration_metadata_csvs:  {help: "Integration metadata CSV, one per GeoMx sample"}
	}
}

task integrate_sample {
	meta {
		description: "Integrate one spatial sample against the scRNA-seq cohort reference"
	}

	input {
		String sample_id
		File   spatial_adata                  # per-sample spatial adata (.h5ad) — raw UMI counts, pre-Harmony
		File   sc_rnaseq_reference_adata      # scRNA-seq cohort_analysis.final_adata_object (.h5ad)
		File   sc_rnaseq_reference_metadata_csv
		File   sc_rnaseq_cell_types_parquet   # cohort_analysis.scanvi_cell_types_parquet

		String batch_key
		Int n_top_genes
		Int n_comps
		Int n_neighbors
		Array[Float] leiden_res

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int disk_size = ceil(size(spatial_adata, "GB") + size(sc_rnaseq_reference_adata, "GB") + 50)

	command <<<
		set -euo pipefail

		echo "TODO: integrate spatial sample ~{sample_id} with scRNA-seq reference"
	>>>

	output {
		File integrated_adata_object  = "~{sample_id}.integrated.h5ad"
		File integration_metadata_csv = "~{sample_id}.integration_metadata.csv"
	}

	runtime {
		docker: "~{container_registry}/integration_tools:1.0.0"
		cpu: 8
		memory: "52 GB"
		disks: "local-disk " + disk_size + " SSD"
		preemptible: 1
		zones: zones
	}
}
