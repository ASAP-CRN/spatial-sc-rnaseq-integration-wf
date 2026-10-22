version 1.0

# One Visium sample: spatial-transcriptomics-wf Visium outputs from spatial-transcriptomics-wf 
# (*cleaned_unfiltered.h5ad, with raw UMIs and all genes--not only HVGs).
# These can be obtained from the dataset buckets:
# gs://asap-raw-<dataset_name>/workflow_execution/spatial_visium/preprocess/counts_to_adata/<version>
# Source: spatial_visium preprocess.qc_adata_objects
struct SpatialVisiumSample {
	String sample_id
	File qc_adata_object
}

# Processed Visium outputs from spatial-transcriptomics-wf used as inputs to this workflow.
# Each project provides a list of per-sample *cleaned_unfiltered.h5ad files to deconvolute independently.
struct SpatialVisiumProject {
	String asap_team_id
	String asap_dataset_id
	String asap_dataset_doi_url

	Array[SpatialVisiumSample] samples

	String raw_data_bucket
}
