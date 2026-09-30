version 1.0

# Processed outputs from sc-rnaseq-wf used as inputs to this workflow.
# Source fields come from the cohort_analysis output block of sc-rnaseq-wf.
struct ScRnaseqProject {
	String asap_team_id
	String asap_dataset_id
	String asap_dataset_doi_url

	File final_adata_object          # cohort_analysis.final_adata_object (.h5ad) — Harmony-integrated, annotated
	File final_metadata_csv          # cohort_analysis.final_metadata_csv
	File scanvi_cell_types_parquet   # cohort_analysis.scanvi_cell_types_parquet — scANVI cell type predictions

	String raw_data_bucket
	Array[String] staging_data_buckets
}

# One Visium sample: per-sample post-QC adata (raw UMI counts) from spatial-transcriptomics-wf.
# Source: spatial_visium preprocess.qc_adata_objects — one file per sample, post-QC.
struct SpatialVisiumSample {
	String sample_id
	File   qc_adata_object           # spatial_visium preprocess.qc_adata_objects (.h5ad) — per-sample, post-QC, raw UMI counts
}

# Processed outputs from spatial-transcriptomics-wf (Visium) used as inputs to this workflow.
# Each project provides a list of per-sample pre-QC adatas to integrate independently.
struct SpatialVisiumProject {
	String asap_team_id
	String asap_dataset_id
	String asap_dataset_doi_url

	Array[SpatialVisiumSample] samples

	String raw_data_bucket
	Array[String] staging_data_buckets
}
