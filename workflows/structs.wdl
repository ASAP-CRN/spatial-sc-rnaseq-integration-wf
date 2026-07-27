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

# One Visium sample: per-sample post-QC adata (raw UMI counts, pre-Harmony) from spatial-transcriptomics-wf.
# Source: spatial_visium preprocess.qc_adata_objects — one file per sample, before Harmony batch correction.
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

# One GeoMx sample: per-sample Q3-normalized adata (pre-Harmony) from spatial-transcriptomics-wf.
# Source: spatial_geomx process_to_adata.processed_adata_objects — one file per sample.
# Q3 normalization is a within-sample/within-segment step; cross-sample Harmony correction
# runs later (in spatial_geomx cohort_analysis) and is NOT present in this object.
struct SpatialGeomxSample {
	String sample_id
	File   processed_adata_object    # spatial_geomx process_to_adata.processed_adata_objects (.h5ad) — per-sample, Q3-normalized, pre-Harmony
}

# Processed outputs from spatial-transcriptomics-wf (GeoMx) used as inputs to this workflow.
# Each project provides a list of per-sample processed adatas to integrate independently.
struct SpatialGeomxProject {
	String asap_team_id
	String asap_dataset_id
	String asap_dataset_doi_url

	Array[SpatialGeomxSample] samples

	String raw_data_bucket
	Array[String] staging_data_buckets
}
