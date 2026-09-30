### Human brain 10Xv3 single nucleus transcriptomes
https://alleninstitute.github.io/abc_atlas_access/descriptions/WHB-10Xv3.html

Redirects to a couple of resources:
* Clustering analysis of over 3 million single nucleus transcriptomes spanning the whole adult human brain [(Cell and gene metadata)](https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/index.html#metadata/WHB-10Xv3/20241115/)
* Whole human brain taxonomy of cell types [(WHB-taxonomy)](https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/index.html#metadata/WHB-taxonomy/20240330/)

They include Allen Institute cell type annotations of the Siletti et al 2023 dataset, including `DA VGLUT2` neurons (siletti_allen_cluster `Splat_395`), which were not included in the Siletti deposition into CELLxGENE. In particular, these two files:
1. WHB-10Xv3/20241115/cell_metadata.csv.gz
2. WHB-taxonomy/20240330/cluster_to_cluster_annotation_membership.csv.gz

The data flow for the Silleti CELLxGENE h5ad `barcode` to be mapped to `siletti_allen_cluster` is as follows:


| Step | File | Role |
|---|---|---|
| 1 | `reference_data/siletti/{dissection_*}.h5ad` | Raw Siletti dissections H5ADs from CELLxGENE — `obs_names` are the barcodes everything keys off |
| 2 | `reference_data/allen-brain-cell-atlas/metadata/WHB-10Xv3/20241115/cell_metadata.csv.gz` | Barcode → `cluster_alias` (Allen's per-cell join key), columns `cell_label`, `cluster_alias` |
| 3 | `reference_data/allen-brain-cell-atlas/metadata/WHB-taxonomy/20240330/cluster_to_cluster_annotation_membership.csv.gz` | `cluster_alias` → annotation terms at 4 levels (`neurotransmitter`, `subcluster`, `cluster`, `supercluster`), keyed by `cluster_annotation_term_set_name` |
