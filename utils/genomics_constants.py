# Human mitochondrial protein-coding genes (GRCh38 Ensembl IDs).
MT_GENE_ENSEMBL_IDS = frozenset({
    "ENSG00000198888",  # MT-ND1
    "ENSG00000198763",  # MT-ND2
    "ENSG00000198804",  # MT-CO1
    "ENSG00000198712",  # MT-CO2
    "ENSG00000228253",  # MT-ATP8
    "ENSG00000198899",  # MT-ATP6
    "ENSG00000198938",  # MT-CO3
    "ENSG00000198840",  # MT-ND3
    "ENSG00000212907",  # MT-ND4L
    "ENSG00000198886",  # MT-ND4
    "ENSG00000198786",  # MT-ND5
    "ENSG00000198695",  # MT-ND6
    "ENSG00000198727",  # MT-CYB
})

# Kamath et al. (Nat Neurosci 2022) DA marker gene symbols; resolved to Ensembl IDs at
# runtime via resolve_marker_ensembl_ids (utils/h5ad_utils.py).
KAMATH_DA_MARKER_SYMBOLS = (
    # Pan-DA identity
    "TH", "SLC6A3", "SLC18A2", "DDC",
    # Branch markers (the two main DA subdivisions)
    "SOX6", "CALB1",
    # Subtype-defining markers (e.g. SOX6_AGTR1, CALB1_GEM)
    "AGTR1", "DDT", "GFRA2", "PART1", "CALCR", "CRYM", "CCDC68", "GEM", "PPP1R17", "RBP4", "TRHR",
)
