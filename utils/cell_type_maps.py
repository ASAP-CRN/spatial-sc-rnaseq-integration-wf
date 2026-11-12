# Harmonizes Kamath et al. (2022) cell-type nomenclature (68 labels) against the
# Siletti/Allen Brain Cell Atlas nomenclature (41 labels).
#
# Derived empirically: Siletti cells were annotated with scANVI independently in two builds:
# Build    Cell-type labels source  repo, infile source                                                                        Column
# -------  -----------------------  -----------------------------------------------------------------------------------------  --------------
# Build 1  Kamath                   spatial-sc-rnaseq-integration-wf, deconvolution_reference_300k_c2lUda_Kamath200Rand.h5ad  cell_type_pred
# Build 2  Siletti/Allen            sc-rnaseq-wf, 5.0.0/2026-10-19T14-59-48Z/asap-cohort.final_metadata.csv.gz                C_scANVI
#
# The majority vote per Kamath-trained label, across Siletti cells where both predictions are available, defines the mapping.
# The 10 Kamath DA subtypes (SOX6_*/CALB1_*) are mapped to themselves.
# Unless indicated otherwise, the mappings are considered high-confidence based on the majority vote (frac ≥ 0.90)
#
# Note: several Ex_* and Inh_* labels map to "Splatter". Siletti's excitatory/inhibitory categories are defined
# by other brain regions (cortex, thalamus, amygdala, cerebellum, MGE/CGE lineage), not substantia nigra, so most 
# of Kamath's SN-specific Ex_*/Inh_* subtypes have no real anatomical match -- except Inh_OTX2_CASR, which correctly 
# lands on "Midbrain-derived inhibitory".
KAMATH_TO_SILETTI_ALLEN = {
    "Astro_CYP4F12": "Astrocyte",
    "Astro_GBP2_SPOCD1": "Astrocyte",
    "Astro_GJB6_OXTR": "Astrocyte",
    "Astro_GLYATL2": "Astrocyte",
    "Astro_GUCY1A2": "Astrocyte",
    "Astro_SERPINA3": "Astrocyte",
    "Astro_SIDT1": "Astrocyte",
    "Astro_VIM_LHX2": "Astrocyte",
    "Astro_VIM_TNFSRF12A": "Astrocyte",
    "CALB1_CALCR": "CALB1_CALCR",
    "CALB1_CRYM_CCDC68": "CALB1_CRYM_CCDC68",
    "CALB1_GEM": "CALB1_GEM",
    "CALB1_PPP1R17": "CALB1_PPP1R17",
    "CALB1_RBP4": "CALB1_RBP4",
    "CALB1_TRHR": "CALB1_TRHR",
    "Endo_COL6A3": "Fibroblast", # HIGH-CONFIDENCE (frac=0.994, n=310; second best [Oligodendrocyte precursor, frac=0.006, n=2])
    "Endo_DCN_ABCC9": "Vascular",
    "Endo_IL27RA": "Vascular",
    "Endo_MET": "Vascular",
    "Endo_NOTCH3_PLK2": "Vascular",
    "Endo_SLIT3": "Vascular",
    "Endo_SNTG2": "Vascular",
    "Ependyma_ZBBX": "Ependymal",  # LOW-CONFIDENCE (frac=0.870, n=23; second best [Astrocyte, frac=0.087, n=2])
    "Ex_CYP2J2": "Splatter",  # LOW-CONFIDENCE (frac=0.736, n=1076; second best [Thalamic excitatory, frac=0.261, n=281])
    "Ex_EBF2_CTC-552D5.1": "Splatter",
    "Ex_LAMP5_BAIAP3": "Splatter",
    "Ex_LAMP5_NTNG2": "Splatter",  # LOW-CONFIDENCE (frac=0.778, n=9; second best [Lower rhombic lip, frac=0.111, n=1])
    "Ex_MYO5B": "Splatter",
    "Ex_OPRD1": "Lower rhombic lip",  # LOW-CONFIDENCE (frac=0.644, n=236; second best [Splatter, frac=0.352, n=83])
    "Ex_POSTN": "Lower rhombic lip",  # LOW-CONFIDENCE (frac=0.492, n=61; second best [Splatter, frac=0.492, n=30])
    "Ex_PPP1R1C": "Splatter",
    "Ex_SATB2": "Medium spiny neuron",  # LOW-CONFIDENCE (frac=0.690, n=1499; second best [Hippocampal dentate gyrus, frac=0.241, n=362])
    "Ex_VWA5B1_CALB1": "Splatter",
    "Inh_IGFBP5": "Splatter",  # LOW-CONFIDENCE (frac=0.749, n=994; second best [CGE interneuron, frac=0.136, n=135])
    "Inh_INHBA": "Splatter",
    "Inh_OTX2_CASR": "Midbrain-derived inhibitory",
    "Inh_PAX5_CCBE1": "Splatter",
    "Inh_PAX5_VCAN": "Splatter",
    "Inh_PRLR_RP11-384J4.2": "Splatter",
    "Inh_SIX3": "Splatter",
    "MG_CCL3": "Microglia",  # LOW-CONFIDENCE (frac=0.833, n=18; second best [Miscellaneous, frac=0.167, n=3])
    "MG_CECR2_FGL1": "Microglia",
    "MG_FOSL2": "Microglia",
    "MG_GPNMB_LPL": "Microglia",
    "MG_GPNMB_SULT1C2": "Microglia",
    "MG_GPNMB_SUSD1": "Microglia",
    "MG_MGAM": "Miscellaneous",  # HIGH-CONFIDENCE (frac=0.972, n=353; second best [Oligodendrocyte, frac=0.017, n=6])
    "MG_MKI67": "Microglia",
    "MG_OPRM1": "Microglia",
    "MG_SPON1": "Microglia",
    "MG_TSPO_VIM": "Microglia",
    "Macro_CD200R1": "Microglia",
    "OPC_ADM": "Oligodendrocyte precursor",
    "OPC_CACNG4": "Oligodendrocyte precursor",
    "OPC_HOXD3": "Oligodendrocyte precursor",
    "OPC_KIAA0040": "Oligodendrocyte precursor",
    "OPC_MDFI": "Committed oligodendrocyte precursor",
    "Olig_ENPP6_ACTN2": "Oligodendrocyte",
    "Olig_ENPP6_EMILIN2": "Oligodendrocyte",
    "Olig_ENPP6_LUCAT1": "Oligodendrocyte",
    "Olig_PLXDC2": "Oligodendrocyte",
    "Olig_PLXDC2_KCNAB1": "Oligodendrocyte",
    "Olig_PLXDC2_KCNK10": "Oligodendrocyte",
    "Olig_PLXDC2_SFRP1": "Oligodendrocyte",
    "SOX6_AGTR1": "SOX6_AGTR1",
    "SOX6_DDT": "SOX6_DDT",
    "SOX6_GFRA2": "SOX6_GFRA2",
    "SOX6_PART1": "SOX6_PART1",
}
