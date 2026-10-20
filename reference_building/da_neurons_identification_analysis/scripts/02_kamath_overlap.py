"""All-vs-all barcode overlap analysis for Kamath H5AD files."""
import sys
import anndata as ad
import pandas as pd
from pathlib import Path
from itertools import combinations

repo_root = next(p for p in Path(__file__).resolve().parents if (p / ".git").exists())
ROOT = repo_root / "reference_building" / "da_neurons_identification_analysis" / "reference_data" / "kamath"

FILES = [
    "human_astrocytes_10x_scrna_seq",
    "human_da_neurons_10x_scrna_seq",
    "human_endothelial_cells_10x_scrna_seq",
    "human_microglia_10x_scrna_seq",
    "human_non_da_neurons_10x_scrna_seq",
    "human_nurr_negative_nuclei_10x_scrna_seq",
    "human_nurr_positive_nuclei_10x_scrna_seq",
    "human_oligodendrocytes_10x_scrna_seq",
    "human_opc_cells_10x_scrna_seq",
]

LABELS = {
    "human_astrocytes_10x_scrna_seq":           "astrocytes",
    "human_da_neurons_10x_scrna_seq":           "da_neurons",
    "human_endothelial_cells_10x_scrna_seq":    "endothelial",
    "human_microglia_10x_scrna_seq":            "microglia",
    "human_non_da_neurons_10x_scrna_seq":       "non_da_neurons",
    "human_nurr_negative_nuclei_10x_scrna_seq": "nurr_neg",
    "human_nurr_positive_nuclei_10x_scrna_seq": "nurr_pos",
    "human_oligodendrocytes_10x_scrna_seq":     "oligodendro",
    "human_opc_cells_10x_scrna_seq":            "opc",
}

print("Loading barcodes ...")
barcodes = {}
for f in FILES:
    path = ROOT / (f + ".h5ad")
    if not path.exists():
        print(f"  ERROR: file not found: {path}", file=sys.stderr)
        sys.exit(1)
    obs = ad.read_h5ad(path, backed='r').obs_names
    barcodes[f] = set(obs)
    print(f"  {len(obs):>7,}  {LABELS[f]}")

labels = [LABELS[f] for f in FILES]
matrix = pd.DataFrame(index=labels, columns=labels, data=0)

for f in FILES:
    matrix.loc[LABELS[f], LABELS[f]] = len(barcodes[f])

for f1, f2 in combinations(FILES, 2):
    n_overlap = len(barcodes[f1] & barcodes[f2])
    matrix.loc[LABELS[f1], LABELS[f2]] = n_overlap
    matrix.loc[LABELS[f2], LABELS[f1]] = n_overlap

print("\n=== All-vs-all barcode overlap (diagonal = file cell count) ===\n")
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 200)
print(matrix.to_string())

print("\n=== Overlapping pairs (n > 0) ===")
found = False
for f1, f2 in combinations(FILES, 2):
    n = matrix.loc[LABELS[f1], LABELS[f2]]
    if n > 0:
        print(f"  {LABELS[f1]:25s}  ∩  {LABELS[f2]:25s}  =  {n:,}")
        found = True
if not found:
    print("  None — all files have disjoint barcodes.")
