#!/usr/bin/env python3
"""
Collapse cell2location's 68 Kamath cell subtype labels to match
Ma et al. (2025)'s 13 BayesSpace/RCTD cell type labels via `kamath_collapse_dic`,
then compare the two as per-spot proportion vectors using PCC, RMSE, and
Jensen-Shannon divergence -- following Li et al. 2022 Nat Methods' deconvolution
benchmark metrics.

Note: in the workflow, we DON'T assign the per-spot winner across the Ma (2025)
13 cell-type columns. We instead use the continuous proportion values for
downstream analysis. The reason for this is because Dopaminergic_Neuron and
Fibrous_Astrocyte proportions are large enough in nearly every cluster to
swamp the rest, so they would always "win" if we took the argmax, which
would lead to flawed results.
"""

import os
import sys
import argparse
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import pearsonr
from scipy.spatial.distance import jensenshannon

_DATA_DIR = "reference_building/benchmark/reference_data/ma/data"
_CELL2LOCATION_DIR = f"{_DATA_DIR}/reformatted/cell2location_Strategy1"
_MA_BENCHMARK_CSV = f"{_DATA_DIR}/reformatted/per_spot_benchmark_al_dalahmah_sn.csv.gz"
_OUTPUT_DIR = _CELL2LOCATION_DIR
_H5AD_SUFFIX = ".deconvoluted_clustered.h5ad"
_PER_SPOT_OUT = "asap_agreement_per_spot.csv.gz"
_METRICS_OUT = "asap_agreement_metrics.csv"
_ALL_VS_ALL_OUT = "asap_agreement_all_vs_all.csv"
_MARGINALS_OUT = "asap_agreement_marginals.csv"

# Kamath sub-cell type harmonization
kamath_collapse_dic = {
    "Olig_PLXDC2_SFRP1": "Oligodendrocyte",
    "Olig_ENPP6_LUCAT1": "Oligodendrocyte",
    "Olig_PLXDC2": "Oligodendrocyte",
    "Olig_PLXDC2_KCNAB1": "Oligodendrocyte",
    "Olig_ENPP6_EMILIN2": "Oligodendrocyte",
    "Olig_PLXDC2_KCNK10": "Oligodendrocyte",
    "Olig_ENPP6_ACTN2": "Oligodendrocyte",

    "OPC_CACNG4": "Oligod. Precursor",
    "OPC_KIAA0040": "Oligod. Precursor",
    "OPC_ADM": "Oligod. Precursor",
    "OPC_HOXD3": "Oligod. Precursor",
    "OPC_MDFI": "Oligod. Precursor",

    # Astrocytes -- fibrous/protoplasmic not distinguishable from these markers
    "Astro_GJB6_OXTR": "Astrocyte",
    "Astro_SIDT1": "Astrocyte",
    "Astro_CYP4F12": "Astrocyte",
    "Astro_GUCY1A2": "Astrocyte",
    "Astro_SERPINA3": "Astrocyte",
    "Astro_GLYATL2": "Astrocyte",
    "Astro_VIM_LHX2": "Astrocyte",
    "Astro_VIM_TNFSRF12A": "Astrocyte",
    "Astro_GBP2_SPOCD1": "Astrocyte",

    "MG_CECR2_FGL1": "Microglia",
    "MG_OPRM1": "Microglia",
    "MG_GPNMB_SUSD1": "Microglia",
    "MG_TSPO_VIM": "Microglia",
    "MG_GPNMB_SULT1C2": "Microglia",
    "MG_SPON1": "Microglia",
    "MG_GPNMB_LPL": "Microglia",
    "MG_MGAM": "Microglia",
    "MG_FOSL2": "Microglia",
    "MG_MKI67": "Microglia",
    "MG_CCL3": "Microglia",

    "Inh_PRLR_RP11-384J4.2": "Inhibitory neuron",
    "Inh_OTX2_CASR": "Inhibitory neuron",
    "Inh_SIX3": "Inhibitory neuron",
    "Inh_PAX5_CCBE1": "Inhibitory neuron",
    "Inh_IGFBP5": "Inhibitory neuron",
    "Inh_PAX5_VCAN": "Inhibitory neuron",
    "Inh_INHBA": "Inhibitory neuron",

    "CALB1_PPP1R17": "DA neuron",
    "CALB1_CALCR": "DA neuron",
    "CALB1_GEM": "DA neuron",
    "CALB1_CRYM_CCDC68": "DA neuron",
    "CALB1_TRHR": "DA neuron",
    "CALB1_RBP4": "DA neuron",
    "SOX6_AGTR1": "DA neuron",
    "SOX6_PART1": "DA neuron",
    "SOX6_GFRA2": "DA neuron",
    "SOX6_DDT": "DA neuron",

    "Endo_IL27RA": "Endothelial",
    "Endo_DCN_ABCC9": "Endothelial",
    "Endo_COL6A3": "Endothelial",
    "Endo_NOTCH3_PLK2": "Endothelial",
    "Endo_SLIT3": "Endothelial",
    "Endo_SNTG2": "Endothelial",
    "Endo_MET": "Endothelial",

    # Border-associated macrophage
    "Macro_CD200R1": "BA macrophages",

    # Excitatory neurons -- no direct Ma-side counterpart in SN tissue
    "Ex_SATB2": "Excitatory neuron",
    "Ex_POSTN": "Excitatory neuron",
    "Ex_VWA5B1_CALB1": "Excitatory neuron",
    "Ex_EBF2_CTC-552D5.1": "Excitatory neuron",
    "Ex_MYO5B": "Excitatory neuron",
    "Ex_OPRD1": "Excitatory neuron",
    "Ex_PPP1R1C": "Excitatory neuron",
    "Ex_LAMP5_NTNG2": "Excitatory neuron",
    "Ex_CYP2J2": "Excitatory neuron",
    "Ex_LAMP5_BAIAP3": "Excitatory neuron",

    # No Ma-side counterpart in SN tissue
    "Ependyma_ZBBX": "Ependyma",

}

# Ma cell type harmonization
MA_CELL_TYPE_COLS = [
    "Activated_microglia", "Quiescent_microglia", "Dopaminergic_Neuron", "Endothelial",
    "Fibrous_Astrocyte", "GABAergic_Neuron", "MEG3_Interneurons", "BAM", "Oligodendrocyte",
    "OPC", "Protoplasmic_Astrocyte", "Tcells", "VascSmMs",
]
ma_pool_dic = {
    "Activated_microglia": "Microglia",
    "Quiescent_microglia": "Microglia",
    "Dopaminergic_Neuron": "DA neuron",
    "Endothelial": "Endothelial",
    "Fibrous_Astrocyte": "Astrocyte",
    "GABAergic_Neuron": "Inhibitory neuron",
    "MEG3_Interneurons": "Interneuron",
    "BAM": "BA macrophages",
    "Oligodendrocyte": "Oligodendrocyte",
    "OPC": "Oligod. Precursor",
    "Protoplasmic_Astrocyte": "Astrocyte",
    "Tcells": "Tcells",
    "VascSmMs": "VascSmMs",
}


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

class _Tee:
    """Mirror writes to both the original stream and a log file (used for stdout/stderr)."""

    def __init__(self, original, log_file) -> None:
        self._original = original
        self._log = log_file

    def write(self, data: str) -> int:
        self._original.write(data)
        self._log.write(data)
        return len(data)

    def flush(self) -> None:
        self._original.flush()
        self._log.flush()

    def __getattr__(self, name):
        return getattr(self._original, name)


def log_run_command(args: argparse.Namespace) -> None:
    """Log the invocation to {_DATA_DIR}/{script_stem}.log and tee all subsequent output there."""
    log_file = Path(_DATA_DIR) / f"{Path(__file__).stem}.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    f = open(log_file, "w", encoding="utf-8")  # noqa: SIM115 — intentionally kept open
    f.write(f"[{timestamp}] {' '.join(sys.argv)}\n\n")
    sys.stdout = _Tee(sys.stdout, f)
    sys.stderr = _Tee(sys.stderr, f)


# ---------------------------------------------------------------------------
# Collapsing
# ---------------------------------------------------------------------------

def _strip_ct_prefix(name: str) -> str:
    """Strip cell2location's `q05_cell_abundance_w_sf` obsm/obs column prefix."""
    for prefix in ("frac_q05cell_abundance_w_sf_", "q05cell_abundance_w_sf_"):
        if name.startswith(prefix):
            return name[len(prefix):]
    return name


def load_sample_fractions(h5ad_path: Path) -> pd.DataFrame:
    """
    Load a deconvoluted sample's per-spot fine-grained cell-type fractions.

    Parameters
    ----------
    h5ad_path : pathlib.Path
        Path to a <sample_id>.deconvoluted_clustered.h5ad from cell2location_deconv.py.

    Returns
    -------
    pandas.DataFrame
        Spots x fine cell_type_pred labels; each row sums to 1.
    """
    adata = ad.read_h5ad(h5ad_path)
    cell_types = list(adata.obsm["q05_cell_abundance_w_sf"].columns)
    frac_df = adata.obsm["cell_type_fractions"][[f"frac_{ct}" for ct in cell_types]].copy()
    frac_df.columns = [_strip_ct_prefix(ct) for ct in cell_types]
    frac_df.index = adata.obs_names
    return frac_df


def collapse_columns(df: pd.DataFrame, collapse_dic: dict, unmapped_label: str = "Unmapped") -> pd.DataFrame:
    """
    Sum `df`'s columns into pooled categories per `collapse_dic`.

    Parameters
    ----------
    df : pandas.DataFrame
        Spots x fine-grained labels.
    collapse_dic : dict
        Maps each column name in `df` to a pooled category name, or to None
        to route it into `unmapped_label` instead of discarding it silently.
    unmapped_label : str
        Column collecting the summed mass of all None-mapped columns.

    Returns
    -------
    pandas.DataFrame
        Spots x pooled categories (plus `unmapped_label` if any column
        mapped to None).
    """
    targets = {col: (collapse_dic.get(col, unmapped_label) or unmapped_label) for col in df.columns}
    return df.T.groupby(targets).sum().T


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

def build_comparison(cell2location_dir: str, ma_benchmark_csv: str, sample_ids=None) -> tuple:
    """
    Load, collapse, and join every sample's pipeline output against Ma's benchmark.

    Parameters
    ----------
    cell2location_dir : str
        Directory of <sample_id>.deconvoluted_clustered.h5ad files.
    ma_benchmark_csv : str
        Path to per_spot_benchmark_al_dalahmah_sn.csv.gz (03_build_benchmark.py).
    sample_ids : list of str or None
        Restrict to these sample IDs (e.g. "T-0081"); None processes all found.

    Returns
    -------
    tuple of pandas.DataFrame
        (our_all, ma_all): spots x [pooled categories + "sample"], row-aligned
        1:1 between the two frames.
    """
    ma_bench = pd.read_csv(ma_benchmark_csv)
    our_rows, ma_rows = [], []

    h5ad_paths = sorted(Path(cell2location_dir).glob(f"*{_H5AD_SUFFIX}"))
    if sample_ids:
        h5ad_paths = [p for p in h5ad_paths if p.name.replace(_H5AD_SUFFIX, "") in sample_ids]

    for h5ad_path in h5ad_paths:
        sample_name = h5ad_path.name.replace(_H5AD_SUFFIX, "")
        our_frac = load_sample_fractions(h5ad_path)
        our_pooled = collapse_columns(our_frac, kamath_collapse_dic)

        ma_sample_id = sample_name.replace("-", "")
        ma_sample = ma_bench[ma_bench["sample"] == ma_sample_id].set_index("barcode")
        ma_sample = ma_sample.reindex(our_pooled.index)
        ma_pooled = collapse_columns(ma_sample[MA_CELL_TYPE_COLS], ma_pool_dic)

        joined = our_pooled.notna().all(axis=1) & ma_pooled.notna().all(axis=1)
        our_pooled, ma_pooled = our_pooled[joined].copy(), ma_pooled[joined].copy()
        our_pooled["sample"], ma_pooled["sample"] = sample_name, sample_name

        print(f"  {sample_name}: {joined.sum()}/{len(joined)} spots joined")
        our_rows.append(our_pooled)
        ma_rows.append(ma_pooled)

    return pd.concat(our_rows), pd.concat(ma_rows)


def compute_metrics(ours: np.ndarray, ma: np.ndarray) -> dict:
    """
    PCC, RMSE, and Jensen-Shannon divergence between two proportion vectors.

    Parameters
    ----------
    ours : numpy.ndarray
        This pipeline's per-spot proportion for one pooled category.
    ma : numpy.ndarray
        Ma's per-spot (cluster-stamped) proportion for the same category.

    Returns
    -------
    dict
        pcc : float
            Pearson correlation across spots; NaN if either side is constant.
        rmse : float
            Root mean squared error across spots, on the raw 0-1 scale.
        js : float
            Jensen-Shannon divergence (base 2) between the two spot-wise
            spatial distributions, each normalized to sum to 1.
    """
    pcc = pearsonr(ours, ma)[0] if ours.std() > 0 and ma.std() > 0 else float("nan")
    rmse = float(np.sqrt(np.mean((ours - ma) ** 2)))
    ours_norm = ours / ours.sum() if ours.sum() > 0 else ours
    ma_norm = ma / ma.sum() if ma.sum() > 0 else ma
    js = float(jensenshannon(ours_norm, ma_norm, base=2))
    return {"pcc": pcc, "rmse": rmse, "js": js}


def compute_all_metrics(our_all: pd.DataFrame, ma_all: pd.DataFrame) -> pd.DataFrame:
    """Per pooled category present on both sides, compute `compute_metrics` across all spots."""
    categories = sorted((set(our_all.columns) & set(ma_all.columns)) - {"sample"})
    rows = []
    for cat in categories:
        row = compute_metrics(our_all[cat].values, ma_all[cat].values)
        row["category"] = cat
        rows.append(row)
    return pd.DataFrame(rows).set_index("category")[["pcc", "rmse", "js"]].sort_values("pcc", ascending=False)


def compute_all_vs_all(our_all: pd.DataFrame, ma_all: pd.DataFrame) -> pd.DataFrame:
    """
    PCC between every our-pooled-category and every Ma-pooled-category, not
    just matched pairs. Off-diagonal cells reveal confusions `compute_all_metrics`
    (matched pairs only) hides.

    Parameters
    ----------
    our_all : pandas.DataFrame
        Spots x our pooled categories (plus "sample"), from `build_comparison`.
    ma_all : pandas.DataFrame
        Spots x Ma pooled categories (plus "sample"), from `build_comparison`.

    Returns
    -------
    pandas.DataFrame
        Our pooled categories (rows) x Ma pooled categories (columns), PCC values.
    """
    our_cats = [c for c in our_all.columns if c != "sample"]
    ma_cats = [c for c in ma_all.columns if c != "sample"]
    return pd.DataFrame(
        {
            our_cat: {
                ma_cat: compute_metrics(our_all[our_cat].values, ma_all[ma_cat].values)["pcc"]
                for ma_cat in ma_cats
            }
            for our_cat in our_cats
        }
    ).T  # rows = our pooled category, columns = Ma pooled category


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(args: argparse.Namespace) -> None:
    """Build the collapsed-category comparison and report (or write) per-category agreement metrics."""
    print(f"Loading cell2location outputs from {args.cell2location_dir}")
    our_all, ma_all = build_comparison(args.cell2location_dir, args.ma_benchmark_csv, args.samples)
    print(f"\nPooled {len(our_all)} spots across {our_all['sample'].nunique()} sample(s)")

    metrics = compute_all_metrics(our_all, ma_all)
    print("\nPer-category agreement (continuous proportions, not argmax labels):")
    print(metrics.to_string(float_format=lambda x: f"{x:.3f}"))

    all_vs_all = compute_all_vs_all(our_all, ma_all)
    print("\nAll-vs-all agreement (PCC, our pooled category x Ma pooled category):")
    print(all_vs_all.to_string(float_format=lambda x: f"{x:.2f}"))

    os.makedirs(args.output_dir, exist_ok=True)

    # Join on (sample, barcode) explicitly -- barcode alone repeats across
    # samples and would cross-join unrelated spots that share one.
    our_reset = our_all.reset_index(names="barcode").rename(
        columns={c: f"our_{c}" for c in our_all.columns if c != "sample"}
    )
    ma_reset = ma_all.reset_index(names="barcode").rename(
        columns={c: f"ma_{c}" for c in ma_all.columns if c != "sample"}
    )
    merged = our_reset.merge(ma_reset, on=["sample", "barcode"])

    # Summed proportion mass per category -- a continuous "cell count"
    # analogue, since no spot is hard-assigned. Feeds Section 7b's heatmap marginals.
    our_totals = our_all.drop(columns="sample").sum().rename("total")
    ma_totals = ma_all.drop(columns="sample").sum().rename("total")
    marginals = pd.concat([
        our_totals.rename_axis("category").reset_index().assign(side="our"),
        ma_totals.rename_axis("category").reset_index().assign(side="ma"),
    ], ignore_index=True)

    merged.to_csv(os.path.join(args.output_dir, _PER_SPOT_OUT), compression="gzip", index=False)
    metrics.to_csv(os.path.join(args.output_dir, _METRICS_OUT))
    all_vs_all.to_csv(os.path.join(args.output_dir, _ALL_VS_ALL_OUT))
    marginals.to_csv(os.path.join(args.output_dir, _MARGINALS_OUT), index=False)
    print(f"\nWritten {os.path.join(args.output_dir, _PER_SPOT_OUT)}")
    print(f"Written {os.path.join(args.output_dir, _METRICS_OUT)}")
    print(f"Written {os.path.join(args.output_dir, _ALL_VS_ALL_OUT)}")
    print(f"Written {os.path.join(args.output_dir, _MARGINALS_OUT)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--cell2location-dir",
        default=None,
        help="Directory of <sample_id>.deconvoluted_clustered.h5ad files. A relative path is\n"
             "resolved against the directory this command is run from.\n"
             f"[Default: {_CELL2LOCATION_DIR}, relative to the repo root]",
    )
    parser.add_argument(
        "--ma-benchmark-csv",
        default=None,
        help="per_spot_benchmark_al_dalahmah_sn.csv.gz from 03_build_benchmark.py. A relative\n"
             "path is resolved against the directory this command is run from.\n"
             f"[Default: {_MA_BENCHMARK_CSV}, relative to the repo root]",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory for the per-spot and per-category comparison tables. A relative\n"
             "path is resolved against the directory this command is run from.\n"
             f"[Default: {_OUTPUT_DIR}, relative to the repo root]",
    )
    parser.add_argument(
        "--samples",
        nargs="+",
        default=None,
        help="Restrict to these sample IDs (e.g. T-0081 T-2465). Default: all samples found.",
    )
    args = parser.parse_args()

    invocation_cwd = Path.cwd()
    for attr in ("cell2location_dir", "ma_benchmark_csv", "output_dir"):
        value = getattr(args, attr)
        if value is not None:
            setattr(args, attr, str((invocation_cwd / value).resolve()))

    repo_root = next(p for p in Path(__file__).resolve().parents if (p / ".git").exists())
    if Path.cwd() != repo_root:
        print(f"[INFO] Changing working directory to repo root: {repo_root}")
        os.chdir(repo_root)

    if args.cell2location_dir is None:
        args.cell2location_dir = _CELL2LOCATION_DIR
    if args.ma_benchmark_csv is None:
        args.ma_benchmark_csv = _MA_BENCHMARK_CSV
    if args.output_dir is None:
        args.output_dir = _OUTPUT_DIR

    log_run_command(args)
    main(args)
