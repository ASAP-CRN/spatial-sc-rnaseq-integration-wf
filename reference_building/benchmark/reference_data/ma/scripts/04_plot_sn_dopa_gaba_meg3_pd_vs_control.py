#!/usr/bin/env python3
"""
Regenerate the Substantia Nigra dopaminergic/GABAergic/MEG3 interneuron
composition report (PD vs Control) from Ma et al. (2025) Supplementary Data 1
tabs A and C (41467_2025_62478_MOESM3_ESM.xlsx): Mann-Whitney/Friedman/Wilcoxon
significance tests, rendered as a self-contained HTML report (chart, 6x6
p-value matrix, per-sample table).

Usage
-----
  python3 04_plot_sn_dopa_gaba_meg3_pd_vs_control.py
"""

import os
import sys
import json
import argparse
from pathlib import Path
from datetime import datetime

import openpyxl
import pandas as pd
from scipy import stats

_DATA_DIR = "reference_building/benchmark/reference_data/ma/data"
_OUTDIR = f"{_DATA_DIR}/analysis"
_SUPP_XLSX = f"{_DATA_DIR}/original/41467_2025_62478_MOESM3_ESM.xlsx"
_OUTPUT_HTML = f"{_OUTDIR}/da_sn_dopa_gaba_meg3_pd_vs_control.html"
os.makedirs(_OUTDIR, exist_ok=True)

_SN_SHEET_A = "A"
_SN_SHEET_C = "C"
_BRAIN_REGION = "Substantia Nigra"
_NEURON_COLS = ["Dopaminergic", "GABAergic", "MEG3 Interneurons"]

# tab A columns selectable as the chart's point-shape encoding (applied identically to both facets).
_SHAPE_CATEGORY_COLS = ["Neuropath_Diagnosis", "Age", "Sex", "path_braak_nft", "path_braak_asyn", "path_cerad"]

# Marker shapes available; overflow columns fold into a shared "Other" bucket.
_N_SHAPES_AVAILABLE = 8


# ---------------------------------------------------------------------------
# Logging (mirrors 03_build_benchmark.py's convention in this directory)
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
# Data loading
# ---------------------------------------------------------------------------

def load_sn_metadata(xlsx_path: str) -> pd.DataFrame:
    """
    Load per-sample Substantia Nigra metadata from tab A.

    Parameters
    ----------
    xlsx_path : str
        Path to 41467_2025_62478_MOESM3_ESM.xlsx.

    Returns
    -------
    pandas.DataFrame
        One row per Substantia Nigra sample, columns `TB #`, `Condition`,
        and the six shape-category columns in `_SHAPE_CATEGORY_COLS`
        (Neuropath_Diagnosis, Age, Sex, path_braak_nft, path_braak_asyn,
        path_cerad).
    """
    wb = openpyxl.load_workbook(xlsx_path, read_only=True, data_only=True)
    rows = list(wb[_SN_SHEET_A].iter_rows(values_only=True))
    header = rows[0]
    meta = pd.DataFrame(rows[1:], columns=header)
    keep_cols = ["TB #", "brain_region", "Condition"] + _SHAPE_CATEGORY_COLS
    meta = meta[keep_cols].dropna(subset=["TB #"])
    meta_sn = meta[meta["brain_region"] == _BRAIN_REGION].drop_duplicates(subset=["TB #"])
    return meta_sn.drop(columns=["brain_region"])


def load_sn_celltype_counts(xlsx_path: str) -> pd.DataFrame:
    """
    Load per-sample nucleus counts for the Substantia Nigra block of tab C.

    Parameters
    ----------
    xlsx_path : str
        Path to 41467_2025_62478_MOESM3_ESM.xlsx.

    Returns
    -------
    pandas.DataFrame
        One row per Substantia Nigra sample, columns `TB #` plus one column
        per cell type (Fibrous-like astrocytes ... Endothelial).
    """
    wb = openpyxl.load_workbook(xlsx_path, read_only=True, data_only=True)
    rows = list(wb[_SN_SHEET_C].iter_rows(values_only=True))
    header = [c for c in rows[1] if c is not None]

    sample_rows = []
    for row in rows[2:]:
        if row[0] == "sum":
            break
        if row[0] is None:
            continue
        sample_rows.append(row[: len(header)])

    return pd.DataFrame(sample_rows, columns=header).rename(columns={"Sample": "TB #"})


def build_merged_table(xlsx_path: str) -> pd.DataFrame:
    """
    Join Substantia Nigra metadata and cell-type counts, and derive the
    per-sample neuron-subtype share columns used throughout this report.

    Parameters
    ----------
    xlsx_path : str
        Path to 41467_2025_62478_MOESM3_ESM.xlsx.

    Returns
    -------
    pandas.DataFrame
        One row per Substantia Nigra sample with the raw cell-type counts,
        `total_cells`, `total_neurons`, `other_nonneuron_sum`, and
        `<subtype>_pct_of_neurons` / `<subtype>_pct_of_total` for each of
        Dopaminergic, GABAergic, MEG3 Interneurons.
    """
    meta_sn = load_sn_metadata(xlsx_path)
    counts = load_sn_celltype_counts(xlsx_path)
    celltype_cols = [c for c in counts.columns if c != "TB #"]

    merged = counts.merge(meta_sn, on="TB #", how="left")
    unmatched = merged.loc[merged["Condition"].isna(), "TB #"].tolist()
    if unmatched:
        print(f"  WARNING: {len(unmatched)} Substantia Nigra sample(s) had no Condition match in tab A: {unmatched}")

    merged["total_cells"] = merged[celltype_cols].sum(axis=1)
    merged["total_neurons"] = merged[_NEURON_COLS].sum(axis=1)

    other_cols = [c for c in celltype_cols if c not in _NEURON_COLS]
    merged["other_nonneuron_sum"] = merged[other_cols].sum(axis=1)

    for c in _NEURON_COLS:
        merged[f"{c}_pct_of_total"] = merged[c] / merged["total_cells"] * 100
        merged[f"{c}_pct_of_neurons"] = merged[c] / merged["total_neurons"] * 100

    for col in _SHAPE_CATEGORY_COLS:
        n_distinct = merged[col].nunique(dropna=False)
        if n_distinct > _N_SHAPES_AVAILABLE:
            print(
                f"  WARNING: {col} has {n_distinct} distinct values across these samples, "
                f"more than the {_N_SHAPES_AVAILABLE} shapes the chart can draw -- the JS "
                "will fold the overflow into a shared 'Other' marker."
            )

    return merged


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def format_p(p: float) -> str:
    """
    Format a p-value the way the report's prose quotes it: four decimal
    places once it drops below the 0.01 threshold, two decimal places
    otherwise.

    Parameters
    ----------
    p : float
        A p-value.

    Returns
    -------
    str
        The formatted p-value, with no leading "p = ".
    """
    return f"{p:.4f}" if p < 0.01 else f"{p:.2f}"


def format_range(a: float, b: float) -> str:
    """
    Format two values as "lo–hi%" (or a single "v%" if they round the same).

    Parameters
    ----------
    a : float
        First value, in percent.
    b : float
        Second value, in percent.

    Returns
    -------
    str
        The formatted range string, without a trailing "%" on the low end
        when the two differ (e.g. "34–37%"), or a single "58%" when equal.
    """
    lo, hi = sorted((round(a), round(b)))
    return f"{lo}%" if lo == hi else f"{lo}–{hi}%"


def compute_stats(merged: pd.DataFrame) -> dict:
    """
    Compute every number the HTML report's narrative text and charts need.

    Parameters
    ----------
    merged : pandas.DataFrame
        Output of `build_merged_table`.

    Returns
    -------
    dict
        n_control, n_pd : int
            Sample counts, Control / PD.
        mean_da_pctneurons_control, mean_da_pctneurons_pd : float
            Mean Dopaminergic share of neurons.
        p_mw_da_pctneurons : float
            Mann-Whitney p, Dopaminergic share of neurons, PD vs Control.
        mean_da_raw_control, mean_da_raw_pd : float
            Mean Dopaminergic share of all nuclei (all samples).
        p_mw_da_raw : float
            Mann-Whitney p for the row above.
        fans_sorted_samples : list of str
            TB # of samples with zero non-neuronal nuclei recorded (likely FANS-sorted).
        mean_da_excl_control, mean_da_excl_pd : float
            Mean Dopaminergic share of all nuclei, FANS-sorted samples excluded.
        p_mw_da_excl : float
            Mann-Whitney p for the row above.
        p_friedman_control : float
            Friedman test p across the three neuron subtypes, Control.
        p_wilcoxon_da_gaba, p_wilcoxon_gaba_meg3, p_wilcoxon_da_meg3 : float
            Wilcoxon signed-rank p between each subtype pair, Control.
        meg3_range, da_range, gaba_range : str
            Formatted share-of-neurons range across Control and PD, per subtype.
        top_da_samples : list of str
            TB # of the 3 samples with the highest Dopaminergic share of neurons.
    """
    ctrl = merged[merged["Condition"] == "Control"]
    pd_ = merged[merged["Condition"] == "PD"]

    da_pctneurons_ctrl = ctrl["Dopaminergic_pct_of_neurons"].dropna()
    da_pctneurons_pd = pd_["Dopaminergic_pct_of_neurons"].dropna()
    _, p_mw_da_pctneurons = stats.mannwhitneyu(da_pctneurons_pd, da_pctneurons_ctrl, alternative="two-sided")

    da_raw_ctrl = ctrl["Dopaminergic_pct_of_total"].dropna()
    da_raw_pd = pd_["Dopaminergic_pct_of_total"].dropna()
    _, p_mw_da_raw = stats.mannwhitneyu(da_raw_pd, da_raw_ctrl, alternative="two-sided")

    fans_sorted = merged.loc[merged["other_nonneuron_sum"] == 0, "TB #"].tolist()
    no_fans = merged[merged["other_nonneuron_sum"] > 0]
    da_excl_ctrl = no_fans.loc[no_fans["Condition"] == "Control", "Dopaminergic_pct_of_total"].dropna()
    da_excl_pd = no_fans.loc[no_fans["Condition"] == "PD", "Dopaminergic_pct_of_total"].dropna()
    _, p_mw_da_excl = stats.mannwhitneyu(da_excl_pd, da_excl_ctrl, alternative="two-sided")

    _, p_friedman = stats.friedmanchisquare(
        ctrl["Dopaminergic_pct_of_neurons"], ctrl["GABAergic_pct_of_neurons"], ctrl["MEG3 Interneurons_pct_of_neurons"]
    )
    _, p_wilcoxon_da_gaba = stats.wilcoxon(ctrl["Dopaminergic_pct_of_neurons"], ctrl["GABAergic_pct_of_neurons"])
    _, p_wilcoxon_gaba_meg3 = stats.wilcoxon(ctrl["GABAergic_pct_of_neurons"], ctrl["MEG3 Interneurons_pct_of_neurons"])
    _, p_wilcoxon_da_meg3 = stats.wilcoxon(ctrl["MEG3 Interneurons_pct_of_neurons"], ctrl["Dopaminergic_pct_of_neurons"])

    top_da_samples = merged.nlargest(3, "Dopaminergic_pct_of_neurons")["TB #"].tolist()

    return {
        "n_control": len(ctrl),
        "n_pd": len(pd_),
        "mean_da_pctneurons_control": da_pctneurons_ctrl.mean(),
        "mean_da_pctneurons_pd": da_pctneurons_pd.mean(),
        "p_mw_da_pctneurons": p_mw_da_pctneurons,
        "mean_da_raw_control": da_raw_ctrl.mean(),
        "mean_da_raw_pd": da_raw_pd.mean(),
        "p_mw_da_raw": p_mw_da_raw,
        "fans_sorted_samples": fans_sorted,
        "mean_da_excl_control": da_excl_ctrl.mean(),
        "mean_da_excl_pd": da_excl_pd.mean(),
        "p_mw_da_excl": p_mw_da_excl,
        "p_friedman_control": p_friedman,
        "p_wilcoxon_da_gaba": p_wilcoxon_da_gaba,
        "p_wilcoxon_gaba_meg3": p_wilcoxon_gaba_meg3,
        "p_wilcoxon_da_meg3": p_wilcoxon_da_meg3,
        "meg3_range": format_range(
            ctrl["MEG3 Interneurons_pct_of_neurons"].mean(), pd_["MEG3 Interneurons_pct_of_neurons"].mean()
        ),
        "da_range": format_range(da_pctneurons_ctrl.mean(), da_pctneurons_pd.mean()),
        "gaba_range": format_range(ctrl["GABAergic_pct_of_neurons"].mean(), pd_["GABAergic_pct_of_neurons"].mean()),
        "top_da_samples": top_da_samples,
    }


# ---------------------------------------------------------------------------
# HTML rendering
# ---------------------------------------------------------------------------

def build_data_json(merged: pd.DataFrame) -> str:
    """
    Serialize the per-sample records the chart's client-side JS needs.

    Parameters
    ----------
    merged : pandas.DataFrame
        Output of `build_merged_table`.

    Returns
    -------
    str
        Compact JSON array, one object per Substantia Nigra sample.
    """
    cols = (
        ["TB #", "Condition"] + _SHAPE_CATEGORY_COLS + [
            "Dopaminergic", "GABAergic", "MEG3 Interneurons",
            "total_cells", "total_neurons",
            "Dopaminergic_pct_of_neurons", "GABAergic_pct_of_neurons", "MEG3 Interneurons_pct_of_neurons",
        ]
    )
    records = merged[cols].round(3).to_dict(orient="records")
    return json.dumps(records, separators=(",", ":"))


def build_pvalue_matrix_table(merged: pd.DataFrame) -> str:
    """
    Build the 6x6 HTML table of pairwise p-values across the three neuron
    subtypes crossed with Control/PD: paired Wilcoxon within a condition
    (same samples, 3 compositions each), Mann-Whitney U across conditions
    (independent groups).

    Parameters
    ----------
    merged : pandas.DataFrame
        Output of `build_merged_table`.

    Returns
    -------
    str
        A complete `<table>...</table>` HTML fragment.
    """
    subtype_cols = [
        ("Dopaminergic", "Dopaminergic_pct_of_neurons"),
        ("GABAergic", "GABAergic_pct_of_neurons"),
        ("MEG3", "MEG3 Interneurons_pct_of_neurons"),
    ]
    groups = [
        {"label": f"{cond} {name}", "cond": cond, "values": merged.loc[merged["Condition"] == cond, col]}
        for cond in ("Control", "PD")
        for name, col in subtype_cols
    ]

    n = len(groups)
    pmat = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            gi, gj = groups[i], groups[j]
            if gi["cond"] == gj["cond"]:
                _, p = stats.wilcoxon(gi["values"], gj["values"])
            else:
                _, p = stats.mannwhitneyu(gi["values"], gj["values"], alternative="two-sided")
            pmat[i][j] = pmat[j][i] = p

    header = "".join(f"<th>{g['label']}</th>" for g in groups)
    body_rows = []
    for i, gi in enumerate(groups):
        cells = []
        for j in range(n):
            if i == j:
                cells.append('<td class="diag">—</td>')
            else:
                p = pmat[i][j]
                cls = ' class="sig"' if p < 0.05 else ""
                cells.append(f"<td{cls}>{format_p(p)}</td>")
        body_rows.append(f"<tr><th>{gi['label']}</th>{''.join(cells)}</tr>")

    return (
        '<table class="pmatrix">\n'
        f"  <thead><tr><th></th>{header}</tr></thead>\n"
        "  <tbody>\n    " + "\n    ".join(body_rows) + "\n  </tbody>\n"
        "</table>"
    )


def render_html(merged: pd.DataFrame, computed: dict) -> str:
    """
    Fill the report template with the computed statistics and per-sample data.

    Parameters
    ----------
    merged : pandas.DataFrame
        Output of `build_merged_table`.
    computed : dict
        Output of `compute_stats`.

    Returns
    -------
    str
        The complete, self-contained HTML document.
    """
    html = _HTML_TEMPLATE
    tokens = {
        "@@DATA_JSON@@": build_data_json(merged),
        "@@N_CONTROL@@": str(computed["n_control"]),
        "@@N_PD@@": str(computed["n_pd"]),
        "@@MEAN_DA_PCTNEURONS_CONTROL@@": f"{computed['mean_da_pctneurons_control']:.1f}",
        "@@MEAN_DA_PCTNEURONS_PD@@": f"{computed['mean_da_pctneurons_pd']:.1f}",
        "@@P_MW_DA_PCTNEURONS@@": format_p(computed["p_mw_da_pctneurons"]),
        "@@MEAN_DA_RAW_PD@@": f"{computed['mean_da_raw_pd']:.1f}",
        "@@MEAN_DA_RAW_CONTROL@@": f"{computed['mean_da_raw_control']:.1f}",
        "@@P_MW_DA_RAW@@": format_p(computed["p_mw_da_raw"]),
        "@@FANS_SAMPLES@@": ", ".join(computed["fans_sorted_samples"]),
        "@@MEAN_DA_EXCL_PD@@": f"{computed['mean_da_excl_pd']:.1f}",
        "@@MEAN_DA_EXCL_CONTROL@@": f"{computed['mean_da_excl_control']:.1f}",
        "@@P_MW_DA_EXCL@@": format_p(computed["p_mw_da_excl"]),
        "@@MEG3_RANGE@@": computed["meg3_range"],
        "@@DA_RANGE@@": computed["da_range"],
        "@@GABA_RANGE@@": computed["gaba_range"],
        "@@P_FRIEDMAN@@": format_p(computed["p_friedman_control"]),
        "@@P_WILCOXON_DA_GABA@@": format_p(computed["p_wilcoxon_da_gaba"]),
        "@@P_WILCOXON_GABA_MEG3@@": format_p(computed["p_wilcoxon_gaba_meg3"]),
        "@@P_WILCOXON_DA_MEG3@@": format_p(computed["p_wilcoxon_da_meg3"]),
        "@@TOP_DA_SAMPLES@@": ", ".join(computed["top_da_samples"]),
        "@@PVALUE_MATRIX_TABLE@@": build_pvalue_matrix_table(merged),
    }
    for token, value in tokens.items():
        html = html.replace(token, value)
    return html


_HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Dopaminergic, GABAergic &amp; MEG3 neurons — Substantia Nigra, PD vs Control</title>
<style>
.viz-root {
  color-scheme: light;
  --surface-1:      #fcfcfb;
  --page:           #f9f9f7;
  --text-primary:   #0b0b0b;
  --text-secondary: #52514e;
  --text-muted:     #898781;
  --grid:           #e1e0d9;
  --axis:           #c3c2b7;
  --border:         rgba(11,11,11,0.10);
  --series-1:       #2a78d6; /* blue */
  --series-2:       #eb6834; /* orange */
  --series-3:       #1baf7a; /* aqua */
}
@media (prefers-color-scheme: dark) {
  :root:where(:not([data-theme="light"])) .viz-root {
    color-scheme: dark;
    --surface-1:      #1a1a19;
    --page:           #0d0d0d;
    --text-primary:   #ffffff;
    --text-secondary: #c3c2b7;
    --text-muted:     #898781;
    --grid:           #2c2c2a;
    --axis:           #383835;
    --border:         rgba(255,255,255,0.10);
    --series-1:       #3987e5;
    --series-2:       #d95926;
    --series-3:       #199e70;
  }
}
:root[data-theme="dark"] .viz-root {
  color-scheme: dark;
  --surface-1:      #1a1a19;
  --page:           #0d0d0d;
  --text-primary:   #ffffff;
  --text-secondary: #c3c2b7;
  --text-muted:     #898781;
  --grid:           #2c2c2a;
  --axis:           #383835;
  --border:         rgba(255,255,255,0.10);
  --series-1:       #3987e5;
  --series-2:       #d95926;
  --series-3:       #199e70;
}

* { box-sizing: border-box; }
body { margin: 0; }
.viz-root {
  font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
  background: var(--page);
  color: var(--text-primary);
  padding: 32px 20px 64px;
}
.wrap { max-width: 980px; margin: 0 auto; }
h1 { font-size: 21px; font-weight: 650; margin: 0 0 4px; }
.subtitle { color: var(--text-secondary); font-size: 13.5px; margin: 0 0 28px; line-height: 1.5; }
.card {
  background: var(--surface-1);
  border: 1px solid var(--border);
  border-radius: 14px;
  padding: 22px 24px 18px;
  margin-bottom: 22px;
  overflow-x: auto;
}
.card h2 { font-size: 15.5px; font-weight: 650; margin: 0 0 2px; }
.card .finding { font-size: 13px; color: var(--text-secondary); margin: 0 0 14px; line-height: 1.5; max-width: 760px; }
.finding b { color: var(--text-primary); }
.legend { display: flex; gap: 18px; flex-wrap: wrap; font-size: 12.5px; color: var(--text-secondary); margin: 2px 0 10px; }
.legend .item { display: flex; align-items: center; gap: 6px; }
.legend .swatch { width: 11px; height: 11px; border-radius: 3px; flex: none; }
.legend .swatch.circle { border-radius: 50%; }
.legend .swatch.square { border-radius: 2px; }
.legend .swatch.triangle { border-radius: 0; clip-path: polygon(50% 0%, 0% 100%, 100% 100%); }
.legend .swatch.diamond { border-radius: 0; clip-path: polygon(50% 0%, 100% 50%, 50% 100%, 0% 50%); }
.legend .swatch.triangle-down { border-radius: 0; clip-path: polygon(0% 0%, 100% 0%, 50% 100%); }
.legend .swatch.pentagon { border-radius: 0; clip-path: polygon(50% 0%, 100% 38%, 82% 100%, 18% 100%, 0% 38%); }
.legend .swatch.star { border-radius: 0; clip-path: polygon(50% 0%, 61% 35%, 98% 35%, 68% 57%, 79% 91%, 50% 70%, 21% 91%, 32% 57%, 2% 35%, 39% 35%); }
.legend .swatch.cross { border-radius: 0; clip-path: polygon(35% 0%, 65% 0%, 65% 35%, 100% 35%, 100% 65%, 65% 65%, 65% 100%, 35% 100%, 35% 65%, 0% 65%, 0% 35%, 35% 35%); }
.legend.shape-legend .swatch { background: var(--text-secondary) !important; }
.shape-control { display: flex; align-items: center; gap: 8px; font-size: 12.5px; color: var(--text-secondary); margin: -4px 0 8px; }
.shape-control select {
  font: inherit;
  color: var(--text-primary);
  background: var(--surface-1);
  border: 1px solid var(--border);
  border-radius: 8px;
  padding: 4px 8px;
}
.shape-note { font-size: 11.5px; color: var(--text-muted); font-style: italic; margin: -2px 0 10px; }
svg text { fill: var(--text-muted); font-size: 11.5px; }
svg .axis-line { stroke: var(--axis); stroke-width: 1; }
svg .grid-line { stroke: var(--grid); stroke-width: 1; }
svg .median-tick { stroke-width: 2.5; stroke-linecap: round; }
svg .panel-label { fill: var(--text-primary); font-size: 12.5px; font-weight: 650; }
svg .n-label { fill: var(--text-muted); font-size: 11px; }
.pt { cursor: pointer; stroke: var(--surface-1); stroke-width: 2; }
.tooltip {
  position: fixed;
  pointer-events: none;
  background: var(--text-primary);
  color: var(--surface-1);
  font-size: 12px;
  line-height: 1.45;
  padding: 7px 10px;
  border-radius: 8px;
  opacity: 0;
  transform: translate(-50%, -100%);
  transition: opacity 0.08s ease;
  z-index: 10;
  white-space: nowrap;
}
.tooltip.show { opacity: 1; }
.caveat {
  font-size: 12.5px;
  color: var(--text-secondary);
  line-height: 1.6;
  background: var(--surface-1);
  border: 1px solid var(--border);
  border-radius: 14px;
  padding: 16px 20px;
}
.caveat h2 { font-size: 14px; margin: 0 0 8px; color: var(--text-primary); }
.caveat ul { margin: 6px 0 0; padding-left: 18px; }
.caveat li { margin-bottom: 6px; }
details { margin-top: 10px; }
summary { cursor: pointer; font-size: 12.5px; color: var(--text-secondary); }
table { border-collapse: collapse; width: 100%; font-size: 11.5px; margin-top: 10px; }
th, td { text-align: right; padding: 5px 8px; border-bottom: 1px solid var(--border); white-space: nowrap; }
th:first-child, td:first-child { text-align: left; }
th { color: var(--text-muted); font-weight: 600; }
td { color: var(--text-secondary); font-variant-numeric: tabular-nums; }
tr.pd td { color: var(--text-primary); }
table.pmatrix { margin-top: 14px; }
table.pmatrix td.diag { text-align: center; color: var(--text-muted); }
table.pmatrix td.sig { color: var(--text-primary); font-weight: 650; }
.pmatrix-note { font-size: 11.5px; color: var(--text-muted); margin: 8px 0 0; }
</style>
</head>
<body>

<div class="viz-root">
<div class="wrap">
  <h1>Dopaminergic, GABAergic &amp; MEG3 neurons — Substantia Nigra</h1>
  <p class="subtitle">Per-sample snRNA-seq nucleus counts, Table S1 tab C (Ma et&nbsp;al. 2025, <i>Nat Commun</i>, 41467_2025_62478_MOESM3_ESM.xlsx). Condition assigned from tab A metadata by TB&nbsp;#. Substantia Nigra samples: @@N_CONTROL@@ Control, @@N_PD@@ PD.</p>

  <div class="card">
    <h2>1. Is the dopaminergic fraction overrepresented in PD vs Control?</h2>
    <p class="finding">No. Dopaminergic nuclei as a share of all captured neurons (DA / [DA+GABA+MEG3]) are statistically indistinguishable between groups — <b>Control @@MEAN_DA_PCTNEURONS_CONTROL@@% vs PD @@MEAN_DA_PCTNEURONS_PD@@%</b> (mean), Mann–Whitney <b>p = @@P_MW_DA_PCTNEURONS@@</b>. Raw DA share of <i>all</i> nuclei looks nominally higher in PD (@@MEAN_DA_RAW_PD@@% vs @@MEAN_DA_RAW_CONTROL@@%, p = @@P_MW_DA_RAW@@), but that's an artifact of samples (@@FANS_SAMPLES@@) that were apparently FANS/neuron-sorted (zero non-neuronal nuclei) — excluding them gives @@MEAN_DA_EXCL_PD@@% PD vs @@MEAN_DA_EXCL_CONTROL@@% Control, p = @@P_MW_DA_EXCL@@.</p>

    <h2>2. Within a group, is one neuron subtype far more abundant than the others?</h2>
    <p class="finding">Yes, in both groups. <b>MEG3 interneurons dominate</b> (~@@MEG3_RANGE@@ of neuronal nuclei), dopaminergic neurons are second (~@@DA_RANGE@@), and <b>GABAergic neurons are consistently the minority</b> (~@@GABA_RANGE@@). In Control, a Friedman test across all three subtypes is significant (p = @@P_FRIEDMAN@@); pairwise, GABAergic is significantly lower than both Dopaminergic (Wilcoxon p = @@P_WILCOXON_DA_GABA@@) and MEG3 (p = @@P_WILCOXON_GABA_MEG3@@), while MEG3 vs Dopaminergic doesn't reach significance (p = @@P_WILCOXON_DA_MEG3@@).</p>
  </div>

  <div class="card">
    <div class="legend">
      <span class="item"><span class="swatch" style="background:var(--series-1)"></span>Dopaminergic</span>
      <span class="item"><span class="swatch" style="background:var(--series-2)"></span>GABAergic</span>
      <span class="item"><span class="swatch" style="background:var(--series-3)"></span>MEG3 Interneurons</span>
    </div>
    <p style="font-size:12px;color:var(--text-muted);margin:-4px 0 6px">Color = neuron subtype. Point <b>shape</b> additionally encodes a tab A column of your choice, applied identically to both facets:</p>
    <div class="shape-control">
      <label for="shapeCategory">Shape =</label>
      <select id="shapeCategory">
        <option value="none">None (circles)</option>
        <option value="Neuropath_Diagnosis" selected>Neuropath_Diagnosis</option>
        <option value="Age">Age</option>
        <option value="Sex">Sex</option>
        <option value="path_braak_nft">path_braak_nft (Braak NFT stage)</option>
        <option value="path_braak_asyn">path_braak_asyn (Braak &alpha;-syn stage)</option>
        <option value="path_cerad">path_cerad (CERAD score)</option>
      </select>
    </div>
    <div class="legend shape-legend" id="shapeLegend"></div>
    <p class="shape-note" id="shapeNote" style="display:none"></p>
    <svg id="chartB" viewBox="0 0 900 320" width="100%" style="max-width:900px;display:block"></svg>
  </div>

  <div class="card">
    <h2>Pairwise p-values — 3 subtypes &times; Control/PD</h2>
    <div style="overflow-x:auto">@@PVALUE_MATRIX_TABLE@@</div>
    <p class="pmatrix-note">Same-condition pairs: paired Wilcoxon signed-rank (same samples, 3 compositions each). Cross-condition pairs: Mann&ndash;Whitney U (independent samples). Bold = p &lt; 0.05.</p>
  </div>

  <div class="caveat">
    <h2>Notes &amp; caveats</h2>
    <ul>
      <li>Values are shares of <b>captured nuclei</b> in snRNA-seq, not true tissue density or cell survival — dissociation/capture efficiency can differ by sample and doesn't equal in-situ neuron counts (the classic PD finding is nigral DA neuron <i>loss</i>, assessed histologically, not by nucleus-capture fraction).</li>
      <li>Sample(s) @@FANS_SAMPLES@@ have zero non-neuronal nuclei recorded — consistent with neuron-enriched (FANS) sorting rather than whole-tissue snRNA-seq — and are flagged/excluded in the sensitivity check above.</li>
      <li>Group sizes are small (n=@@N_CONTROL@@ Control, n=@@N_PD@@ PD); individual high-DA samples (e.g. @@TOP_DA_SAMPLES@@) pull the PD/Control means but the group medians and rank-based tests tell the same story.</li>
    </ul>
    <details>
      <summary>Show per-sample data (all Substantia Nigra samples)</summary>
      <div style="overflow-x:auto"><table id="dataTable"></table></div>
    </details>
  </div>
</div>
</div>

<div class="tooltip" id="tip"></div>

<script>
const DATA = @@DATA_JSON@@;

const css = getComputedStyle(document.querySelector('.viz-root'));
const col = (name) => css.getPropertyValue(name).trim();

const tip = document.getElementById('tip');
function showTip(evt, html) {
  tip.innerHTML = html;
  tip.style.left = evt.clientX + 'px';
  tip.style.top = (evt.clientY - 10) + 'px';
  tip.classList.add('show');
}
function hideTip() { tip.classList.remove('show'); }

function median(arr) {
  const s = [...arr].sort((a,b)=>a-b);
  const n = s.length;
  const mid = Math.floor(n/2);
  return n % 2 ? s[mid] : (s[mid-1]+s[mid])/2;
}

function svgEl(tag, attrs) {
  const el = document.createElementNS('http://www.w3.org/2000/svg', tag);
  for (const k in attrs) el.setAttribute(k, attrs[k]);
  return el;
}

// jitter, seeded per index so repeat renders are stable within a session
function jitter(seed) {
  const x = Math.sin(seed * 999.7) * 43758.5453;
  return x - Math.floor(x);
}

// Secondary, color-independent encoding; a 9th+ distinct value shares the last ('star') as "Other".
const SHAPE_PALETTE = ['circle', 'square', 'triangle', 'diamond', 'triangle-down', 'cross', 'pentagon', 'star'];

function regularPolygonPoints(cx, cy, r, sides, rotationDeg) {
  const pts = [];
  for (let i = 0; i < sides; i++) {
    const a = (Math.PI / 180) * (rotationDeg + (360 / sides) * i);
    pts.push(`${cx + r * Math.cos(a)},${cy + r * Math.sin(a)}`);
  }
  return pts.join(' ');
}

function starPoints(cx, cy, rOuter, rInner, points, rotationDeg) {
  const pts = [];
  for (let i = 0; i < points * 2; i++) {
    const r = i % 2 === 0 ? rOuter : rInner;
    const a = (Math.PI / 180) * (rotationDeg + (180 / points) * i);
    pts.push(`${cx + r * Math.cos(a)},${cy + r * Math.sin(a)}`);
  }
  return pts.join(' ');
}

function crossPoints(cx, cy, armLen, armHalfWidth) {
  const o = armLen, w = armHalfWidth;
  return [
    [w, -o], [w, -w], [o, -w], [o, w], [w, w], [w, o],
    [-w, o], [-w, w], [-o, w], [-o, -w], [-w, -w], [-w, -o],
  ].map(([dx, dy]) => `${cx + dx},${cy + dy}`).join(' ');
}

function drawMark(shape, cx, cy, color) {
  const s = 5; // same visual weight as a r=5 circle
  if (shape === 'square') {
    return svgEl('rect', {class:'pt', x:cx-s, y:cy-s, width:s*2, height:s*2, rx:1.5, fill:color});
  }
  if (shape === 'triangle') {
    const h = s*2.1;
    const pts = `${cx},${cy-h*0.62} ${cx-h*0.58},${cy+h*0.42} ${cx+h*0.58},${cy+h*0.42}`;
    return svgEl('polygon', {class:'pt', points:pts, fill:color});
  }
  if (shape === 'triangle-down') {
    const h = s*2.1;
    const pts = `${cx-h*0.58},${cy-h*0.42} ${cx+h*0.58},${cy-h*0.42} ${cx},${cy+h*0.62}`;
    return svgEl('polygon', {class:'pt', points:pts, fill:color});
  }
  if (shape === 'diamond') {
    const h = s*1.5;
    const pts = `${cx},${cy-h} ${cx+h},${cy} ${cx},${cy+h} ${cx-h},${cy}`;
    return svgEl('polygon', {class:'pt', points:pts, fill:color});
  }
  if (shape === 'pentagon') {
    return svgEl('polygon', {class:'pt', points: regularPolygonPoints(cx, cy, s*1.3, 5, -90), fill:color});
  }
  if (shape === 'star') {
    return svgEl('polygon', {class:'pt', points: starPoints(cx, cy, s*1.5, s*0.58, 5, -90), fill:color});
  }
  if (shape === 'cross') {
    return svgEl('polygon', {class:'pt', points: crossPoints(cx, cy, s*1.3, s*0.42), fill:color});
  }
  return svgEl('circle', {class:'pt', cx, cy, r:s, fill:color});
}

function drawStrip(svg, groups, opts) {
  // groups: [{label, color, values: [{v, tbNum}], x}]
  const { x0, x1, yTop, yBot, yMax, title } = opts;
  const colW = (x1 - x0) / groups.length;

  // y gridlines + axis
  const ticks = yMax <= 80 ? [0,20,40,60,80] : [0,25,50,75,100];
  ticks.filter(t=>t<=yMax+1).forEach(t => {
    const y = yBot - (t/yMax)*(yBot-yTop);
    svg.appendChild(svgEl('line', {class:'grid-line', x1: x0-6, x2: x1, y1:y, y2:y}));
    const lbl = svgEl('text', {x: x0-10, y: y+4, 'text-anchor':'end'});
    lbl.textContent = t + '%';
    svg.appendChild(lbl);
  });
  svg.appendChild(svgEl('line', {class:'axis-line', x1:x0, x2:x1, y1:yBot, y2:yBot}));

  groups.forEach((g, gi) => {
    const cx = x0 + colW*gi + colW/2;
    const vals = g.values.map(d=>d.v);
    const med = median(vals);
    const medY = yBot - (med/yMax)*(yBot-yTop);

    // points
    g.values.forEach((d, i) => {
      const jx = (jitter(gi*97 + i*13 + 1) - 0.5) * (colW*0.55);
      const py = yBot - (d.v/yMax)*(yBot-yTop);
      const c = drawMark(d.shape || g.shape || 'circle', cx+jx, py, g.color);
      const tipExtra = d.diag ? `<br><span style="opacity:0.75">${d.diag}</span>` : '';
      c.addEventListener('mouseenter', (e)=> showTip(e, `<b>${d.tbNum}</b><br>${g.label}: ${d.v.toFixed(1)}%${tipExtra}`));
      c.addEventListener('mousemove', (e)=> { tip.style.left = e.clientX+'px'; tip.style.top=(e.clientY-10)+'px'; });
      c.addEventListener('mouseleave', hideTip);
      svg.appendChild(c);
    });

    // median tick
    svg.appendChild(svgEl('line', {class:'median-tick', x1:cx-colW*0.32, x2:cx+colW*0.32, y1:medY, y2:medY, stroke:g.color}));

    // x label
    const lbl = svgEl('text', {x:cx, y:yBot+18, 'text-anchor':'middle'});
    lbl.textContent = g.label;
    svg.appendChild(lbl);
    if (g.n !== undefined) {
      const nl = svgEl('text', {class:'n-label', x:cx, y:yBot+32, 'text-anchor':'middle'});
      nl.textContent = 'n=' + g.n;
      svg.appendChild(nl);
    }
  });

  if (title) {
    const t = svgEl('text', {class:'panel-label', x:(x0+x1)/2, y: yTop-26, 'text-anchor':'middle'});
    t.textContent = title;
    svg.appendChild(t);
  }
}

// ---- Chart B: subtype composition, faceted by Control/PD, shape category applied identically to both ----
(function(){
  const svg = document.getElementById('chartB');
  const legendEl = document.getElementById('shapeLegend');
  const noteEl = document.getElementById('shapeNote');
  const selectEl = document.getElementById('shapeCategory');

  const subtypes = [
    {key:'Dopaminergic_pct_of_neurons', label:'Dopaminergic', color: col('--series-1')},
    {key:'GABAergic_pct_of_neurons', label:'GABAergic', color: col('--series-2')},
    {key:'MEG3 Interneurons_pct_of_neurons', label:'MEG3', color: col('--series-3')},
  ];

  // Fixed display order for columns where sorted() would be wrong (Roman
  // numerals, "Control, ..." should lead, age bins aren't alphabetical).
  const FIXED_ORDER = {
    path_braak_nft: ['0', 'I', 'II', 'III', 'IV', 'V', 'VI', 'NA'],
    path_cerad: ['A', 'B', 'C', 'Not eligible', 'NA'],
    Sex: ['F', 'M'],
  };

  function categoryValues(key) {
    const vals = Array.from(new Set(DATA.map(d => String(d[key]))));
    const order = FIXED_ORDER[key];
    if (order) {
      return vals.sort((a, b) => {
        const ia = order.indexOf(a), ib = order.indexOf(b);
        return (ia < 0 ? 999 : ia) - (ib < 0 ? 999 : ib);
      });
    }
    if (key === 'Age') {
      const lo = (s) => s.startsWith('>') ? 999 : parseInt(s, 10);
      return vals.sort((a, b) => lo(a) - lo(b));
    }
    if (key === 'path_braak_asyn') {
      const n = (s) => isNaN(parseFloat(s)) ? 999 : parseFloat(s);
      return vals.sort((a, b) => n(a) - n(b));
    }
    if (key === 'Neuropath_Diagnosis') {
      return vals.sort((a, b) => {
        const ca = a.startsWith('Control') ? 0 : 1, cb = b.startsWith('Control') ? 0 : 1;
        return ca !== cb ? ca - cb : a.localeCompare(b);
      });
    }
    return vals.sort();
  }

  function buildShapeMap(key) {
    if (key === 'none') return { map: {}, vals: [], overflow: 0, directLimit: 0 };
    const vals = categoryValues(key);
    const hasOverflow = vals.length > SHAPE_PALETTE.length;
    const directLimit = hasOverflow ? SHAPE_PALETTE.length - 1 : vals.length;
    const map = {};
    vals.forEach((v, i) => {
      map[v] = i < directLimit ? SHAPE_PALETTE[i] : SHAPE_PALETTE[SHAPE_PALETTE.length - 1];
    });
    return { map, vals, overflow: hasOverflow ? vals.length - directLimit : 0, directLimit };
  }

  function facetGroups(cond, shapeKey, shapeInfo) {
    const rows = DATA.filter(d => d.Condition === cond);
    return subtypes.map(s => ({
      label: s.label, color: s.color, n: rows.length,
      values: rows.map(d => {
        const raw = shapeKey === 'none' ? null : String(d[shapeKey]);
        return {
          v: d[s.key], tbNum: d['TB #'],
          diag: raw !== null ? `${shapeKey}: ${raw}` : undefined,
          shape: raw !== null ? (shapeInfo.map[raw] || 'circle') : 'circle',
        };
      })
    }));
  }

  function renderShapeLegend(shapeKey, shapeInfo) {
    legendEl.innerHTML = '';
    if (shapeKey === 'none') { noteEl.style.display = 'none'; return; }

    function addItem(shape, label) {
      const item = document.createElement('span');
      item.className = 'item';
      const sw = document.createElement('span');
      sw.className = `swatch ${shape}`;
      item.appendChild(sw);
      item.appendChild(document.createTextNode(label));
      legendEl.appendChild(item);
    }

    shapeInfo.vals.slice(0, shapeInfo.directLimit).forEach((v, i) => addItem(SHAPE_PALETTE[i], v));

    if (shapeInfo.overflow > 0) {
      addItem(SHAPE_PALETTE[SHAPE_PALETTE.length - 1], 'Other');
      noteEl.style.display = 'block';
      noteEl.textContent = `${shapeInfo.overflow} additional ${shapeKey} value(s) beyond the ${shapeInfo.directLimit} shown share the "Other" marker.`;
    } else {
      noteEl.style.display = 'none';
    }
  }

  function render(shapeKey) {
    svg.innerHTML = '';
    const shapeInfo = buildShapeMap(shapeKey);
    renderShapeLegend(shapeKey, shapeInfo);

    drawStrip(svg, facetGroups('Control', shapeKey, shapeInfo), {
      x0: 60, x1: 420, yTop: 56, yBot: 260, yMax: 100,
      title: 'Control (n=@@N_CONTROL@@)',
    });
    drawStrip(svg, facetGroups('PD', shapeKey, shapeInfo), {
      x0: 500, x1: 860, yTop: 56, yBot: 260, yMax: 100,
      title: 'PD (n=@@N_PD@@)',
    });
    svg.appendChild(svgEl('line', {x1:460, x2:460, y1:40, y2:276, stroke: col('--border'), 'stroke-width':1}));
  }

  selectEl.addEventListener('change', () => render(selectEl.value));
  render(selectEl.value);
})();

// ---- data table ----
(function(){
  const t = document.getElementById('dataTable');
  const head = ['TB #','Condition','Neuropath_Diagnosis','Age','Sex','path_braak_nft','path_braak_asyn','path_cerad','Dopaminergic','GABAergic','MEG3 Interneurons','Total neurons','DA % of neurons','GABA % of neurons','MEG3 % of neurons'];
  const thead = document.createElement('thead');
  const htr = document.createElement('tr');
  head.forEach(h=>{ const th=document.createElement('th'); th.textContent=h; htr.appendChild(th); });
  thead.appendChild(htr);
  t.appendChild(thead);
  const tbody = document.createElement('tbody');
  [...DATA].sort((a,b)=> a.Condition===b.Condition ? a['TB #'].localeCompare(b['TB #']) : (a.Condition==='PD'?1:-1))
    .forEach(d=>{
      const tr = document.createElement('tr');
      if (d.Condition==='PD') tr.className='pd';
      const cells = [d['TB #'], d.Condition, d.Neuropath_Diagnosis, d.Age, d.Sex, d.path_braak_nft, d.path_braak_asyn, d.path_cerad,
        d.Dopaminergic, d.GABAergic, d['MEG3 Interneurons'],
        d.total_neurons, d.Dopaminergic_pct_of_neurons.toFixed(1)+'%', d.GABAergic_pct_of_neurons.toFixed(1)+'%', d['MEG3 Interneurons_pct_of_neurons'].toFixed(1)+'%'];
      cells.forEach(v=>{ const td=document.createElement('td'); td.textContent=v; tr.appendChild(td); });
      tbody.appendChild(tr);
    });
  t.appendChild(tbody);
})();
</script>

</body>
</html>
"""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(args: argparse.Namespace) -> None:
    """Build the Substantia Nigra composition report and write it to `args.html_out`."""
    merged = build_merged_table(args.xlsx_path)
    print(f"Loaded {len(merged)} Substantia Nigra samples from {args.xlsx_path}")
    print(merged["Condition"].value_counts().to_string())

    computed = compute_stats(merged)
    print("\nComputed statistics:")
    for key, value in computed.items():
        print(f"  {key}: {value}")

    html = render_html(merged, computed)

    os.makedirs(os.path.dirname(args.html_out), exist_ok=True)
    with open(args.html_out, "w", encoding="utf-8") as fh:
        fh.write(html)
    print(f"\nWritten {args.html_out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--xlsx-path",
        default=None,
        help="Ma et al. 2025 supplementary workbook. A relative path is resolved against\n"
             "the directory this command is run from (not the repo root).\n"
             f"[Default: {_SUPP_XLSX}, relative to the repo root]",
    )
    parser.add_argument(
        "--html-out",
        default=None,
        help="Output path for the rendered HTML report. A relative path is resolved\n"
             "against the directory this command is run from (not the repo root).\n"
             f"[Default: {_OUTPUT_HTML}, relative to the repo root]",
    )
    args = parser.parse_args()

    # Resolve any user-supplied path against the invocation cwd *before* the
    # chdir below -- otherwise a relative --xlsx-path/--html-out silently
    # gets reinterpreted as relative to the repo root instead.
    invocation_cwd = Path.cwd()
    if args.xlsx_path is not None:
        args.xlsx_path = str((invocation_cwd / args.xlsx_path).resolve())
    if args.html_out is not None:
        args.html_out = str((invocation_cwd / args.html_out).resolve())

    repo_root = next(p for p in Path(__file__).resolve().parents if (p / ".git").exists())
    if Path.cwd() != repo_root:
        print(f"[INFO] Changing working directory to repo root: {repo_root}")
        os.chdir(repo_root)

    # The defaults are repo-root-relative and only resolve correctly once
    # the chdir above has happened, so fill them in afterwards.
    if args.xlsx_path is None:
        args.xlsx_path = _SUPP_XLSX
    if args.html_out is None:
        args.html_out = _OUTPUT_HTML

    log_run_command(args)
    main(args)
