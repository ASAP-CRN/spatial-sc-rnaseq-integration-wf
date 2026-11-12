"""
Shared Jupyter-notebook utilities for spatial-sc-rnaseq-integration-wf.

Extracted from notebooks/explore_integration_kamath_random.ipynb and
notebooks/visualize_deconvolution.ipynb, where near-identical versions of
these had been duplicated. Where the two notebooks disagreed,
visualize_deconvolution.ipynb's version was kept as the default (see
`hierarchical_leaf_order`'s `normalize_rows` for the one behavioural
difference that was preserved as an opt-in rather than silently dropped).

NOTE on bootstrapping: a notebook can't `import notebook_utils` before it
knows where this file is, so each notebook still needs a few lines of inline
code to locate and `sys.path.insert` the `utils/` directory. But that search
only needs to find `utils/notebook_utils.py` itself -- a plain ancestor walk,
since this file is mirrored into both the live repo and a claude_outputs
copy alike. It does NOT need the live-repo-specific `.git`/sibling-name logic
that `get_repo_root` below uses, which is why `get_repo_root` can live here:
by the time a notebook calls it, this module is already imported.
"""

import colorsys
import importlib.util
import subprocess

import ipywidgets as widgets
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import scanpy as sc
import squidpy as sq
from IPython.display import HTML, display
from pathlib import Path
from plotly.subplots import make_subplots
from scipy.cluster.hierarchy import leaves_list, linkage

GOLDEN_ANGLE = 0.6180339887498949  # golden ratio conjugate


def get_repo_root(repo_name):
    """
    Find a repo's root directory by walking up from the current working
    directory, checking both each ancestor and `ancestor / repo_name` for a
    `.git` directory.

    Parameters
    ----------
    repo_name : str
        Name of the repo directory to resolve (e.g. "spatial-sc-rnaseq-integration-wf").
        Needed because this can be called from a non-git mirror (e.g. a
        claude_outputs copy), which has no `.git` of its own but typically
        sits under the same parent directory as a real checkout of `repo_name`.

    Returns
    -------
    pathlib.Path
        Absolute path to the repo root (the directory containing `.git`).
    """
    return next(
        p for base in (Path.cwd(), *Path.cwd().parents)
        for p in (base, base / repo_name)
        if (p / ".git").exists()
    )


def ensure_cached(gcs_dir, local_dir, filename):
    """
    Return the local cache path for `filename`, downloading it from `gcs_dir`
    first if not already cached under `local_dir`.

    Parameters
    ----------
    gcs_dir : str
        GCS directory `filename` lives under.
    local_dir : str
        Local cache directory to check/write `filename` under.
    filename : str
        Basename of the file under `gcs_dir` (e.g. a `*.h5ad`).

    Returns
    -------
    pathlib.Path
        Local path to the (now-guaranteed-present) cached file.
    """
    local_path = Path(local_dir, filename)
    if not local_path.exists():
        local_path.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            ["gcloud", "storage", "cp", f"{gcs_dir}/{filename}", str(local_path)],
            check=True,
        )
    return local_path


def list_available_and_cached(gcs_dir, local_dir, suffix):
    """
    List sample names found under a GCS directory, partitioned by whether
    each is already cached locally.

    Parameters
    ----------
    gcs_dir : str
        GCS directory to list (e.g. a `*_DATASET_GCS_DIR`).
    local_dir : str
        Local cache directory to check each sample against.
    suffix : str
        Filename suffix shared by every sample under `gcs_dir` (e.g.
        ".cleaned_unfiltered.h5ad"), stripped to recover each sample's name.

    Returns
    -------
    tuple of list of str
        available : every sample name found under `gcs_dir`, sorted.
        cached : the subset of `available` already present under `local_dir`.
        missing : the rest of `available`.
    """
    ls = subprocess.run(
        ["gcloud", "storage", "ls", f"{gcs_dir}/*{suffix}"],
        capture_output=True, text=True,
    )
    available = sorted(
        Path(line.strip()).name.removesuffix(suffix)
        for line in ls.stdout.splitlines() if line.strip()
    )
    cached = [s for s in available if Path(local_dir, s + suffix).exists()]
    missing = [s for s in available if s not in set(cached)]
    return available, cached, missing

BIG_DROPDOWN_CSS = """
<style>
.big-dropdown select {
    font-size: 26px !important;
    color: #1a1a1a !important;
    background-color: #fff3b0 !important;
    border: 2px solid #4a90d9 !important;
}
.big-dropdown .widget-label {
    font-size: 26px !important;
}
</style>
"""


def inject_big_dropdown_css():
    """Display the CSS for the "big-dropdown" widget class; call once per notebook."""
    display(HTML(BIG_DROPDOWN_CSS))


def dropdown(big=False, **kwargs):
    """
    An `ipywidgets.Dropdown`, optionally pre-sized and styled via the "big-dropdown" CSS class.

    Every dropdown in these notebooks sets `style={"description_width":
    "initial"}` so the description label isn't truncated -- this applies
    that default either way, so there's no need for a bare `widgets.Dropdown(...)`.

    Parameters
    ----------
    big : bool
        Apply the bigger, higher-contrast "big-dropdown" look (wider layout
        + the CSS class). Needs `inject_big_dropdown_css()` called once
        somewhere earlier in the notebook for the class to have any effect.
    **kwargs
        Forwarded to `ipywidgets.Dropdown` (e.g. `options`, `description`,
        `layout`). `style` (and `layout`, if `big`) default if not given.

    Returns
    -------
    ipywidgets.Dropdown
    """
    kwargs.setdefault("style", {"description_width": "initial"})
    if big:
        kwargs.setdefault("layout", widgets.Layout(width="700px", height="56px"))
    widget = widgets.Dropdown(**kwargs)
    if big:
        widget.add_class("big-dropdown")
    return widget


def golden_angle_palette(categories, saturation=0.85, value=0.75):
    """
    Assign each category a maximally-distinct hue, stepped by the golden angle.

    Stepping hues by the golden angle (rather than evenly by 1/n) means
    alphabetically (and therefore often visually/positionally) adjacent
    categories never land on similar hues.

    Parameters
    ----------
    categories : iterable of str
        Category labels to assign colours to; order determines hue assignment.
    saturation : float
        HSV saturation, passed to `colorsys.hsv_to_rgb`.
    value : float
        HSV value/brightness, passed to `colorsys.hsv_to_rgb`.

    Returns
    -------
    dict
        category : tuple
            (r, g, b), each channel in [0, 1].
    """
    return {
        category: colorsys.hsv_to_rgb((i * GOLDEN_ANGLE) % 1, saturation, value)
        for i, category in enumerate(categories)
    }


def hierarchical_leaf_order(matrix_df, normalize_rows=False):
    """
    Leaf order for a matrix's rows via average-linkage hierarchical clustering.

    Parameters
    ----------
    matrix_df : pandas.DataFrame
        Rows to order; clustering runs on these rows' values (or their
        row-sum-normalized profiles, if `normalize_rows`).
    normalize_rows : bool
        Normalize each row to sum to 1 before clustering, so clustering
        reflects each row's value distribution rather than its magnitude.
        Needed for e.g. raw count crosstabs, where row totals can differ by
        orders of magnitude and would otherwise dominate the clustering.

    Returns
    -------
    pandas.Index
        `matrix_df.index`, reordered by clustering leaf order (unchanged,
        in original order, if `matrix_df` has fewer than 3 rows).
    """
    if matrix_df.shape[0] < 3:
        return matrix_df.index
    values = matrix_df.values
    if normalize_rows:
        row_sums = matrix_df.sum(axis=1).replace(0, 1)
        values = matrix_df.div(row_sums, axis=0).values
    link = linkage(values, method="average", metric="euclidean")
    return matrix_df.index[leaves_list(link)]


def crosstab_sort_order(matrix, row_totals, col_totals, sort_by, normalize_rows=False):
    """
    Row/column order for a crosstab-style matrix, by `sort_by`.

    Parameters
    ----------
    matrix : pandas.DataFrame
        Used only for `sort_by="cluster"` (hierarchical leaf order on rows,
        and on `matrix.T` for columns).
    row_totals : pandas.Series
        Per-row totals, used for `sort_by="counts (max to min)"`.
    col_totals : pandas.Series
        Per-column totals, used for `sort_by="counts (max to min)"`.
    sort_by : str
        "alphanumeric", "cluster", or "counts (max to min)" (default).
    normalize_rows : bool
        Forwarded to `hierarchical_leaf_order` for `sort_by="cluster"` --
        True for raw-count crosstabs (row totals vary by orders of
        magnitude), False for a matrix already on a comparable scale (e.g.
        a correlation matrix).

    Returns
    -------
    tuple of pandas.Index
        (row_order, col_order)
    """
    if sort_by == "alphanumeric":
        return matrix.index.sort_values(), matrix.columns.sort_values()
    if sort_by == "cluster":
        return (
            hierarchical_leaf_order(matrix, normalize_rows=normalize_rows),
            hierarchical_leaf_order(matrix.T, normalize_rows=normalize_rows),
        )
    return (
        row_totals.sort_values(ascending=False).index,
        col_totals.sort_values(ascending=False).index,
    )


def crosstab_heatmap_figure(
    z,
    row_order,
    col_order,
    row_marginal,
    col_marginal,
    height,
    width,
    heatmap_kwargs=None,
    row_ticktext=None,
    col_ticktext=None,
    marginal_frac=0.10,
    marginal_color="#6B6470",
    row_axis_reversed=False,
    title=None,
    margin=None,
):
    """
    Heatmap-with-marginals figure: a matrix heatmap plus log-scale row/column
    total bars on shared axes, so rows/columns never drift out of sync.

    Parameters
    ----------
    z : pandas.DataFrame
        Values to plot as the heatmap, already reindexed to `row_order`/`col_order`.
    row_order : pandas.Index
        Row category order (top to bottom, unless `row_axis_reversed`).
    col_order : pandas.Index
        Column category order (left to right).
    row_marginal : pandas.Series
        Per-row totals, indexed like `row_order`, plotted as the left bar panel.
    col_marginal : pandas.Series
        Per-column totals, indexed like `col_order`, plotted as the bottom bar panel.
    height : int
        Figure height in pixels.
    width : int
        Figure width in pixels.
    heatmap_kwargs : dict
        Forwarded to `go.Heatmap` (e.g. `colorscale`, `zmin`, `zmid`, `zmax`,
        `colorbar`, `text`, `texttemplate`, `hovertemplate`).
    row_ticktext : list of str
        Row tick labels (e.g. colourized HTML); defaults to `row_order` as-is.
    col_ticktext : list of str
        Column tick labels; defaults to `col_order` as-is.
    marginal_frac : float
        Fraction of the figure given to each marginal bar panel.
    marginal_color : str
        Bar colour for both marginals.
    row_axis_reversed : bool
        Reverse the shared row axis so the first entry in `row_order` renders
        at the top; otherwise it renders at the bottom (Plotly's default for
        a categorical axis).
    title : str
        Figure title.
    margin : dict
        Forwarded to `fig.update_layout(margin=...)`.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    heatmap_kwargs = heatmap_kwargs or {}
    row_ticktext = row_ticktext if row_ticktext is not None else list(row_order)
    col_ticktext = col_ticktext if col_ticktext is not None else list(col_order)

    fig = make_subplots(
        rows=2, cols=2,
        shared_xaxes=True, shared_yaxes=True,
        row_heights=[1 - marginal_frac, marginal_frac],
        column_widths=[marginal_frac, 1 - marginal_frac],
        horizontal_spacing=0.01, vertical_spacing=0.02,
    )

    fig.add_trace(
        go.Bar(
            x=row_marginal.values, y=row_marginal.index, orientation="h",
            marker_color=marginal_color, customdata=row_marginal.values,
            hovertemplate="%{y}: %{customdata:,.0f}<extra></extra>",
        ),
        row=1, col=1,
    )
    fig.add_trace(
        go.Bar(
            x=col_marginal.index, y=col_marginal.values,
            marker_color=marginal_color, customdata=col_marginal.values,
            hovertemplate="%{x}: %{customdata:,.0f}<extra></extra>",
        ),
        row=2, col=2,
    )
    fig.add_trace(
        go.Heatmap(
            z=z.values, x=list(z.columns), y=list(z.index),
            xgap=1, ygap=1,
            **heatmap_kwargs,
        ),
        row=1, col=2,
    )

    fig.update_xaxes(type="log", autorange="reversed", showticklabels=False, row=1, col=1)
    row_yaxis_kwargs = dict(
        showticklabels=True, tickmode="array",
        tickvals=list(row_order), ticktext=row_ticktext,
    )
    if row_axis_reversed:
        row_yaxis_kwargs["autorange"] = "reversed"
    fig.update_yaxes(row=1, col=1, **row_yaxis_kwargs)
    fig.update_xaxes(showticklabels=False, row=1, col=2)
    fig.update_yaxes(showticklabels=False, row=1, col=2)
    fig.update_yaxes(type="log", autorange="reversed", showticklabels=False, row=2, col=2)
    fig.update_xaxes(
        showticklabels=True, tickangle=45, tickmode="array",
        tickvals=list(col_order), ticktext=col_ticktext,
        row=2, col=2,
    )

    fig.update_layout(
        height=height, width=width, showlegend=False,
        title=title, margin=margin,
    )
    return fig


def load_kamath_collapse_and_palette(repo_root, strategy=1):
    """
    Load `kamath_collapse_dic` from a `07_agreement_matrix_reference_RefStrategyN.py`
    script and build a matching uber-category colour palette.

    Parameters
    ----------
    repo_root : pathlib.Path
        Repository root (e.g. from `get_repo_root`).
    strategy : int
        Which RefStrategyN script to load. `kamath_collapse_dic` and
        `MA_CELL_TYPE_COLS` are identical across all 3 strategies, so the
        default (1) is fine regardless of which strategy's data is being
        displayed.

    Returns
    -------
    tuple
        compare_mod : module
            The imported `07_agreement_matrix_reference_RefStrategyN` module
            (exposes `kamath_collapse_dic`, `MA_CELL_TYPE_COLS`, etc.).
        uber_categories : list of str
            Sorted distinct uber categories (via `coarse_category`) across
            `kamath_collapse_dic`'s keys.
        uber_palette : dict
            uber_category : (r, g, b) tuple, from `golden_angle_palette`.
    """
    script_path = repo_root / (
        "reference_building/benchmark/reference_data/ma/scripts/"
        f"07_agreement_matrix_reference_RefStrategy{strategy}.py"
    )
    spec = importlib.util.spec_from_file_location("compare_uber_module", script_path)
    compare_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(compare_mod)

    uber_categories = sorted(set(coarse_category(ct) for ct in compare_mod.kamath_collapse_dic.keys()))
    uber_palette = golden_angle_palette(uber_categories)
    return compare_mod, uber_categories, uber_palette


def load_all_vs_all(
    csv_dir,
    all_vs_all_filename="asap_agreement_all_vs_all.csv",
    marginals_filename="asap_agreement_marginals.csv",
):
    """
    Load an all-vs-all PCC matrix plus its row/column marginal totals.

    Parameters
    ----------
    csv_dir : pathlib.Path
        Directory containing the two CSVs (one strategy's cell2location output dir).
    all_vs_all_filename : str
        Filename of the all-vs-all PCC matrix CSV.
    marginals_filename : str
        Filename of the row/column marginal totals CSV.

    Returns
    -------
    tuple
        all_vs_all : pandas.DataFrame
            Our pooled category (rows) x Ma pooled category (columns), PCC values.
        row_counts : pandas.Series
            Per-row (our-side) total proportion mass, reindexed to all_vs_all's rows.
        col_counts : pandas.Series
            Per-column (Ma-side) total proportion mass, reindexed to all_vs_all's columns.
    """
    all_vs_all = pd.read_csv(Path(csv_dir, all_vs_all_filename), index_col=0)
    marginals_ = pd.read_csv(Path(csv_dir, marginals_filename))
    row_counts = marginals_[marginals_["side"] == "our"].set_index("category")["total"].reindex(all_vs_all.index)
    col_counts = marginals_[marginals_["side"] == "ma"].set_index("category")["total"].reindex(all_vs_all.columns)
    return all_vs_all, row_counts, col_counts


def build_all_vs_all_figure(
    all_vs_all, row_counts, col_counts, sort_by, title,
    wide_factor=1.0, height_factor=1.0,
    y_annotation_offset=-0.30, x_annotation_offset=-0.30,
):
    """
    Heatmap-with-marginals figure comparing our pooled categories vs. Ma's.

    Parameters
    ----------
    all_vs_all : pandas.DataFrame
        Our pooled category (rows) x Ma pooled category (columns), PCC values.
    row_counts : pandas.Series
        Per-row (our-side) total proportion mass.
    col_counts : pandas.Series
        Per-column (Ma-side) total proportion mass.
    sort_by : str
        One of `crosstab_sort_order`'s supported sort modes.
    title : str
        Figure title.
    wide_factor : float
        Multiplier on the computed figure width -- shrink this when several
        of these figures are shown side by side.
    height_factor : float
        Multiplier on the computed figure height.
    y_annotation_offset : float
        Paper-coordinate y position for the "Ma (2025) — RCTD" x-axis label
        annotation; more negative pushes it further below the plot.
    x_annotation_offset : float
        Paper-coordinate x position for the "ASAP — Cell2Location" y-axis
        label annotation; more negative pushes it further left of the plot.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    row_order, col_order = crosstab_sort_order(all_vs_all, row_counts, col_counts, sort_by)

    data = all_vs_all.loc[row_order, col_order]
    row_marg = row_counts.loc[row_order]
    col_marg = col_counts.loc[col_order]

    height = max(500, len(data) * 35 + 220) * height_factor
    width = max(800, len(data.columns) * 40 + 280) * wide_factor

    fig = crosstab_heatmap_figure(
        data, row_order, col_order, row_marg, col_marg, height, width,
        heatmap_kwargs=dict(
            colorscale="RdBu", reversescale=True, zmin=-1, zmid=0, zmax=1,
            colorbar=dict(title="PCC"),
            hovertemplate="ASAP=%{y}<br>Ma=%{x}<br>PCC=%{z:.2f}<extra></extra>",
        ),
        marginal_frac=0.15,
        row_axis_reversed=True,  # preserves this plot's existing top=highest-count row ordering
        title=title,
        margin=dict(l=150, b=140, t=60, r=60),
    )
    fig.add_annotation(
        text="<b>Ma (2025) — RCTD</b>", xref="paper", yref="paper",
        x=0.5, y=y_annotation_offset, xanchor="center", yanchor="top", showarrow=False,
    )
    fig.add_annotation(
        text="<b>ASAP — Cell2Location</b>", xref="paper", yref="paper",
        x=x_annotation_offset, y=0.5, xanchor="center", yanchor="middle", textangle=-90, showarrow=False,
    )
    return fig


def strip_known_prefix(name, prefixes):
    """
    Strip the first of `prefixes` that `name` starts with.

    Parameters
    ----------
    name : str
        String to strip a prefix from.
    prefixes : tuple of str
        Candidate prefixes, checked in order; the first match is stripped.

    Returns
    -------
    str
        `name` with its matching prefix removed, or `name` unchanged if none match.
    """
    for prefix in prefixes:
        if name.startswith(prefix):
            return name[len(prefix):]
    return name


def coarse_category(label, sep="_"):
    """Coarse category = the prefix of `label` before the first `sep`, e.g. "Astro_GJB6_OXTR" -> "Astro"."""
    return label.split(sep, 1)[0]


def strip_q05_abundance_prefix(name):
    """
    Strip cell2location's `(frac_)q05cell_abundance_w_sf_` obsm/obs column prefix.

    Parameters
    ----------
    name : str
        Column name, e.g. "frac_q05cell_abundance_w_sf_CALB1_CALCR" or
        "q05cell_abundance_w_sf_CALB1_CALCR".

    Returns
    -------
    str
        `name` with the matching prefix removed, or unchanged if neither matches.
    """
    return strip_known_prefix(name, prefixes=("frac_q05cell_abundance_w_sf_", "q05cell_abundance_w_sf_"))


def selected_cell_types(cell_types, cell_types_, view, da_uber_categories):
    """
    Filter (cell_types, cell_types_) down to DA-only cell types when requested.

    Parameters
    ----------
    cell_types : list of str
        Long (`.obsm`-style) names, e.g. from one strategy's adata.
    cell_types_ : list of str
        Short names (prefix stripped), same length/order as `cell_types`.
    view : str
        Current "Show:" dropdown value; "DAs only" restricts to
        `da_uber_categories`, anything else keeps every cell type.
    da_uber_categories : set of str
        Uber categories (from `coarse_category`) considered DA, e.g.
        {"CALB1", "SOX6"}.

    Returns
    -------
    tuple
        cell_types_sel : list
            Long (`.obsm`-style) names for the selected cell types.
        cell_types_sel_ : list
            Short names (prefix stripped) for the selected cell types.
    """
    pairs = [
        (ct, ct_) for ct, ct_ in zip(cell_types, cell_types_)
        if view != "DAs only" or coarse_category(ct_) in da_uber_categories
    ]
    cell_types_sel, cell_types_sel_ = zip(*pairs) if pairs else ((), ())
    return list(cell_types_sel), list(cell_types_sel_)


def ensure_uber_colors(cols, uber_palette):
    """
    Extend `uber_palette` in place with any uber category in `cols` not
    already covered.

    Parameters
    ----------
    cols : iterable of str
        Fine-grained cell type names (e.g. a load's `cell_types_`).
    uber_palette : dict
        Uber category (from `coarse_category`) to colour, extended in place.

    Returns
    -------
    None
    """
    missing = sorted({coarse_category(c) for c in cols} - set(uber_palette))
    if missing:
        uber_palette.update(golden_angle_palette(list(uber_palette) + missing))


def colorize_label_html(label, prefixes, color="red"):
    """Wrap `label` in a coloured HTML `<span>` if it starts with any of `prefixes`."""
    if str(label).startswith(tuple(prefixes)):
        return f'<span style="color:{color}">{label}</span>'
    return label


def plot_umaps_categ(
    adata,
    panel_colors=("source", "cell_type_pred", "da_subtype"),
    plot_w=400, plot_h=400,
    margin_l=60, margin_t=50, margin_b=50,
    legend_w=260, cbar_w=100,
    point_size=3, opacity=0.6,
    priority_cols=("source", "cell_type_pred", "cell_type_scanvi", "da_subtype"),
    native_cols_by_source=None,
    display_label_remap=None,
    remap_cols=("cell_type_pred", "cell_type_scanvi", "da_subtype"),
):
    """
    N-panel interactive UMAP, colouring by obs columns.

    `number_of_panels` is `len(panel_colors)` -- one column per panel, 1-3.
    Lasso/box-select on any panel updates that panel's entry in the returned
    dict live -- read it any time after selecting (no need to re-run the cell).

    Parameters
    ----------
    adata : anndata.AnnData
        Must have `obsm['X_umap']` (2D UMAP coordinates).
    panel_colors : tuple of str
        One obs column per panel; `len(panel_colors)` sets the panel count (1-3).
    plot_w : int
        Per-panel plot width in pixels.
    plot_h : int
        Per-panel plot height in pixels.
    margin_l : int
        Left margin in pixels.
    margin_t : int
        Top margin in pixels.
    margin_b : int
        Bottom margin in pixels.
    legend_w : int
        Right-side legend width in pixels, used when the colour column is categorical.
    cbar_w : int
        Right-side colourbar width in pixels, used when the colour column is continuous.
    point_size : int
        Initial marker size (each panel's own slider can change it after render).
    opacity : float
        Initial marker opacity (each panel's own slider can change it after render).
    priority_cols : tuple of str
        Columns listed first (if present) in the "Colour by" dropdown and in hover text.
    native_cols_by_source : dict or None
        Maps a `source` value to a tuple of its own native label columns,
        appended to `priority_cols` (e.g. {"kamath": ("cell_type_scanvi",)}).
        None/empty skips this extension.
    display_label_remap : dict or None
        Maps raw category values to display labels, applied only to
        `remap_cols` for presentation (the underlying h5ad labels are unchanged).
        None/empty skips remapping.
    remap_cols : tuple of str
        Columns `display_label_remap` applies to, if given.

    Returns
    -------
    dict
        panel_color : list
            obs_names (barcodes) currently lasso/box-selected in that panel.
            Empty until a selection is made.
    """
    number_of_panels = len(panel_colors)
    if not 1 <= number_of_panels <= 3:
        raise ValueError(f"panel_colors must have 1-3 entries, got {number_of_panels}")

    native_cols_by_source = native_cols_by_source or {}
    display_label_remap = display_label_remap or {}

    umap_df = pd.DataFrame(
        adata.obsm["X_umap"], columns=["UMAP_1", "UMAP_2"], index=adata.obs_names,
    )
    umap_df = pd.concat([umap_df, adata.obs], axis=1)
    umap_df["_barcode"] = umap_df.index  # customdata needs a column, not the index itself

    for col in remap_cols:
        if col not in umap_df.columns:
            continue
        if isinstance(umap_df[col].dtype, pd.CategoricalDtype):
            rename = {k: v for k, v in display_label_remap.items() if k in umap_df[col].cat.categories}
            umap_df[col] = umap_df[col].cat.rename_categories(rename)
        else:
            umap_df[col] = umap_df[col].replace(display_label_remap)

    priority = list(priority_cols)
    for cols in native_cols_by_source.values():
        priority.extend(c for c in cols if c not in priority)
    all_obs = adata.obs.columns.tolist()
    color_options = (
        [c for c in priority if c in all_obs]
        + sorted(c for c in all_obs if c not in priority)
    )
    hover_cols = [c for c in priority if c in umap_df.columns]

    def is_categorical(col):
        s = umap_df[col]
        return (
            pd.api.types.is_object_dtype(s)
            or pd.api.types.is_bool_dtype(s)
            or hasattr(s, "cat")
            or (pd.api.types.is_integer_dtype(s) and s.nunique() < 30)
        )

    def build_figure(color_by, size, alpha, selected_barcodes, out):
        is_cat = is_categorical(color_by)
        r_margin = legend_w if is_cat else cbar_w
        total_w = margin_l + plot_w + r_margin

        hover = {c: True for c in hover_cols if c != color_by}
        common = dict(
            data_frame=umap_df, x="UMAP_1", y="UMAP_2",
            color=color_by, hover_data=hover, opacity=alpha,
            custom_data=["_barcode"], title=f"<b>{color_by}</b>",
        )
        fig = px.scatter(**common, color_continuous_scale=None if is_cat else "Viridis")
        if is_cat:
            fig.update_layout(legend=dict(
                title=None, itemsizing="constant",
                x=1.01, xanchor="left", y=0.5, yanchor="middle",
            ))

        fig.update_traces(marker_size=size)
        fig.update_layout(
            autosize=False, width=total_w, height=plot_h + margin_t + margin_b,
            xaxis_title="UMAP 1", yaxis_title="UMAP 2",
            margin=dict(autoexpand=False, l=margin_l, r=r_margin, t=margin_t, b=margin_b),
            dragmode="lasso",  # click-drag lasso-selects by default; box-select still in the toolbar
        )

        fig = go.FigureWidget(fig)

        # A categorical `color` builds one trace per category, so point_inds are
        # trace-local -- customdata (the barcode) stays point-aligned regardless,
        # so read barcodes from it directly instead of re-deriving row positions.
        per_trace_selection = [[] for _ in fig.data]

        def make_on_select(trace_idx):
            def on_select(trace, points, selector):
                per_trace_selection[trace_idx] = [
                    trace.customdata[i][0] for i in points.point_inds
                ]
                selected_barcodes[:] = [b for sel in per_trace_selection for b in sel]
                with out:
                    print(f"  [{color_by}] {len(selected_barcodes):,} cells selected")
            return on_select

        for i, trace in enumerate(fig.data):
            trace.on_selection(make_on_select(i))

        return fig

    def make_panel(default_color):
        if default_color not in color_options:
            print(f"  Warning: '{default_color}' not in adata.obs — defaulting to '{color_options[0]}'")
            default_color = color_options[0]

        color_w = widgets.Dropdown(
            options=color_options, value=default_color,
            description="Colour by:", style={"description_width": "initial"},
            layout=widgets.Layout(width="340px"),
        )
        size_w = widgets.IntSlider(
            value=point_size, min=1, max=10, step=1, description="Point size:",
            layout=widgets.Layout(width="320px"),
        )
        opacity_w = widgets.FloatSlider(
            value=opacity, min=0.1, max=1.0, step=0.05, description="Opacity:",
            layout=widgets.Layout(width="320px"),
        )
        out = widgets.Output()
        selected_barcodes = []

        def redraw(_=None):
            out.clear_output(wait=True)
            with out:
                display(build_figure(color_w.value, size_w.value, opacity_w.value, selected_barcodes, out))

        color_w.observe(redraw, names="value")
        size_w.observe(redraw, names="value")
        opacity_w.observe(redraw, names="value")

        vbox = widgets.VBox([color_w, size_w, opacity_w, out])

        def init():
            with out:
                display(build_figure(color_w.value, size_w.value, opacity_w.value, selected_barcodes, out))

        return vbox, init, default_color, selected_barcodes

    panels = [make_panel(c) for c in panel_colors]
    for _, init, _, _ in panels:
        init()
    display(widgets.HBox([box for box, _, _, _ in panels]))

    return {color: selected_barcodes for _, _, color, selected_barcodes in panels}


def spatial_scatter_grid(
    adata, columns, titles=None,
    cmap="magma", vmin=0, vmax_percentile=99.2, vmax_floor=None,
    ncols=4, size=1.3, img_res_key="hires",
    title_color_fn=None,
    shared_vmax=False,
):
    """
    Grid of spatial scatter panels, cropped to the tissue's own spot extent.

    One panel per entry in `columns` -- each is both the `color=` argument to
    `squidpy.pl.spatial_scatter` (an `.obs` column or a gene/var name --
    squidpy resolves either) and (unless `shared_vmax`) the source of that
    panel's own vmax. An entry of `None` renders a plain tissue-only panel
    instead (no color overlay, no vmax computed) -- pass an explicit string in
    the matching position of `titles` for it (e.g. "Tissue only"), since
    there's no column/gene name to fall back to.

    Parameters
    ----------
    adata : anndata.AnnData
        Must have `obsm['spatial']` and `uns['spatial']` (Visium scale factors).
    columns : list of str or None
        One panel per entry; `None` renders a tissue-only panel.
    titles : list of str or None
        Panel titles, same length/order as `columns`. Defaults to `columns`
        (pass a separate list if a panel's display title differs from its
        underlying column/gene name, e.g. a column copied into `.obs` under
        a prefixed key, or any `None` entries in `columns`).
    cmap : str
        Colormap.
    vmin : float or None
        Colour-scale floor, forwarded to `squidpy.pl.spatial_scatter`. Explicit
        0 by default to avoid squidpy falling back to each column's own
        minimum. Pass None to restore squidpy's own default.
    vmax_percentile : float
        Percentile (of non-NaN values) used as the vmax -- each panel's own
        values if `shared_vmax` is False, the pooled values of all panels if True.
    vmax_floor : float or None
        Minimum vmax -- sparse columns can have a near-zero percentile that
        saturates the panel to a single flat colour; floor it if given.
    ncols : int
        Grid columns.
    size : float
        Spot size, forwarded to `squidpy.pl.spatial_scatter`.
    img_res_key : str
        Image resolution key, forwarded to `squidpy.pl.spatial_scatter`.
    title_color_fn : callable or None
        `title_color_fn(column) -> colour`, applied to each panel's title.
        None (default) skips title colouring.
    shared_vmax : bool
        False (default) scales each panel independently from its own values --
        appropriate when panels are on genuinely different scales (e.g.
        different genes). True computes ONE vmax from all panels' pooled
        values and applies it everywhere -- needed when panels should stay
        visually comparable in absolute magnitude (e.g. q05 cell-type
        abundance across cell types): independent per-panel scaling can make
        a cell type whose true values are tiny (e.g. max 0.05) look just as
        "hot" as one whose true values are much larger (e.g. max 8.5).

    Returns
    -------
    matplotlib.figure.Figure or None
        Avoids `plt.show()` unreliability when called more than once per cell.
    """
    if not columns:
        print("spatial_scatter_grid: empty `columns` -- nothing to plot.")
        return None

    titles = list(titles) if titles is not None else list(columns)

    library_id = list(adata.uns["spatial"].keys())[0]
    scalef = adata.uns["spatial"][library_id]["scalefactors"][f"tissue_{img_res_key}_scalef"]
    coords = adata.obsm["spatial"] * scalef
    # nanmin/nanmax: a single spot with a NaN coordinate (seen in real data) would
    # otherwise poison the whole crop computation via plain min/max.
    (x_min, y_min), (x_max, y_max) = np.nanmin(coords, axis=0), np.nanmax(coords, axis=0)
    x_pad, y_pad = (x_max - x_min) * 0.05, (y_max - y_min) * 0.05

    # Close any figures left open by a previous call -- squidpy's colorbar/legend
    # creation can touch pyplot's global current-figure state even when given an
    # explicit `ax=`, which otherwise bleeds across repeated calls in one cell
    # (e.g. one `widgets.Output()` per strategy in a loop).
    plt.close("all")

    nrows = int(np.ceil(len(columns) / ncols))
    fig, axs = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows))
    axs = np.atleast_1d(axs).ravel()

    column_values = {}
    for column in columns:
        if column is None:
            continue
        if column in adata.obs.columns:
            values = adata.obs[column].values
        else:
            values = adata[:, column].X
            values = values.toarray() if hasattr(values, "toarray") else np.asarray(values)
        column_values[column] = np.asarray(values).flatten()

    shared_vmax_value = None
    if shared_vmax and column_values:
        shared_vmax_value = np.nanpercentile(np.concatenate(list(column_values.values())), vmax_percentile)
        if vmax_floor is not None:
            shared_vmax_value = max(shared_vmax_value, vmax_floor)

    for ax, column, title in zip(axs, columns, titles):
        if column is None:
            sq.pl.spatial_scatter(adata, color=None, size=size, img_res_key=img_res_key, title=title, ax=ax)
        else:
            if shared_vmax:
                vmax = shared_vmax_value
            else:
                vmax = np.nanpercentile(column_values[column], vmax_percentile)
                if vmax_floor is not None:
                    vmax = max(vmax, vmax_floor)
            sq.pl.spatial_scatter(
                adata, color=column, cmap=cmap, vmin=vmin, vmax=vmax, size=size,
                img_res_key=img_res_key, title=title, ax=ax,
            )
        ax.set_xlim(x_min - x_pad, x_max + x_pad)
        ax.set_ylim(y_max + y_pad, y_min - y_pad)
        if title_color_fn is not None and column is not None:
            ax.title.set_color(title_color_fn(column))

    for ax in axs[len(columns):]:
        ax.axis("off")

    plt.tight_layout()
    return fig


def spatial_plot_genes(
    adata, genes, sample_name,
    mode="raw_and_norm",
    cmap="viridis",
    normalize_target_sum=1e4,
    raw_vmax_floor=1,
    norm_vmax_floor=0.01,
    ncols=6,
):
    """
    Render spatial scatter figure(s) for a gene list -- raw UMI counts and/or
    log1p-normalized expression -- via `spatial_scatter_grid`, raw first.

    Mutates `adata` in place: each plotted gene's values are copied into
    `.obs` under a value-type-suffixed key (e.g. "TH (raw UMI)"), and if
    `mode` includes "norm", `.X` is overwritten via `sc.pp.normalize_total`/
    `sc.pp.log1p` -- this affects ALL genes in `adata`, not just those plotted.

    Parameters
    ----------
    adata : anndata.AnnData
        Loaded spatial sample with raw counts in `.X` and `obsm['spatial']`/
        `uns['spatial']` (Visium scale factors).
    genes : list of str
        Gene symbols to plot.
    sample_name : str
        Sample name, used only for the missing-gene warning message.
    mode : str
        One of "raw", "norm", "raw_and_norm" (default) -- which figure(s) to render.
    cmap : str
        Colormap, forwarded to `spatial_scatter_grid`.
    normalize_target_sum : float
        `target_sum` forwarded to `sc.pp.normalize_total`. Only used if `mode`
        includes "norm".
    raw_vmax_floor : float
        `vmax_floor` for the raw-UMI figure. Raw counts are integers, so a
        floor around 1 avoids a degenerate near-zero colour scale on sparse
        samples.
    norm_vmax_floor : float
        `vmax_floor` for the log1p-normalized figure. Normalized values are
        small fractions, so a floor around 0.01 is appropriate instead --
        sharing one floor between both value types collapses them to the
        same degenerate scale on sparse samples.
    ncols : int
        Panels per row (including the leading "Tissue only" panel), forwarded
        to `spatial_scatter_grid`. Default 6 fits the common case -- Tissue
        only + the 4 Pan-DA markers + one more (e.g. a highly-expressed
        control gene) -- in a single row; longer gene lists wrap to extra
        rows instead of squeezing every panel into one row.

    Returns
    -------
    None
    """
    if mode not in ("raw", "norm", "raw_and_norm"):
        raise ValueError(f'mode must be one of "raw", "norm", "raw_and_norm" -- got {mode!r}')

    adata.var_names_make_unique()

    genes_present = [g for g in genes if g in adata.var_names]
    genes_missing = [g for g in genes if g not in adata.var_names]
    if genes_missing:
        print(f"WARNING: {len(genes_missing)} gene(s) not in {sample_name}'s gene panel: {genes_missing}")

    need_raw = mode in ("raw", "raw_and_norm")
    need_norm = mode in ("norm", "raw_and_norm")

    raw_columns = []
    if need_raw:
        for gene in genes_present:
            raw_vals = adata[:, gene].X
            raw_vals = raw_vals.toarray().flatten() if hasattr(raw_vals, "toarray") else raw_vals.flatten()
            raw_col = f"{gene} (raw UMI)"
            adata.obs[raw_col] = raw_vals
            raw_columns.append(raw_col)

        fig_raw = spatial_scatter_grid(
            adata, [None] + raw_columns, titles=["Tissue only"] + raw_columns,
            cmap=cmap, vmax_floor=raw_vmax_floor, ncols=ncols,
        )
        if fig_raw is not None:
            display(fig_raw)
            plt.close(fig_raw)

    if need_norm:
        sc.pp.normalize_total(adata, target_sum=normalize_target_sum)
        sc.pp.log1p(adata)

        norm_columns = []
        for gene in genes_present:
            norm_vals = adata[:, gene].X
            norm_vals = norm_vals.toarray().flatten() if hasattr(norm_vals, "toarray") else norm_vals.flatten()
            norm_col = f"{gene} (norm log1p)"
            adata.obs[norm_col] = norm_vals
            norm_columns.append(norm_col)

        fig_norm = spatial_scatter_grid(
            adata, [None] + norm_columns, titles=["Tissue only"] + norm_columns,
            cmap=cmap, vmax_floor=norm_vmax_floor, ncols=ncols,
        )
        if fig_norm is not None:
            display(fig_norm)
            plt.close(fig_norm)


def plot_umaps_by_gene(
    adata,
    panel_genes=None,
    n_panels=4,
    suggested_genes=(),
    plot_w=320, plot_h=320,
    margin_l=55, margin_t=50, margin_b=50, cbar_w=90,
    point_size=3, opacity=0.6,
):
    """
    N-panel UMAP with a free-text gene box per panel; colours by log-normalized expression of any gene.

    Parameters
    ----------
    adata : anndata.AnnData
        Must have `obsm['X_umap']` and `var['gene_symbol']` (a symbol lookup
        for a `var_names` index that isn't already gene symbols, e.g. Ensembl IDs).
    panel_genes : tuple of str or None
        Initial gene (or Ensembl ID) per panel; `n_panels` empty panels if None.
    n_panels : int
        Number of panels, used only when `panel_genes` is None.
    suggested_genes : tuple of str
        Autocomplete suggestions offered in each panel's gene combobox.
    plot_w : int
        Per-panel plot width in pixels.
    plot_h : int
        Per-panel plot height in pixels.
    margin_l : int
        Left margin in pixels.
    margin_t : int
        Top margin in pixels.
    margin_b : int
        Bottom margin in pixels.
    cbar_w : int
        Right-side colourbar width in pixels.
    point_size : int
        Initial marker size (each panel's own slider can change it after render).
    opacity : float
        Initial marker opacity (each panel's own slider can change it after render).

    Returns
    -------
    None
    """
    if panel_genes is None:
        panel_genes = ("",) * n_panels

    umap_coords = pd.DataFrame(
        adata.obsm["X_umap"], columns=["UMAP_1", "UMAP_2"], index=adata.obs_names,
    )

    def resolve_expression(query):
        """Return (expression Series, display label) or (None, None) if not found."""
        query = query.strip()
        if not query:
            return None, None
        try:
            df = sc.get.obs_df(adata, keys=[query.upper()], gene_symbols="gene_symbol")
            return df[query.upper()], query.upper()
        except KeyError:
            pass
        if query in adata.var_names:
            symbol = adata.var.loc[query, "gene_symbol"]
            label = symbol if pd.notna(symbol) else query
            return sc.get.obs_df(adata, keys=[query])[query], label
        return None, None

    def build_figure(gene_query, size, alpha):
        total_w = margin_l + plot_w + cbar_w
        expr, label = resolve_expression(gene_query)
        if expr is None:
            fig = go.Figure()
            if gene_query.strip():
                text = f"'{gene_query}' not found in this h5ad's gene panel ({adata.n_vars:,} genes)"
                color = "red"
            else:
                text = "Type a gene symbol or Ensembl ID above"
                color = "gray"
            fig.add_annotation(
                text=text,
                xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False,
                font=dict(size=13, color=color),
            )
            fig.update_layout(width=total_w, height=plot_h + margin_t + margin_b)
            return fig

        df = umap_coords.copy()
        df["expression"] = expr.values
        # sort_order equivalent (as in sc.pl.umap): draw highest-expressing
        # cells last (on top) so they don't get buried under low/zero-
        # expression cells -- plotly renders scatter points in row order,
        # later rows on top. NaNs sorted first so they never cover real values.
        df = df.sort_values("expression", ascending=True, na_position="first")

        fig = px.scatter(
            df, x="UMAP_1", y="UMAP_2", color="expression",
            color_continuous_scale="Viridis", opacity=alpha,
            title=f"<b>{label}</b>",
        )
        fig.update_traces(marker_size=size)
        fig.update_layout(
            autosize=False, width=total_w, height=plot_h + margin_t + margin_b,
            xaxis_title="UMAP 1", yaxis_title="UMAP 2",
            coloraxis_colorbar=dict(title="log-norm expr"),
            margin=dict(autoexpand=False, l=margin_l, r=cbar_w, t=margin_t, b=margin_b),
        )
        return fig

    def make_panel(default_gene):
        gene_w = widgets.Combobox(
            value=default_gene, options=list(suggested_genes), ensure_option=False,
            description="Gene:", style={"description_width": "initial"},
            layout=widgets.Layout(width="280px"), continuous_update=False,
        )
        size_w = widgets.IntSlider(
            value=point_size, min=1, max=10, step=1, description="Point size:",
            layout=widgets.Layout(width="260px"),
        )
        opacity_w = widgets.FloatSlider(
            value=opacity, min=0.1, max=1.0, step=0.05, description="Opacity:",
            layout=widgets.Layout(width="260px"),
        )
        out = widgets.Output()

        def redraw(_=None):
            out.clear_output(wait=True)
            with out:
                display(build_figure(gene_w.value, size_w.value, opacity_w.value))

        gene_w.observe(redraw, names="value")
        size_w.observe(redraw, names="value")
        opacity_w.observe(redraw, names="value")

        vbox = widgets.VBox([gene_w, size_w, opacity_w, out])

        def init():
            with out:
                display(build_figure(gene_w.value, size_w.value, opacity_w.value))

        return vbox, init

    panels = [make_panel(g) for g in panel_genes]
    for _, init in panels:
        init()
    display(widgets.HBox([box for box, _ in panels]))
