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

import ipywidgets as widgets
from IPython.display import HTML, display
from pathlib import Path
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


def big_dropdown(**kwargs):
    """
    An `ipywidgets.Dropdown` pre-sized and styled via the "big-dropdown" CSS class.

    Parameters
    ----------
    **kwargs
        Forwarded to `ipywidgets.Dropdown` (e.g. `options`, `description`).
        `style`/`layout` default to the standard big-dropdown look if not given.

    Returns
    -------
    ipywidgets.Dropdown
        Styled dropdown; still needs `inject_big_dropdown_css()` called once
        somewhere earlier in the notebook for the class to have any effect.
    """
    kwargs.setdefault("style", {"description_width": "initial"})
    kwargs.setdefault("layout", widgets.Layout(width="700px", height="56px"))
    dropdown = widgets.Dropdown(**kwargs)
    dropdown.add_class("big-dropdown")
    return dropdown


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


def colorize_label_html(label, prefixes, color="red"):
    """Wrap `label` in a coloured HTML `<span>` if it starts with any of `prefixes`."""
    if str(label).startswith(tuple(prefixes)):
        return f'<span style="color:{color}">{label}</span>'
    return label
