import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, issparse


def extract_raw_counts(adata, source_name: str = "") -> csr_matrix:
    """
    Extract raw integer count matrix from an AnnData object.

    Checks layers["counts"], then adata.raw.X, then adata.X in order.
    Validates that the matrix is integer-valued (whole-number floats such
    as 1.0 pass; fractional values such as 1.5 raise) and casts to int32.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData to extract counts from.
    source_name : str, optional
        Label used in error/info messages (e.g. file basename).

    Returns
    -------
    csr_matrix
        Sparse int32 count matrix (cells × genes).

    Raises
    ------
    ValueError
        If no raw count layer is found, or the matrix contains non-integer
        values.
    """
    if "counts" in adata.layers:
        mat = adata.layers["counts"]
    elif adata.raw is not None:
        mat = adata.raw.X
    else:
        pfx = f"{source_name}: " if source_name else ""
        print(
            f"  [INFO] {pfx}no layers['counts'] or adata.raw — "
            "using adata.X (integer-valued check happens below)"
        )
        mat = adata.X

    mat = csr_matrix(mat) if not issparse(mat) else mat.tocsr()

    if not np.issubdtype(mat.dtype, np.integer) and mat.data.size > 0:
        non_int_mask = mat.data % 1 != 0
        if non_int_mask.any():
            examples = mat.data[non_int_mask][:5].tolist()
            pfx = f"{source_name}: " if source_name else ""
            raise ValueError(
                f"{pfx}count matrix contains non-integer values "
                f"(e.g. {examples}). Raw counts must be integer-valued. "
                "Use a raw H5AD download or check --counts-layer."
            )

    if mat.dtype != np.int32:
        mat = mat.astype(np.int32)

    return mat


def resolve_marker_ensembl_ids(
    feature_name: pd.Series,
    marker_symbols: tuple,
    source_label: str,
) -> list[str]:
    """
    Resolve gene symbols to Ensembl IDs via a source's var['feature_name'].

    Parameters
    ----------
    feature_name : pd.Series
        A source's `adata.var['feature_name']`, indexed by Ensembl ID.
    marker_symbols : tuple
        Gene symbols to resolve.
    source_label : str
        Label for the warning message (e.g. 'Menon').

    Returns
    -------
    list[str]
        Ensembl IDs found. Unmatched symbols are warned about, not raised on.
    """
    symbols_str = feature_name.astype(str)
    symbol_to_ensembl = pd.Series(symbols_str.index, index=symbols_str.values)
    symbol_to_ensembl = symbol_to_ensembl[~symbol_to_ensembl.index.duplicated(keep="first")]
    found = symbol_to_ensembl.reindex(marker_symbols)
    missing = [s for s, e in zip(marker_symbols, found) if pd.isna(e)]
    if missing:
        print(f"  Warning: {len(missing)}/{len(marker_symbols)} marker symbols "
              f"not found in {source_label}'s var['feature_name']: {missing}")
    return found.dropna().tolist()


def apply_cell2location_gene_filter(
    adata,
    filter_genes_mode: str,
    da_marker_ensembl_ids: list = None,
):
    """
    Reduce `adata` to cell2location's data-driven recommended gene set.

    Parameters
    ----------
    adata : anndata.AnnData
        Raw-count AnnData (cells x genes), var_names = Ensembl gene IDs.
    filter_genes_mode : str
        "cell2location_default" : cell2location.utils.filtering.filter_genes
            with the tutorial's defaults (cell_count_cutoff=5,
            cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12).
        "cell2location_u_da" : cell2location_default, unioned with
            `da_marker_ensembl_ids`.
    da_marker_ensembl_ids : list, optional
        Ensembl IDs to union back in when filter_genes_mode is
        "cell2location_u_da" (ignored otherwise); required in that mode.

    Returns
    -------
    anndata.AnnData
        `adata` subset to the selected genes, original gene order preserved.

    Raises
    ------
    ValueError
        filter_genes_mode == "cell2location_u_da" but da_marker_ensembl_ids
        is None.
    """
    import matplotlib
    matplotlib.use("Agg")
    from cell2location.utils.filtering import filter_genes as c2l_filter_genes

    n_before = adata.n_vars
    selected = set(c2l_filter_genes(
        adata, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12
    ))
    print(f"  cell2location filter_genes: {n_before:,} → {len(selected):,} genes")

    if filter_genes_mode == "cell2location_u_da":
        if da_marker_ensembl_ids is None:
            raise ValueError(
                "filter_genes_mode='cell2location_u_da' requires da_marker_ensembl_ids."
            )
        added = [g for g in da_marker_ensembl_ids if g in adata.var_names and g not in selected]
        selected.update(added)
        print(f"  + {len(added)} Kamath DA marker gene(s) added back via union "
              f"(present in the gene set but dropped by the cell2location filter)")

    keep = [g for g in adata.var_names if g in selected]
    print(f"  Final gene set: {len(keep):,} genes")
    return adata[:, keep].copy()
