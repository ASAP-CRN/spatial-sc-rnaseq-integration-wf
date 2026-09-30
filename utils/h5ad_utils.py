"""
Generic h5ad / AnnData utilities shared across repo-root scripts.
"""

import numpy as np
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
