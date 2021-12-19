import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc


def LogNormalize(data: ad.AnnData, scale_factor: int):
    """Normalizes data. Reimplemented Seurat LogNormalize() function:
    https://satijalab.org/seurat/reference/normalizedata

    Args:
        data (ad.AnnData): Dataset
        scale_factor (int): Sets the scale factor for cell-level normalization.
    Returns:
        ad.AnnData: Normalized dataset with normalized counts in data.X
        and non-normalized counts in data.layer["counts"].
    """

    # Recalculate count sum for each cell as gene QC filters could be applied
    count_sum = np.sum(data.X, axis=1)
    data.obs["n_counts_recalc"] = np.asarray(count_sum).reshape(-1)

    # Save non-normalized counts in separate layer
    data.layers["counts"] = data.X.copy()

    # Normalize & save normalized counts in data.X
    data.X /= data.obs["n_counts_recalc"].values[:, None]
    data.X *= scale_factor
    data.X = sc.pp.log1p(data.X)
    return data
