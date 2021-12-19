import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc


def LogNormalize(data: ad.AnnData, scale_factor: int):
    """Normalizes data. Reimplemented Seurat LogNormalize() function:
    https://satijalab.org/seurat/reference/normalizedata

    Args:
        data (ad.AnnData): Dataset
    Returns:
        ad.AnnData: Normalized dataset with normalized counts in data.X
        and non-normalized counts in data.layer["counts"].
    """

    # Save non-normalized counts in separate layer
    data.layers["counts"] = data.X.copy()

    # Normalize & Log-transform
    data.X *= data.obs['size_factors'].values[:,None]
    data.X = sc.pp.log1p(data.X)
    return data
