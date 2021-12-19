Scripts related to the single cell RNAseq data processing.

1. Input raw data required for running this pipeline:
    * Sparse matrix in AnnData (h5ad) format.
<br>

2. Create config.yaml file with QC settings. Example file:
    ```python3
    mito_max_fraction: 20

    cell_min_counts: 500

    cell_max_counts: 40000

    cell_min_genes: 200

    cell_max_genes: 6500

    gene_min_cells : 20
    ```
<br>

3. For QC calculation & filtering, run *scrna_qc.standard_qc* function with arguments:
    * data: initial AnnData dataset.
    * config: dictionary with QC metrics.
<br>

4. Get filtered file in AnnData (h5ad) format.
<br>

5. For normalization, run *scrna_normalization.LogNormalize* function with arguments:
    * data: filtered AnnData dataset with data.obs['size_factors'] calculated option

This is a LogNormalize() Seurat function reimplemented in python, you can learn more about original function [here](https://satijalab.org/seurat/reference/normalizedata).

6. Get filtered and normalized file in AnnData (h5ad) format.
