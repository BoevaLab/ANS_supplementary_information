import scanpy as sc
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.sparse import diags
from scipy.stats import median_abs_deviation


def shifted_transformation(adata, y0=1):
    """
    From Twitter post https://twitter.com/Sanbomics/status/1647654042749874177?s=20
    Refering to publication by Ahlmann-Eltze & Huber.

    Ahlmann-Eltze, C., Huber, W. Comparison of transformations for single-cell RNA-seq data.
    Nat Methods (2023). https://doi.org/10.1038/s41592-023-01814-1
    """
    target_sum = np.mean(adata.X.sum(axis=1))
    print(
        f"Mean shift logarithm normalization with normalization target count {target_sum}"
    )
    size_factors = adata.X.sum(axis=1) / target_sum

    adata.X = diags(1 / size_factors.A1).dot(adata.X)
    adata.X.data = np.log(adata.X.data + y0)
    adata.uns["log1p"] = {"base": None}
    return adata


def is_outlier(adata, metric: str, nmads: int):
    """
    The method is taken from the tutorial https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html#filtering-low-quality-reads

    Heumos, L., Schaar, A.C., Lance, C. et al. Best practices for single-cell analysis across modalities. Nat Rev Genet (2023). https://doi.org/10.1038/s41576-023-00586-w
    """
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


def filter_low_quality_reads(
    adata,
    mad_tot_cnt=5,
    mad_ngenes_cnt=5,
    nr_top_genes=20,
    mad_pct_cnt_top_genes=5,
    mad_pct_mt=3,
    min_pct_mt=8,
):
    """
    The method is taken from the tutorial https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html#filtering-low-quality-reads

    Heumos, L., Schaar, A.C., Lance, C. et al. Best practices for single-cell analysis across modalities.
    Nat Rev Genet (2023). https://doi.org/10.1038/s41576-023-00586-w
    """
    if "mt" not in adata.var:
        # get mitochondrial genes
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
    if "ribo" not in adata.var:
        # get ribosomal genes
        adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    if "hb" not in adata.var:
        # get hemoglobin genes.
        adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))

    # compute the quality control metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        inplace=True,
        percent_top=[nr_top_genes],
        log1p=True,
    )

    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", mad_tot_cnt)
        | is_outlier(adata, "log1p_n_genes_by_counts", mad_ngenes_cnt)
        | is_outlier(
            adata, f"pct_counts_in_top_{nr_top_genes}_genes", mad_pct_cnt_top_genes
        )
    )
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", mad_pct_mt) | (
        adata.obs["pct_counts_mt"] > min_pct_mt
    )

    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
    return adata


def filtercells(adata, sample_col="orig.ident", params_cell_filtering={}):
    """
    Filter loww quality reads per sample as suggested by Heumos et al.

    Heumos, L., Schaar, A.C., Lance, C. et al. Best practices for single-cell analysis across modalities.
    Nat Rev Genet (2023). https://doi.org/10.1038/s41576-023-00586-w
    """
    nr_cells_orig = adata.shape[0]
    adatas = {}
    for sid, sample_data in adata.obs.groupby(sample_col):
        adatas[sid] = adata[sample_data.index,].copy()
    for key, curr_adata in adatas.items():
        adatas[key] = filter_low_quality_reads(curr_adata, **params_cell_filtering)
    adata = sc.concat(list(adatas.values()), join="outer", merge="same")
    nr_cells_filt = adata.shape[0]
    print(
        f"Filtering {(nr_cells_orig-nr_cells_filt)} of {nr_cells_orig} low quality cells"
        f"({np.round((nr_cells_orig-nr_cells_filt)/nr_cells_orig *100, decimals=2)}%)."
    )
    return adata


def filtergenes(adata, pct=0.01):
    """
    Remove genes that are not present in at least 1% of all cells. We do the same as it was done in CanSig.

    CanSig: Discovering de novo shared transcriptional programs in single cancer cells
    Josephine Yates, Florian Barkmann, Paweł Czyż, Marc Glettig, Frederieke Lohmann,
    Richard von der Horst, Elia Saquand, Nicolas Volken, Agnieszka Kraft, Valentina Boeva,
    bioRxiv 2022.04.14.488324; doi: https://doi.org/10.1101/2022.04.14.488324
    """
    nr_cells, nr_genes = adata.shape
    gene_expr_in_cells_cnts = adata.X.getnnz(axis=0)
    enough_genes = gene_expr_in_cells_cnts - nr_cells * pct
    print(
        f"Filtering {np.sum(enough_genes < 0)} of {nr_genes} genes"
        f"({np.round((np.sum(enough_genes < 0))/nr_genes *100, decimals=2)}%)."
    )
    adata = adata[:, enough_genes >= 0].copy()
    return adata


def preprocess_dataset(
    adata,
    filter_cells=True,
    filter_genes=True,
    shift_method="mean",
    params_cell_filtering={},
    show=False,
):
    if show:
        # CREATE PLOT
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
        )
        sns.histplot(adata.obs["total_counts"], bins=100, kde=False)
        # sc.pl.violin(adata, 'total_counts')
        sc.pl.violin(adata, "pct_counts_mt")
        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
        plt.show()

    # FILTER CELLS
    if filter_cells:
        adata = filtercells(adata, params_cell_filtering=params_cell_filtering)
    
    if show:
        # CREATE PLOT
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
        )
        sns.histplot(adata.obs["total_counts"], bins=100, kde=False)
        # sc.pl.violin(adata, 'total_counts')
        sc.pl.violin(adata, "pct_counts_mt")
        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
        plt.show()

    # FILTER GENES
    if filter_genes:
        adata = filtergenes(adata)

    adata.layers["counts"] = adata.X

    # FILTER NORMALIZE
    if shift_method == "median":
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        adata.uns["log1p"]["base"] = None
    elif shift_method == "mean":
        adata = shifted_transformation(adata)
    elif shift_method == "CP10k":
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.uns["log1p"]["base"] = None
    else:
        raise ValueError(
            "Unknown shift transformation method! Can choose between mean, median and CP10k."
        )

    return adata
