import argparse
import json
import os
import sys
from datetime import datetime

import pandas as pd
import scanpy as sc
from signaturescoring import score_signature
from signaturescoring.utils.utils import (
    get_bins_wrt_avg_gene_expression, get_mean_and_variance_gene_expression)

sys.path.append("../..")
from data.constants import DATASETS, BASE_PATH_DATA, BASE_PATH_EXPERIMENTS
from data.load_data import load_datasets
from data.preprocessing.preprocess_pbmc_helper import preprocess_dataset
from experiments.experiment_utils import (AttributeDict,
                                          get_scoring_method_params)

ALLOWED_CELL_SUBTYPES = ['B memory kappa', 'B naive kappa', 'B naive lambda', 'CD14 Mono', 'CD16 Mono', 'CD4 CTL', 'CD4 Naive', 'CD4 TCM_1', 'CD4 TCM_3', 'CD4 TEM_1', 'CD4 TEM_3', 'CD8 Naive', 'CD8 TEM_1', 'CD8 TEM_2', 'CD8 TEM_4', 'CD8 TEM_5', 'MAIT', 'NK_1', 'NK_2', 'NK_3', 'Platelet', 'cDC2_2']
SC_METHODS = ["adjusted_neighborhood_scoring", "seurat_scoring", "seurat_ag_scoring", "seurat_lvg_scoring", "scanpy_scoring"]

def _generate_storing_path(base_path, subtype):
    st_path_dec = os.path.join(base_path, subtype)

    if not os.path.exists(st_path_dec):
        os.makedirs(st_path_dec)
        sc.logging.info(f"> Created new directory with path {st_path_dec}")

    return st_path_dec


def _get_genes_of_last_two_expression_bins(adata, n_bins=25):
    df_mean_var = get_mean_and_variance_gene_expression(adata)
    df_mean_var = df_mean_var.sort_values(by="mean", ascending=True)
    gene_bins = get_bins_wrt_avg_gene_expression(df_mean_var["mean"], n_bins)
    print("last bin: ", (gene_bins == (n_bins - 1)).sum())
    print("2nd last bin: ", (gene_bins == (n_bins - 2)).sum())
    genes_of_lat_2_bins = gene_bins[(gene_bins == (n_bins - 1)) | (gene_bins == (n_bins - 2))].index.tolist()
    print("last two bins: ", len(genes_of_lat_2_bins))
    return genes_of_lat_2_bins


def _compute_score_matrix(adata, sig_genes, sc_method_full, sc_method, sc_params):
    scores = pd.DataFrame(
        0, index=adata.obs_names, columns=(["scoring_method"] + sig_genes)
    )
    scores.loc[:, "scoring_method"] = sc_method_full

    score_name = sc_params["score_name"]
    for gene in sig_genes:
        score_signature(method=sc_method, adata=adata, gene_list=[gene], **sc_params)
        scores.loc[:, gene] = adata.obs[score_name].copy()
    print((scores != 0).sum().sum())
    print((scores == 0).sum().sum())
    return scores

def load_pbmc_data():
    fn_data = os.path.join(BASE_PATH_DATA, 'raw_data/pbmc_citeseq.h5ad')
    adata = sc.read_h5ad(fn_data)

    adata = adata.raw.to_adata()
    adata.var_names = adata.var['_index']
    adata.var_names.name = None
    adata.var.columns = ['gene_names']

    if 'mt' not in adata.var:
        # get mitochondrial genes
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
    if 'ribo' not in adata.var:         
        # get ribosomal genes
        adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    if 'hb' not in adata.var:
        # get hemoglobin genes.
        adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))

    return adata


def main(config):
    # create storing path
    print('Create scoring path.')
    storing_path = _generate_storing_path(config.base_storing_path, config.subtype)

    # load dataset
    print('Load entire PBMC dataset.')
    adata = load_pbmc_data()
    adata = adata[adata.obs['celltype.l3']==config.subtype,:].copy()
    
    # preprocess dataset
    print(f'Preprocessing data with cell-type {config.subtype}.')
    adata = preprocess_dataset(adata,
                               params_cell_filtering=dict(mad_tot_cnt=5, 
                                                          mad_ngenes_cnt=5, 
                                                          nr_top_genes=20,
                                                          mad_pct_cnt_top_genes=5,
                                                          mad_pct_mt=5,
                                                          min_pct_mt=9),
                                )

    # get genes of last expression bin
    sig_genes = _get_genes_of_last_two_expression_bins(adata)

    # iterate overscoring methods and compute
    scoring_methods_params = get_scoring_method_params(SC_METHODS)

    for sc_method_full, (sc_method, sc_params) in scoring_methods_params.items():
        print(f'Scoring with {sc_method_full}.')
        # if we are scoring with ANS we need to score for all genes in the expression bin (incl. c/2 largest genes)
        # plus we need to define th gene
        if sc_method_full == "adjusted_neighborhood_scoring":
            sc_params["remove_genes_with_invalid_control_set"] = False

        scores = _compute_score_matrix(
            adata, sig_genes, sc_method_full, sc_method, sc_params
        )
        scores.to_csv(os.path.join(storing_path, f"scores_{sc_method_full}.csv"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run dataset composition experiment.")

    # cell sub CONFIGURATION
    parser.add_argument(
        "--subtype",
        choices=ALLOWED_CELL_SUBTYPES,
        default="Platelet",
        help="Indicate on which cell subtype to run the experiment.",
    )
    # STORING CONFIGURATION
    parser.add_argument(
        "--base_storing_path",
        default=os.path.join(BASE_PATH_EXPERIMENTS, "control_genes_selection/mean_var_per_gene_scores"),
        help="Path where to store the results.",
    )

    args = AttributeDict(vars(parser.parse_args()))
    
    main(args)