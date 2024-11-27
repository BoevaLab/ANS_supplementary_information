#!/usr/bin/env python
# coding: utf-8
import argparse
import json
import os.path
import sys
from datetime import datetime

import numpy as np
import scanpy as sc
from scipy.sparse import diags

sys.path.append('..')
from experiments.experiment_utils import AttributeDict
from constants import DATASETS, NORM_METHODS, BASE_PATH_RAW_CANCER, BASE_PATH_CANSIG_PP_CANCER, BASE_PATH_PREPROCESSED

sc.settings.verbosity = 2  # verbosity: errors (0), warnings (1), info (2), hints (3)


def shifted_transformation(adata, y0=1):
    """
    From Twitter post https://twitter.com/Sanbomics/status/1647654042749874177?s=20
    Refering to publication by Ahlmann-Eltze & Huber.
    
    Ahlmann-Eltze, C., Huber, W. Comparison of transformations for single-cell RNA-seq data. 
    Nat Methods (2023). https://doi.org/10.1038/s41592-023-01814-1
    """
    target_sum = np.mean(adata.X.sum(axis=1))
    sc.logging.info(f'Mean shift logarithm normalization with normalization target count {target_sum}')
    size_factors = adata.X.sum(axis=1) / target_sum

    adata.X = diags(1 / size_factors.A1).dot(adata.X)
    adata.X.data = np.log(adata.X.data + y0)
    adata.uns["log1p"] = {"base": None}

    return adata


def filtergenes(adata, pct=0.01):
    # remove genes that are not present in at least 1% of all cells
    nr_cells, nr_genes = adata.shape
    gene_expr_in_cells_cnts = adata.X.getnnz(axis=0)
    enough_genes = gene_expr_in_cells_cnts - nr_cells * pct
    sc.logging.info(f'Filtering {np.sum(enough_genes < 0)} of {nr_genes} genes'
                    f'({np.round((np.sum(enough_genes < 0)) / nr_genes * 100, decimals=2)}%).')
    adata = adata[:, enough_genes >= 0].copy()
    return adata


def preprocess_data(adata, filter_genes=True, shift_method='mean'):
    if filter_genes:
        # Since we might have removed cells we need to refilter the genes, as they are filtered based on the
        # percentage of available cells in the data
        adata = filtergenes(adata)

    if shift_method == 'median':
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        adata.uns['log1p']['base'] = None
    elif shift_method == 'mean':
        adata = shifted_transformation(adata)
    elif shift_method == 'CP10k':
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.uns['log1p']['base'] = None
    else:
        raise ValueError('Unknown shift transformation method! Can choose between mean, median, CP10k.')

    return adata


def preprocess_dataset(adata, filter_genes=True, shift_method='mean', sample_based=False, sample_col='sample_id'):
    """
    Preprocessing pipeline for  a given dataset pre-preprocessed with [1], filter genes as in [1] and further
    normalize with [2] implemented in [3]

    [1] CanSig: Discovering de novo shared transcriptional programs in single cancer cells Josephine Yates,
    Florian Barkmann, Paweł Czyż, Marc Glettig, Frederieke Lohmann, Richard von der Horst, Elia Saquand,
    Nicolas Volken, Agnieszka Kraft, Valentina Boeva, bioRxiv 2022.04.14.488324;
    doi:https://doi.org/10.1101/2022.04.14.488324

    [2] Ahlmann-Eltze, C., Huber, W. Comparison of transformations for single-cell RNA-seq data. Nat Methods (
    2023). https://doi.org/10.1038/s41592-023-01814-1
    [3] https://twitter.com/Sanbomics/status/1647654042749874177?s=20

    Args:
        adata: loaded CanSig preprocesed dataset
        filter_genes: indicates whether genes need to be filtered
        shift_method: indicates which method should be used for the shift-algorithm
        sample_based: indicates if the dataset should be preprocessed per sample (True) or on the entire dataset
        sample_col: indicates the column, where sample identifiers are stored

    Returns:
        Preprocessed adata.
    """

    # sc.logging.info input configuration
    sc.logging.info(f'filter_genes={filter_genes}, shift_method={shift_method}, sample_based={sample_based}')
    # remove cells that were undecided in malignancy from CanSig pipeline
    if 'malignant_key' in adata.obs:
        adata = adata[adata.obs.malignant_key != 'undecided', :].copy()
        adata.obs.malignant_key = adata.obs.malignant_key.astype('category')

    adata.layers["counts"] = adata.X

    if sample_based:
        adatas = {}
        for group in adata.obs.groupby(sample_col):
            adatas[group[0]] = adata[group[1].index,].copy()
        del adata
        for key, curr_adata in adatas.items():
            adatas[key] = preprocess_data(curr_adata, filter_genes, shift_method)

        adata = sc.concat(list(adatas.values()), join='outer', merge='first')
        del adatas
    else:
        adata = preprocess_data(adata, filter_genes, shift_method)
    if 'mt' in adata.var:
        adata.var.mt = adata.var.mt.astype(bool)
    if 'cnv_called' in adata.var:
        adata.var.cnv_called = adata.var.cnv_called.astype(bool)

    return adata


def get_appendix(normalization_method, pp_sample_based):
    if normalization_method == 'median':
        appendix = '_med_per_sid' if pp_sample_based else '_med'
    elif normalization_method == 'CP10k':
        appendix = '_cp10k_per_sid' if pp_sample_based else '_cp10k'
    else:
        appendix = '_per_sid' if pp_sample_based else ''
    sc.logging.info(appendix)
    return appendix


def main(config):
    appendix = get_appendix(config.norm_method, config.sample_based)
    if config.dataset in ['ovarian_malignant', 'skin_malignant']:
        fn_data = os.path.join(BASE_PATH_RAW_CANCER, f'{config.dataset}.h5ad')
    else:
        fn_data = os.path.join(BASE_PATH_CANSIG_PP_CANCER, f'{config.dataset}.h5ad')
    fn_output = os.path.join(BASE_PATH_PREPROCESSED, f'pp_{config.dataset}{appendix}.h5ad')

    adata = sc.read_h5ad(fn_data)
    #if config.dataset == 'luad_xing':
    #    adata = adata[adata.obs.sample_id.str.startswith('SSN')].copy()

    adata = preprocess_dataset(adata, shift_method=config.norm_method, sample_based=config.sample_based)

    adata.write(fn_output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Proprocess a cancer dataset.')

    parser.add_argument('--dataset', choices=DATASETS,
                        default='crc',
                        help="""Indicate which cancer dataset to preprocess. We have the following: 
                        
                        """)
    parser.add_argument('--norm_method', choices=NORM_METHODS,
                        default='mean', help='Indicate normalization method used in the preprocessing.')
    parser.add_argument('--sample_based', action='store_true', help='Preprocess each sample individually, '
                                                                    'else the entire dataset.')

    args = AttributeDict(vars(parser.parse_args()))

    start = datetime.now()
    sc.logging.info(f'Preprocessing with following configuration:\n{json.dumps(args, indent=4)}')

    main(args)

    sc.logging.info(f'Finished experiment in total {datetime.now() - start} time.')
