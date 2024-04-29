import os
import sys
import pandas as pd
import scanpy as sc

sys.path.append('..')
sys.path.append('.')

from .constants import (
    BASE_PATH_RAW_CANCER,
    BASE_PATH_RAW_PBMC,
    BASE_PATH_PREPROCESSED,
    BASE_PATH_DGEX_CANCER,
    DATASETS,
    CANCER_DATASETS,
    NORM_METHODS,
    DGEX_GRANULARITY,
    PCTGS
)


def get_appendix(norm_method, sample_based):
    if norm_method not in NORM_METHODS:
        raise KeyError(f'Unknown normalization method {norm_method}. Allowed values [mean, median, CP10k] ')

    if norm_method == 'median':
        appendix = '_med_per_sid' if sample_based else '_med'
    elif norm_method == 'CP10k':
        appendix = '_cp10k_per_sid' if sample_based else '_cp10k'
    else:
        appendix = '_per_sid' if sample_based else ''
    return appendix


def get_fn_dge(norm_method, sample_based, dge_on_all, intersect_pct, min_log2fc, pval, ranked_means=True):
    if dge_on_all == 'all':
        fn = os.path.join(f'{norm_method}_norm',
                          f'dgex_on_all_sid',
                          f'min_log2fc_{min_log2fc}_pval_{pval}',
                          f'dgex_genes.csv')
    elif dge_on_all == 'pseudobulk':
        fn = os.path.join(f'{norm_method}_norm',
                          f'on_pseudobulk',
                          f'dgex_genes.csv')
    else:
        fn = os.path.join(f'{norm_method}_norm',
                          f'dgex_on_each_sid',
                          f'min_log2fc_{min_log2fc}_pval_{pval}',
                          'dgex_genes_ranked_means.csv' if ranked_means else f'dgex_genes_intersec'
                                                                             f'_{int(round(intersect_pct * 100))}_psid.csv'
                          )
    return fn


def load_datasets(dataset_name, preprocessed=True, norm_method='mean', sample_based=False):
    assert dataset_name in DATASETS, f'Unknown dataset_name {dataset_name}. Allowed datasets {DATASETS}.'

    if preprocessed:
        appendix = get_appendix(norm_method, sample_based)
        if appendix is None:
            appendix = ''

    if dataset_name == 'pbmc_b_mono_nk':
        fn = os.path.join(BASE_PATH_PREPROCESSED,
                          f'pp_pbmc_b_mono_nk{appendix}.h5ad') if preprocessed else BASE_PATH_RAW_PBMC
    elif dataset_name == 'pbmc_b_subtypes':
        fn = os.path.join(BASE_PATH_PREPROCESSED,
                          f'pp_pbmc_b_subtypes{appendix}.h5ad') if preprocessed else BASE_PATH_RAW_PBMC
    elif dataset_name == 'pbmc_cd4_subtypes':
        fn = os.path.join(BASE_PATH_PREPROCESSED,
                          f'pp_pbmc_cd4_subtypes{appendix}.h5ad') if preprocessed else BASE_PATH_RAW_PBMC
    elif dataset_name == 'pbmc_cd8_subtypes':
        fn = os.path.join(BASE_PATH_PREPROCESSED,
                          f'pp_pbmc_cd8_subtypes{appendix}.h5ad') if preprocessed else BASE_PATH_RAW_PBMC
    elif dataset_name == 'pbmc':
        assert not preprocessed, 'No preprocessed version available for pbmc dataset.'
        fn = BASE_PATH_RAW_PBMC
    else:
        fn = os.path.join(BASE_PATH_PREPROCESSED,
                          f'pp_{dataset_name}{appendix}.h5ad') if preprocessed else os.path.join(BASE_PATH_RAW_CANCER,
                                                                                                 f'{dataset_name}.h5ad')
    sc.logging.info(f'Load datasets with path {fn}')
    return sc.read_h5ad(fn)


def load_dgex_genes_for_mal_cells(dataset_name,
                                  norm_method='mean',
                                  sample_based=False,
                                  dge_on_all='pseudobulk',
                                  intersect_pct=0.9,
                                  min_log2fc=1,
                                  pval=0.05,
                                  ranked_means=True):
    assert dataset_name in CANCER_DATASETS, f'Unknown dataset_name {dataset_name}. Allowed datasets {CANCER_DATASETS}.'
    assert dge_on_all in DGEX_GRANULARITY, f'The variable \'dge_on_all\' needs to be in {DGEX_GRANULARITY}.'
    if (dge_on_all == 'individual') and (intersect_pct not in PCTGS):
        raise KeyError(f'The argument \'intersect_pct\' ({intersect_pct}) requires to be one in {PCTGS}.')

    filename = get_fn_dge(norm_method, sample_based, dge_on_all, intersect_pct, min_log2fc, pval, ranked_means)
    fn = os.path.join(BASE_PATH_DGEX_CANCER, dataset_name, filename)
    sc.logging.info(f'Load dgex_genes with path {fn}')
    return pd.read_csv(fn)


def load_non_rel_genes_for_mal_cells(dataset_name,
                                     norm_method='mean',
                                     sample_based=False,
                                     dge_on_all='pseudobulk'):
    assert dataset_name in CANCER_DATASETS, f'Unknown dataset_name {dataset_name}. Allowed datasets {CANCER_DATASETS}.'
    filename = os.path.join(f'{norm_method}_norm', f'on_pseudobulk', f'all_genes.csv')
    fn = os.path.join(BASE_PATH_DGEX_CANCER, dataset_name, filename)
    sc.logging.info(f'Load dgex_genes with path {fn}')
    return pd.read_csv(fn)
