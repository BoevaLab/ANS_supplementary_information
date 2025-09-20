import os
import sys
from pathlib import Path
from typing import Dict, List

import decoupler as dc
import pandas as pd
import scanpy as sc

sys.path.append('.')
PARENT_PATH = Path(__file__).parent
sys.path.append(str(PARENT_PATH))


from constants import (
# from constants import (
    # Paths
    BASE_PATH_RAW_CANCER,
    BASE_PATH_RAW_PBMC,
    BASE_PATH_PREPROCESSED,
    BASE_PATH_DGEX_CANCER,
    BASE_PATH_ANNOT_CANCER,
    BASE_PATH_CANSIG_PP_CANCER,
    PBMC_DEXG,

    # Dataset definitions
    DATASETS,
    CANCER_DATASETS,
    PBMC_DATASETS,
    DATASETS_WITH_ANNOTATIONS,

    # Configuration
    CANCER_DATASET_SIGS_CONFIGS,
    PBMC_DATASET_SIGS_CONFIGS,

    # Parameters
    NORM_METHODS,
    DGEX_GRANULARITY,
    PCTGS,
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


def load_datasets(dataset_name, preprocessed=True, cansig_pp=True, norm_method='mean', sample_based=False):
    assert dataset_name in DATASETS, f'Unknown dataset_name {dataset_name}. Allowed datasets {DATASETS}.'

    if preprocessed:
        appendix = get_appendix(norm_method, sample_based)
        if appendix is None:
            appendix = ''

    if dataset_name in PBMC_DATASETS:
        fn = os.path.join(BASE_PATH_PREPROCESSED,
                          f'pp_{dataset_name}{appendix}.h5ad') if preprocessed else BASE_PATH_RAW_PBMC
    elif dataset_name in CANCER_DATASETS:
        if preprocessed:
            fn = os.path.join(BASE_PATH_PREPROCESSED, f'pp_{dataset_name}{appendix}.h5ad')
        else:
            if cansig_pp:
                fn = os.path.join(BASE_PATH_CANSIG_PP_CANCER, f'{dataset_name}.h5ad')
            else:
                fn = os.path.join(BASE_PATH_RAW_CANCER, f'{dataset_name}.h5ad')
    else:
        raise KeyError(f'Unknown dataset_name {dataset_name}. Allowed datasets {DATASETS}.')

    sc.logging.info(f'Load datasets with path {fn}')
    adata = sc.read_h5ad(fn)

    return adata


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


def _remove_overlapping_signature_genes(sig_genes):
    sig_genes_old = sig_genes.copy()
    genes_removed = set()
    for cell_type, genes in sig_genes.items():
        for cell_type_other, genes_other in sig_genes_old.items():
            if cell_type == cell_type_other:
                continue
            intersecting_genes = set(genes) & set(genes_other)
            if len(intersecting_genes) > 0:
                genes = list(set(genes) - intersecting_genes)
                genes_removed.update(intersecting_genes)
        sig_genes[cell_type] = genes

    print(f"Removed {genes_removed} overlapping genes.")
    for key, value in sig_genes.items():
        print(f"Signature {key} length: before={len(sig_genes_old[key])}, after={len(value)}")

    return sig_genes


def _load_pbmc_dex_genes(
        filter_genes: bool,
        threshold_pval: float = 0.01,
        threshold_log2fc: float = 0.5
) -> pd.DataFrame:
    dge_genes = pd.read_csv(PBMC_DEXG)

    if filter_genes:
        print('Shape DEX genes BEFORE filtering', dge_genes.shape)
        dge_genes = dge_genes[
            (dge_genes['P Value'] <= threshold_pval) &
            (dge_genes['Average Log Fold Change'] >= threshold_log2fc)
            ].copy().reset_index(drop=True)
        print('Shape DEX genes AFTER filtering', dge_genes.shape)

    return dge_genes


def _load_pbmc_signatures(
        dataset: str,
        filter_genes: bool = True,
        threshold_pval: float = 0.01,
        threshold_log2fc: float = 0.5
):
    DE_of_celltypes = _load_pbmc_dex_genes(
        filter_genes=filter_genes,
        threshold_pval=threshold_pval,
        threshold_log2fc=threshold_log2fc
    )
    high_to_low_lbl_map = PBMC_DATASET_SIGS_CONFIGS[dataset].high_level_low_level_mapping

    signatures = {}
    for row in high_to_low_lbl_map.items():
        cell_type, subtypes = row
        all_genes_for_cell_type = DE_of_celltypes[DE_of_celltypes['Cell Type'].isin(subtypes)]['Gene'].tolist()
        signatures[cell_type] = sorted(list(set(all_genes_for_cell_type)))

    return signatures


def _process_csv_signatures(file_path: Path):
    """Process CSV signature files into a standardized format."""
    signatures = pd.read_csv(file_path)
    signatures = signatures.to_dict('series')  # Convert to series first
    return {k: sorted(v.dropna().tolist()) for k, v in signatures.items()}


def _process_gmt_signatures(file_path: Path) -> Dict[str, List[str]]:
    """Process GMT signature files into a standardized format."""
    signatures = dc.read_gmt(file_path)
    signatures = signatures.groupby('source').target.apply(lambda x: sorted(x.unique()))
    return signatures.to_dict()


def _load_cancer_signatures(dataset: str) -> Dict[str, List[str]]:
    if dataset not in CANCER_DATASET_SIGS_CONFIGS:
        raise ValueError(f"Unknown dataset. Please choose from: {', '.join(list(CANCER_DATASET_SIGS_CONFIGS.keys()))}")

    config = CANCER_DATASET_SIGS_CONFIGS[dataset]
    file_path = Path(BASE_PATH_ANNOT_CANCER) / config.file_path

    if not file_path.exists():
        raise FileNotFoundError(f"Signature file not found: {file_path}")

    pp_fun = _process_gmt_signatures if config.file_type == 'gmt' else _process_csv_signatures
    signatures = pp_fun(file_path)

    if config.name_transform:
        signatures = {config.name_transform(k): v for k, v in signatures.items()}

    return signatures


def load_signatures(
        dataset: str,
        remove_overlapping_genes: bool = False,
        **pbmc_kwargs,
):
    if dataset not in DATASETS_WITH_ANNOTATIONS:
        raise ValueError(f"Unknown datasets. Please provide a datasets from {DATASETS_WITH_ANNOTATIONS}")

    if dataset.startswith('pbmc'):
        sigs = _load_pbmc_signatures(dataset, **pbmc_kwargs)
    else:
        sigs = _load_cancer_signatures(dataset)

    if remove_overlapping_genes:
        sigs = _remove_overlapping_signature_genes(sigs)
    return sigs
