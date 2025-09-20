import argparse
import json
import os
import sys

import decoupler as dc
import matplotlib.pyplot as plt
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

sys.path.append('..')
sys.path.append('../..')

from constants import CANCER_DATASETS, NORM_METHODS, BASE_PATH_DGEX_CANCER
from load_data import load_datasets
from experiments.experiment_utils import AttributeDict

# Global settings
sc.settings.verbosity = 2
plt.rcParams.update({'pdf.fonttype': 42, 'font.family': 'sans-serif', 'font.sans-serif': 'Arial', 'font.size': 14})


def get_storing_path(base_storing_path, dataset, norm_method):
    storing_path = os.path.join(base_storing_path,
                                dataset,
                                f'{norm_method}_norm',
                                'on_pseudobulk')
    if not os.path.exists(storing_path):
        os.makedirs(storing_path)
    sc.logging.info(f'> Selected storing path {storing_path}')
    return storing_path


def main(config):
    # get storing path
    sc.logging.info(f'Get storing path.')
    storing_path = get_storing_path(BASE_PATH_DGEX_CANCER, config.dataset, config.norm_method)

    sc.logging.info(f'Storing experiment configuration...')
    with open(os.path.join(storing_path, 'non_rel_genes_config.txt'), "w") as fp:
        json.dump(config, fp, indent=4)

        # load data
    sc.logging.info(f'Load dataset.')
    adata = load_datasets(config.dataset, preprocessed=True, norm_method=config.norm_method)
    adata.uns['log1p']['base'] = None

    adata.layers['normalized'] = adata.X.copy()

    # to avoid bug in python deseq2 impementation rename column of interest 
    sc.logging.info(
        f'Map \'malignant_key\' to \'condition\' feature with \'malignant\'=\'B\' and \'non-malignant\'=\'A\'')
    adata.obs['condition'] = adata.obs['malignant_key'].astype(str).copy()
    adata.obs['condition'] = adata.obs['condition'].map({
        'malignant': 'B',
        'non-malignant': 'A',
    })
    adata.obs['condition'] = adata.obs['condition'].astype('category')

    # Get pseudo-bulk profile
    sc.logging.info(f'Get pseudo-bulk profiles per sample and malignancy.')
    pdata = dc.get_pseudobulk(adata,
                              sample_col='sample_id',
                              groups_col='condition',
                              layer='counts',
                              mode='sum',
                              min_cells=10,
                              min_counts=1000,
                              )

    sc.logging.info(f'Filter genes by expression: min_count=10, min_total_count=15.')
    nr_genes_before_filtering = pdata.shape[1]
    genes = dc.filter_by_expr(pdata, group='condition', min_count=10, min_total_count=15)
    nr_genes_after_filtering = len(genes)

    # Filter by these genes
    sc.logging.info(f'> Filtering {nr_genes_before_filtering - nr_genes_after_filtering} genes (before='
                    f'{nr_genes_before_filtering}, after={nr_genes_after_filtering})')
    pdata = pdata[:, genes].copy()

    # Build DESeq2 object
    sc.logging.info(f'Build DESeq2 object.')
    dds = DeseqDataSet(
        adata=pdata,
        design_factors='condition',
        refit_cooks=True,
        n_cpus=config.ncpu,
    )

    # Compute LFCs
    sc.logging.info(f'Compute LFCs (natural scale).')
    dds.deseq2()

    sc.logging.info(f'Compute DeseqStats.')
    # stat_res = DeseqStats(dds, n_cpus=config.ncpu, joblib_verbosity=0)
    stat_res = DeseqStats(dds)

    # Compute Wald test
    sc.logging.info(f'Compute Wald test.')
    stat_res.summary()

    # Shrink LFCs
    sc.logging.info(f'Shrink LFCs.')
    stat_res.lfc_shrink(coeff='condition_B_vs_A')

    sc.logging.info(f'LogFC are now shrunken {stat_res.shrunk_LFCs}')

    # Extract results
    sc.logging.info(f'Extract results.')
    results_df = stat_res.results_df

    # Getting important genes
    sc.logging.info(f'Select genes with |log2FoldChange|<={config.max_abs_log2fc} and padj > {config.min_padj}.')
    dgex_genes = results_df[
        (results_df.padj > config.min_padj) & (results_df.log2FoldChange.abs() <= config.max_abs_log2fc)]
    dgex_genes = dgex_genes.sort_values(by=['padj', 'log2FoldChange'], ascending=[True, False])
    dgex_genes = dgex_genes.reset_index(names='genes')
    dgex_genes.to_csv(os.path.join(storing_path, f'non_rel_genes.csv'), index=False)
    results_df.sort_values(by=['padj', 'log2FoldChange'], ascending=[True, False]).reset_index(names='genes').to_csv(
        os.path.join(storing_path, f'all_genes.csv'), index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run DGEX on given dataset to non-relevant genes for malignant vs. '
                                                 'non-malignant distinction.')

    # DATASET CONFIGURATION
    parser.add_argument('--dataset', choices=CANCER_DATASETS,
                        default='crc', help='indicate dataset on which to create signature.')
    parser.add_argument('--norm_method', choices=NORM_METHODS,
                        default='mean', help='Indicate normalization method used in the preprocessing.')
    # DGEX CONFIGURATION
    parser.add_argument('--max_abs_log2fc', default=0.5, type=float, help='Minimum log2FC for DGEX.')
    parser.add_argument('--min_padj', default=0.01, type=float, help='Maximum adjusted P-value for DGEX.')

    # COMPUTATION CONFIGURATION
    parser.add_argument('--ncpu', default=8, type=int, help='Number of CPUs.')

    # STORING CONFIGURATION
    parser.add_argument('--compute_umaps', action='store_true', help='Whether to compute and store the dataset UMAP.')

    args = AttributeDict(vars(parser.parse_args()))

    sc.logging.info(f'Creating malignant signatures with pseudobulks with the following '
                    f'configuration:\n{json.dumps(args, indent=4)}')

    main(args)
