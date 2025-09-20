import argparse
import json
import os
import sys
import textwrap
import warnings
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from matplotlib.patches import Patch
from scipy.stats import mannwhitneyu, shapiro, normaltest, levene, ttest_ind
from signaturescoring import score_signature
from signaturescoring.utils.utils import get_mean_and_variance_gene_expression
from statannotations.Annotator import Annotator

sys.path.append('../..')
from data.load_data import load_datasets, load_dgex_genes_for_mal_cells
from data.constants import BASE_PATH_EXPERIMENTS, SCORING_METHODS, CANCER_DATASETS, NORM_METHODS, METHOD_WO_MEAN, DGEX_GRANULARITY
from experiments.experiment_utils import AttributeDict

# Global settings
sc.settings.verbosity = 2
plt.rcParams.update({'pdf.fonttype': 42, 'font.family': 'sans-serif', 'font.sans-serif': 'Arial', 'font.size': 14})
warnings.filterwarnings("ignore")


def wrap_labels(ax, width, break_long_words=False):
    """
    Method to wrap ticklabels to a certain length.
    Args:
        ax: Figure axis
        width: Desired max width of a label before breaking.
        break_long_words: Indicate whether long words should be broken.
    """
    labels = []
    for label in ax.get_xticklabels():
        text = label.get_text()
        labels.append(textwrap.fill(text, width=width,
                                    break_long_words=break_long_words))
    ax.set_xticklabels(labels, rotation=0)


def get_storing_path(base_storing_path, dataset, norm_method="mean", dge_on_all='pseudobulk', min_log2fc=2,
                     pval=0.01):
    """
    Method to create an experiment specific path.
    Args:
        base_storing_path: Base-experiment path.
        dataset: Datasetname --> see constants.py
        norm_method: Method with which data was normalized.
        dge_on_all: Indicates on which level the malignant signature was created (on 'all', 'individual' or
        'bulkified' samples)
        min_log2fc: Minimum log2FC required for the malignant signature.
        pval: Maximum log2FC required for the malignant signature

    Returns:
        Path where the experiment's data and figures will be stored.
    """
    if dge_on_all == 'all':
        storing_path = os.path.join(base_storing_path,
                                    dataset,
                                    f'{norm_method}_norm',
                                    f'dgex_on_all_sid',
                                    f'min_log2fc_{min_log2fc}_pval_{pval}')
    elif dge_on_all == 'pseudobulk':
        storing_path = os.path.join(base_storing_path,
                                    dataset,
                                    f'{norm_method}_norm',
                                    f'dgex_on_pseudobulk')
    else:
        storing_path = os.path.join(base_storing_path,
                                    dataset,
                                    f'{norm_method}_norm',
                                    f'dgex_on_each_sid',
                                    f'min_log2fc_{min_log2fc}_pval_{pval}')
    if not os.path.exists(storing_path):
        os.makedirs(storing_path)
    return storing_path


def get_malignant_signature(dataset, norm_method, dge_on_all, intersect_pct, min_log2fc, pval, ranked_means,
                            sort_values_by, nr_sig_genes, return_random=False):
    """
    Method to retrieve and prepare desired malignant signature.
    Args:
        dataset: Datasetname --> see constants.py
        norm_method: Method with which data was normalized.
        dge_on_all: Indicates on which level the malignant signature was created (on 'all', 'individual' or
        'bulkified' samples)
        intersect_pct (float, [0,1]): Used if dge_on_all=individual and ranked_means=False. Then a signature is
        selected, where a gene appears in at least intersect_pct of the samples.
        min_log2fc: Minimum log2FC required for the malignant signature.
        pval:  Maximum log2FC required for the malignant signature
        ranked_means: Used if dge_on_all=individual. Indicates whether the on individual samples computed signature
        is sorted by the mean log2FC rank (True) or only genes are selected that appear in at least intersect_pct of
        the samples.
        sort_values_by: Used if dge_on_all=individual and ranked_means=False. Then sort_values_by indicates if the
        signature should be sorted by mean or median log2FC
        nr_sig_genes: Desired signature length.
        return_random: Indicates if a random ordering of the differentially expressed genes should be returned
        additionally.

    Returns:
        If return_random=True, two signatures one with decreasing log2FC and with random gene ordering are returned.
        Otherwise, only one with decreasing log2FC is returned.

    """
    wc = load_dgex_genes_for_mal_cells(dataset,
                                       norm_method=norm_method,
                                       sample_based=False,
                                       dge_on_all=dge_on_all,
                                       intersect_pct=intersect_pct,
                                       min_log2fc=min_log2fc,
                                       pval=pval,
                                       ranked_means=ranked_means)
    if dge_on_all == 'all':
        wc = wc.sort_values(by='logfoldchanges', ascending=False)
        col = 'names'
    elif dge_on_all == 'pseudobulk':
        wc = wc.sort_values(by=['padj', 'log2FoldChange'], ascending=[True, False])
        col = 'genes'
    elif ranked_means:
        wc = wc.sort_values(by='mean_ranked_log2FC', ascending=True)
        col = 'names'
    else:
        wc = wc.sort_values(by=sort_values_by, ascending=False)
        col = 'names'

    gene_list_most_de = wc[0:nr_sig_genes][col].tolist()
    if return_random:
        gene_list_random_de = wc.sample(nr_sig_genes)[col].tolist()
        return gene_list_most_de, gene_list_random_de
    else:
        return gene_list_most_de


def score_independently(dataset, norm_method, sample_based, sig_most_dge, sig_random_dge=None):
    if sample_based:
        pp_appendix = '_ppsi'
    else:
        pp_appendix = '_ppas'

    adata = load_datasets(dataset, norm_method=norm_method, sample_based=sample_based)
    if 'log1p' not in adata.uns.keys():
        adata.uns["log1p"] = {"base": None}
    else:
        adata.uns['log1p']['base'] = None

    # We separate the dataset to score on each sample individually.
    adatas = {}
    for group in adata.obs.groupby('sample_id'):
        adatas[group[0]] = adata[group[1].index,].copy()
    del adata

    sc_names_most_dge_si = []
    if sig_random_dge is not None:
        sc_names_random_dge_si = []
    for key, adata in adatas.items():
        start = datetime.now()
        sc.logging.info(f'Running scoring on sample {key} ...')
        for sc_method in SCORING_METHODS:
            scoring_method = sc_method['scoring_method']
            sc_params = sc_method['sc_params'].copy()
            sc_params['score_name'] = sc_params['score_name'] + '_most_dge_si' + pp_appendix
            sc_names_most_dge_si.append(sc_params['score_name'])
            sc.settings.verbosity = 0

            score_signature(method=scoring_method,
                            adata=adata,
                            gene_list=sig_most_dge,
                            **sc_params)

            if sig_random_dge is not None:
                sc_params = sc_method['sc_params'].copy()
                # si_ppsi = score_individually_preprocess_each_sample_individually
                sc_params['score_name'] = sc_params['score_name'] + '_random_dge_si' + pp_appendix
                sc_names_random_dge_si.append(sc_params['score_name'])

                score_signature(adata=adata,
                                gene_list=sig_random_dge,
                                method=scoring_method,
                                **sc_params)

            sc.settings.verbosity = 2

        sc.logging.info(f'Duration scoring all methods {datetime.now() - start}.')

    sc_names_most_dge_si = list(set(sc_names_most_dge_si))
    if sig_random_dge is not None:
        sc_names_random_dge_si = list(set(sc_names_random_dge_si))

    adata = sc.concat(list(adatas.values()), join='outer', merge='same')

    scores_most_dge_si = adata.obs[sc_names_most_dge_si].copy()
    if sig_random_dge is not None:
        scores_random_dge_si = adata.obs[sc_names_random_dge_si].copy()
        return scores_most_dge_si, scores_random_dge_si
    else:
        return scores_most_dge_si


def score_together(dataset, norm_method, sample_based, sig_most_dge, sig_random_dge=None):
    adata = load_datasets(dataset, norm_method=norm_method, sample_based=sample_based)
    if 'log1p' not in adata.uns.keys():
        adata.uns["log1p"] = {"base": None}
    else:
        adata.uns['log1p']['base'] = None

    df_mean_var = get_mean_and_variance_gene_expression(adata, estim_var=False, show_plots=False, store_path=None)

    sc_names_most_dge_as = []
    if sig_random_dge is not None:
        sc_names_random_dge_as = []
    for sc_method in SCORING_METHODS:
        adata.X = adata.X.tocsc()
        start = datetime.now()
        scoring_method = sc_method['scoring_method']
        sc.logging.info(f'Running scoring with scoring method {scoring_method}')

        sc_params = sc_method['sc_params'].copy()
        sc_params['score_name'] = sc_params['score_name'] + '_most_dge_all_samples'
        sc_names_most_dge_as.append(sc_params['score_name'])
        sc.settings.verbosity = 0
        if scoring_method in METHOD_WO_MEAN:
            score_signature(method=scoring_method,
                            adata=adata,
                            gene_list=sig_most_dge,
                            **sc_params)
        else:
            score_signature(method=scoring_method,
                            adata=adata,
                            gene_list=sig_most_dge,
                            df_mean_var=df_mean_var,
                            **sc_params)

        if sig_random_dge is not None:
            sc_params = sc_method['sc_params'].copy()
            sc_params['score_name'] = sc_params['score_name'] + '_random_dge_all_samples'
            sc_names_random_dge_as.append(sc_params['score_name'])

            if scoring_method in METHOD_WO_MEAN:
                score_signature(method=scoring_method,
                                adata=adata,
                                gene_list=sig_random_dge,
                                **sc_params)
            else:
                score_signature(method=scoring_method,
                                adata=adata,
                                gene_list=sig_random_dge,
                                df_mean_var=df_mean_var,
                                **sc_params)
        sc.settings.verbosity = 2
        sc.logging.info(f'Duration scoring {scoring_method} {datetime.now() - start}.')

    return adata


def get_pvals_chemistry(adata):
    subnames = ['Scanpy', 'Seurat', 'ANS', 'Jasmine', 'UCell']
    cols = [col for col in adata.obs.columns if any(map(col.__contains__, subnames))]
    res = adata.obs.groupby(by=['sample_id', 'malignant_key', 'SINGLECELL_TYPE'])[cols].mean().reset_index()
    res = res.dropna(axis=0, how='all', subset=res.columns.tolist()[3:])
    res = res.melt(id_vars=['sample_id', 'SINGLECELL_TYPE', 'malignant_key'],
                   var_name='scoring_method',
                   value_name='scores')

    name_mapping = {'all_samples': 'Scoring all samples together',
                    'si_ppas': 'Scoring each sample individually (preprocessed together)',
                    'si_ppsi': 'Scoring each sample individually (preprocessed independently)',
                    }

    res['scoring_mode'] = res.scoring_method.apply(lambda x: name_mapping['_'.join(x.split('_')[-2:])])
    res['scoring_method'] = res.scoring_method.apply(lambda x: '_'.join(x.split('_')[0:-4]))

    pvals = []
    score_means = []
    alpha = 0.05
    for (key, df) in res.groupby(['scoring_method', 'scoring_mode']):
        mal_p2 = df.scores[(df.SINGLECELL_TYPE == 'SC3Pv2') & (df.malignant_key == 'malignant')].values
        mal_p3 = df.scores[(df.SINGLECELL_TYPE == 'SC3Pv3') & (df.malignant_key == 'malignant')].values

        are_normal_shapiro = shapiro(mal_p3).pvalue > alpha and shapiro(mal_p2).pvalue > alpha
        are_normal_normaltest = normaltest(mal_p3).pvalue > alpha and normaltest(mal_p2).pvalue > alpha
        variance_test = levene(mal_p3, mal_p2).pvalue > alpha
        pvals.append({'scoring_method': key[0],
                      'scoring_mode': key[1],
                      'MannWhitneyU p-val': mannwhitneyu(mal_p3, mal_p2).pvalue,
                      'MannWhitneyU statistic': mannwhitneyu(mal_p3, mal_p2).statistic,
                      'are_normal_shapiro': are_normal_shapiro,
                      'are_normal_normaltest': are_normal_normaltest,
                      'equal_variance': variance_test,
                      'ttest p-val': ttest_ind(mal_p3, mal_p2, equal_var=variance_test).pvalue,
                      'ttest statistic': ttest_ind(mal_p3, mal_p2, equal_var=variance_test).statistic,
                      })
        score_means.append({
            'scoring_method': key[0],
            'scoring_mode': key[1],
            'mal_p2_mean': np.mean(mal_p2),
            'mal_p2_var': np.var(mal_p2),
            'mal_p2_min': np.min(mal_p2),
            'mal_p2_max': np.max(mal_p2),
            'mal_p3_mean': np.mean(mal_p3),
            'mal_p3_var': np.var(mal_p3),
            'mal_p3_min': np.min(mal_p3),
            'mal_p3_max': np.max(mal_p3),
        })

    pvals = pd.DataFrame(pvals)
    score_means = pd.DataFrame(score_means)
    return pvals, score_means


def get_stats_over_modes(adata):
    subnames = ['Scanpy', 'Seurat', 'ANS', 'Jasmine', 'UCell']
    cols = [col for col in adata.obs.columns if any(map(col.__contains__, subnames))]
    stat_methods = ['mean', 'var', 'std', 'min', 'max', 'median']
    name_mapping = {'all_samples': 'Scoring all samples together',
                    'si_ppas': 'Scoring each sample individually (preprocessed together)',
                    'si_ppsi': 'Scoring each sample individually (preprocessed independently)',
                    }
    aggregated_data = adata.obs.groupby('malignant_key')[cols].agg({col: stat_methods for col in cols})
    aggregated_data = aggregated_data.T
    aggregated_data = aggregated_data.reset_index()
    aggregated_data.columns = ['scoring_method', 'stat_method', 'malignant', 'non-malignant']
    aggregated_data['scoring_mode'] = aggregated_data.scoring_method.apply(lambda x: name_mapping['_'.join(x.split('_')[-2:])])
    aggregated_data['scoring_method'] = aggregated_data.scoring_method.apply(lambda x: '_'.join(x.split('_')[0:-4]))
    return aggregated_data


def plot_pval_heatmaps(pvals, storing_path):
    col_order = ['ANS', 'Seurat', 'Seurat_AG', 'Seurat_LVG', 'Scanpy', 'Jasmine_LH', 'Jasmine_OR', 'UCell']
    plt.figure(figsize=(6, 10))
    sns.heatmap(pvals[pvals.scoring_mode == 'Scoring all samples together'].set_index('scoring_method')[
                    ['MannWhitneyU p-val', 'ttest p-val']].loc[col_order], annot=True)
    plt.yticks(rotation=0);
    plt.title('Scoring all samples together');
    plt.savefig(os.path.join(storing_path, f'pvals_mal_means_chem_not_corrected_all.svg'), format='svg')

    plt.figure(figsize=(6, 10))
    sns.heatmap(pvals[pvals.scoring_mode == 'Scoring each sample individually (preprocessed together)'].set_index(
        'scoring_method')[['MannWhitneyU p-val', 'ttest p-val']].loc[col_order], annot=True)
    plt.yticks(rotation=0);
    plt.title('Scoring each sample individually (preprocessed together)');
    plt.savefig(os.path.join(storing_path, f'pvals_mal_means_chem_not_corrected_ppas.svg'), format='svg', bbox_inches='tight')

    plt.figure(figsize=(6, 10))
    sns.heatmap(pvals[pvals.scoring_mode == 'Scoring each sample individually (preprocessed independently)'].set_index(
        'scoring_method')[['MannWhitneyU p-val', 'ttest p-val']].loc[col_order], annot=True)
    plt.yticks(rotation=0);
    plt.title('Scoring each sample individually (preprocessed independently)');
    plt.savefig(os.path.join(storing_path, f'pvals_mal_means_chem_not_corrected_ppsi.svg'), format='svg', bbox_inches='tight')


def create_strip_plot_crc(adata, exclude_dge_si_ppsi=False, exclude_jasmine=False):
    plt.rcParams.update({'pdf.fonttype': 42, 'font.family': 'sans-serif', 'font.sans-serif': 'Arial', 'font.size': 14})
    means_std_per_sample = []
    for (key, adata_obs) in adata.obs.groupby('sample_id'):
        curr_means = adata_obs.groupby(by=['SINGLECELL_TYPE', 'malignant_key']).mean()
        curr_means['sample_id'] = key
        curr_means['scoring_mode'] = 'Scoring on all samples' if exclude_dge_si_ppsi else 'Scoring on all samples (' \
                                                                                          'preprocessed together)'
        means_std_per_sample.append(curr_means)
    scmeans_per_chem_mal = pd.concat(means_std_per_sample).reset_index()
    scmeans_per_chem_mal = scmeans_per_chem_mal.drop(columns=['TumorSize', 'SizeQuantile', 'Age',
                                                              'n_counts', 'n_genes_by_counts',
                                                              'total_counts', 'total_counts_mt',
                                                              'pct_counts_mt', 'S_score', 'G2M_score',
                                                              'iCMS2_GT', 'iCMS3_GT',
                                                              'log_counts', ])
    scmeans_per_chem_mal.dropna(inplace=True)
    # change structure of data
    scmeans_per_chem_mal = scmeans_per_chem_mal.melt(
        id_vars=['sample_id', 'SINGLECELL_TYPE', 'malignant_key', 'scoring_mode'],
        var_name='scoring_method',
        value_name='scores'
    )
    if exclude_dge_si_ppsi:
        scmeans_per_chem_mal = scmeans_per_chem_mal[~scmeans_per_chem_mal.scoring_method.str.contains('si_ppsi')].copy()
        scmeans_per_chem_mal.loc[scmeans_per_chem_mal.scoring_method.str.contains(
            'si_ppas'), 'scoring_mode'] = 'Scoring samples individually'
    else:
        scmeans_per_chem_mal.loc[scmeans_per_chem_mal.scoring_method.str.contains(
            'si_ppsi'), 'scoring_mode'] = 'Scoring samples individually (preprocessed individually)'
        scmeans_per_chem_mal.loc[scmeans_per_chem_mal.scoring_method.str.contains(
            'si_ppas'), 'scoring_mode'] = 'Scoring samples individually (preprocessed together)'

    scmeans_per_chem_mal['mal_nd_chem'] = scmeans_per_chem_mal['malignant_key'].astype(str) + \
                                          ' with ' + scmeans_per_chem_mal['SINGLECELL_TYPE'].astype(str)
    scmeans_per_chem_mal['scmethod_nd_chem'] = scmeans_per_chem_mal['scoring_mode'].astype(
        str) + ' with ' + scmeans_per_chem_mal['SINGLECELL_TYPE'].astype(str)

    most_means = scmeans_per_chem_mal[
        scmeans_per_chem_mal.scoring_method.str.contains('most_dge')].copy()
    # random_means = sample_means_per_chem_nd_mal[
    #     sample_means_per_chem_nd_mal.scoring_method.str.contains('random_dge')].copy()

    most_means.scoring_method = most_means.scoring_method.apply(lambda x: '_'.join(x.split('_')[0:-4]))
    # random_means.scoring_method = random_means.scoring_method.apply(lambda x: '_'.join(x.split('_')[0:-4]))

    most_means = most_means.sort_values(by='scoring_method')
    # random_means = random_means.sort_values(by='scoring_method')

    mean_of_sid_means = most_means.groupby(
        by=['scoring_method', 'scmethod_nd_chem', 'malignant_key']).mean().reset_index()

    sns.set_style("ticks")
    if exclude_dge_si_ppsi:
        order_x = ['Scoring samples individually with SC3Pv2',
                   'Scoring samples individually with SC3Pv3',
                   'Scoring on all samples with SC3Pv2',
                   'Scoring on all samples with SC3Pv3']
        height = 5
        aspect = 0.6
        pairs = [(('Scoring samples individually with SC3Pv2', 'malignant'),
                  ('Scoring samples individually with SC3Pv3', 'malignant')),
                 (('Scoring on all samples with SC3Pv2', 'malignant'),
                  ('Scoring on all samples with SC3Pv3', 'malignant')), ]
    else:
        order_x = ['Scoring samples individually (preprocessed individually) with SC3Pv2',
                   'Scoring samples individually (preprocessed individually) with SC3Pv3',
                   'Scoring samples individually (preprocessed together) with SC3Pv2',
                   'Scoring samples individually (preprocessed together) with SC3Pv3',
                   'Scoring on all samples (preprocessed together) with SC3Pv2',
                   'Scoring on all samples (preprocessed together) with SC3Pv3']
        height = 6
        aspect = 0.75
        pairs = [(('Scoring samples individually (preprocessed individually) with SC3Pv2', 'malignant'),
                  ('Scoring samples individually (preprocessed individually) with SC3Pv3', 'malignant')),
                 (('Scoring samples individually (preprocessed together) with SC3Pv2', 'malignant'),
                  ('Scoring samples individually (preprocessed together) with SC3Pv3', 'malignant')),
                 (('Scoring on all samples (preprocessed together) with SC3Pv2', 'malignant'),
                  ('Scoring on all samples (preprocessed together) with SC3Pv3', 'malignant'))]

    if exclude_jasmine:
        col_order = ['ANS', 'Seurat', 'Seurat_AG', 'Seurat_LVG', 'Scanpy', 'UCell']
        col_wrap = 3
    else:
        col_order = ['ANS', 'Seurat', 'Seurat_AG', 'Seurat_LVG', 'Scanpy', 'Jasmine_LH', 'Jasmine_OR', 'UCell']
        col_wrap = 4

    args = dict(
        x='scmethod_nd_chem',
        order=order_x,
        y='scores',
        hue='malignant_key',
        col='scoring_method',
        col_order=col_order,
        col_wrap=col_wrap,
        kind='strip',
        height=height,
        aspect=aspect,
        legend=True
    )

    g = sns.catplot(
        data=most_means,
        dodge=False,
        **args
    )

    for name, ax in g.axes_dict.items():
        # Subset the data based on the 'scoring_method' column
        subset_data = most_means.loc[most_means['scoring_method'] == name, :].copy()

        annot = Annotator(ax, pairs, **args, data=subset_data)
        annot.configure(test='Mann-Whitney', loc='inside', verbose=0)
        annot.apply_test().annotate()

    g.add_legend(fontsize=16, loc='upper right', bbox_to_anchor=(1, 1))
    g.legend.set_title('Celltype', prop={'size': 18})

    palette = {'malignant': 'black', 'non-malignant': 'black'}
    for i, ax in enumerate(g.axes):
        ax.axvspan(xmin=-0.35, xmax=0.35, facecolor='grey', alpha=0.5, zorder=0)
        ax.axvspan(xmin=0.65, xmax=1.35, facecolor='tan', alpha=0.5, zorder=0)
        if not exclude_dge_si_ppsi:
            ax.axvline(1.5, c='black', ls=':')

        ax.axvspan(xmin=1.65, xmax=2.35, facecolor='grey', alpha=0.5, zorder=0)
        ax.axvspan(xmin=2.65, xmax=3.35, facecolor='tan', alpha=0.5, zorder=0)
        if not exclude_dge_si_ppsi:
            ax.axvline(3.5, c='black', ls=':')
            ax.axvspan(xmin=3.65, xmax=4.35, facecolor='grey', alpha=0.5, zorder=0)
            ax.axvspan(xmin=4.65, xmax=5.35, facecolor='tan', alpha=0.5, zorder=0)
        sns.scatterplot(y="scores", x="scmethod_nd_chem", hue='malignant_key',
                        data=mean_of_sid_means[
                            mean_of_sid_means.scoring_method == ax.title.get_text().split(' = ')[-1]], marker='_',
                        s=1000,
                        color='k', ax=ax, legend=False, palette=palette, zorder=2)

    g.set_axis_labels("", "Mean scores per sample and celltype", size=18)
    if exclude_dge_si_ppsi:
        g.set(xticks=([0.5, 2.5]))
        g.set_xticklabels(["Scoring samples individually", "Scoring on all samples"], size=16)
    else:
        g.set(xticks=([0.5, 2.5, 4.5]))
        g.set_xticklabels(["Scoring samples individually (preprocessed individually)",
                           "Scoring samples individually (preprocessed together)",
                           "Scoring on all samples (prepro-cessed together)"], size=14)
    for ax in g.axes[col_wrap:]:
        wrap_labels(ax, 11, break_long_words=False)

    g.set_yticklabels(np.round(g.axes[0].get_yticks(), decimals=2), size=16)
    g.set_titles("{col_name}", size=18, weight='normal')

    legend_elements = [Patch(facecolor='grey', edgecolor=None,
                             label='SC3Pv2'),
                       Patch(facecolor='tan', edgecolor=None,
                             label='SC3Pv3')]
    plt.legend(handles=legend_elements, title="Chemistry", fontsize=16, frameon=False, facecolor=None,
               title_fontsize=18, loc='upper right', bbox_to_anchor=(1.75, 1.1))
    return plt.gcf()


def create_strip_plot_other(dataset, adata, exclude_dge_si_ppsi=False, exclude_jasmine=False):
    plt.rcParams.update({'pdf.fonttype': 42, 'font.family': 'sans-serif', 'font.sans-serif': 'Arial', 'font.size': 14})
    means_std_per_sample = []
    for (key, adata_obs) in adata.obs.groupby('sample_id'):
        curr_means = adata_obs.groupby(by='malignant_key').mean()
        curr_means['sample_id'] = key
        curr_means['scoring_mode'] = 'Scoring on all samples' if exclude_dge_si_ppsi else 'Scoring on all samples (' \
                                                                                          'preprocessed together)'
        means_std_per_sample.append(curr_means)
    scmeans_per_mal = pd.concat(means_std_per_sample).reset_index()

    if dataset == 'escc':
        cols = ['n_counts', 'n_genes_by_counts', 'total_counts',
                'total_counts_mt', 'pct_counts_mt', 'S_score', 'G2M_score', 'AP_GT',
                'Cycling_GT', 'Epi1_GT', 'Epi2_GT', 'Mes_GT', 'Mucosal_GT', 'Oxd_GT',
                'Stress_GT', 'log_counts', ]
    elif dataset == 'luad':
        cols = ['age', 'n_genes_by_counts', 'total_counts',
                'total_counts_mito', 'pct_counts_mito', 'is_primary_data', 'n_genes',
                'total_counts_mt', 'pct_counts_mt', 'S_score', 'G2M_score',
                'log_counts', ]
    else:
        raise KeyError(f'dataset {dataset} unknown.')

    scmeans_per_mal = scmeans_per_mal.drop(columns=cols)
    scmeans_per_mal.dropna(inplace=True)
    scmeans_per_mal = scmeans_per_mal.melt(
        id_vars=['sample_id', 'malignant_key', 'scoring_mode'],
        var_name='scoring_method',
        value_name='scores'
    )
    if exclude_dge_si_ppsi:
        scmeans_per_mal = scmeans_per_mal[~scmeans_per_mal.scoring_method.str.contains('si_ppsi')].copy()
        scmeans_per_mal.loc[scmeans_per_mal.scoring_method.str.contains(
            'si_ppas'), 'scoring_mode'] = 'Scoring samples individually'
    else:
        scmeans_per_mal.loc[scmeans_per_mal.scoring_method.str.contains(
            'si_ppsi'), 'scoring_mode'] = 'Scoring samples individually (preprocessed individually)'
        scmeans_per_mal.loc[scmeans_per_mal.scoring_method.str.contains(
            'si_ppas'), 'scoring_mode'] = 'Scoring samples individually (preprocessed together)'

    most_means = scmeans_per_mal[scmeans_per_mal.scoring_method.str.contains('most_dge')].copy()
    # random_means = sample_means_per_chem_nd_mal[
    #     sample_means_per_chem_nd_mal.scoring_method.str.contains('random_dge')].copy()

    most_means.scoring_method = most_means.scoring_method.apply(lambda x: '_'.join(x.split('_')[0:-4]))
    # random_means.scoring_method = random_means.scoring_method.apply(lambda x: '_'.join(x.split('_')[0:-4]))

    most_means = most_means.sort_values(by='scoring_method')
    # random_means = random_means.sort_values(by='scoring_method')

    mean_of_sid_means = most_means.groupby(by=['scoring_method', 'scoring_mode', 'malignant_key']).mean().reset_index()

    sns.set_style("ticks")
    if exclude_dge_si_ppsi:
        order_x = ['Scoring samples individually',
                   'Scoring on all samples']
        height = 5
        aspect = 0.6
    else:
        order_x = ['Scoring samples individually (preprocessed individually)',
                   'Scoring samples individually (preprocessed together)',
                   'Scoring on all samples (preprocessed together)']
        height = 6
        aspect = 0.5

    if exclude_jasmine:
        col_order = ['ANS', 'Seurat', 'Seurat_AG', 'Seurat_LVG', 'Scanpy', 'UCell']
        col_wrap = 3
    else:
        col_order = ['ANS', 'Seurat', 'Seurat_AG', 'Seurat_LVG', 'Scanpy', 'Jasmine_LH', 'Jasmine_OR', 'UCell']
        col_wrap = 4

    g = sns.catplot(
        data=most_means,
        x='scoring_mode',
        order=order_x,
        y='scores',
        hue='malignant_key',
        col='scoring_method',
        col_order=col_order,
        col_wrap=col_wrap,
        kind='strip',
        dodge=False,
        height=height,
        aspect=aspect,
        legend=True,
    )
    g.add_legend(fontsize=16, loc='upper right', bbox_to_anchor=(1, 1))
    g.legend.set_title('Celltype', prop={'size': 18})

    palette = {'malignant': 'black', 'non-malignant': 'black'}
    for i, ax in enumerate(g.axes):
        sns.scatterplot(y="scores", x="scoring_mode", hue='malignant_key',
                        data=mean_of_sid_means[
                            mean_of_sid_means.scoring_method == ax.title.get_text().split(' = ')[-1]], marker='_',
                        s=1000,
                        color='k', ax=ax, legend=False, palette=palette, zorder=2)

    g.set_axis_labels("", "Mean scores per sample and celltype", size=18)
    if exclude_dge_si_ppsi:
        g.set(xticks=([0, 1]))
        g.set_xticklabels(["Scoring samples individually", "Scoring on all samples"], size=16)
    else:
        g.set(xticks=([0, 1, 2]))
        g.set_xticklabels(["Scoring samples individually (preprocessed individually)",
                           "Scoring samples individually (preprocessed together)",
                           "Scoring on all samples (prepro-cessed together)"], size=14)
    for ax in g.axes[col_wrap:]:
        wrap_labels(ax, 11, break_long_words=False)

    g.set_yticklabels(np.round(g.axes[0].get_yticks(), decimals=2), size=16)
    g.set_titles("{col_name}", size=18, weight='normal')

    return plt.gcf()


def create_strip_plot(dataset, adata, exclude_dge_si_ppsi=False, exclude_jasmine=False):
    if dataset == 'crc':
        return create_strip_plot_crc(adata, exclude_dge_si_ppsi, exclude_jasmine)
    elif dataset in ['escc', 'luad']:
        return create_strip_plot_other(dataset, adata, exclude_dge_si_ppsi, exclude_jasmine)
    else:
        raise KeyError(f'Cannot create strip plot for dataset {dataset}. ')


def main(config):
    # get storing path
    storing_path = get_storing_path(config.base_storing_path,
                                    dataset=config.dataset,
                                    norm_method=config.norm_method,
                                    dge_on_all=config.dge_on_all,
                                    min_log2fc=config.min_log2fc,
                                    pval=config.adj_pval)
    config['storing_path'] = storing_path
    sc.logging.info(f'-----------------------------------------------------------------')
    sc.logging.info(f'Got the storing path {storing_path}. Storing experiment configuration ...')

    sc.logging.info(f'-----------------------------------------------------------------')
    sc.logging.info(f'Get the signatures for malignant cells in decreasing and random log2FC order.')
    start = datetime.now()
    gene_list_most_de = get_malignant_signature(
        dataset=config.dataset,
        norm_method=config.norm_method,
        dge_on_all=config.dge_on_all,
        intersect_pct=config.intersect_pct,
        min_log2fc=config.min_log2fc,
        pval=config.adj_pval,
        ranked_means=config.ranked_means,
        sort_values_by=config.sort_log2FC_by,
        nr_sig_genes=config.sig_length
    )
    sc.logging.info(f'    > Got gene_list_most_de of length {len(gene_list_most_de)} '
                    f'(Duration {datetime.now() - start}).')

    sc.logging.info(f'-----------------------------------------------------------------')
    sc.logging.info(f'Score each sample individually on sample-based preprocessed data.')
    start = datetime.now()
    # scores_most_dge_si_ppsi, scores_random_dge_si_ppsi = score_ind_pp_ind(
    scores_most_dge_si_ppsi = score_independently(
        dataset=config.dataset,
        norm_method=config.norm_method,
        sample_based=True,
        sig_most_dge=gene_list_most_de,
        # sig_random_dge=gene_list_random_de

    )
    sc.logging.info(f'> Duration {datetime.now() - start}.')

    sc.logging.info(f'-----------------------------------------------------------------')
    sc.logging.info(f'Score each sample individually on data preprocessed over all samples.')
    start = datetime.now()
    # scores_most_dge_si_ppas, scores_random_dge_si_ppas = score_ind_pp_ind(
    scores_most_dge_si_ppas = score_independently(
        dataset=config.dataset,
        norm_method=config.norm_method,
        sample_based=False,
        sig_most_dge=gene_list_most_de,
        # sig_random_dge=gene_list_random_de
    )
    sc.logging.info(f'> Duration {datetime.now() - start}.')

    sc.logging.info(f'-----------------------------------------------------------------')
    sc.logging.info(f'Score all samples together.')
    start = datetime.now()
    adata = score_together(
        dataset=config.dataset,
        norm_method=config.norm_method,
        sample_based=False,
        sig_most_dge=gene_list_most_de,
        # sig_random_dge=gene_list_random_de
    )

    # Concatenate results
    adata.obs = pd.concat([adata.obs,
                           scores_most_dge_si_ppsi,
                           scores_most_dge_si_ppas,
                           # scores_random_dge_si_ppsi,
                           # scores_random_dge_si_ppas
                           ], axis=1)

    sc.logging.info(f'> Duration {datetime.now() - start}.')

    sc.logging.info(f'-----------------------------------------------------------------')
    fn = os.path.join(storing_path, f'{config.dataset}_adata_obs.csv')
    sc.logging.info(f'Storing current adata.obs. at {fn}')
    if config.save:
        adata.obs.to_csv(fn)

    sc.logging.info(f'-----------------------------------------------------------------')
    sc.logging.info(f'Get distribution information per malignancy and scoring_mode group.')
    start = datetime.now()
    stat_data = get_stats_over_modes(adata)
    if config.save:
        stat_data.to_csv(os.path.join(storing_path, f'distribution_info_scores.csv'))
    sc.logging.info(f'> Duration {datetime.now() - start}.')

    
    if config.dataset == 'crc':
        sc.logging.info(f'-----------------------------------------------------------------')
        sc.logging.info(f'Compare sample mean scores of malignant cells between the two chemistries.')
        start = datetime.now()
    
        pvals, score_means = get_pvals_chemistry(adata)
    
        if config.save:
            plt.rcParams.update(
                {'pdf.fonttype': 42, 'font.family': 'sans-serif', 'font.sans-serif': 'Arial', 'font.size': 16})
            pvals.to_csv(os.path.join(storing_path, f'mal_mean_scores_chemistry_pvals.csv'))
            score_means.to_csv(os.path.join(storing_path, f'mal_mean_and_var_scores.csv'))
            plot_pval_heatmaps(pvals, storing_path)
    
        sc.logging.info(f'> Duration {datetime.now() - start}.')
    
    sc.logging.info(f'-----------------------------------------------------------------')
    sc.logging.info(f'Create strip plot of per sample mean scores grouped by chemistry and malignancy.')
    start = datetime.now()
    ext_storing_path = os.path.join(storing_path, 'strip_plots')
    if config.save and not os.path.exists(ext_storing_path):
        os.makedirs(ext_storing_path)
    plt.rcParams.update({'pdf.fonttype': 42, 'font.family': 'sans-serif', 'font.sans-serif': 'Arial', 'font.size': 16})
    fig = create_strip_plot(config.dataset, adata, exclude_dge_si_ppsi=False, exclude_jasmine=False)
    if config.save:
        fig.savefig(os.path.join(ext_storing_path, f'strip_sample_mean_scores_wppsi_wjas.svg'), format='svg', bbox_inches='tight')
    
    fig = create_strip_plot(config.dataset, adata, exclude_dge_si_ppsi=True, exclude_jasmine=False)
    if config.save:
        fig.savefig(os.path.join(ext_storing_path, f'strip_sample_mean_scores_woppsi_wjas.svg'), format='svg', bbox_inches='tight')
    
    fig = create_strip_plot(config.dataset, adata, exclude_dge_si_ppsi=False, exclude_jasmine=True)
    if config.save:
        fig.savefig(os.path.join(ext_storing_path, f'strip_sample_mean_scores_wppsi_wojas.svg'), format='svg', bbox_inches='tight')
    
    fig = create_strip_plot(config.dataset, adata, exclude_dge_si_ppsi=True, exclude_jasmine=True)
    if config.save:
        fig.savefig(os.path.join(ext_storing_path, f'strip_sample_mean_scores_woppsi_wojas.svg'), format='svg', bbox_inches='tight')

    sc.logging.info(f'> Duration {datetime.now() - start}.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run dataset composition experiment.')

    # DATASET CONFIGURATION
    parser.add_argument('--dataset', choices=CANCER_DATASETS, default='crc', help='Indicate dataset on which '
                                                                                  'to create signature.')
    parser.add_argument('--norm_method', choices=NORM_METHODS,
                        default='mean', help='Indicate normalization method used in the preprocessing.')

    # EXPERIMENT CONFIGURATION
    parser.add_argument('--dge_on_all', default='pseudobulk', choices=DGEX_GRANULARITY,
                        help='Use signature for malignant cells computed on all samples. Default use signature '
                             'computed on each sample individually and take genes that appear in intersect_pct % of '
                             'the samples.')
    parser.add_argument('--ranked_means', action='store_true', help='Use dgex genes sorted by mean rank log2FC. Only '
                                                                    'used if dge_on_all=individual.')
    parser.add_argument('--intersect_pct', default=0.9, type=float, help="Defines the percentage of samples a gene "
                                                                         "needs to appear in the per sample DGEX to "
                                                                         "be considered as overall DGEX. Only "
                                                                         "important for dge_on_all=individual with "
                                                                         "ranked_means=False.")
    parser.add_argument('--min_log2fc', default=2, type=float, help='Log2FC threshold for DGEX gene selection.')
    parser.add_argument('--adj_pval', default=0.01, type=float, help='Adjusted p-value threshold for DGEX gene '
                                                                     'selection.')
    parser.add_argument('--sort_log2FC_by', default='median_log2FC', choices=['mean_log2FC', 'median_log2FC'],
                        help="Defines the order in which to sort the DGEX genes found in several samples, "
                             "when dge_on_all=individual.")

    parser.add_argument('--sig_length', default=100, type=int, help="Define the length of the signature use for the "
                                                                    "experiment.")

    # STORING CONFIGURATION
    parser.add_argument('--save', action='store_true', help='Figures and data only stored at storing path it flag is '
                                                            'set.')
    parser.add_argument('--base_storing_path',
                        default=os.path.join(BASE_PATH_EXPERIMENTS,'data_composition_experiments'), help='Path where to store the results.')

    args = AttributeDict(vars(parser.parse_args()))

    start = datetime.now()
    sc.logging.info(f'Data-composition experiment with the following configuration:\n{json.dumps(args, indent=4)}')

    main(args)

    sc.logging.info(f'Finished experiment in total {datetime.now() - start} time.')
