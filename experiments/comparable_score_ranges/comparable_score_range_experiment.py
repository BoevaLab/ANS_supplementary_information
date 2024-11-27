import argparse
import os
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, balanced_accuracy_score, f1_score, jaccard_score
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

import signaturescoring as ssc
from signaturescoring.scoring_methods.gmm_postprocessing import GMMPostprocessor

sys.path.append('../..')

from ANS_supplementary_information.data.constants import (DATASETS_WITH_ANNOTATIONS, BASE_PATH_RESULTS, SCORING_METHODS,
                                                          METHOD_WITH_GENE_POOL)
from ANS_supplementary_information.data.load_data import load_datasets, load_signatures

from helper_methods import get_violin_all_methods


def get_storing_path(base_storing_path, dataset, remove_overlapping_genes, verbose=False):
    storing_path = Path(base_storing_path) / dataset
    storing_path = storing_path / (
        'without_overlapping_genes' if remove_overlapping_genes else 'with_overlapping_genes')

    if not storing_path.exists():
        storing_path.mkdir(parents=True)
        if verbose:
            print(f"Created storing path {storing_path}")
    return storing_path


def _get_gene_pool(adata, signatures):
    all_sig_genes = set()
    for key, val in signatures.items():
        all_sig_genes.update(val)
    gene_pool = list(set(adata.var_names).difference(all_sig_genes))
    return gene_pool


def score_signatures_with_all_methods(adata, signatures, use_gene_pool: bool, verbose=False):
    if use_gene_pool:
        gene_pool = _get_gene_pool(adata, signatures)
    else:
        gene_pool = None

    added_cols = defaultdict(list)
    for sig_name, sig in signatures.items():
        for sc_method in SCORING_METHODS:
            scoring_method = sc_method['scoring_method']
            sc_params = sc_method['sc_params'].copy()

            if verbose:
                print(f"Scoring {sig_name} with {scoring_method}")

            col_name = sig_name.replace(' ', '_')
            prev_score_name = sc_params['score_name']
            sc_params['score_name'] = f"{col_name}_{sc_params['score_name']}_scores"
            if use_gene_pool and scoring_method in METHOD_WITH_GENE_POOL:
                sc_params['gene_pool'] = gene_pool

            ssc.score_signature(
                method=scoring_method,
                adata=adata,
                gene_list=sig,
                **sc_params
            )

            added_cols[prev_score_name].append(sc_params['score_name'])
    return added_cols, adata


def lbls_from_method_scores(adata, method_name, method_scores, include_undefined=False):
    lbl_col = f'{method_name}_label'
    max_score_col = f'{method_name}_max_score'

    adata.obs[lbl_col] = adata.obs.loc[:, method_scores].idxmax(axis=1)
    adata.obs[max_score_col] = adata.obs.loc[:, method_scores].max(axis=1)
    if include_undefined:
        adata.obs.loc[adata.obs[max_score_col] < 0, lbl_col] = 'Undefined'

    adata.obs[lbl_col] = adata.obs[lbl_col].apply(lambda x: x.split(f'_{method_name}_')[0])
    adata.obs[lbl_col] = adata.obs[lbl_col].apply(lambda x: x.replace('_', ' '))
    return adata, lbl_col


def assign_labels_from_scores(adata, score_cols, include_undefined=False):
    all_cols = []

    label_cols = {}
    for method_name, method_scores in score_cols.items():
        adata, new_lbl_col = lbls_from_method_scores(adata, method_name, method_scores,
                                                     include_undefined=include_undefined)
        label_cols[method_name] = new_lbl_col
        all_cols += method_scores + [new_lbl_col]

    return all_cols, label_cols, adata


def get_lbl_assignment_performance(adata, y_true_col, y_pred_col, label_names, avg='weighted'):
    conf_mat = confusion_matrix(adata.obs[y_true_col],
                                adata.obs[y_pred_col],
                                normalize='true',
                                labels=label_names)

    bal_acc = balanced_accuracy_score(adata.obs[y_true_col], adata.obs[y_pred_col])
    f1_val = f1_score(adata.obs[y_true_col], adata.obs[y_pred_col], average=avg)
    jacc_score = jaccard_score(adata.obs[y_true_col], adata.obs[y_pred_col], average=avg)
    return conf_mat, bal_acc, f1_val, jacc_score


def get_gmm_perf(adata, y_true_col, score_cols, K=3, avg='weighted'):
    gmm_post = GMMPostprocessor(n_components=K)

    store_name_pred, store_names_proba, _ = gmm_post.fit_and_predict(adata, score_cols)

    assignments = gmm_post.assign_clusters_to_signatures(adata, score_cols, store_names_proba, plot=False)

    gmm_cols = []
    for key, val in assignments.items():
        if key == "rest":
            warnings.warn(
                f"Couldn't could find correct GMM cluster and signature mapping. Skipping {score_cols} is advised"
            )
            continue
        adata.obs[key + "_gmm"] = adata.obs[val].copy()
        gmm_cols.append(key + "_gmm")

    y = adata.obs[y_true_col].values
    y_pred = adata.obs[gmm_cols].idxmax(axis=1)
    y_pred = (y_pred.apply(lambda x: x.rsplit("_", 2)[1])).values

    return (
        f1_score(y, y_pred, average=avg),
        balanced_accuracy_score(y, y_pred),
        jaccard_score(y, y_pred, average=avg),
    )


def get_information_from_scores(adata, y_true_col, score_cols, nfold=10, max_iter=1000):
    X = adata.obs.loc[:, score_cols].values
    y = adata.obs[y_true_col].values

    model = Pipeline([
        ('scaler', StandardScaler()),
        ('logreg', LogisticRegression(max_iter=max_iter, penalty=None))
    ])

    cv = StratifiedKFold(n_splits=nfold)

    perfs = {}
    for met in ['balanced_accuracy', 'f1_weighted', 'jaccard_weighted']:
        scores = cross_val_score(model, X, y, cv=cv, scoring=met)
        perfs[met] = scores
    agg_perfs = {k: (np.mean(v), np.std(v)) for k, v in perfs.items()}
    return agg_perfs


def get_all_performances(score_cols, label_cols, adata, y_true_col, signature_order, nfold=10):
    metrics = defaultdict(dict)
    for method_name, method_scores in score_cols.items():
        method_metrics = {}

        lbl_col = label_cols[method_name]
        # Hard label assignment
        conf_mat, bal_acc, f1_val, jacc_score = get_lbl_assignment_performance(adata,
                                                                               y_true_col=y_true_col,
                                                                               y_pred_col=lbl_col,
                                                                               label_names=signature_order)

        # GMM
        f1_gmm, bal_acc_gmm, jacc_score_gmm = get_gmm_perf(adata, y_true_col, method_scores, K=len(method_scores))

        method_metrics.update({
            'conf_mat': conf_mat,
            'balanced_accuracy': bal_acc,
            'f1_score': f1_val,
            'jaccard_score': jacc_score,
            'gmm_balanced_accuracy': bal_acc_gmm,
            'gmm_f1_score': f1_gmm,
            'gmm_jaccard_score': jacc_score_gmm
        })

        # Logreg
        logreg_metrics = get_information_from_scores(adata, y_true_col, method_scores, nfold=nfold)

        method_metrics.update({f'logreg_{key}_{nfold}cv_mean': val[0] for key, val in logreg_metrics.items()})
        method_metrics.update({f'logreg_{key}_{nfold}cv_std': val[1] for key, val in logreg_metrics.items()})

        metrics[method_name] = method_metrics

    return metrics


def main(args):
    # prepare storing path
    storing_path = get_storing_path(
        args.base_storing_path,
        args.dataset,
        args.remove_overlapping_genes,
        args.verbose
    )

    # Load data
    adata = load_datasets(args.dataset, preprocessed=True, norm_method='mean', sample_based=False)
    if args.verbose:
        print(adata)

    # Test if ground truth annotation column is available
    if args.gt_annotation_col not in adata.obs.columns:
        raise KeyError(
            f'Ground truth annotation column {args.gt_annotation_col} not in adata columns {adata.obs.columns}')

    # Load signatures
    signatures = load_signatures(args.dataset, args.remove_overlapping_genes)
    signature_order = sorted(list(signatures.keys()))
    if args.verbose:
        for cell_type, genes in signatures.items():
            print(cell_type, ' with signature length: ', len(genes))

    # Score signatures with all methods
    score_cols, adata = score_signatures_with_all_methods(adata, signatures, args.use_gene_pool, verbose=args.verbose)

    # Label assignment from scores
    all_cols, label_cols, adata = assign_labels_from_scores(adata, score_cols, include_undefined=False)

    # Get performance: hard label assignment, logreg, and label assignment with GMMs
    metrics = get_all_performances(score_cols, label_cols, adata, args.gt_annotation_col, signature_order, args.nfolds)
    metrics_df = pd.DataFrame(metrics)

    # Save data
    metrics_df.to_csv(storing_path / 'metrics.csv')
    if args.verbose:
        print(f"""Saved metrics to {storing_path / 'metrics.csv'}""")

    if args.save_adata:
        adata.write(storing_path / 'adata.h5ad')
        if args.verbose:
            print(f"""Saved adata to {storing_path / 'adata.h5ad'}""")


    # Create plots
    ## Umap plots
    ## Violin plots
    ## Confusion matrices


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run comparable score range experiment.")

    parser.add_argument("--dataset", type=str, choices=DATASETS_WITH_ANNOTATIONS)
    parser.add_argument("--gt_annotation_col", type=str, help="Ground truth annotation column.")
    parser.add_argument("--remove_overlapping_genes", action="store_true")
    parser.add_argument("--use_gene_pool", action="store_true",
                        help="Use gene pool for methods that support")
    parser.add_argument("--nfolds", type=int, default=10, help="Number of folds for cross-validation.")
    parser.add_argument("--save_adata", action="store_true",
                        help="Figures and data only stored at storing path it flag is set.")
    parser.add_argument("--base_storing_path", default=os.path.join(BASE_PATH_RESULTS, "comparable_score_ranges"),
                        help="Path where to store the results.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing files.")
    parser.add_argument("--verbose", action="store_true", help="Print additional information.")

    args = parser.parse_args()

    main(args)
