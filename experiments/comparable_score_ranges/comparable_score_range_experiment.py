import argparse
import json
import os
import sys
from collections import defaultdict
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, balanced_accuracy_score, f1_score, jaccard_score
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

import signaturescoring as ssc

sys.path.append('../../data')

from constants import (DATASETS_WITH_ANNOTATIONS, BASE_PATH_RESULTS, SCORING_METHODS,
                       METHOD_WITH_GENE_POOL, VIOLIN_PLOT_CONFIG)
from load_data import load_datasets, load_signatures

from helper_methods import get_violin_all_methods, prepare_data_for_violin_plot, plot_confusion_matrix


def get_storing_path(base_storing_path, dataset, remove_overlapping_genes, use_gene_pool, verbose=False):
    storing_path = Path(base_storing_path) / dataset
    if use_gene_pool:
        storing_path = storing_path / 'with_gene_pool'
    else:
        storing_path = storing_path / 'without_gene_pool'

    if remove_overlapping_genes:
        storing_path = storing_path / 'without_overlapping_genes'
    else:
        storing_path = storing_path / 'with_overlapping_genes'

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

        method_metrics.update({
            'conf_mat': conf_mat,
            'balanced_accuracy': bal_acc,
            'f1_score': f1_val,
            'jaccard_score': jacc_score
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
        args.use_gene_pool,
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
    if args.dataset == 'breast_malignant':
        adata.obs[args.gt_annotation_col] = adata.obs[args.gt_annotation_col].map(
            {str(i): f'GM{i}' for i in range(1, 8)})
    elif args.dataset.startswith('skin_malignant'):
        y_true_col = args.gt_annotation_col
        adata.obs[y_true_col] = adata.obs[y_true_col].astype(str)
        adata.obs.loc[adata.obs[y_true_col] == 'Tumor_KC_Cyc', y_true_col] = 'Tumor KC Cycling'
        adata.obs[y_true_col] = adata.obs[y_true_col].astype('category')
        adata.obs[y_true_col] = adata.obs[y_true_col].map(
            {val: val.replace('_', ' ') for val in adata.obs[y_true_col].unique()})

    # Load signatures
    signatures = load_signatures(args.dataset, args.remove_overlapping_genes)
    signature_order = sorted(list(signatures.keys()))
    if args.verbose:
        for cell_type, genes in signatures.items():
            print(cell_type, ' with signature length: ', len(genes))
    if args.save_signatures:
        try:
            fn = storing_path / f'{args.dataset}_signautres.json'
            with open(fn, 'w') as f:
                json.dump(signatures, f, indent=4)
            print(f"Dictionary successfully saved to {fn}")
        except Exception as e:
            print(f"Error saving dictionary: {e}")

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
    ## Violin plots
    df_melted = prepare_data_for_violin_plot(adata, args.gt_annotation_col, score_cols)
    fig = get_violin_all_methods(df_melted,
                                 args.gt_annotation_col,
                                 hue_order=signature_order,
                                 **asdict(VIOLIN_PLOT_CONFIG[args.dataset])
                                 )
    fig.savefig(storing_path / "violin_all_methods.pdf", bbox_inches='tight')
    fig.savefig(storing_path / "violin_all_methods.svg", bbox_inches='tight')
    if args.verbose:
        print(f"""Saved violin plot to {storing_path / "violin_all_methods.pdf"}""")

    ## Confusion matrices
    for key, val in metrics.items():
        conf_mat = val['conf_mat']
        fig = plot_confusion_matrix(conf_mat, signature_order, key)
        fig.savefig(storing_path / f"confusion_matrix_{key}.pdf", bbox_inches='tight')
        fig.savefig(storing_path / f"confusion_matrix_{key}.svg", bbox_inches='tight')
        if args.verbose:
            print(f"""Saved confusion matrix to {storing_path / f"confusion_matrix_{key}.pdf"}""")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run comparable score range experiment.")

    parser.add_argument("--dataset", type=str, choices=DATASETS_WITH_ANNOTATIONS)
    parser.add_argument("--sample_col", type=str, help="Patient/sample column.")
    parser.add_argument("--gt_annotation_col", type=str, help="Ground truth annotation column.")
    parser.add_argument("--remove_overlapping_genes", action="store_true")
    parser.add_argument("--use_gene_pool", action="store_true",
                        help="Use gene pool for methods that support")
    parser.add_argument("--nfolds", type=int, default=10, help="Number of folds for cross-validation.")
    parser.add_argument("--save_signatures", action="store_true",
                        help="Whether to store the signatures used for scoring.")
    parser.add_argument("--save_adata", action="store_true",
                        help="Figures and data only stored at storing path it flag is set.")
    parser.add_argument("--base_storing_path", default=os.path.join(BASE_PATH_RESULTS, "comparable_score_ranges"),
                        help="Path where to store the results.")
    parser.add_argument("--verbose", action="store_true", help="Print additional information.")

    args = parser.parse_args()

    print(args)

    main(args)
