import argparse
import itertools
import json
import multiprocessing
import os
import sys
from collections import Counter
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
from signaturescoring import score_signature
from signaturescoring.scoring_methods.gmm_postprocessing import GMMPostprocessor
from signaturescoring.utils.utils import get_mean_and_variance_gene_expression, check_signature_genes
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import (
    balanced_accuracy_score,
    f1_score,
    jaccard_score,
    make_scorer,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

sys.path.append("..")
from experiment_utils import AttributeDict
sys.path.append("../..")
from data.constants import NORM_METHODS, BASE_PATH_DATA, BASE_PATH_EXPERIMENTS
from data.load_data import load_datasets

sc.settings.verbosity = 2

## GLOBAL VARIABLES
SC_NAMES = [
    "ANS",
    "Seurat",
    "Seurat_AG",
    "Seurat_LVG",
    "Scanpy",
    "Jasmine_LH",
    "Jasmine_OR",
    "UCell",
]

PERFORMANCE_COLS = [
    "Scoring method",
    "Information quantity (cross-validated F1-score)",
    "Hard labeling on scores (F1-score weighted)",
    "Hard labeling on scaled scores (F1-score weighted)",
    "Rediscovery score (F1-score weighted for unsupervised clustering)",
    "Information quantity (cross-validated Balanced accuracy)",
    "Hard labeling on scores (Balanced accuracy)",
    "Hard labeling on scaled scores (Balanced accuracy)",
    "Rediscovery score (Balanced accuracy for unsupervised clustering)",
    "Information quantity (cross-validated Jaccard-score)",
    "Hard labeling on scores (Jaccard-score weighted)",
    "Hard labeling on scaled scores (Jaccard-score weighted)",
    "Rediscovery score (Jaccard-score weighted for unsupervised clustering)",
]


## PERFORMANCE ESTIMATION METHODS
def _estimate_with_log_reg(X, y):
    # Load your data into X and y
    # Split your data into training and test sets if necessary

    # Define the number of outer and inner folds
    outer_cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
    inner_cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

    # Define the range of C values to test
    Cs = [0.001, 0.01, 0.1, 1, 10]

    # Initialize lists to store results
    f1_weighted_scores = []
    balanced_accuracy_scores = []
    jaccard_weighted_scores = []

    # Outer cross-validation loop
    for train_index, test_index in outer_cv.split(X, y):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

        # Initialize the StandardScaler
        scaler = StandardScaler()

        # Fit and transform the scaler on the training data
        X_train_scaled = scaler.fit_transform(X_train)

        # Create a logistic regression model with cross-validation
        logreg_cv = LogisticRegressionCV(
            Cs=Cs,
            cv=inner_cv,
            scoring=make_scorer(f1_score, average="weighted"),  # 'f1_weighted' scoring
            random_state=42,
            max_iter=10000,
            class_weight="balanced",
            refit=True,
        )

        # Fit the model on the scaled training data
        logreg_cv.fit(X_train_scaled, y_train)

        # Use the best model to predict on the test data
        X_test_scaled = scaler.transform(X_test)
        y_pred = logreg_cv.predict(X_test_scaled)

        # Calculate scores on the test data
        # f1_weighted
        f1_weighted_scores.append(f1_score(y_test, y_pred, average="weighted"))
        # balanced accuracy
        balanced_accuracy_scores.append(balanced_accuracy_score(y_test, y_pred))
        # jaccard weighted
        jaccard_weighted_scores.append(jaccard_score(y_test, y_pred, average="weighted"))

    # Calculate the average F1-weighted score across folds
    return (
        np.mean(f1_weighted_scores),
        np.mean(balanced_accuracy_scores),
        np.mean(jaccard_weighted_scores),
    )


def _get_hard_labeling_perf(df, y):
    y_pred = df.idxmax(axis=1)
    print(y_pred.unique())
    y_pred = (y_pred.apply(lambda x: x.rsplit("_", 1)[1])).values
    print(np.unique(y_pred))
    return (
        f1_score(y, y_pred, average="weighted"),
        balanced_accuracy_score(y, y_pred),
        jaccard_score(y, y_pred, average="weighted"),
    )


def _get_params_to_test(n_params, min_sc, max_sc, step=0.025):
    a = np.arange(min_sc, max_sc, step)
    params_list = [params for params in _generate_params(a, n_params)]
    return params_list


def _generate_params(a, n):
    if n == 1:
        return [(x,) for x in a]
    else:
        return [params + (x,) for x in a for params in _generate_params(a, n - 1)]


def _get_scaled_hard_labeling_one_item(df, y, scalars):
    df = df.copy()
    for i, sc_val in enumerate(scalars):
        df.iloc[:, i] = df.iloc[:, i] + sc_val

    y_pred = df.idxmax(axis=1)
    y_pred = (y_pred.apply(lambda x: x.rsplit("_", 1)[1])).values
    return (
        f1_score(y, y_pred, average="weighted"),
        balanced_accuracy_score(y, y_pred),
        jaccard_score(y, y_pred, average="weighted"),
    )


def _get_scaled_hard_lableling_perf(df, y, stride=0.025):
    min_score = df.values.min()
    max_score = df.values.max()
    n_cols = df.shape[1]
    params_list = _get_params_to_test(n_cols, min_score, max_score, stride)

    args = [(df, y, params) for params in params_list]

    with multiprocessing.Pool(processes=6) as pool:
        results = pool.starmap(_get_scaled_hard_labeling_one_item, args)

    results = np.array(results)
    return tuple(np.max(results, axis=0))


def _get_gmm_perf(adata, curr_score_names, y, K=3):
    sc.settings.verbosity = 0
    gmm_post = GMMPostprocessor(n_components=K)

    store_name_pred, store_names_proba, _ = gmm_post.fit_and_predict(
        adata, curr_score_names
    )
    assignments = gmm_post.assign_clusters_to_signatures(
        adata, curr_score_names, store_names_proba, plot=False
    )
    gmm_cols = []
    for key, val in assignments.items():
        if key == "rest":
            sc.logging.warning(
                f"Couldn't could find correct GMM cluster and signature mapping. Skipping {curr_score_names} is adviced"
            )
            continue
        adata.obs[key + "_gmm"] = adata.obs[val].copy()
        gmm_cols.append(key + "_gmm")

    y_pred = adata.obs[gmm_cols].idxmax(axis=1)
    y_pred = (y_pred.apply(lambda x: x.rsplit("_", 2)[1])).values
    sc.settings.verbosity = 2

    return (
        f1_score(y, y_pred, average="weighted"),
        balanced_accuracy_score(y, y_pred),
        jaccard_score(y, y_pred, average="weighted"),
    )


## DATASET LOADING METHODS
def _only_desired_keys(SG_subtypes, desired_keys):
    return {k: v for k, v in SG_subtypes.items() if k in desired_keys}


def _get_celltype_sigs_b_mono_nk(DE_of_celltypes, adata, desired_keys):
    celltype_sets_l2 = {}
    for group in adata.obs[["celltype.l2", "celltype.l3"]].groupby(by="celltype.l2"):
        celltype_sets_l2[group[0]] = list(np.unique(group[1]["celltype.l3"]))
    celltype_sets_l1 = {}
    for group in adata.obs[["celltype.l1", "celltype.l2"]].groupby(by="celltype.l1"):
        new_celltype_list = []
        for celltype in list(np.unique(group[1]["celltype.l2"])):
            new_celltype_list.append(celltype_sets_l2[celltype])
        celltype_sets_l1[group[0]] = list(itertools.chain(*new_celltype_list))

    grouped_DE_of_celltypes = DE_of_celltypes.groupby(by="Cell Type")

    SG_subtypes = {}
    for key, subtypes in celltype_sets_l1.items():
        sig_genes = set()
        for subtype in subtypes:
            if subtype not in grouped_DE_of_celltypes.groups.keys():
                continue
            group = grouped_DE_of_celltypes.get_group(subtype)
            sig_genes.update(
                group.sort_values(by="Average Log Fold Change", ascending=False)[
                    "Gene"
                ].iloc[0:300]
            )
        SG_subtypes[key] = sig_genes
    return _only_desired_keys(SG_subtypes, desired_keys)


def _get_celltype_sigs_b_subtypes(DE_of_celltypes, desired_keys):
    subtypes_B = np.unique(
        DE_of_celltypes[DE_of_celltypes["Cell Type"].str.contains("B ")]["Cell Type"]
    )
    SG_subtypes_B = {}
    for subtype in subtypes_B:
        SG_subtypes_B[subtype] = list(
            DE_of_celltypes[DE_of_celltypes["Cell Type"] == subtype]["Gene"]
        )
    SG_subtypes_B["B intermediate"] = set(SG_subtypes_B["B intermediate kappa"]).union(
        set(SG_subtypes_B["B intermediate lambda"])
    )
    SG_subtypes_B["B memory"] = set(SG_subtypes_B["B memory kappa"]).union(
        set(SG_subtypes_B["B memory lambda"])
    )
    SG_subtypes_B["B naive"] = set(SG_subtypes_B["B naive kappa"]).union(
        set(SG_subtypes_B["B naive lambda"])
    )
    for subtype in subtypes_B:
        SG_subtypes_B.pop(subtype, None)

    return _only_desired_keys(SG_subtypes_B, desired_keys)


def _remove_overlapping_genes(SG_subtypes):
    # Step 1: Create a frequency dictionary
    frequency_dict = Counter(
        item for sublist in SG_subtypes.values() for item in sublist
    )

    # Step 2: Remove strings appearing in more than one list
    return {
        key: [item for item in value if frequency_dict[item] == 1]
        for key, value in SG_subtypes.items()
    }


def _remove_genes_w_too_high_expr(adata, SG_subtypes, ctrl_size=100):
    for key, val in SG_subtypes.items():
        SG_subtypes[key] = list(val)
        SG_subtypes[key]  = check_signature_genes(adata.var_names, val)
    df_mean_var = get_mean_and_variance_gene_expression(adata)
    gene_means = df_mean_var["mean"].copy()
    sorted_gene_means = gene_means.sort_values()
    not_allowed = set(
        sorted_gene_means.iloc[-(ctrl_size // 2) :].index.tolist()
        + sorted_gene_means.iloc[: (ctrl_size // 2)].index.tolist()
    )
    return {
        key: [x for x in gene_list if not (x in not_allowed)]
        for key, gene_list in SG_subtypes.items()
    }


def _get_celltype_signatures(config, adata):
    DE_of_celltypes = pd.read_csv(config.DE_of_celltypes_fn)

    if config.dataset == "pbmc_b_subtypes":
        desired_keys = (
            ["B naive", "B memory"]
            if config.nr_sigs == 2
            else ["B naive", "B memory", "B intermediate"]
        )
        SG_subtypes = _get_celltype_sigs_b_subtypes(DE_of_celltypes, desired_keys)

    else:
        desired_keys = ["B", "Mono"] if config.nr_sigs == 2 else ["B", "Mono", "NK"]
        SG_subtypes = _get_celltype_sigs_b_mono_nk(DE_of_celltypes, adata, desired_keys)

    if not config.overlapping_sigs:
        SG_subtypes = _remove_overlapping_genes(SG_subtypes)
    
    SG_subtypes = _remove_genes_w_too_high_expr(adata, SG_subtypes, config.ctrl_size)
    
    for k,v in SG_subtypes.items():
        print(f'Signature for subtype {k} contains {len(v)} genes.')
    
    return SG_subtypes


def _get_sc_method_params(gene_pool, ctrl_size, n_bins):
    return [
        {
            "scoring_method": "scanpy_scoring",
            "sc_params": {
                "ctrl_size": ctrl_size,
                "n_bins": n_bins,
                "score_name": "Scanpy",
            },
        },
        {
            "scoring_method": "seurat_scoring",
            "sc_params": {
                "ctrl_size": ctrl_size,
                "n_bins": n_bins,
                "score_name": "Seurat",
                "gene_pool": gene_pool,
            },
        },
        {
            "scoring_method": "adjusted_neighborhood_scoring",
            "sc_params": {
                "ctrl_size": ctrl_size,
                "score_name": "ANS",
                "gene_pool": gene_pool,
            },
        },
        {
            "scoring_method": "seurat_ag_scoring",
            "sc_params": {
                "n_bins": n_bins,
                "score_name": "Seurat_AG",
                "gene_pool": gene_pool,
            },
        },
        {
            "scoring_method": "seurat_lvg_scoring",
            "sc_params": {
                "ctrl_size": ctrl_size,
                "n_bins": n_bins,
                "lvg_computation_version": "v1",
                "lvg_computation_method": "seurat",
                "score_name": "Seurat_LVG",
                "gene_pool": gene_pool,
            },
        },
        {
            "scoring_method": "ucell_scoring",
            "sc_params": {
                "score_name": "UCell",
            },
        },
        {
            "scoring_method": "jasmine_scoring",
            "sc_params": {
                "score_method": "likelihood",
                "score_name": "Jasmine_LH",
            },
        },
        {
            "scoring_method": "jasmine_scoring",
            "sc_params": {
                "score_method": "oddsratio",
                "score_name": "Jasmine_OR",
            },
        },
    ]


## MAIN COMPUTATION CODE
def main(config):
    print('### load dataset')
    adata = load_datasets(config.dataset, norm_method=config.norm_method)
    if "log1p" not in adata.uns.keys():
        adata.uns["log1p"] = {"base": None}
    else:
        adata.uns["log1p"]["base"] = None

    print('### get signatures to score')
    SG_subtypes = _get_celltype_signatures(config, adata)

    print('### create gene pool')
    all_sig_genes = set()
    for key, val in SG_subtypes.items():
        all_sig_genes.update(val)
    gene_pool = list(set(adata.var_names).difference(all_sig_genes))

    print('### get the available scoring methods')
    scoring_methods = _get_sc_method_params(gene_pool, config.ctrl_size, config.n_bins)

    print('### score all scoring methods')
    sc.settings.verbosity = 0
    scoring_names = []
    for sc_method in scoring_methods:
        scoring_method = sc_method["scoring_method"]
        sc_params = sc_method["sc_params"]

        for k1, v1 in SG_subtypes.items():
            curr_sc_params = sc_params.copy()
            curr_sc_params["score_name"] = curr_sc_params["score_name"] + "_" + k1

            score_signature(
                method=scoring_method, adata=adata, gene_list=v1, **curr_sc_params
            )

            scoring_names.append(curr_sc_params["score_name"])
    sc.settings.verbosity = 2
    scoring_names = [
        x
        for x in adata.obs.columns
        if any([y in SC_NAMES or y == "Jasmine" for y in x.split("_")])
    ]
    scoring_names = sorted(scoring_names)
    
    print('### save adata if demanded ')
    if config.save_adata:
        adata.write_h5ad(
            os.path.join(
                config.base_storing_path,
                f"adata_{config.dataset}_{config.nr_sigs}_3_{'over' if config.overlapping_sigs else 'non_over'}.h5ad",
            )
        )

    # terminate computation if we are only scoring for two signatures
    if config.nr_sigs == 2:
        sc.logging.warning(f'Hard labeling for 2/3 not implemented yet. Thus terminate experiment now ...')
        return

    print('### evaluate all performances')
    performances = []
    for i in range(0, len(scoring_names), config.nr_sigs):
        sc_method_name = scoring_names[i].rsplit("_", 1)[0]
        sc.logging.info(f"Evaluate performances for {sc_method_name}")

        curr_sc_names = scoring_names[i : (i + config.nr_sigs)]
        
        X = adata.obs[curr_sc_names].values
        y = adata.obs["celltype.l2"].values if config.dataset == 'pbmc_b_subtypes' else adata.obs["celltype.l1"].values 
        print("> Estimate performance with (supervised) LogReg")
        f1_logreg, balacc_logreg, jac_logreg = _estimate_with_log_reg(X, y)
        print("> Get hard labeling performance")
        f1_hard_lbl, balacc_hard_lbl, jac_hard_lbl = _get_hard_labeling_perf(
            adata.obs[curr_sc_names], y
        )
        print("> Get scaled hard labeling performance")
        (
            f1_sc_hard_lbl,
            balacc_sc_hard_lbl,
            jac_sc_hard_lbl,
        ) = _get_scaled_hard_lableling_perf(adata.obs[curr_sc_names], y, stride=config.scale_stride)
        print("> Get performance with GMM probabilities")
        f1_gmm, balacc_gmm, jac_gmm = _get_gmm_perf(adata, curr_sc_names, y)
        performances.append(
            (
                sc_method_name,
                f1_logreg,
                f1_hard_lbl,
                f1_sc_hard_lbl,
                f1_gmm,
                balacc_logreg,
                balacc_hard_lbl,
                balacc_sc_hard_lbl,
                balacc_gmm,
                jac_logreg,
                jac_hard_lbl,
                jac_sc_hard_lbl,
                jac_gmm,
            )
        )

    df = pd.DataFrame(performances, columns=PERFORMANCE_COLS)

    df.to_csv(
        os.path.join(
            config.base_storing_path,
            f"performances_{config.dataset}_{config.nr_sigs}_3_{'over' if config.overlapping_sigs else 'non_over'}.csv",
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run dataset composition experiment.")

    # DATASET CONFIGURATION
    parser.add_argument(
        "--dataset",
        choices=["pbmc_b_subtypes", "pbmc_b_mono_nk"],
        default="pbmc_b_subtypes",
    )
    parser.add_argument("--norm_method", choices=NORM_METHODS, default="mean")

    parser.add_argument("--n_bins", type=int, default=25)

    parser.add_argument("--ctrl_size", type=int, default=100)

    parser.add_argument("--nr_sigs", choices=[2, 3], type=int, default=3)

    parser.add_argument("--overlapping_sigs", action="store_true")

    parser.add_argument(
        "--DE_of_celltypes_fn",
        default=os.path.join(BASE_PATH_DATA, "annotations/citeseq_pbmc/DE_by_celltype.csv"),
    )

    parser.add_argument("--scale_stride", type=float, default=0.025)

    # STORING CONFIGURATION
    parser.add_argument(
        "--save_adata",
        action="store_true",
        help="Figures and data only stored at storing path it flag is " "set.",
    )
    parser.add_argument(
        "--base_storing_path",
        default=os.path.join(BASE_PATH_EXPERIMENTS, "comparable_score_ranges"),
        help="Path where to store the results.",
    )

    args = AttributeDict(vars(parser.parse_args()))

    start = datetime.now()
    sc.logging.info(
        f"Comparable score ranges experiment with the following configuration:\n{json.dumps(args, indent=4)}"
    )

    main(args)
