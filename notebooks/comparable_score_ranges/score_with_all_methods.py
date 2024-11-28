import sys
import textwrap
from collections import defaultdict

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, confusion_matrix, f1_score
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

import signaturescoring as ssc

sys.path.append('/Users/lciernik/Documents/TUB/projects/ans_scoring/ANS_supplementary_information')
from data.constants import SCORING_METHODS


def score_signatures_with_all_methods(adata, sig_genes, verbose=False):
    added_cols = defaultdict(list)
    for sig_name, sig in sig_genes.items():
        for sc_method in SCORING_METHODS:
            scoring_method = sc_method['scoring_method']
            sc_params = sc_method['sc_params'].copy()

            if verbose:
                print(f"Scoring {sig_name} with {scoring_method}")

            col_name = sig_name.replace(' ', '_')
            prev_score_name = sc_params['score_name']
            sc_params['score_name'] = f"{col_name}_{sc_params['score_name']}_scores"

            ssc.score_signature(
                method=scoring_method,
                adata=adata,
                gene_list=sig,
                **sc_params
            )

            added_cols[prev_score_name].append(sc_params['score_name'])
    return added_cols, adata


def label_assignment_from_scores(adata, method_name, method_scores, include_undefined=False):
    lbl_col = f'{method_name}_label'
    max_score_col = f'{method_name}_max_score'

    adata.obs[lbl_col] = adata.obs.loc[:, method_scores].idxmax(axis=1)
    adata.obs[max_score_col] = adata.obs.loc[:, method_scores].max(axis=1)
    if include_undefined:
        adata.obs.loc[adata.obs[max_score_col] < 0, lbl_col] = 'Undefined'

    adata.obs[lbl_col] = adata.obs[lbl_col].apply(lambda x: x.split(f'_{method_name}_')[0])
    adata.obs[lbl_col] = adata.obs[lbl_col].apply(lambda x: x.replace('_', ' '))
    return adata, lbl_col


def get_lbl_assignment_performance(adata, y_true_col, y_pred_col, label_names, avg='weighted'):
    conf_mat = confusion_matrix(adata.obs[y_true_col],
                                adata.obs[y_pred_col],
                                normalize='true',
                                labels=label_names)

    bal_acc = balanced_accuracy_score(adata.obs[y_true_col], adata.obs[y_pred_col])
    f1_val = f1_score(adata.obs[y_true_col], adata.obs[y_pred_col], average=avg)

    return conf_mat, bal_acc, f1_val


def plot_confusion_matrix(
        conf_mat,
        label_names,
        method_name,
        # figsize=(6, 5),
        base_size=0.8,
        textwrap_width=8,
        xrotation=0,
        cbar=False,
        vmin=None,
        vmax=None,
        fontsizes={'title': 12, 'labels': 11, 'ticks': 11, 'legend': 11}

):
    figsize = (base_size * len(label_names), base_size * len(label_names))
    fig = plt.figure(figsize=figsize)

    g = sns.heatmap(
        conf_mat * 100,
        annot=True,
        fmt=".2f",
        cmap='coolwarm',
        annot_kws={"fontsize": fontsizes['labels']},
        cbar=cbar,
        vmin=vmin,
        vmax=vmax
    )

    new_labels = [textwrap.fill(label, width=textwrap_width) for label in label_names]

    g.set_title(f'{method_name}', fontsize=fontsizes['title'])
    g.set_ylabel('True', fontsize=fontsizes['labels'])
    g.set_xlabel('Predicted', fontsize=fontsizes['labels'])
    g.set_xticklabels(new_labels, rotation=xrotation, fontsize=fontsizes['ticks'])
    g.set_yticklabels(new_labels, rotation=0, fontsize=fontsizes['ticks'])
    g.tick_params(axis='x', labelsize=fontsizes['ticks'], width=0.85)
    g.tick_params(axis='y', labelsize=fontsizes['ticks'], width=0.85)
    return fig


def get_information_from_scores(adata, y_true_col, scores, nfold=10, metric='balanced_accuracy', max_iter=1000):
    X = adata.obs.loc[:, scores].values
    y = adata.obs[y_true_col].values

    model = Pipeline([
        ('scaler', StandardScaler()),
        ('logreg', LogisticRegression(max_iter=max_iter))
    ])

    cv = StratifiedKFold(n_splits=nfold)

    scores = cross_val_score(model, X, y, cv=cv, scoring=metric)
    return scores


def get_violin(df, scores_cols, y_true_col, figsize=(12, 6)):
    sns.set(style='ticks')
    melted_df = pd.melt(
        df,
        id_vars=y_true_col,
        value_vars=scores_cols,
        var_name='Signature',
        value_name='Score'
    )
    fig = plt.figure(figsize=figsize)
    sns.violinplot(
        data=melted_df,
        x=y_true_col,
        y='Score',
        hue='Signature',
    )
    return fig


def remove_overlapping_signature_genes(sig_genes):
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


def prepare_data_for_violin_plot(adata, y_true_col, score_cols):
    tmp = adata.obs.copy()
    tmp = tmp.reset_index(names=['old_index'])
    dfs = []
    for method_name, method_scores in score_cols.items():
        new_col_names = [x.split("_" + method_name + "_")[0].replace('_', ' ') for x in method_scores]
        df = tmp.loc[:, [y_true_col] + method_scores].copy()
        df.columns = [y_true_col] + new_col_names
        df['Scoring method'] = method_name
        dfs.append(df)
    df = pd.concat(dfs, ignore_index=True).reset_index(drop=True)
    df_melted = df.melt(id_vars=[y_true_col, 'Scoring method'], var_name='Signature', value_name='Scores')
    return df_melted


def get_violin_all_methods(
        df_melted,
        y_true_col,
        hue_order,
        height=1.95,
        aspect=0.925,
        textwrap_width=8,
        sharey=False,
        wspace=0.05,
        col_wrap=4,
        legend_bbox_anchor=(1.15, 1),
        fontsizes=dict(title=12, labels=11, ticks=11, legend=11),
):
    sns.set_style("ticks")
    g = sns.catplot(
        data=df_melted,
        x=y_true_col,
        order=hue_order,
        y="Scores",
        hue="Signature",
        hue_order=hue_order,
        kind="violin",
        col='Scoring method',
        col_wrap=col_wrap,
        height=height,
        aspect=aspect,
        density_norm='width',
        inner=None,
        linewidth=0.5,
        sharey=sharey
    )

    g.set_titles("{col_name}", size=fontsizes['title'])
    g.set_xlabels("")
    g.set_ylabels("Scores", size=fontsizes['labels'])
    for ax in g.axes.flat:
        ax.tick_params(axis='y', labelsize=fontsizes['ticks'], length=2.5, width=0.85, pad=-0.01)
        ax.tick_params(axis='x', labelsize=fontsizes['ticks'], length=2.5, width=0.85)

        list_xticks_lbls = [label.get_text() for label in ax.get_xticklabels()]
        if list_xticks_lbls:
            labels = [textwrap.fill(label, width=textwrap_width) for label in list_xticks_lbls]
            ax.set_xticks(ax.get_xticks(), labels=labels)

    #### SPAN and Separator ####
    unique_values = hue_order
    n_unique = len(unique_values)

    x_positions = np.arange(len(unique_values) + 1) - 0.5

    colors = sns.color_palette("tab10", n_unique)  # Get a color palette

    for k, ax in enumerate(g.axes.flat):
        # Color ground truth
        for i, color in enumerate(colors):
            ax.axvspan(x_positions[i], x_positions[i + 1],
                       color=color, alpha=0.3, ymin=0,
                       ymax=0.05)  # Adjust ymin and ymax to control the height of the span

        # Add vertical lines to separate the spans
        for pos in x_positions[1:-1]:
            ax.axvline(pos, color='grey', lw=0.85, ls=':', alpha=0.5, zorder=-1)

    #### LEGEND ####
    handles1, labels1 = g.legend.legendHandles, [text.get_text() for text in g.legend.get_texts()]

    colors = sns.color_palette("tab10", len(unique_values))
    handles2 = [mpatches.Patch(color=colors[i], alpha=0.3, label=unique_values[i]) for i in range(len(unique_values))]
    labels2 = unique_values

    title1 = mpatches.Patch(color='none', label="Signature Legend")
    title2 = mpatches.Patch(color='none', label="y_true_col Legend")

    merged_handles = [title1] + handles1 + [title2] + handles2
    merged_labels = ["Signature"] + labels1 + ["Ground truth cell type"] + list(labels2)

    g.legend.remove()
    g.fig.legend(handles=merged_handles, labels=merged_labels, frameon=False, bbox_to_anchor=legend_bbox_anchor,
                 borderaxespad=0.,
                 fontsize=fontsizes['legend'])

    # plt.tight_layout()

    #### Axis line width ####
    for ax in g.axes.flat:
        for spine in ax.spines.values():
            spine.set_linewidth(0.85)  # Set axis line width

    g.fig.subplots_adjust(hspace=0.2, wspace=wspace)

    return g.fig


def save_close_or_show(fig, save, save_path):
    if save:
        fig.savefig(save_path, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved figure at {save_path}.")
    else:
        plt.show(fig)
