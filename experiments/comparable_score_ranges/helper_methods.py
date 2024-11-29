import textwrap

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def plot_confusion_matrix(
        conf_mat,
        label_names,
        method_name,
        base_size=0.9,
        textwrap_width=7,
        xrotation=0,
        cbar=False,
        vmin=None,
        vmax=None,
        fontsizes={'title': 12, 'labels': 11, 'ticks': 10, 'legend': 11}

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
