import textwrap

import matplotlib.patches as mpatches
import numpy as np
import seaborn as sns


def get_violin_all_methods(
        df_melted,
        y_true_col,
        hue_order,
        height=1.95,
        aspect=0.925,
        textwrap_width=8,
        sharey=False,
        wspace=0.05,
        legend_bbox_anchor=(1.15, 1),
        fontsizes=dict(title=12, labels=11, ticks=9, legend=11),
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
        col_wrap=4,
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
