import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import pandas as pd
import time

sys.path.append('../..')
from data.constants import BASE_PATH_EXPERIMENTS

plt.rcParams.update(
    {
        "pdf.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": "Arial",
        "font.size": 10,
    }
)

# ALLOWED_CELL_SUBTYPES = ['B memory kappa', 'B naive kappa', 'B naive lambda', 'CD14 Mono', 'CD16 Mono', 'CD4 CTL', 'CD4 Naive', 'CD4 TCM_1', 'CD4 TCM_3', 'CD4 TEM_1', 'CD4 TEM_3', 'CD8 Naive', 'CD8 TEM_1', 'CD8 TEM_2', 'CD8 TEM_4', 'CD8 TEM_5', 'MAIT', 'NK_1', 'NK_2', 'NK_3', 'Platelet', 'cDC2_2']
ALLOWED_CELL_SUBTYPES = ['B memory kappa', 'CD14 Mono', 'CD8 TEM_2', 'NK_3']


def _plot_lines(melted_all_scores, bin_cuttoff, fifty_last_genes, hue_order, with_eb=False):
    cm = 1 / 2.54  # centimeters in inches
    plt.figure(figsize=(11 * cm, 6.5 * cm))
    ax = sns.lineplot(
        data=melted_all_scores,
        x="gene",
        y="score",
        style="scoring_method",
        style_order=hue_order,
        hue="scoring_method",
        hue_order=hue_order,
        errorbar="sd" if with_eb else None,
    )

    ax.set_title(f"Control genes selection bias for subtype {folder_path.name}", fontsize=10)
    plt.xlabel("Genes ordered according to average expression values (top 8%).", fontsize=10)
    plt.ylabel("Average score of a 1-gene signature", fontsize=10)

    plt.xticks([])
    plt.yticks(fontsize=8)

    plt.axvline(bin_cuttoff, ls=":", c="grey", label="Cutoff bins")
    plt.axvline(fifty_last_genes, ls=":", c="grey", label="Cutoff bins")
    #     all_xticks = ax.get_xticks()
    #     subset_x_ticks = all_xticks[::50]
    #     ax.set_xticks(subset_x_ticks)

    return plt.gcf()


def _create_and_store_plot(folder_path, storing_path):
    sc_method_name_mapping = {
        "adjusted_neighborhood_scoring": "ANS",
        "seurat_scoring": "Seurat",
        "seurat_ag_scoring": "Seurat_AG",
        "seurat_lvg_scoring": "Seurat_LVG",
        "scanpy_scoring": "Scanpy",
    }
    hue_order = list(sc_method_name_mapping.values())

    print("> read all .csv files containing cell scores")
    dfs = []
    for file in folder_path.glob("*.csv"):
        df = pd.read_csv(file)
        dfs.append(df.copy())
    all_scores = pd.concat(dfs, axis=0)
    all_scores.columns = ["sample_id"] + list(all_scores.columns)[1:]
    bin_cuttoff = (len(all_scores.columns) - 2) // 2
    fifty_last_genes = (len(all_scores.columns) - 2) - 50

    print("> melt all scores to make it ready for seaborn lineplot")
    melted_all_scores = pd.melt(
        all_scores,
        id_vars=["sample_id", "scoring_method"],
        var_name="gene",
        value_name="score",
    )

    melted_all_scores.scoring_method = melted_all_scores.scoring_method.map(
        sc_method_name_mapping
    )

    print("> plot score lines WITHOUT errorbounds")
    fig_wo_eb = _plot_lines(melted_all_scores, bin_cuttoff, fifty_last_genes, hue_order, with_eb=False)
    print("> plot score lines WITH errorbounds")
    fig_w_eb = _plot_lines(melted_all_scores, bin_cuttoff, fifty_last_genes, hue_order, with_eb=True)

    print("> store figure")
    fig_wo_eb.tight_layout()
    fig_wo_eb.savefig(
        storing_path / f"bias_{folder_path.name.replace(' ', '_')}_wo_eb.pdf",
        format="pdf",
    )
    fig_w_eb.tight_layout()
    fig_w_eb.savefig(
        storing_path / f"bias_{folder_path.name.replace(' ', '_')}_w_eb.pdf",
        format="pdf",
    )


if __name__ == "__main__":
    root_exp_dir = Path(
        os.path.join(BASE_PATH_EXPERIMENTS,"control_genes_selection/mean_var_per_gene_scores")
        )
    storing_path = root_exp_dir / "plots"

    # Create the directory if it doesn't exist
    storing_path.mkdir(parents=True, exist_ok=True)

    for folder_path in root_exp_dir.rglob("*"):
        if folder_path.is_dir() and folder_path.name in ALLOWED_CELL_SUBTYPES:
            if "plots" in folder_path.name:
                continue
            print(f"Processing subtype {folder_path}")
            start_time = time.time()
            _create_and_store_plot(folder_path, storing_path)
            elapsed_time = time.time() - start_time
            print(f"> Time elapsed for {folder_path.name}: {elapsed_time} seconds")

    print("Finished figure creation for control selection bias.")
