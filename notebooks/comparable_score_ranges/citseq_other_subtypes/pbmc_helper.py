import pandas as pd


def load_dex_genes(filter_genes=False, threshold_pval=0.01, threshold_log2fc=0.5):
    dge_genes_fn = ('/Users/lciernik/Documents/TUB/projects/ans_scoring/data/annotations/citeseq_pbmc/DE_by_celltype.csv')
    dge_genes = pd.read_csv(dge_genes_fn)

    if filter_genes:
        print('Shape DEX genes BEFORE filtering', dge_genes.shape)
        dge_genes = dge_genes[
            (dge_genes['P Value'] <= threshold_pval) &
            (dge_genes['Average Log Fold Change'] >= threshold_log2fc)
            ].copy().reset_index(drop=True)
        print('Shape DEX genes AFTER filtering', dge_genes.shape)
    return dge_genes
