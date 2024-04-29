import sys
import warnings

import scanpy as sc

sys.path.append('..')
sys.path.append('../..')
from ..data.load_data import load_dgex_genes_for_mal_cells


def get_scoring_method_params(sc_method_name):
    scoring_methods = {
        "adjusted_neighborhood_scoring": ("adjusted_neighborhood_scoring", {"ctrl_size": 100, "score_name": "ANS"}),
        "seurat_scoring": ("seurat_scoring", {"ctrl_size": 100, "n_bins": 25, "score_name": "Seurat"}),
        "seurat_ag_scoring": ("seurat_ag_scoring", {"n_bins": 25, "score_name": "Seurat_AG", }),
        "seurat_lvg_scoring": ("seurat_lvg_scoring", {"ctrl_size": 100, "n_bins": 25,
                                                      "lvg_computation_version": "v1",
                                                      "lvg_computation_method": "seurat",
                                                      "score_name": "Seurat_LVG", }),
        "scanpy_scoring": ("scanpy_scoring", {"ctrl_size": 100, "n_bins": 25, "score_name": "Scanpy", }),
        "jasmine_scoring_lh": ("jasmine_scoring", {"score_method": 'likelihood', "score_name": "Jasmine_LH", }),
        "jasmine_scoring_or": ("jasmine_scoring", {"score_method": 'oddsratio', "score_name": "Jasmine_OR", }),
        "ucell_scoring": ("ucell_scoring", {"score_name": "UCell", "maxRank": 1500, }),
    }
    if sc_method_name=='all':
        return scoring_methods
    elif isinstance(sc_method_name, list) and all([x in scoring_methods.keys() for x in sc_method_name]):
        return dict((key, scoring_methods[key]) for key in sc_method_name)
    elif sc_method_name not in scoring_methods.keys():
        raise KeyError(f'Unknown scoring method {sc_method_name}. Allowed methods are {list(scoring_methods.keys())}')
    return scoring_methods[sc_method_name]


def get_malignant_signature(dataset, norm_method, sample_based, dge_on_all, intersect_pct, min_log2fc, pval,
                            ranked_means, sort_values_by, sig_length=None, most_dge=True):
    wc = load_dgex_genes_for_mal_cells(dataset,
                                       norm_method=norm_method,
                                       sample_based=sample_based,
                                       dge_on_all=dge_on_all,
                                       intersect_pct=intersect_pct,
                                       min_log2fc=min_log2fc,
                                       pval=pval,
                                       ranked_means=ranked_means)
    sc.logging.info(f'> Got {len(wc)} available DGEX genes for malignant cells.')
    if (sig_length is not None) and (sig_length > len(wc)):
        warnings.warn(f'Passed sig_length {sig_length} is longer than the number of available DGEX genes. Running '
                      f'experiment with max sig_length len(gene_list).')
        sig_length = len(wc)

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

    if (sig_length is None) or (not isinstance(sig_length, int)):
        return wc[col].tolist()
    elif most_dge:
        return wc[0:sig_length][col].tolist()
    else:
        return wc.sample(sig_length)[col].tolist()


class AttributeDict(dict):
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
