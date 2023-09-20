import glob
import os
import sys

import pandas as pd
import pyreadr
import scanpy as sc
from signaturescoring import score_signature

sys.path.append('../..')
from data.load_data import load_datasets, load_dgex_genes_for_mal_cells
from data.constants import BASE_PATH_EXPERIMENTS

sc.settings.verbosity = 2


############################################################################################################
# HELPER METHODS
def load_adata_gene_list(dataset, sig_len=100):
    adata = load_datasets(dataset)
    wc = load_dgex_genes_for_mal_cells(dataset)
    wc = wc.sort_values(by=['padj', 'log2FoldChange'], ascending=[True, False])
    gene_list = wc[0:sig_len].genes.tolist()
    return adata, gene_list


def get_ranked_data(dataset):
    path = os.path.join(BASE_PATH_EXPERIMENTS, f'construction_scoring_methods/{dataset}_debug_data_ranked.rds')
    result = pyreadr.read_r(path)  # also works for RData
    df = result[None]
    df = df.set_index('rn')
    return df.T


def get_ucell_scores_R(dataset):
    path = os.path.join(BASE_PATH_EXPERIMENTS, f'construction_scoring_methods/{dataset}_debug_ucell_scores.csv')
    df = pd.read_csv(path)
    df = df.set_index('Unnamed: 0')
    df.index.name = None
    return df


def score_data(adata, gene_list):
    scoring_method = "ucell_scoring"
    sc_params = {"score_name": "UCell_Python",
                 "maxRank": 1500, }
    score_signature(method=scoring_method,
                    adata=adata,
                    gene_list=gene_list,
                    **sc_params)
    return adata


############################################################################################################
# "DEFINE DATASETS
# datasets = ['luad']
datasets = ['crc', 'escc', 'luad']

############################################################################################################
sc.logging.info("LOAD LARGE DATASETS")
adatas = {}
gene_lists = {}
for ds in datasets:
    adatas[ds], gene_lists[ds] = load_adata_gene_list(ds)

############################################################################################################
sc.logging.info("REDUCE DATASETS")
base_path = os.path.join(BASE_PATH_EXPERIMENTS, 'construction_scoring_methods/')

fns = glob.glob(os.path.join(base_path, '*_sample_cells.csv'))
for fn in fns:
    ds = (os.path.basename(fn)).split('_')[0]
    if ds not in adatas.keys():
        continue
    samples = pd.read_csv(fn).columns.tolist()
    adatas[ds] = adatas[ds][samples, :].copy()

for ds in datasets:
    sc.logging.info(f'{ds}, {adatas[ds].shape}')

############################################################################################################
sc.logging.info("GET DATASETS WITH RANKED GENES FROM R VERSION")
ranked_data = {}
for ds in datasets:
    ranked_df = get_ranked_data(ds)
    ranked_df.columns.name = None
    ranked_data[ds] = ranked_df

############################################################################################################
sc.logging.info("GET UCELL SCORES FROM R VERSION AND ADD TO THE ANNDATA.OBS DATAFRAMES")
for ds in datasets:
    adatas[ds].obs['UCell_R'] = get_ucell_scores_R(ds)

############################################################################################################
sc.logging.info("SCORE DATASETS WITH UCELL PYTHON VERSION")
for ds in datasets:
    adatas[ds] = score_data(adatas[ds], gene_lists[ds])

