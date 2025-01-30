## DEFINED CONSTANTS
####################################################################################
### PATHS
## !!!!!! TODO: CHANGE path to folder where data and experiments should go !!!!!!!
from typing import Dict, List, Tuple

BASE_PATH_DRIVE = '/Users/lciernik/Documents/TUB/projects/ans_scoring'

BASE_PATH_EXPERIMENTS = BASE_PATH_DRIVE + '/experiments'
BASE_PATH_RESULTS = BASE_PATH_DRIVE + '/results_2'
BASE_PATH_DATA = BASE_PATH_DRIVE + '/data'
BASE_PATH_RAW_CANCER = BASE_PATH_DATA + '/raw_data'
BASE_PATH_CANSIG_PP_CANCER = BASE_PATH_DATA + '/cansig_processed'
BASE_PATH_RAW_PBMC = BASE_PATH_DATA + '/raw_data/pbmc_citeseq.h5ad'
BASE_PATH_PREPROCESSED = BASE_PATH_DATA + '/preprocessed'
BASE_PATH_DGEX_CANCER = BASE_PATH_DATA + '/dgex_genes'
BASE_PATH_ANNOT_CANCER = BASE_PATH_DATA + '/annotations/cancer'
BASE_PATH_ANNOT_PBMC = BASE_PATH_DATA + '/annotations/citeseq_pbmc'

####################################################################################
### AVAILABLE DATASETS
# 'crc' : Karin Pelka, Matan Hofree, Jonathan H. Chen, Siranush Sarkizova, Joshua D. Pirl, Vjola
# Jorgji, Alborz Bejnood, Danielle Dionne, William H. Ge, Katherine H. Xu, Sherry X. Chao,
# Daniel R. Zollinger, David J. Lieb, Jason W. Reeves, Christopher A. Fuhrman, Margaret L.
# Hoang, Toni Delorey, Lan T. Nguyen, Julia Waldman, Max Klapholz, Isaac Wakiro, Ofir Cohen,
# Julian Albers, Christopher S. Smillie, Michael S. Cuoco, Jingyi Wu, Mei-ju Su, Jason Yeung,
# Brinda Vijaykumar, Angela M. Magnuson, Natasha Asinovski, Tabea Moll, Max N. Goder-Reiser,
# Anise S. Applebaum, Lauren K. Brais, Laura K. DelloStritto, Sarah L. Denning, Susannah T.
# Phillips, Emma K. Hill, Julia K. Meehan, Dennie T. Frederick, Tatyana Sharova, Abhay Kanodia,
# Ellen Z. Todres, Judit Jané-Valbuena, Moshe Biton, Benjamin Izar, Conner D. Lambden,
# Thomas E. Clancy, Ronald Bleday, Nelya Melnitchouk, Jennifer Irani, Hiroko Kunitake,
# David L. Berger, Amitabh Srivastava, Jason L. Hornick, Shuji Ogino, Asaf Rotem,
# Sébastien Vigneau, Bruce E. Johnson, Ryan B. Corcoran, Arlene H. Sharpe, Vijay K. Kuchroo,
# Kimmie Ng, Marios Giannakis, Linda T. Nieman, Genevieve M. Boland, Andrew J. Aguirre,
# Ana C. Anderson, Orit Rozenblatt-Rosen, Aviv Regev, Nir Hacohen, # Spatially organized
# multicellular immune hubs in human colorectal cancer, Cell, Volume 184, Issue 18, 2021,
# Pages 4734-4752.e20, ISSN 0092-8674, https://doi.org/10.1016/j.cell.2021.08.003 (
# https://www.sciencedirect.com/science/article/pii/S0092867421009454)

# 'escc' : Zhang, X., Peng, L., Luo, Y. et al. Dissecting esophageal squamous-cell carcinoma
# ecosystem by single-cell transcriptomic analysis. Nat Commun 12, 5291 (2021).
# https://doi.org/10.1038/s41467-021-25539-x

# 'luad_xing' : Xudong Xing et al., Decoding the multicellular ecosystem of lung adenocarcinoma
# manifested as pulmonary subsolid nodules by single-cell RNA sequencing.Sci. Adv.7,
# eabd9738(2021).DOI:10.1126/sciadv.abd9738

# 'luad_atlas' : Stefan Salcher, Gregor Sturm, Lena Horvath, Gerold Untergasser, Christiane
# Kuempers, Georgios Fotakis, Elisa Panizzolo, Agnieszka Martowicz, Manuel Trebo, Georg Pall,
# Gabriele Gamerith, Martina Sykora, Florian Augustin, Katja Schmitz, Francesca Finotello,
# Dietmar Rieder, Sven Perner, Sieghart Sopper, Dominik Wolf, Andreas Pircher,
# Zlatko Trajanoski, High-resolution single-cell atlas reveals diversity and plasticity of
# tissue-resident neutrophils in non-small cell  lung cancer, Cancer Cell, Volume 40, Issue 12,
# 2022, Pages 1503-1520.e8, ISSN 1535-6108, https://doi.org/10.1016/j.ccell.2022.10.008.

# 'breast' and subsets of it: Wu, S.Z., Al-Eryani, G., Roden, D.L. et al. A single-cell and
# spatially resolved atlas of human breast cancers. Nat Genet 53, 1334–1347 (2021).
# https://doi.org/10.1038/s41588-021-00911-1

# 'melanoma' : Tirosh I, Izar B, Prakadan SM, Wadsworth MH 2nd, Treacy D, Trombetta JJ,
# Rotem A, Rodman C, Lian C, Murphy G, Fallahi-Sichani M, Dutton-Regester K, Lin JR, Cohen O,
# Shah P, Lu D, Genshaft AS, Hughes TK, Ziegler CG, Kazer SW, Gaillard A, Kolb KE, Villani AC,
# Johannessen CM, Andreev AY, Van Allen EM, Bertagnolli M, Sorger PK, Sullivan RJ, Flaherty KT,
# Frederick DT, Jané-Valbuena J, Yoon CH, Rozenblatt-Rosen O, Shalek AK, Regev A, Garraway LA.
# Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. Science.
# 2016 Apr 8;352(6282):189-96. doi: 10.1126/science.aad0501. PMID: 27124452; PMCID: PMC4944528.

# 'skin': Ji et al. 2020 skin dataset Ji AL, Rubin AJ, Thrane K, Jiang S, Reynolds DL, Meyers RM,
# Guo MG, George BM, Mollbrink A, Bergenstråhle J, Larsson L, Bai Y, Zhu B, Bhaduri A, Meyers JM,
# Rovira-Clavé X, Hollmig ST, Aasi SZ, Nolan GP, Lundeberg J, Khavari PA.
# Multimodal Analysis of Composition and Spatial Architecture in Human Squamous Cell Carcinoma.
# Cell. 2020 Jul 23;182(2):497-514.e22. doi: 10.1016/j.cell.2020.05.039. Epub 2020 Jun 23. Erratum in: Cell.
# 2020 Sep 17;182(6):1661-1662. doi: 10.1016/j.cell.2020.08.043. PMID: 32579974; PMCID: PMC7391009.

# 'ovarian': # Ovarian cancer dataset
# Vázquez-García, I., Uhlitz, F., Ceglia, N. et al. Ovarian cancer mutational processes drive
# site-specific immune evasion. Nature 612, 778–786 (2022). https://doi.org/10.1038/s41586-022-05496-1

# 'pbmc': #TODO: add citation

from dataclasses import dataclass

DATASETS = [
    'crc', 'escc', 'luad_xing', 'breast', 'breast_large', 'breast_small', 'breast_malignant',
    'melanoma', 'luad_kim_malignant', 'luad_kim_malignant_2', 'skin_malignant', 'skin_malignant_2',
    'ovarian_malignant', 'ovarian_malignant_2',
    'pbmc_b_mono_nk', 'pbmc_b_subtypes', 'pbmc_cd4_subtypes', 'pbmc_cd8_subtypes', 'pbmc'
]

CANCER_DATASETS = DATASETS[0:-5]
PBMC_DATASETS = DATASETS[-5:]

DATASETS_WITH_ANNOTATIONS = [
    'breast_malignant', 'luad_kim_malignant', 'luad_kim_malignant_2', 'skin_malignant', 'skin_malignant_2',
    'ovarian_malignant', 'ovarian_malignant_2', 'pbmc_b_mono_nk', 'pbmc_b_subtypes', 'pbmc_cd4_subtypes',
    'pbmc_cd8_subtypes',
]


@dataclass
class CancerSignatureConfig:
    """Configuration for signature file loading."""
    file_path: str
    file_type: str = 'csv'
    name_transform: callable = None


CANCER_DATASET_SIGS_CONFIGS = {
    'breast_malignant': CancerSignatureConfig(
        file_path='breast_wu/wu_6.csv'
    ),
    'luad_kim_malignant': CancerSignatureConfig(
        file_path='luad_kim/kim_3.csv'
    ),
    'luad_kim_malignant_2': CancerSignatureConfig(
        file_path='luad_kim/kim_3.csv'
    ),
    'skin_malignant': CancerSignatureConfig(
        file_path='skin_ji/gene_sets.gmt',
        file_type='gmt',
        name_transform=lambda k: k.replace('_', ' ')
    ),
    'skin_malignant_2': CancerSignatureConfig(
        file_path='skin_ji/gene_sets.gmt',
        file_type='gmt',
        name_transform=lambda k: k.replace('_', ' ')
    ),
    'ovarian_malignant': CancerSignatureConfig(
        file_path='ovarian_vazquez/ovarian_states.csv'
    ),
    'ovarian_malignant_2': CancerSignatureConfig(
        file_path='ovarian_vazquez/ovarian_states.csv'
    )
}

PBMC_DEXG = BASE_PATH_ANNOT_PBMC + '/DE_by_celltype.csv'


@dataclass
class PBMCSignatureConfig:
    """Configuration for signature file loading."""
    high_level_low_level_mapping: Dict[str, List[str]]


PBMC_DATASET_SIGS_CONFIGS = {
    'pbmc_b_mono_nk': PBMCSignatureConfig(
        high_level_low_level_mapping={
            'B': ['B naive kappa', 'B memory kappa', 'B naive lambda',
                  'B memory lambda', 'B intermediate kappa',
                  'B intermediate lambda', 'Plasma', 'Plasmablast'],
            'Mono': ['CD14 Mono', 'CD16 Mono'],
            'NK': ['NK_2', 'NK_4', 'NK_1', 'NK_3', 'NK_CD56bright',
                   'NK Proliferating']
        }),
    'pbmc_b_subtypes': PBMCSignatureConfig(
        high_level_low_level_mapping={
            'B intermediate': ['B intermediate kappa', 'B intermediate lambda'],
            'B memory': ['B memory kappa', 'B memory lambda'],
            'B naive': ['B naive kappa', 'B naive lambda']
        }),
    'pbmc_cd4_subtypes': PBMCSignatureConfig(
        high_level_low_level_mapping={
            'CD4 CTL': ['CD4 CTL'],
            'CD4 Naive': ['CD4 Naive'],
            'CD4 Proliferating': ['CD4 Proliferating'],
            'CD4 TCM': ['CD4 TCM_1', 'CD4 TCM_3', 'CD4 TCM_2'],
            'CD4 TEM': ['CD4 TEM_3', 'CD4 TEM_1', 'CD4 TEM_2', 'CD4 TEM_4'],
            'Treg': ['Treg Naive', 'Treg Memory']
        }),
    'pbmc_cd8_subtypes': PBMCSignatureConfig(
        high_level_low_level_mapping={
            'CD8 Naive': ['CD8 Naive', 'CD8 Naive_2'],
            'CD8 Proliferating': ['CD8 Proliferating'],
            'CD8 TCM': ['CD8 TCM_1', 'CD8 TCM_3', 'CD8 TCM_2'],
            'CD8 TEM': ['CD8 TEM_2', 'CD8 TEM_1', 'CD8 TEM_4',
                        'CD8 TEM_5', 'CD8 TEM_6', 'CD8 TEM_3'],
        }),
}


@dataclass
class ViolinPlotConfig:
    """Configuration for signature file loading."""
    textwrap_width: int = 8
    height: float = 1.95
    aspect: float = 0.925
    sharey: bool = False
    wspace: float = 0.05
    col_wrap: int = 4
    legend_bbox_anchor: Tuple[float, float] = (1.15, 1)
    fontsizes: Dict[str, int] = None

    def __post_init__(self):
        if self.fontsizes is None:
            self.fontsizes = {'title': 12, 'labels': 11, 'ticks': 11, 'legend': 11}


VIOLIN_PLOT_CONFIG = {
    'breast_malignant': ViolinPlotConfig(
        aspect=2.5,
        wspace=0.075,
        col_wrap=2,
        legend_bbox_anchor=(1.125, 1),
    ),
    'luad_kim_malignant': ViolinPlotConfig(
        aspect=1.15,
        wspace=0.15,
        legend_bbox_anchor=(1.13, 1),
    ),
    'luad_kim_malignant_2': ViolinPlotConfig(
        aspect=1.15,
        wspace=0.15,
        legend_bbox_anchor=(1.13, 1),
    ),
    'skin_malignant': ViolinPlotConfig(
        textwrap_width=9,
        aspect=1.85,
        wspace=0.1,
        legend_bbox_anchor=(1.05, 1),
        col_wrap=2,
        fontsizes={'title': 12, 'labels': 11, 'ticks': 10, 'legend': 11}
    ),
    'skin_malignant_2': ViolinPlotConfig(
        textwrap_width=9,
        aspect=1.85,
        wspace=0.1,
        legend_bbox_anchor=(1.05, 1),
        col_wrap=2,
        fontsizes={'title': 12, 'labels': 11, 'ticks': 10, 'legend': 11}
    ),
    'ovarian_malignant': ViolinPlotConfig(
        textwrap_width=7,
        height=2,
        aspect=3,
        legend_bbox_anchor=(1, 1),
        col_wrap=2,
        fontsizes={'title': 12, 'labels': 10, 'ticks': 10, 'legend': 10}
    ),
    'ovarian_malignant_2': ViolinPlotConfig(
        textwrap_width=7,
        height=2,
        aspect=3,
        legend_bbox_anchor=(1, 1),
        col_wrap=2,
        fontsizes={'title': 12, 'labels': 10, 'ticks': 10, 'legend': 10}
    ),
    'pbmc_b_mono_nk': ViolinPlotConfig(
        textwrap_width=7,
        sharey=True,
        aspect=1.05,
        legend_bbox_anchor=(1.125, 1),
        fontsizes={'title': 12, 'labels': 11, 'ticks': 10, 'legend': 11}
    ),
    'pbmc_b_subtypes': ViolinPlotConfig(
        textwrap_width=7,
        sharey=True,
        aspect=1.05,
        legend_bbox_anchor=(1.075, 1),
        fontsizes={'title': 12, 'labels': 11, 'ticks': 10, 'legend': 11}
    ),
    'pbmc_cd4_subtypes': ViolinPlotConfig(
        textwrap_width=6,
        height=2.5,
        aspect=1.75,
        col_wrap=2,
        wspace=0.075,
        legend_bbox_anchor=(1.05, 1),
    ),
    'pbmc_cd8_subtypes': ViolinPlotConfig(
        textwrap_width=5,
        col_wrap=2,
        aspect=1.5,
        wspace=0.2,
        legend_bbox_anchor=(1.05, 1),
    ),
}

####################################################################################
# AVAILABLE NORMALIZATION METHODS
NORM_METHODS = ['mean', 'median', 'CP10k']

####################################################################################
# APPLICATION LEVEL OF DGEX
DGEX_GRANULARITY = ['all', 'individual', 'pseudobulk']

# FOR SAMPLE BASED DGEX PERCENTAGES OF OVERLAP
PCTGS = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.925, 0.95, 0.975, 0.99, 1]

####################################################################################
# DEFINED SCORING METHODS AND SETTINGS
SCORING_METHODS = [
    {
        "scoring_method": "adjusted_neighborhood_scoring",
        "sc_params": {
            "ctrl_size": 100,
            "score_name": "ANS",
        },
    },
    {
        "scoring_method": "seurat_scoring",
        "sc_params": {
            "ctrl_size": 100,
            "n_bins": 25,
            "score_name": "Seurat",
        },
    },
    {
        "scoring_method": "seurat_ag_scoring",
        "sc_params": {
            "n_bins": 25,
            "score_name": "Seurat_AG",
        },
    },
    {
        "scoring_method": "seurat_lvg_scoring",
        "sc_params": {
            "ctrl_size": 100,
            "n_bins": 25,
            "lvg_computation_version": "v1",
            "lvg_computation_method": "seurat",
            "score_name": "Seurat_LVG",
        },
    },
    {
        "scoring_method": "scanpy_scoring",
        "sc_params": {
            "ctrl_size": 100,
            "n_bins": 25,
            "score_name": "Scanpy",
        },
    },

    {
        "scoring_method": "jasmine_scoring",
        "sc_params": {
            "score_method": 'likelihood',
            "score_name": "Jasmine_LH",
        },
    },
    {
        "scoring_method": "jasmine_scoring",
        "sc_params": {
            "score_method": 'oddsratio',
            "score_name": "Jasmine_OR",
        },
    },
    {
        "scoring_method": "ucell_scoring",
        "sc_params": {
            "score_name": "UCell",
            "maxRank": 1500,
        },
    },
]
METHOD_WO_MEAN = ['scanpy_scoring', 'ucell_scoring', 'jasmine_scoring']
METHOD_WITH_GENE_POOL = ['adjusted_neighborhood_scoring', 'seurat_scoring', 'seurat_ag_scoring', 'seurat_lvg_scoring']
