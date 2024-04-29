## DEFINED CONSTANTS
####################################################################################
### PATHS
## !!!!!! TODO: CHANGE path to folder where data and experiments should go !!!!!!!
BASE_PATH_DRIVE = '/Users/lciernik/Documents/TUB/projects/ans_scoring'

BASE_PATH_EXPERIMENTS = BASE_PATH_DRIVE + '/experiments'
BASE_PATH_DATA = BASE_PATH_DRIVE + '/data'
BASE_PATH_RAW_CANCER = BASE_PATH_DATA +'/cansig_processed'
BASE_PATH_RAW_PBMC = BASE_PATH_DATA +'/raw_data/pbmc_citeseq.h5ad'
BASE_PATH_PREPROCESSED = BASE_PATH_DATA +'/preprocessed'
BASE_PATH_DGEX_CANCER = BASE_PATH_DATA +'/dgex_genes'

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
# 'melanoma' : Tirosh I, Izar B, Prakadan SM, Wadsworth MH 2nd, Treacy D, Trombetta JJ,
# Rotem A, Rodman C, Lian C, Murphy G, Fallahi-Sichani M, Dutton-Regester K, Lin JR, Cohen O,
# Shah P, Lu D, Genshaft AS, Hughes TK, Ziegler CG, Kazer SW, Gaillard A, Kolb KE, Villani AC,
# Johannessen CM, Andreev AY, Van Allen EM, Bertagnolli M, Sorger PK, Sullivan RJ, Flaherty KT,
# Frederick DT, Jané-Valbuena J, Yoon CH, Rozenblatt-Rosen O, Shalek AK, Regev A, Garraway LA.
# Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. Science.
# 2016 Apr 8;352(6282):189-96. doi: 10.1126/science.aad0501. PMID: 27124452; PMCID: PMC4944528.

DATASETS = ['crc', 'escc', 'luad_xing', 'breast_large', 'breast_small', 'melanoma', 'pbmc_b_mono_nk', 'pbmc_b_subtypes', 'pbmc_cd4_subtypes', 'pbmc_cd8_subtypes', 'pbmc']
CANCER_DATASETS = DATASETS[0:-4]

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
