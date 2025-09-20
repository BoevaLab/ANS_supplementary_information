import sys

sys.path.append('..')
sys.path.append('../..')

from constants import BASE_PATH_RAW_CANCER, BASE_PATH_CANSIG_PP_CANCER

############################################################################################################################
#### Ovarian cancer ####
# Vázquez-García, I., Uhlitz, F., Ceglia, N. et al. Ovarian cancer mutational 
# processes drive site-specific immune evasion. Nature 612, 778–786 (2022). 
# https://doi.org/10.1038/s41586-022-05496-1
print("Subsetting ovarian cancer dataset for only malignant cells...")
ovarian_adata = sc.read_h5ad(os.path.join(BASE_PATH_RAW_CANCER, 'OV_Vazquez_10X.h5ad'))
print(f">> {ovarian_adata.shape=}")
ovarian_adata = ovarian_adata[ovarian_adata.obs['cell_type'] == 'Ovarian.cancer.cell'].copy()
ovarian_adata = ovarian_adata[~ovarian_adata.obs['cluster_label'].str.startswith('Ciliated.cell', na=False)].copy()
ovarian_adata = ovarian_adata[~ovarian_adata.obs['cluster_label'].isna()].copy()
print(f">> {ovarian_adata.shape=}")
print(f">> {ovarian_adata.obs['cluster_label'].value_counts(dropna=False).sort_index()}")
ovarian_adata.write(os.path.join(BASE_PATH_RAW_CANCER, 'ovarian_malignant.h5ad'))


############################################################################################################################
#### Skin cancer ####
# Ji AL, Rubin AJ, Thrane K, Jiang S, Reynolds DL, Meyers RM, Guo MG, George BM, Mollbrink A, Bergenstråhle J, Larsson L, 
# Bai Y, Zhu B, Bhaduri A, Meyers JM, Rovira-Clavé X, Hollmig ST, Aasi SZ, Nolan GP, Lundeberg J, Khavari PA. Multimodal 
# Analysis of Composition and Spatial Architecture in Human Squamous Cell Carcinoma. Cell. 2020 Jul 23;182(2):497-514.e22. 
# doi: 10.1016/j.cell.2020.05.039. Epub 2020 Jun 23. Erratum in: Cell. 2020 Sep 17;182(6):1661-1662. doi: 10.1016/j.cell.2020.08.043. 
# PMID: 32579974; PMCID: PMC7391009.

print("Subsetting skin cancer dataset for only malignant cells...")
skin_adata = sc.read_h5ad(os.path.join(BASE_PATH_CANSIG_PP_CANCER, 'skin.h5ad'))
print(f">> {skin_adata.shape=}")
skin_adata = skin_adata[skin_adata.obs['level2_celltype'].isin(['TSK', 'Tumor_KC_Basal', 'Tumor_KC_Cyc', 'Tumor_KC_Diff'])].copy()
print(f">> {skin_adata.shape=}")
print(f">> {skin_adata.obs['level2_celltype'].value_counts(dropna=False).sort_index()}")
skin_adata.write(os.path.join(BASE_PATH_CANSIG_PP_CANCER, 'skin_malignant.h5ad'))


############################################################################################################################
#### Lung cancer ####
print("Subsetting LUAD KIM cancer dataset for only malignant cells...")
luad_kim_adata = sc.read_h5ad(os.path.join(BASE_PATH_CANSIG_PP_CANCER, 'luad_kim.h5ad'))
print(f">> {luad_kim_adata.shape=}")
samples_in_adata = luad_kim_adata.obs.index.to_list()
cell_labels = pd.read_table(Path(BASE_PATH_ANNOT_CANCER) / 'luad_kim' / "GSE131907_Lung_Cancer_cell_annotation.txt")
cell_labels['Index'] = cell_labels['Index'].str.replace('_', '-')
cell_labels = cell_labels.set_index('Index')
cell_labels = cell_labels.loc[samples_in_adata]

y_true_col = 'Cell_subtype'
sample_col = 'sample'

luad_kim_adata.obs[y_true_col] = cell_labels[y_true_col].str.lower()
luad_kim_adata = luad_kim_adata[luad_kim_adata.obs[luad_kim_adata.obs[y_true_col].str.startswith('ts', na=False)].index].copy()

print(f">> {luad_kim_adata.shape=}")
print(f">> {luad_kim_adata.obs[y_true_col].value_counts().sort_index()}")
luad_kim_adata.write(os.path.join(BASE_PATH_CANSIG_PP_CANCER, 'luad_kim_malignant.h5ad'))

print(f"Now we also subset the data from 3CA, to compare the preprocessing steps")
ds_obj = DATA3CA(data_folder=Path(BASE_PATH_RAW_CANCER) / "luad_kim_2020/Data_Kim2020_Lung")
ds_obj.validate_data_consistency()
luad_kim_3ca_adata = sc.AnnData(
    X=ds_obj.counts,
    obs=ds_obj.cells.reset_index().merge(ds_obj.metadata, on='sample').set_index('cell_name'),
    var=ds_obj.genes
)
luad_kim_3ca_adata = luad_kim_3ca_adata[luad_kim_3ca_adata.obs['cell_type'] == 'Malignant'].copy()
samples_in_adata = luad_kim_3ca_adata.obs.index.to_list()
cell_labels = pd.read_table(Path(BASE_PATH_ANNOT_CANCER) / 'luad_kim' /"GSE131907_Lung_Cancer_cell_annotation.txt")
cell_labels = cell_labels.set_index('Index')
cell_labels = cell_labels.loc[samples_in_adata, :]
y_true_col = 'Cell_subtype'
sample_col = 'sample'
luad_kim_3ca_adata.obs[y_true_col] = cell_labels[y_true_col].str.lower()
luad_kim_3ca_adata = luad_kim_3ca_adata[luad_kim_3ca_adata.obs[luad_kim_3ca_adata.obs[y_true_col].str.startswith('ts', na=False)].index].copy()

luad_kim_3ca_adata.write(os.path.join(BASE_PATH_CANSIG_PP_CANCER, 'luad_kim_malignant_3ca.h5ad'))


############################################################################################################################
#### Breast cancer ####
print("Subsetting breast cancer dataset for only malignant cells...")

full_breast_adata = sc.read_h5ad(os.path.join(BASE_PATH_CANSIG_PP_CANCER, 'breast.h5ad'))
print(f">> {full_breast_adata.shape=}")
breast_adata = full_breast_adata[full_breast_adata.obs.malignant_key == 'malignant'].copy()

print(f">> {breast_adata.shape=}")
print(f">> {breast_adata.obs.malignant_key.value_counts().sort_index()}")
breast_adata.write(os.path.join(BASE_PATH_CANSIG_PP_CANCER, 'breast_malignant.h5ad'))


print("Subsetting breast cancer dataset for basal-like samples...")
basal_samples = [
    'CID3586', 'CID3963',
    'CID4465', 'CID4495',
    'CID44971', 'CID4513',
    'CID4515', 'CID4523',
    ]
breast_basal_adata = full_breast_adata[full_breast_adata.obs.sample_id.isin(basal_samples)].copy()

print(f">> {breast_basal_adata.shape=}")
print(f">> {breast_basal_adata.obs.sample_id.value_counts().sort_index()}")
breast_basal_adata.write(os.path.join(BASE_PATH_CANSIG_PP_CANCER, 'breast_basal_like_samples.h5ad'))


