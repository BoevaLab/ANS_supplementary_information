## TODO: set correct working directory
setwd('.....')
library(Seurat)
library(SeuratDisk)
library(SingleCellExperiment)
library(data.table)
## TODO: Path to the 'construction_scoring_methods/helper_functions.R'
source('....')

## TODO: set base stroing path, e.g., 'experiments/construction_scoring_methods/'
base_path <- '....'

dataset <- 'crc'

## TODO: set base path where preprocessed datasets in RDS format lie
(data_rds <- paste0('....','pp_',dataset,'.rds'))
(storing_path <- paste0(base_path, dataset,'/'))
(sample_path <- paste0(base_path, dataset,'_sample_cells.csv'))

#load data as SingleCellExperiment
data_sce <- readRDS(data_rds)
data_seurat <- as.Seurat(data_sce, counts = 'counts', data = 'X')

sample_cells <- colnames(read.csv(sample_path))
sample_cells <- gsub("\\.", "-", sample_cells)
data_seurat <- data_seurat[, sample_cells]

rm(data_sce)
rm(data_rds)
gc()

#load gene signature
## TODO: set path where signatures are stored 'data/dgex_genes'
dgex_genes <- paste0('.....',dataset,'/mean_norm/on_pseudobulk/dgex_genes.csv')
wc <- read.csv(dgex_genes)
gene_list <- wc[order(wc$pvalue, -wc$log2FoldChange),][1:100,]$genes
markers <- list(markers_for_mal_cells = gene_list)
rm(wc)
rm(gene_list)
gc()



#### Score with AddModuleScore_UCell
set.seed(123)

obj <- data_seurat
#obj <- AddModuleScore_UCell(obj, features = markers)

features <- check_signature_names(markers)

assay=NULL

if (is.null(assay)) {
  assay <- Seurat::DefaultAssay(obj)
}

maxRank=1500
chunk.size=1000
BPPARAM=NULL
ncores=1
storeRanks=FALSE
w_neg=1
slot="data"
ties.method="average"
force.gc=FALSE
name="_UCell"

meta.list <- calculate_Uscore(
  Seurat::GetAssayData(obj, slot, assay=assay),
  features=features, maxRank=maxRank,
  chunk.size=chunk.size, w_neg=w_neg,
  ncores=ncores, BPPARAM=BPPARAM, ties.method=ties.method,
  force.gc=force.gc, storeRanks=storeRanks, name=name)

meta.merge <- lapply(meta.list,function(x) rbind(x[["cells_AUC"]]))
meta.merge <- Reduce(rbind, meta.merge)

obj <- Seurat::AddMetaData(obj, as.data.frame(meta.merge))

curr_storing_path <- paste0(base_path, dataset,'_debug_ucell_scores.csv')
write.csv(FetchData(obj, vars='markers_for_mal_cells_UCell'), curr_storing_path)
gc()



