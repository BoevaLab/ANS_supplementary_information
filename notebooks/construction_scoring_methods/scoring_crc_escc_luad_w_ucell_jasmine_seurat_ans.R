################################################################################
################################################################################
################################################################################
################################################################################

## TODO: set correct working directory
setwd('.....')
library(Seurat)
library(UCell)
library(SeuratDisk)
library(SingleCellExperiment)
source('MT/JASMINE-main/JASMINE_V1_11October2021.r')
source('MT/ANS_signature_scoring/src_R/adjusted_neighborhood_scoring.R')

## TODO: set base storing path
base_storing_path <- '......'

## TODO: set base path preprocessed seurat objects 
base_path_preprocessed_seurat <- '....'

## TODO: set base path to gene signatures
base_path_signatures <- '....'


## CRC
data_rds <- paste0(base_path_preprocessed_seurat, 'pp_crc.rds')
storing_path <- paste0(base_storing_path, 'crc/',sep = "")

#load data as SingleCellExperiment
crc_data_sce <- readRDS(data_rds)

#load gene signature
dgex_genes <- paste0(base_path_signatures, 'crc/mean_norm/on_pseudobulk/dgex_genes.csv')
wc <- read.csv(dgex_genes)
gene_list <- wc[order(wc$pvalue, -wc$log2FoldChange),][1:100,]$genes
rm(wc)
gc()

#### Score with Jasmine OR
set.seed(123)

jas_or_scores  =  JASMINE(assay(crc_data_sce, 'X'), gene_list, method=c('oddsratio')) ## calling JASMINE with oddsratio
curr_storing_path <- paste0(storing_path, 'crc_jas_or_scores.csv')
write.csv(jas_or_scores, curr_storing_path, row.names = FALSE)
rm(jas_or_scores)
gc()

#### Score with Jasmine LH
set.seed(123)
jas_lh_scores  =  JASMINE(assay(crc_data_sce, 'X'), gene_list, method = c('likelihood')) ## calling JASMINE with likelihood
curr_storing_path <- paste0(storing_path, 'crc_jas_lh_scores.csv')
write.csv(jas_lh_scores, curr_storing_path, row.names = FALSE)
rm(jas_lh_scores)
gc()

#### Score with AddModuleScore
#load data in Seurat object
crc_data_seurat <- as.Seurat(crc_data_sce, counts = 'counts', data = 'X')
rm(crc_data_sce)
markers <- list(markers_for_mal_cells = gene_list)
set.seed(123)
crc_data_seurat <- AddModuleScore(crc_data_seurat, features = markers, name = "AddModuleScore")
curr_storing_path <- paste0(storing_path, 'crc_addmodulescore_scores.csv')
write.csv(FetchData(crc_data_seurat, vars='AddModuleScore1'), curr_storing_path)
gc()

#### Score with AddModuleScore_UCell
set.seed(123)
crc_data_seurat <- AddModuleScore_UCell(crc_data_seurat, features = markers)
curr_storing_path <- paste0(storing_path, 'crc_ucell_scores.csv')
write.csv(FetchData(crc_data_seurat, vars='markers_for_mal_cells_UCell'), curr_storing_path)
gc()

#### Score with AdjustedNeighborhoodScoring
crc_data_seurat <- as.Seurat(crc_data_sce, counts = 'counts', data = 'X')
rm(crc_data_sce)
markers <- list(markers_for_mal_cells = gene_list)
set.seed(123)
crc_data_seurat <- AdjustedNeighborhoodScoring(crc_data_seurat, features = markers)
curr_storing_path <- paste0(storing_path, 'crc_ans_scores.csv')
write.csv(FetchData(crc_data_seurat, vars='ANS_scores1'), curr_storing_path)
rm(list = ls())
gc()

################################################################################
################################################################################
################################################################################
################################################################################
## TODO: set correct working directory
setwd('.....')
library(Seurat)
library(UCell)
library(SeuratDisk)
library(SingleCellExperiment)
source('MT/JASMINE-main/JASMINE_V1_11October2021.r')
source('MT/ANS_signature_scoring/src_R/adjusted_neighborhood_scoring.R')

## TODO: set base storing path
base_storing_path <- '......'

## TODO: set base path preprocessed seurat objects 
base_path_preprocessed_seurat <- '....'

## TODO: set base path to gene signatures
base_path_signatures <- '....'

## ESCC
data_rds <- paste0(base_path_preprocessed_seurat, 'pp_escc.rds')
storing_path <- paste0(base_storing_path, 'escc/')

#load data as SingleCellExperiment
escc_data_sce <- readRDS(data_rds)

#load gene signature
dgex_genes <- paste0(base_path_signatures, 'escc/mean_norm/on_pseudobulk/dgex_genes.csv')
wc <- read.csv(dgex_genes)
gene_list <- wc[order(wc$pvalue, -wc$log2FoldChange),][1:100,]$genes
rm(wc)
gc()

#### Score with Jasmine OR
set.seed(123)
jas_or_scores  =  JASMINE(assay(escc_data_sce, 'X'), gene_list, method=c('oddsratio')) ## calling JASMINE with oddsratio
gc()
curr_storing_path <- paste0(storing_path, 'escc_jas_or_scores.csv')
write.csv(jas_or_scores, curr_storing_path, row.names = FALSE)
rm(jas_or_scores)
gc()

#### Score with Jasmine LH
set.seed(123)
jas_lh_scores  =  JASMINE(assay(escc_data_sce, 'X'), gene_list, method = c('likelihood')) ## calling JASMINE with likelihood
gc()
curr_storing_path <- paste0(storing_path, 'escc_jas_lh_scores.csv')
write.csv(jas_lh_scores, curr_storing_path, row.names = FALSE)
rm(jas_lh_scores)
gc()

#### Score with AddModuleScore
#load data in Seurat object
escc_data_seurat <- as.Seurat(escc_data_sce, counts = 'counts', data = 'X')
rm(escc_data_sce)
markers <- list(markers_for_mal_cells = gene_list)
set.seed(123)
escc_data_seurat <- AddModuleScore(escc_data_seurat, features = markers, name = "AddModuleScore")
curr_storing_path <- paste0(storing_path, 'escc_addmodulescore_scores.csv')
write.csv(FetchData(escc_data_seurat, vars='AddModuleScore1'), curr_storing_path)
gc()

#### Score with AddModuleScore_UCell
set.seed(123)
escc_data_seurat <- AddModuleScore_UCell(escc_data_seurat, features = markers)
curr_storing_path <- paste0(storing_path, 'escc_ucell_scores.csv')
write.csv(FetchData(escc_data_seurat, vars='markers_for_mal_cells_UCell'), curr_storing_path)
gc()

#### Score with AdjustedNeighborhoodScoring
escc_data_seurat <- as.Seurat(escc_data_sce, counts = 'counts', data = 'X')
rm(escc_data_sce)
markers <- list(markers_for_mal_cells = gene_list)
set.seed(123)
escc_data_seurat <- AdjustedNeighborhoodScoring(escc_data_seurat, features = markers)
curr_storing_path <- paste0(storing_path, 'escc_ans_scores.csv')
write.csv(FetchData(escc_data_seurat, vars='ANS_scores1'), curr_storing_path)
rm(list = ls())
gc()


################################################################################
################################################################################
################################################################################
################################################################################
## LUAD
## TODO: set correct working directory
setwd('.....')
library(Seurat)
library(UCell)
library(SeuratDisk)
library(SingleCellExperiment)
source('MT/JASMINE-main/JASMINE_V1_11October2021.r')
source('MT/ANS_signature_scoring/src_R/adjusted_neighborhood_scoring.R')


## TODO: set base storing path
base_storing_path <- '......'

## TODO: set base path preprocessed seurat objects 
base_path_preprocessed_seurat <- '....'

## TODO: set base path to gene signatures
base_path_signatures <- '....'

data_rds <- paste0(base_path_preprocessed_seurat, 'pp_luad.rds')
storing_path <- paste0(base_storing_path, 'luad/',sep = "")

#load data as SingleCellExperiment
luad_data_sce <- readRDS(data_rds)

#load gene signature
dgex_genes <- paste0(base_path_signatures, 'luad/mean_norm/on_pseudobulk/dgex_genes.csv')
wc <- read.csv(dgex_genes)
gene_list <- wc[order(wc$pvalue, -wc$log2FoldChange),][1:100,]$genes
rm(wc)
gc()

#### Score with Jasmine OR
set.seed(123)

jas_or_scores  =  JASMINE(assay(luad_data_sce, 'X'), gene_list, method=c('oddsratio')) ## calling JASMINE with oddsratio
curr_storing_path <- paste0(storing_path, 'luad_jas_or_scores.csv')
write.csv(jas_or_scores, curr_storing_path, row.names = FALSE)
rm(jas_or_scores)
gc()

#### Score with Jasmine LH
set.seed(123)
jas_lh_scores  =  JASMINE(assay(luad_data_sce, 'X'), gene_list, method = c('likelihood')) ## calling JASMINE with likelihood
curr_storing_path <- paste0(storing_path, 'luad_jas_lh_scores.csv')
write.csv(jas_lh_scores, curr_storing_path, row.names = FALSE)
rm(jas_lh_scores)
gc()

#### Score with AddModuleScore
#load data in Seurat object
luad_data_seurat <- as.Seurat(luad_data_sce, counts = 'counts', data = 'X')
rm(luad_data_sce)
markers <- list(markers_for_mal_cells = gene_list)
set.seed(123)
luad_data_seurat <- AddModuleScore(luad_data_seurat, features = markers, name = "AddModuleScore")
curr_storing_path <- paste0(storing_path, 'luad_addmodulescore_scores.csv')
write.csv(FetchData(luad_data_seurat, vars='AddModuleScore1'), curr_storing_path)
gc()

#### Score with AddModuleScore_UCell
set.seed(123)
luad_data_seurat <- AddModuleScore_UCell(luad_data_seurat, features = markers)
curr_storing_path <- paste0(storing_path, 'luad_ucell_scores.csv')
write.csv(FetchData(luad_data_seurat, vars='markers_for_mal_cells_UCell'), curr_storing_path)
gc()

#### Score with AdjustedNeighborhoodScoring
luad_data_seurat <- as.Seurat(luad_data_sce, counts = 'counts', data = 'X')
rm(luad_data_sce)
markers <- list(markers_for_mal_cells = gene_list)
set.seed(123)
luad_data_seurat <- AdjustedNeighborhoodScoring(luad_data_seurat, features = markers)
curr_storing_path <- paste0(storing_path, 'luad_ans_scores.csv')
write.csv(FetchData(luad_data_seurat, vars='ANS_scores1'), curr_storing_path)
rm(list = ls())
gc()