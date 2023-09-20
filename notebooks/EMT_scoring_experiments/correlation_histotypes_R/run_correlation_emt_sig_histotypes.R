##
## The following notebook computes the correlations of the ESCC- and LUAD- specific signature 
## with histopathological subtypes. 
## To run the following script ensure to download the TCGA dataset as indicated in the paper. 
##

library("stringr")
library("dplyr")
library("GSVA")


pathExpression <- <pathToTCGAexpression>
outputPath <- <outputPath>
pathIn <- "Lung1_escc2.txt"
pathClinical <- "TCGA-CDR-SupplementalTableS1.csv"


#List all cancer types available for analysis
cancers <- gsub(".txt", "", list.files(pathExpression))

#Read the signature from file
signature_table <- read.table(pathIn, header=TRUE, sep='\n')
signature <- list()
signature[["EMT"]] <- signature_table[,1]

#Read in clinical data for the cancer
clinical = read.csv(pathClinical, header=TRUE, sep='\t', fill=TRUE, row.names=1)


#Dat frame to store results
result_kluskal <- data.frame(cancer = character(), feature = character(), pval = numeric())

#Calculate signature score in each cancer and check link with histological subtype:
for (can in cancers){
  #Read in expression
  expression <- read.table(paste0(pathExpression, can, ".txt"), header=TRUE, row.names=1, sep='\t')
  #Calculate signature score
  gsva.score <- gsva(as.matrix(expression), signature, method="ssgsea")
  scores = data.frame(t(gsva.score))
  #Extract clinical data for a particular cancer type
  cancer_clinical <- clinical[clinical$type == can,]
  #Unify sample names
  sample_labels <- gsub("-",".",data.frame(unlist(cancer_clinical$bcr_patient_barcode))[,1])
  rownames(cancer_clinical) <- sample_labels
  #Remove any samples with multiple values of gene expression
  duplicated_expression <- which(duplicated(str_sub(rownames(scores), start=1, end=12)))
  if (length(duplicated_expression)>0){
    scores <- data.frame(scores[-duplicated_expression,])
    rownames(scores) <- rownames(t(gsva.score))[-duplicated_expression]
  }
  rownames(scores) <- str_sub(rownames(scores), start=1, end=12)
  #Clear the data
  cancer_clinical[cancer_clinical == "[Not Applicable]"] <- NA
  cancer_clinical[cancer_clinical == "[Not Available]"] <- NA
  cancer_clinical[cancer_clinical == "#N/A"] <- NA
  cancer_clinical[cancer_clinical==""] <- NA
  #Select only overlaping samples in expression and clinical data
  keep <- intersect(rownames(scores), rownames(cancer_clinical))
  cancer_clinical <- cancer_clinical[keep,]
  scores <- scores[keep,]
  #Combine scores and clinical data together
  cancer_clinical <- cbind(cancer_clinical,scores) 
  check_remove <- c("NA", NA)
  unknown = which(clinical[,"histological_type"] %in% check_remove)
  if (length(unknown) > 0) {
    #Remove samples without histological type specified and perform analysis
    cancer_clinical <- cancer_clinical[-unknown,]
    if (length(unique(cancer_clinical[,"histological_type"]))>1){
      #Runs only if more than 1 histological type available
      an.error.occured <- FALSE
      tryCatch( { result <- kruskal.test(cancer_clinical$scores, cancer_clinical$histological_type);  }
                , error = function(e) {an.error.occured <<- TRUE})
      if (an.error.occured == TRUE){
        print("Error")
        next
      }else{
        result =  kruskal.test(cancer_clinical$scores, cancer_clinical$histological_type)
        out <- data.frame(cancer = can, feature = 'histological_type', pval = result$p.value)
        result_kluskal <- rbind(result_kluskal, out)
      }
    }
  } else {
    if (length(unique(cancer_clinical[,'histological_type']))>1){
      #Runs only if more than 1 histological type available
      an.error.occured <- FALSE
      tryCatch( { result <- kruskal.test(cancer_clinical$scores, cancer_clinical$histological_type);  }
                , error = function(e) {an.error.occured <<- TRUE})
      if (an.error.occured == TRUE){
        print("Error")
        next
      } else{
        result =  kruskal.test(cancer_clinical$scores, cancer_clinical$histological_type)
        out <- data.frame(cancer = can, feature = "histological_type", pval = result$p.value)
        result_kluskal <- rbind(result_kluskal, out)
      }
    }
  }
}

#Calculate FDR corrected p-values:
result_kluskal$padjusted <- p.adjust(result_kluskal$pval, method="fdr")

#Save results
write.table(result_kluskal, file = paste0(outputPath, "EMT_histotype.txt"), quote=FALSE, sep='\t', row.names=FALSE)

# > sessionInfo()
# R version 4.1.2 (2021-11-01)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 22.04.3 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
# 
# locale:
#   [1] LC_CTYPE=C.UTF-8           LC_NUMERIC=C              
# [3] LC_TIME=de_CH.UTF-8        LC_COLLATE=C.UTF-8        
# [5] LC_MONETARY=de_CH.UTF-8    LC_MESSAGES=C.UTF-8       
# [7] LC_PAPER=de_CH.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=de_CH.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] GSVA_1.42.0      dplyr_1.0.0.9000 stringr_1.4.0