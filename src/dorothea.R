## Transcription factor activity analysis using Dorothea coupled to 
## viper

## Dependencies
library(Seurat)
library(dplyr)
library(dorothea)
library(bcellViper)
library(viper)
library(pheatmap)

set.seed(333)

path2project <- '/Users/carlosramirez/sc/tcd8ExscSeq/'
setwd(path2project)


## Loading data
seurat <- readRDS('data/maike2020/hofmann_hcv_seu.rds')
seurat

###########################################################
## Implementing Dorothea
regulons <- dorothea_hs %>%
        filter(confidence %in% c("A", "B"))

seurat <- run_viper(
        seurat, 
        regulons,
        options = list(method = "scale", 
                       minsize = 1, 
                       eset.filter = FALSE,
                       cores = 3, 
                       verbose = FALSE)
)
seurat@assays$dorothea@data[1:5, 1:5]
## Extracting regulons
dorothea.tfa <- seurat@assays$dorothea@data

rm(seurat)

################################################################
## Loading SCENIC results
scenic <- read.csv('analysis/scenic/scenic_hgv_6000/aucell.csv')
head(scenic)

################################################################
## Comparing both results

## Processing Dorothea table
dorothea.tfa <- t(dorothea.tfa)
colnames(dorothea.tfa) <- paste0('dorothea_', 
                                 colnames(dorothea.tfa))

## Processing scenic results
rownames(scenic) <- scenic$Cell
scenic <- select(scenic, -Cell)
colnames(scenic) <- gsub('\\.', '', colnames(scenic))
colnames(scenic) <- paste0('scenic_', colnames(scenic))
scenic[1:5, 1:5]

## Correlating results
cor.tbl <- cor(dorothea.tfa, scenic)

pdf('figures/scenic_vs_dorothea.pdf',
    height = 20, width = 20)
pheatmap(cor.tbl)
dev.off()
