## Difussion map pseudotime analysis
library(Seurat)
library(destiny)
library(dplyr)
library(Biobase)
library(ggplot2)
library(slingshot)
library(scater)

set.seed(333)

path2project  <- '/Users/carlosramirez/sc/tcd8ExscSeq/'
setwd(path2project)

## Loading seurat object with integrated data Miller-Hoffmann
seurat  <- readRDS('analysis/integrated_miller_hoffmann.rds')

## Subsetting 
#seurat <- subset(seurat, 
#                 dataset == 'Hoffmann' &       ## All Hoffmann cell
#                 predicted != 'Proliferating') ## except proliferating
seurat <- subset(seurat,                        ## All cells
                 predicted %in% c('Progenitor Ex',
                                  'Terminally Ex')) ## except proliferating


## Creation of the eset object 
eset <- ExpressionSet(
        assayData = seurat@assays$SCT@scale.data,
        phenoData = AnnotatedDataFrame(seurat@meta.data)
)


## Difussion map
difmap <- DiffusionMap(eset, 
                       n_pcs = 20, 
                       n_eigs = 20, 
                       k=5)

## Plotting
pdf('figures/diffusion_map_integrated_Miller_Hoffmann_exhausted.pdf')
plot(difmap, 1:2, col_by='predicted') +
        theme_bw() +
        theme(panel.grid = element_blank())
dev.off()


