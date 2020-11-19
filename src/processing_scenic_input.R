## Preprocessing scenic input

## Installing Seurat if is not already installed
if ( ! requireNamespace('Seurat', quietly = TRUE) ) {
        devtools::install_github(repo = 'satijalab/seurat', ref = 'develop')
}

## Dependencies
library(tidyverse)
library(Seurat)

path2project <- '/Users/carlosramirez/sc/tcd8ExscSeq/'
setwd(path2project)

file_path <- 'data/maike2020/nina_thimme_raw_counts.csv'

## Input: A seurat object
## Performs normalization, scaling, pca and umap 
## Output: A processed seurat object
st_workflow <- function(
        seurat_object,
        n_features = 3000,
        n_pca_dims = 15
){
        cat('Normalizing and finding variable features\n')
        seurat.p <- NormalizeData(seurat_object) %>%
                FindVariableFeatures(selection.method = 'vst',
                                     nfeatures = n_features)
        cat('Scaling and projection\n')
        seurat.p <- ScaleData(seurat.p, 
                              verbose = FALSE) %>% 
                RunPCA(npcs = n_pca_dims, 
                       verbose = FALSE) %>%
                RunUMAP(reduction = "pca", 
                        dims = 1:n_pca_dims) 
        return(seurat.p)
}


## Maike Hoffman data preprocessing
## gene names must be standardized
himmer <- read.csv(file_path, 
                   header = TRUE)
r_names <- himmer$X
himmer <- select(himmer, -X)
himmer_mtx <- apply(himmer, 2, as.numeric)
rownames(himmer_mtx) <- gsub('__chr.*', '', r_names) 

## Seurat object construction
himmer_seu <- CreateSeuratObject(counts = himmer_mtx, 
                                 project = 'himmer', 
                                 min.cells = 1, 
                                 min.features = 1,
                                 assay = 'RNA')

himmer_seu <- SCTransform(himmer_seu, 
                          variable.features.n = 10000)

dir.create('analysis/')
dir.create('analysis/scenic/')

write.table(t(himmer_seu@assays$SCT@scale.data),
            file = 'analysis/scenic/maike_scenic_input.tsv',
            sep = '\t')

