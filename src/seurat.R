## Installing Seurat if is not already installed
if ( ! requireNamespace('Seurat', quietly = TRUE) ) {
     devtools::install_github(repo = 'satijalab/seurat', ref = 'develop')
}

## Dependencies
library(tidyverse)
library(Seurat)

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
himmer <- read.csv('data/maike2020/nina_thimme_raw_counts.csv', 
                   header = TRUE)
r_names <- himmer$X
himmer <- select(himmer, -X)
himmer_mtx <- apply(himmer, 2, as.numeric)
rownames(himmer_mtx) <- gsub('__chr.*', '', r_names) 

## Seurat object construction
himmer_seu <- CreateSeuratObject(counts = himmer_mtx, project = 'himmer', min.cells = 1, assay = 'rna')
himmer_seu <- st_workflow(himmer_seu, n_features = 1000, n_pca_dims = 100)

## Splitting Hoffman data for liger pipeline
patients <- himmer_seu$orig.ident
himmer_df <- as.data.frame(t(himmer_mtx))
himmer_sp <- split(himmer_df, f = patients) 
himmer_sp <- lapply(himmer_sp, function(x) t(x) )

## Standardizing gene names
## Gene names should be in upper case in order to align with human data.
## Run once in order to create miller2019_upper_Case/
if ( ! dir.exists('data/miller2019_upper_Case')) {
      dir.create('data/miller2019_upper_Case')  
      file.copy('data/miller2019/matrix.mtx', 'data/miller2019_upper_Case/')
      file.copy('data/miller2019/barcodes.tsv', 'data/miller2019_upper_Case/')
      genes_lc <- read.table('data/miller2019/genes.tsv', header = FALSE)
      genes_lc$V2 <- toupper(genes_lc$V2)
      write_tsv(genes_lc, 'data/miller2019_upper_Case/genes.tsv', col_names = FALSE)
}

## Creating seurat objects for every sample
himmer_sp_seu <- lapply(1:length(himmer_sp), 
       function(i) CreateSeuratObject(counts = himmer_sp[i][[1]], 
                                      project = 'hcv', 
                                      assay = 'RNA', 
                                      min.cells = 1, 
                                      min.features = 1) 
)
names(himmer_sp_seu) <- names(himmer_sp)

### Seurat standard preprocessing
miller_mtx <- Read10X('data/miller2019_upper_Case/')
miller <- CreateSeuratObject(counts = miller_mtx, 
                             project = 'lcmv', 
                             assay = 'RNA', 
                             min.cells = 1, min.features = 200)
miller <- NormalizeData(miller)
miller <- FindVariableFeatures(miller, selection.method = 'vst', nfeatures = 3000)
miller <- ScaleData(miller, verbose = FALSE)
miller <- RunPCA(miller, npcs = 30, verbose = FALSE)
miller <- RunUMAP(miller, reduction = "pca", dims = 1:30)
miller <- FindNeighbors(miller, dims = 1:30, verbose = FALSE)
miller <- FindClusters(miller, resolution= 0.1, verbose = FALSE)

## Miller Cluster annotation
miller$orig.ident <- miller$seurat_clusters
miller_ann <- data.frame('cluster' = miller$orig.ident,
                         'cell_type' = plyr::mapvalues(x = miller$orig.ident,
                                                       from = as.factor(c(3,0,4,2,1)), 
                                                       to = c('Proliferating', 'Effector', 'Naive',
                                                              'Progenitor Ex', 'Terminally Ex')
                         )
)
miller <- AddMetaData(miller, metadata = miller_ann)

## Loading Miller annotations
miller$'cell_type' <- miller_ann$cell_type

## merging datasets
merged_seu <- himmer_sp_seu
merged_seu$'miller' <- miller

cat('Saving Miller seurat object and annotations\n')
write_rds(miller_ann, 'data/miller2019/miller_annotations.rds')
write_rds(miller, 'analysis/miller_seu.rds')

## checking no empty intersection of genes
length( intersect( rownames(himmer_sp$DW), rownames(miller_mtx)))

cat('Seurat standard preprocessing')
for (i in 1:length(merged_seu)) {
        merged_seu[[i]] <- NormalizeData(merged_seu[[i]], verbose = FALSE)
        merged_seu[[i]] <- FindVariableFeatures(merged_seu[[i]], selection.method = "vst", 
                                                nfeatures = 2000, verbose = FALSE)
        merged_seu[[i]] <- ScaleData(merged_seu[[i]], verbose = FALSE)
        merged_seu[[i]] <- RunPCA(merged_seu[[i]], npcs = 30, verbose = FALSE)
        merged_seu[[i]] <- RunUMAP(merged_seu[[i]], reduction = "pca", dims = 1:30)
        merged_seu[[i]] <- SCTransform(merged_seu[[i]], verbose = FALSE)
}

cat('Performing integration\n')
features <- SelectIntegrationFeatures(object.list = merged_seu, nfeatures = 3000)
merged_seu <- PrepSCTIntegration(object.list = merged_seu, anchor.features = features)
reference_dataset <- which(names(merged_seu) == "miller")
anchors <- FindIntegrationAnchors(object.list = merged_seu, normalization.method = "SCT", 
                                       anchor.features = features, reference = reference_dataset)
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", k.weight = 50)

integrated <- RunPCA(object = integrated, verbose = FALSE)
integrated <- RunUMAP(object = integrated, dims = 1:30)
integrated <- FindNeighbors(integrated, dims = 1:15, verbose = FALSE)
integrated <- FindClusters(integrated, resolution = 0.1, verbose = FALSE, weights = 100, node.sizes = 10)

integrated$cell_type <- sapply(integrated$cell_type,
                               function(x) ifelse(is.na(x), 'Maike', x) )
cat('Finishing integration\n')

cat('Dataset annotation\n')
integrated$'dataset' <- sapply(integrated$orig.ident, function(x) 
                               ifelse(x=='lcmv', 'Miller', 'Maike'))

cat('Cell type prediction\n')
integrated$'predicted' <- plyr::mapvalues(as.character(integrated$seurat_clusters), 
                                          from = c('0', '1', '2', '3'),
                                          to = c('Effector', 'Terminally Ex',
                                                 'Progenitor Ex', 'Proliferating'))
write_rds(integrated, 'analysis/integrated_miller_maike_seu.rds')

## Subsetting integrated data to Maike imputed dataset  
integrated_sub <- subset(integrated, dataset != 'Miller')
write_rds(integrated_sub, 'analysis/integrated_maike_imputed_seu.rds')
