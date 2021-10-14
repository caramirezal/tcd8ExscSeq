## This scripts contains the analysis of CD8 T cells scRNA-Seq data from
## Nina Thimme et al, 2020 from the Maike Hofmann group:
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150305


## Dependencies
library(Seurat)
library(tidyverse)
library(dplyr)
library(reshape2)
library(viridis)

## Path to working directory
path2project <- '/media/ag-cherrmann/cramirez/tcd8ExscSeq/'



#####################################################################################
##                                                                                 ##
##            Processing scRNA-Seq data containing pre- and post-DAA               ##
##            conditions from Thimme N et al, 2020                                 ##
##                                                                                 ##
#####################################################################################

## Creating temporal folder to store count matrices
dir.create(path = paste0(path2project, 'analysis'))
temp_folder <- paste0(path2project, 'analysis/thimme2020_tmp')
dir.create(path = temp_folder)
tar.path <- paste0(temp_folder, '/GSE150305_RAW.tar')

## Downloading data
download.file(
  url = 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE150305&format=file',
  destfile = tar.path
)
## Untar files
untar(
  tarfile = tar.path,
  exdir = paste0(temp_folder, '/GSE150305_RAW')
)
## Checking downloading
csv.files <- list.files(paste0(temp_folder, '/GSE150305_RAW'),
                        full.names = TRUE)

## Loading files
counts.mtx.list <- lapply(
  csv.files,
  function(file) {
    counts <- read_tsv(file)
    rownames(counts) <- counts$GENEID
    counts <- as.data.frame(counts)
    counts <- dplyr::select(counts, -GENEID)
    return(counts)
})
names(counts.mtx.list) <- gsub('.*\\/|\\..*', '', csv.files)
lapply(counts.mtx.list, dim)

## Reading into a seurat objects
seurat.list <- lapply(names(counts.mtx.list),
       function(sample){
         counts <- counts.mtx.list[sample][[1]]
         colnames(counts) <- paste0(sample, '_', colnames(counts))
         CreateSeuratObject(
           counts = counts,
           project = sample,
           assay = 'RNA',
           min.cells = 1,
           min.features = 1
         )
})
names(seurat.list) <- names(counts.mtx.list)

## Concatenating seurat objects
seurat <- merge(
  x=seurat.list[1][[1]],
  y = seurat.list[2:length(seurat.list)]
)

## Annotation
seurat$'patient' <- gsub('.*_Pat', 'Patient', colnames(seurat))
seurat$patient <- gsub('HCV_.*', 'HCV', seurat$patient)
seurat$'condition' <- gsub('Patient._', '', seurat$patient)
table(seurat$condition)

#####################################################################################
##                                                                                 ##
##.                     Standard seurat pre-processing                             ##
##                                                                                 ##
#####################################################################################

grep('^MT-', rownames(seurat), value = TRUE) %>% head

## QC METRICS
seurat$'perc.mito' <- PercentageFeatureSet(seurat, pattern = '^MT-')
seurat$'perc.KCNQ1OT1' <- PercentageFeatureSet(seurat, pattern = 'KCNQ1OT1')

## Pre-filtering
seurat <- subset(seurat, nCount_RNA > 800 & perc.KCNQ1OT1 < 2)


seurat <- FindVariableFeatures(seurat, nfeatures = 3000)
seurat <- SCTransform(seurat)
seurat <- RunPCA(seurat, npcs = 50)
seurat <- RunUMAP(seurat, dims = 1:20)

umapByPatient <- DimPlot(seurat,
                        group.by = 'patient',
                        label = TRUE,
                        repel = TRUE) +
                    NoLegend() +
                    labs(title='Patient')
umapByCondition <- DimPlot(seurat, group.by = 'condition',
                           label = TRUE,
                           repel = TRUE) +
                       theme_void() +
                      theme(legend.position = 'none') +
                      ggtitle('Condition')
umap.counts <- FeaturePlot(seurat, features = 'nCount_RNA') +
                  scale_color_viridis() +
                  theme_void() +
                  ggtitle('Total UMIs')

p <- umapBySample + umapByCondition + umap.counts
ggsave(paste0(path2project, '/figures/reprocessing_Thimme_2020.pdf'),
       plot = p, width = 16, height = 6)
