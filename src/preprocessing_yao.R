## Preprocessing Yao et al, 2019

library(dplyr)
library(GEOquery)
library(Seurat)
library(ggplot2)

set.seed(333)

path2project <- '/Users/carlosramirez/sc/tcd8ExscSeq/'
setwd(path2project)

## Downloading data
url_data <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE119940&format=file'
dir.create('data/yao2019')
dir.create('data/yao2019/tmp')
download.file(url = url_data, 
              destfile = 'data/yao2019/tmp/GSE119940_RAW.tar')
untar(tarfile = 'data/yao2019/tmp/GSE119940_RAW.tar', 
      exdir = 'data/yao2019/tmp/GSE119940_RAW')
file.remove('data/yao2019/tmp/GSE119940_RAW.tar')
files <- list.files('data/yao2019/tmp/GSE119940_RAW')
files

########################################################################
## Sample annotations
gse <- getGEO('GSE119940', GSEMatrix = TRUE)
anns <- pData(gse[[1]])
dim(anns)

#######################################################################
## Reading files into Seurat
dir.create('data/yao2019/tmp/10x')
seurat.list <- lapply(1:nrow(anns), 
       function(i){
          folder <- paste0('data/yao2019/tmp/10x/', anns$geo_accession[i])
          dir.create(folder)
          file.copy(from = paste0('data/yao2019/tmp/GSE119940_RAW/',
                                  gsub('.*suppl/', '', 
                                       anns$supplementary_file_1[i])),
                    to = paste0(folder, '/barcodes.tsv.gz'))
          file.copy(from = paste0('data/yao2019/tmp/GSE119940_RAW/',
                                  gsub('.*suppl/', '', 
                                       anns$supplementary_file_2[i])),
                    to = paste0(folder, '/features.tsv.gz'))
          file.copy(from = paste0('data/yao2019/tmp/GSE119940_RAW/',
                                  gsub('.*suppl/', '', 
                                       anns$supplementary_file_3[i])),
                    to = paste0(folder, '/matrix.mtx.gz'))
          Read10X(folder)
})
names(seurat.list) <- anns$title

## Creating seurat objects
seurat.list <- lapply(1:length(seurat.list), function(i){ 
        mtx <- seurat.list[i][[1]]
        colnames(mtx) <- paste0(names(seurat.list)[i], '_',
                                colnames(mtx))
        seurat <- CreateSeuratObject(
                counts = mtx,
                project = names(seurat.list)[i],
                assay = 'RNA',
                min.cells = 5,
                min.features = 10
        )
        seurat$'tag' <- anns$title[i]
        return(seurat)
})
names(seurat.list) <- anns$title

########################################################################
## Integration of samples
seurat.list <- lapply(seurat.list, SCTransform)
features <- SelectIntegrationFeatures(object.list = seurat.list)
seurat.list <- PrepSCTIntegration(object.list = seurat.list, 
                                anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = seurat.list, 
                                  normalization.method = "SCT", 
                                  anchor.features = features)
integrated <- IntegrateData(anchorset = anchors, 
                              normalization.method = "SCT")
rm(seurat.list)
rm(anchors)
unlink('')

integrated

## Dimension reduction
integrated <- RunPCA(integrated)
integrated <- RunUMAP(integrated, dims = 1:30)

saveRDS(integrated, 
        'data/yao2019/yao2019_seurat.rds',
        compress = TRUE) 



