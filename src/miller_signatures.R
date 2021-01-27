## Visualization of scarches results

## Dependencies
library(Seurat)
library(dplyr)
library(readxl)
library(rio)
library(dittoSeq)
library(viridis)

set.seed(333)

## Setting working directory
path2project <- '/Users/carlosramirez/sc/tcd8ExscSeq/'
setwd(path2project)

#####################################################################
## Loading Hofmann data

## Data preprocessing
## gene names must be standardized
himmer <- read.csv('data/maike2020/nina_thimme_raw_counts.csv', 
                   header = TRUE)
r_names <- himmer$X
himmer <- select(himmer, -X)
himmer_mtx <- apply(himmer, 2, as.numeric)
rownames(himmer_mtx) <- gsub('__chr.*', '', r_names) 

## Seurat seurat construction
seurat <- CreateSeuratObject(counts = himmer_mtx, 
                                 project = 'himmer', 
                                 min.cells = 3, 
                                 min.features = 3,
                                 assay = 'RNA')
remove(himmer, himmer_mtx)

#####################################################################
## Adding Nina annotations

## Reading Nina annotations
nina.anns <- read_xlsx('data/cell_annotation_nina.xlsx') %>%
        as.data.frame()
rownames(nina.anns) <- nina.anns$CELLID
head(nina.anns)

## Annotating seurat file
seurat$'nina_annotations' <- 'NA'
seurat$nina_annotations[rownames(nina.anns)] <- nina.anns$TYPE
seurat$nina_annotations %>% table()

seurat <- NormalizeData(seurat) %>%
                ScaleData()
seurat <- FindVariableFeatures(seurat, 
                               nfeatures = 3000)
seurat <- RunPCA(seurat, verbose = FALSE)
seurat <- RunUMAP(seurat, dims = 1:30)

## Corrected Hofmann data
pdf('figures/hcv_umap_before_correction.pdf')
dittoDimPlot(subset(seurat, 
                    nina_annotations != 'NA'), 
             "nina_annotations", size = 2.5)
dittoDimPlot(subset(seurat, 
                    nina_annotations != 'NA'), 
             "orig.ident", size = 2.5)
dev.off()

#####################################################################
## Correcting by sample
## Using all genes
seurat.list <- lapply(unique(seurat$orig.ident), 
                      function(x) subset(seurat, orig.ident == x))
names(seurat.list) <- unique(seurat$orig.ident)
seurat.list <- lapply(seurat.list, 
                      function(x) SCTransform(x))

seurat.features <- SelectIntegrationFeatures(object.list = seurat.list, 
                                               nfeatures = 3000)
seurat.list <- PrepSCTIntegration(object.list = seurat.list, 
                                  anchor.features = seurat.features, 
                                  verbose = FALSE)
seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, 
                                           normalization.method = "SCT", 
                                           anchor.features = seurat.features, 
                                           verbose = FALSE, 
                                         reduction = 'cca', 
                                         k.anchor = 30, 
                                         k.score = 30, 
                                         k.filter = 10)
seurat.integrated <- IntegrateData(anchorset = seurat.anchors, 
                                   normalization.method = "SCT", 
                                   verbose = FALSE)

seurat.integrated <- FindVariableFeatures(seurat.integrated, 
                                          nfeatures = 3000)
seurat.integrated <- RunPCA(seurat.integrated, 
                     verbose = FALSE)
seurat.integrated <- RunUMAP(seurat.integrated, 
                      dims = 1:30)

## Corrected Hofmann data
pdf('figures/hofmann_corrected.pdf')
DimPlot(subset(seurat.integrated, 
               nina_annotations != 'NA') , 
        group.by = c('orig.ident', 
                     'nina_annotations'), 
        pt.size = 2.3)
dev.off()

#####################################################################
## Scoring using Progenitor 

## Reading Progenitor-like signature in Miller
miller_deg_path <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41590-019-0312-6/MediaObjects/41590_2019_312_MOESM4_ESM.xlsx'

## LCMV TEx Progenitor vs Terminally DEG
lcmv_deg <- import(miller_deg_path, which = 1)
str(lcmv_deg) 

## TIL TEx Progenitor vs Terminally DEG
til_deg <- import(miller_deg_path, which = 2)

signature.lcmv <- lcmv_deg %>%
        mutate(padj = as.numeric(padj)) %>%
        filter(padj < 0.05) %>%
        arrange(desc(log2FoldChange)) %>%
        head(100) %>%
        select(Gene) %>%
        unlist() %>%
        toupper() %>%
        unname()

seurat.integrated <- AddModuleScore(
        seurat.integrated,
        features = list(lcmv=signature.lcmv),
        nbin = 5, 
        name = 'Progenitor-like'
)
pdf('figures/umap_miller_prog-like_signature.pdf')
FeaturePlot(subset(seurat.integrated, 
                   nina_annotations != 'NA'), 
            reduction = 'umap',
            features = c('Progenitor.like1'), 
            pt.size = 2.5) +
                scale_color_viridis()
dev.off()

## Plot nina annotations
pdf('figures/vlnplot_miller_prog-like_signature.pdf')
seurat.integrated@meta.data$nina_annotations <- factor(seurat.integrated@meta.data$nina_annotations,
                                                       levels = c('mem', 'trans', 'exh'))
dittoPlot(subset(seurat.integrated, 
                 nina_annotations != 'NA'), 
          "Progenitor.like1", 
          group.by = "nina_annotations", 
          plots = c("jitter", "vlnplot", "boxplot"),
          jitter.color = "blue", jitter.size = 0.5,
          boxplot.color = "white", boxplot.width = 0.1,
          boxplot.fill = FALSE) +
        ggtitle('HCV-infected patients samples') +
        ylab('Miller Progenitor-like Signature Score') +
        xlab('')
dev.off()

pdf('figures/umap_corrected_nina_annotations.pdf')
dittoDimPlot(subset(seurat.integrated, 
                    nina_annotations != 'NA'), 
             "nina_annotations", size = 2.5)
dittoDimPlot(subset(seurat.integrated, 
                    nina_annotations != 'NA'), 
             "orig.ident", size = 2.5)
dittoDimPlot(subset(seurat.integrated, 
                    nina_annotations %in% c('exh', 'mem')), 
             "nina_annotations", 
             size = 2.5,
             do.ellipse = TRUE, 
             do.label = TRUE)
dev.off()

##############################################################
## Trying to improve clustering using DEG for dimensional
## reduction

seurat.integrated$nina_annotations %>% table()
seurat.filtered <- subset(seurat.integrated,
                          nina_annotations != 'NA')
Idents(seurat.filtered) <- seurat.filtered$nina_annotations
degs <- FindAllMarkers(seurat.filtered,
                                  only.pos = TRUE)
degs <- arrange(degs, desc(avg_logFC))
dim(degs)

seurat.filtered <- RunPCA(seurat.filtered, 
                            verbose = FALSE, 
                          features = degs$gene, 
                          npcs = 20)
seurat.filtered <- RunUMAP(seurat.filtered, 
                             dims = 1:20)

pdf('figures/improved_umap_only_degs.pdf')
DimPlot(seurat.filtered, 
        group.by = 'nina_annotations', reduction = 'pca')
DimPlot(subset(seurat.filtered, nina_annotations != 'trans'), 
        group.by = 'nina_annotations', reduction = 'umap')
dittoDimPlot(subset(seurat.filtered, 
                    nina_annotations %in% c('exh', 'mem')), 
             "nina_annotations", 
             size = 2.5,
             do.ellipse = TRUE, 
             do.label = TRUE)
FeaturePlot(subset(seurat.filtered, 
                   nina_annotations != 'trans'), 
            reduction = 'umap',
            features = c('Progenitor.like1'), 
            pt.size = 3) +
        scale_color_viridis()
dev.off()

#################################################################
## Calculating DEGs between Memory-like and Terminally exhausted
## cells
Idents(seurat.integrated) <- seurat.integrated$nina_annotations
degs <- FindMarkers(
        subset(seurat.integrated, nina_annotations != 'trans'), 
        ident.1 = 'mem', 
        logfc.threshold = 0, 
        min.pct = 0
)
degs<- degs %>% 
        filter(p_val_adj < 0.05) 
dim(degs)
head(degs)
saveRDS(degs, 
        'analysis/mem_vs_exh_hofmann_degs.rds',
        compress = TRUE)


