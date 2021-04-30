## This script contains the code to reproduce the figures

## Dependencies
library(Seurat)
library(dplyr)
library(readxl)
library(rio)
library(dittoSeq)
library(viridis)
library(randomForest)
library(destiny)
library(Biobase)

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
seurat.integrated$nina_annotations <- factor(seurat$nina_annotations,
                                             c('mem', 'trans', 'exh'))
dittoDimPlot(subset(seurat.integrated, 
                    nina_annotations != 'NA'), 
             "nina_annotations", size = 2.5)
dittoDimPlot(subset(seurat.integrated, 
                    nina_annotations != 'NA'), 
             "orig.ident", size = 2.5)
dev.off()

#####################################################################
## Scoring Progenitor-like signature in Hofmann data 

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

########################################################################
## Loading TF activities from SCENIC
scenic <- read.csv('analysis/scenic/scenic_hgv_6000/aucell.csv')

## Merge Adding scenic results to metadata in seurat file
hofmann <- subset(seurat.integrated, cells = scenic$Cell)
any( ! rownames(hofmann@meta.data) == scenic$Cell )
hofmann@meta.data <- cbind(hofmann@meta.data, 
                           scenic)

## Selection of the TFs using a random forest predicting
## memory vs terminally phenotypes
n <- 4
mtx <- hofmann@meta.data %>%
        select(nina_annotations, `ALX4...`:`ZSCAN26...`) %>%
        filter(nina_annotations %in% c('mem', 'exh')) %>%
        mutate(nina_annotations= as.character(nina_annotations))
rforest <- randomForest(as.factor(nina_annotations)~., 
                        data = mtx) 
varImpPlot(rforest)
tfs <- rforest$importance %>% 
        as.data.frame() %>%
        arrange(desc(MeanDecreaseGini)) %>%
        rownames()
top_n_tfs <- tfs[1:4]
top_n_tfs


## Selection of the genes to perform the diffussion map
## DEG between Memory-like vs Exhausted cells 
degs <- read.table('analysis/degs_mem_vs_exh.tsv.gzip')
dim(degs)

####################################################################
## Diffusion maps

## Using only memory and exhausted cells
hofmann.exh <- subset(hofmann, nina_annotations %in% c('mem', 
                                                       'exh'))

## Creation of the eset object 
eset <- ExpressionSet(
        assayData = hofmann.exh@assays$integrated@scale.data,
        phenoData = AnnotatedDataFrame(select(hofmann.exh@meta.data,
                                              nina_annotations, 
                                              Progenitor.like1, tfs[1:10]))
)
eset


## Difussion map
difmap <- DiffusionMap(eset, 
                       vars = degs$gene,
                       n_pcs = 50, 
                       n_eigs = 100, 
                       k=50)

res <- cbind(difmap@data_env$data@phenoData@data, 
             DC1=difmap$DC1,
             DC2=difmap$DC2,
             DC3=difmap$DC3,
             DC4=difmap$DC4) 

my_theme <- theme_bw() +
        theme(panel.grid = element_blank()) 
## Plotting
#pdf('figures/diffusion_map_TF_activities_exhausted.pdf',
#    width = 15, height = 10)
categories <- res %>% ggplot(aes(x=DC1, y=DC2,
                                 colour=nina_annotations)) +
        geom_point(size=1.5) + my_theme 
prdm1 <- res %>% ggplot(aes(x=DC1, y=DC2,
                            colour=`PRDM1...`)) +
        geom_point(size=2) + scale_color_viridis()  + my_theme
stat1 <- res %>% ggplot(aes(x=DC1, y=DC2,
                            colour=`STAT1...`)) +
        geom_point(size=2) + scale_color_viridis()  + my_theme
fos <- res %>% ggplot(aes(x=DC1, y=DC2,
                          colour=`FOS...`)) +
        geom_point(size=2) + scale_color_viridis()  + my_theme
junb <- res %>% ggplot(aes(x=DC1, y=DC2,
                           colour=`JUNB...`)) +
        geom_point(size=2) + scale_color_viridis()  + my_theme
prog.signature <- res %>% ggplot(aes(x=DC1, y=DC2,
                                     colour=`Progenitor.like1`)) +
        geom_point(size=2) + scale_color_viridis()  + my_theme
gridExtra::grid.arrange(categories,
                        prdm1,
                        stat1,
                        fos,
                        junb,
                        prog.signature,
                        ncol=3)
