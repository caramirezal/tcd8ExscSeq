## Difussion map pseudotime analysis
library(Seurat)
library(destiny)
library(dplyr)
library(Biobase)
library(ggplot2)
library(slingshot)
library(scater)
library(viridis)
library(randomForest)

set.seed(333)

path2project  <- '/Users/carlosramirez/sc/tcd8ExscSeq/'
setwd(path2project)

## Loading seurat object with integrated data Miller-Hoffmann
seurat  <- readRDS('data/maike2020/hofmann_hcv_seu.rds')

## Loading TF activities from SCENIC
scenic <- read.csv('analysis/scenic/scenic_hgv_6000/aucell.csv')

## Merge Adding scenic results to metadata in seurat file
hofmann <- subset(seurat, cells = scenic$Cell)
any( ! rownames(hofmann@meta.data) == scenic$Cell )
hofmann@meta.data <- cbind(hofmann@meta.data, 
                            scenic)

###################################################################
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

## Saving random forest result
rforest$importance %>% 
        as.data.frame() %>%
        arrange(desc(MeanDecreaseGini)) %>%
        add_rownames(var = 'gene') %>%
        mutate(gene=gsub('\\.', '', gene)) %>%
        saveRDS('analysis/tfa_rforest_ranked.rds', 
                compress = TRUE)

## Visualization

## PRDM1
p1 <- FeaturePlot(subset(hofmann, 
                         nina_annotations %in% c('mem', 'exh')),
            features = top_n_tfs[1], 
            pt.size = 2) + scale_color_viridis()
## STAT1
p2 <- FeaturePlot(subset(hofmann, 
                         nina_annotations %in% c('mem', 'exh')),
            features = top_n_tfs[2], 
            pt.size = 2) + scale_color_viridis()
## STAT1
p3 <- FeaturePlot(subset(hofmann, 
                         nina_annotations %in% c('mem', 'exh')),
            features = top_n_tfs[3], 
            pt.size = 2) + scale_color_viridis()
## STAT1
p4 <- FeaturePlot(subset(hofmann, 
                         nina_annotations %in% c('mem', 'exh')),
            features = top_n_tfs[4], 
            pt.size = 2) + scale_color_viridis()
gridExtra::grid.arrange(grobs=list(p1, p2, p3, p4), 
                        ncol=2)

#####################################################################
## Selection of the genes to perform the diffussion map
## DEG between Memory-like vs Exhausted cells 

degs <- read.table('analysis/degs_mem_vs_exh.tsv.gzip')
dim(degs)

####################################################################

#hofmann.exh <- subset(hofmann, nina_annotations %in% c('mem', 
#                                                       'exh', 'trans'))
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
#dev.off()




