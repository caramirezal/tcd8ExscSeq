## Analysis of TF activities in CD8+ T cells

## Dependencies
library(Seurat)
library(dplyr)
library(randomForest)
library(reshape2)
library(ggplot2)
library(viridis)
library(patchwork)

path2project  <- '/Users/carlosramirez/sc/tcd8ExscSeq/'
setwd(path2project)

## Loading seurat object with integrated data Miller-Hoffmann
seurat  <- readRDS('analysis/integrated_miller_hoffmann.rds')

## Loading TF activities from SCENIC
scenic <- read.csv('analysis/scenic/scenic_hgv_6000/aucell.csv')

## Merge Adding scenic results to metadata in seurat file
hoffmann <- subset(seurat, cells = scenic$Cell)
any( ! rownames(hoffmann@meta.data) == scenic$Cell )
hoffmann@meta.data <- cbind(hoffmann@meta.data, 
                         scenic)

## Which TFs are differentially active in progenitor vs
## exhausted CD8 T cells?
## Running random forest to assess variable importance

## Extracting dataframe
dframe <- hoffmann %>% 
             subset(predicted %in% c('Progenitor Ex',
                                     'Terminally Ex')) %>%
                FetchData(vars = c(grep('\\.\\.\\.', 
                                        colnames(hoffmann@meta.data),
                                        value = TRUE),
                                   'predicted'))
dim(dframe)                
head(dframe[, 1:5])        

## Random forest
rf <- randomForest(as.factor(predicted)~., data = dframe)
varimp <- rf$importance %>%
                as.data.frame() %>%
                add_rownames(var = 'transcription_factor') %>%
                arrange(desc(MeanDecreaseGini))

pdf('figures/rf_var_imp_tfa_precursors_vs_terminally_ex.pdf')
varImpPlot(rf)
dev.off()

pdf('figures/vln_top_varimp_tfs_precursors_vs_terminally_ex.pdf')
VlnPlot(subset(hoffmann, predicted %in% c('Progenitor Ex',
                                          'Terminally Ex')),
        features = varimp$transcription_factor[1:6],
        group.by = 'predicted')
dev.off()

## Bar plot visualisation of TF activity in Progenitor vs Terminally
## Exhausted
n <- 20
tf_means <- dframe %>%
        filter(predicted == 'Progenitor Ex') %>%
        select(-predicted) %>%
        apply(2, mean) 
tfs_sorted <- tf_means[varimp$transcription_factor[1:n]] %>%
                sort(decreasing = FALSE) %>% names()
tfs_sorted <- gsub('\\.', '', tfs_sorted)
pdf('figures/boxplot_top_varimp_tfs_precursors_vs_terminally_ex.pdf')
dframe %>% 
        melt() %>%
        filter(variable %in% varimp$transcription_factor[1:n]) %>%
        mutate(variable=gsub('\\.\\.\\.', '', variable)) %>%
        mutate(variable=factor(variable, 
                               levels = tfs_sorted)) %>%
        ggplot(aes(x=variable, y=log(value+1),
                   fill=predicted)) +
                geom_boxplot() + 
                coord_flip() +
                theme_bw() +
                theme(panel.grid = element_blank(),
                      legend.title = element_blank())
dev.off()

################################################################
## Visualize TF activity in clusters

tfs <- c('JUNB...',
         'FOSB...',
         'MYC...',
         'EGR1...',
         'PAX5...',
         'ZNF274...',
         'KLF2...')


tfa <- FeaturePlot(subset(hoffmann,
                   predicted %in% c('Progenitor Ex',
                                    'Terminally Ex')),
            features = , 
            pt.size = 2.2) 
plots <- list()
for (i in 1:length(tfs)) {
      plots[[i]] <- FeaturePlot(subset(hoffmann,
                                     predicted %in% c('Progenitor Ex',
                                                      'Terminally Ex')),
                              features = tfs[i], 
                              pt.size = 2.2) + 
                        scale_color_viridis()
}

pdf('figures/tfa_umap_top_varimp.pdf')
wrap_plots(plots[1:4], 
           guides = 'collect', 
           ncol = 2)
dev.off()

