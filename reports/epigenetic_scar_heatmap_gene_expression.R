## This scripts contains the analysis of CD8 T cells scRNA-Seq data from
## Nina Thimme et al, 2020 from the Maike Hofmann group:
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150305


## Dependencies
library(Seurat)
library(tidyverse)
library(dplyr)
library(reshape2)
library(viridis)
library(stringr)
library(reshape2)
library(ggpubr)
library(heatmaply)
library(RColorBrewer)
library(ggrepel)
library(scales)

## Path to working directory
path2project <- '/media/ag-cherrmann/cramirez/tcd8ExscSeq/'

## Path to mapped regions to gene annotations and TF motifs
## This file was generated using motifMatchr and GRanges
## as described in the network_inference.Rmd file
ranges.ann.path <- paste0(path2project,
                          '/analysis/motifs_CHaRs_preVsPost_treatment.rds')

## Path to the peaks annotated by scarred, gained and reversed
## categories
ranges.categories.path <- paste0(path2project,
                            '/analysis/regions_categories_epigenetic_scar.tsv.gz')

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
#download.file(
#  url = 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE150305&format=file',
#  destfile = tar.path
#)
## Untar files
#untar(
#  tarfile = tar.path,
#  exdir = paste0(temp_folder, '/GSE150305_RAW')
#)
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
seurat <- subset(seurat, nCount_RNA > 800 & perc.KCNQ1OT1 < 2 )


seurat <- FindVariableFeatures(seurat, nfeatures = 3000)
seurat <- SCTransform(seurat)
seurat <- RunPCA(seurat, npcs = 50)
seurat <- RunUMAP(seurat, dims = 1:20)

## Plotting UMAPs by different categories
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


######################################################################################
##                                                                                  ##
##                    Filtering UMAP outliers.                                      ##
##                                                                                  ##
######################################################################################

seurat <- FindNeighbors(seurat, dims = 1:20)
seurat <- FindClusters(seurat, resolution = 0.1)

## before filtering
umapSeuratClusters <- DimPlot(seurat,
                         group.by = 'seurat_clusters',
                         label = TRUE,
                         repel = TRUE) +
  NoLegend() +
  labs(title='Seurat clusters')
umapSeuratClusters

seurat.filtered <- subset(seurat, seurat_clusters != '2')

## after filtering
patient.filtered.umap <- DimPlot(seurat.filtered,
                              group.by = 'patient',
                              label = TRUE,
                              repel = TRUE) +
  NoLegend() +
  labs(title='Seurat clusters')
patient.filtered.umap


######################################################################################
##                                                                                  ##
##                    Analysis of only treated patients                             ##
##                                                                                  ##
######################################################################################

seurat.treatments <- subset(seurat.filtered, patient %in% c('Patient3_cHCV', 'Patient3_cured_HCV',
                                                              'Patient6_cHCV', 'Patient6_cured_HCV'))
seurat.treatments <- RunPCA(seurat.treatments, npcs = 50)
seurat.treatments <- RunUMAP(seurat.treatments, dims = 1:20, min.dist = 0.01)

## Plotting UMAPs by different categories
umapByPatient.treatments <- DimPlot(seurat.treatments,
                         group.by = 'patient',
                         repel = TRUE) +
  labs(title='Patient')
umapByCondition.treatments <- DimPlot(seurat.treatments, group.by = 'condition',
                           repel = TRUE) +
  theme_void() +
  ggtitle('Condition - Treated patients')
umapByPatient.treatments +
  umapByCondition.treatments


######################################################################################
##                                                                                  ##
##                 Getting mapped genes to regions                                  ##
##                                                                                  ##
######################################################################################

## Regions mapped to gene annotations and motifs
mapped.regions <- readRDS(ranges.ann.path)

######################################################################################
##                                                                                  ##
##                 Categorizing regions                                            ##
##                                                                                  ##
######################################################################################
## Next we categorize genes by scarred, reversed and gained after DAA treatment

## The data was taken from Kathleen B. Y. et al, 2021
## https://www.nature.com/articles/s41590-021-00979-1?proof=t#additional-information
peaks_url <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41590-021-00979-1/MediaObjects/41590_2021_979_MOESM3_ESM.xlsx'

download.file(peaks_url, '/media/ag-cherrmann/cramirez/tcd8ExscSeq/data/41590_2021_979_MOESM3_ESM.xlsx')
peaks.df <- readxl::read_xlsx('/media/ag-cherrmann/cramirez/tcd8ExscSeq/data/41590_2021_979_MOESM3_ESM.xlsx',
                              sheet = 'Supplementary Table 2')
dim(peaks.df)
head(peaks.df)

peaks.df <- peaks.df %>%
  mutate(region_reg_type=ifelse(`log2FC: Group HCV_pre vs HCV_post`> 0,
                                'positive', 'negative')) %>%
  mutate(region_deg=(as.numeric(`q-value: Group HCV_pre vs HCV_post`)) < 0.05 )
is.na(peaks.df$region_deg) %>% sum
peaks.df$'category' <- sapply(1:nrow(peaks.df),
                              function(i){
                                if ( is.na( peaks.df$region_deg[i] ) ) {
                                  return('NA')
                                } else {
                                  if ( peaks.df$region_deg[i] == FALSE ) {
                                    return('Scarred')
                                  }

                                  if ( peaks.df$region_deg[i] == TRUE) {
                                    if ( peaks.df$region_reg_type[i] == 'positive') {
                                      return('Reversed')
                                    }
                                    if ( peaks.df$region_reg_type[i] == 'negative') {
                                      return('Gained')
                                    }
                                  }
                                }

                              })

print('')
## Categorization of regions by the change in epigenetic landscape
dpeaks.categories <- peaks.df %>%
  ggplot(aes(x=`P1_HCV_pre`,
             y=`P1_HCV_post`,
             colour=category)) +
  geom_point() +
  labs(x='HCV - DAA Pre-treatment Peak Score',
       y='HCV - DAA Post-treatment Peak Score') +
  theme_classic()
#dpeaks.categories


## Checking order
any(mapped.regions.df$start != peaks.df$Start & mapped.regions.df$end != peaks.df$End)
mapped.regions.df <- as.data.frame(mapped.regions@rowRanges)
mapped.regions.df$'categories' <- peaks.df$category


######################################################################################
##                                                                                  ##
##                   Analysis of motifs in CHaRs                                    ##
##                                                                                  ##
######################################################################################

freqs.tbl <- mapped.regions.df %>%
  pull(TF) %>%
  str_split(pattern = ';') %>%
  unlist() %>%
  toupper() %>%
  table()
freqs.df <- data.frame(TF = names(freqs.tbl),
                       frequency=freqs.tbl)
freqs.df <- freqs.df %>%
  arrange(desc(`frequency.Freq`)) %>%
  mutate(rank=n():1) %>%
  mutate(gene_label=ifelse((n()-30)<rank, TF, ''))
freqs.df %>%
  ggplot(aes(x=rank, y=`frequency.Freq`,
             label=gene_label)) +
  geom_point() +
  xlim(-10, 2000) +
  geom_text_repel(max.overlaps = 100,
                  force = 20,
                  colour='Red') +
  theme_bw() +
  labs(x='Rank', y='Frequency')



######################################################################################
##                                                                                  ##
##                      Gene expression by categories                               ##
##                                                                                  ##
######################################################################################

## Splitting genes by categories
categories.df <- filter(mapped.regions.df, genes != '')
categories.df <- mutate(categories.df, gene_chr=paste0(genes, '--', seqnames))
genesByCategory <- split(categories.df$gene_chr, categories.df$categories)
genesByCategory <- genesByCategory[c('Gained', 'Reversed', 'Scarred')]
genesByCategory <- lapply(genesByCategory, unique)
genesByCategory <- lapply(genesByCategory,
                          function(x) intersect(x, rownames(seurat.filtered)))
genesByCategory$'Random genes' <- sample(rownames(seurat), size = 735)
seurat.patient3 <- subset(seurat.filtered, patient %in% c('Patient3_cHCV', 'Patient3_cured_HCV'))

mtxByCategory <- lapply(genesByCategory,
                        function(genes) FetchData(seurat.patient3,
                                                  vars = c(genes,
                                                           'condition')))
## Mean values
mtxByCategory <- lapply(mtxByCategory, function(mtx){
  df <- group_by(mtx, condition)
  summ.df <- summarise_all(df, mean)
  summ.df
})
## Plotting
boxplots.list <- lapply(names(mtxByCategory), function(category){
  mtxByCategory[category][[1]] %>%
    melt %>%
    ggplot(aes(x=condition, y=log(value))) +
    geom_jitter() +
    geom_boxplot() +
    ggtitle(category) +
    stat_compare_means(label.y = 1.5, label = "p.signif",
                       colour='red', size=5) +
    theme_bw() +
    theme(axis.text.x = element_text(hjust = 1, angle=45)) +
    labs(x='', y='Log(mean(gene expression))')
})
boxplot.pat3 <- gridExtra::grid.arrange(grobs=boxplots.list, ncol=4)
plot(boxplot.pat3)

####################################################################

heatmap.mtx <- lapply(mtxByCategory,
       function(mtx){
         df <- select(mtx, -condition)
         mean.vals <- apply(df, 1, mean)
         mean.vals
       }) %>%
       unlist() %>%
       matrix(ncol=2, byrow = TRUE)
rownames(heatmap.mtx) <- names(mtxByCategory)
colnames(heatmap.mtx) <- c('cHCV', 'Cured')

heatmap.ggplot <- heatmap.mtx %>%
  scale() %>%
  as.data.frame() %>%
  add_rownames(var = 'condition') %>%
  melt() %>%
  ggplot(aes(x=variable,
             y=factor(condition,
                      levels=c('Random genes',
                               'Scarred',
                               'Gained',
                               'Reversed')),
             fill=value)) +
         geom_tile(width=1,
                   height=1,
                   size=1,
                   colour='white') +
         scale_fill_gradient2(low = 'chartreuse4',
                              mid = 'Yellow',
                              high = 'Red') +
  coord_equal() +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank()) +
  labs(x='', y='')
heatmap.ggplot

heatmap.plot <-pheatmap::pheatmap(
  heatmap.mtx,
  color = rev(brewer.pal(10, 'RdYlGn')),
  cellwidth = 50,
  cellheight = 50
)

####################################################################

## A bias for low gene expression can be observed in the patient 6
## in the cells of the cured patients
seurat.filtered <- subset(seurat, patient %in% c('Patient6_cHCV',
                                                 'Patient6_cured_HCV'))
mtxByCategory <- lapply(genesByCategory,
                        function(genes) FetchData(seurat.filtered,
                                                  vars = c(genes,
                                                           'condition')))
## Mean values
mtxByCategory <- lapply(mtxByCategory, function(mtx){
  df <- group_by(mtx, condition)
  summ.df <- summarise_all(df, mean)
  summ.df
})
## Plotting
boxplots.list <- lapply(names(mtxByCategory), function(category){
  mtxByCategory[category][[1]] %>%
    melt %>%
    ggplot(aes(x=condition, y=log(value))) +
    geom_jitter() +
    geom_boxplot() +
    ggtitle(category) +
    stat_compare_means(label.y = 1.5) +
    theme_bw() +
    labs(x='', y='Log(mean(gene expression))')
})
boxplot.pat6 <- gridExtra::grid.arrange(grobs=boxplots.list, ncol=4)
boxplot.pat6





######################################################################################
##                                                                                  ##
##                          Plotting                                                ##
##                                                                                  ##
######################################################################################

p <-  ( umapByPatient + umapByCondition + umap.counts ) /
           ( dpeaks.categories + boxplot.pat3 + heatmap.ggplot )
ggsave(paste0(path2project, '/figures/reprocessing_Thimme_2020.pdf'),
       plot = p, width = 16, height = 8)
