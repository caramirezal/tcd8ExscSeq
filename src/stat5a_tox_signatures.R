## This script contains an analysis of TOX and STAT5a signatures in HCV specific
## comparing treated vs un treated cells

library(Seurat)
library(tidyverse)
library(rio)
library(ggrepel)
library(reshape2)
library(ggpmisc)
library(ggpubr)

path2seurat <- '/media/ag-cherrmann/cramirez/tcd8ExscSeq/data/maike2020/hofmann_hcv_seu.rds'
seurat <- readRDS(path2seurat)

DimPlot(seurat)

seurat[[]] %>% tail()


##----------------
## Loading annotations of the cells from cured patients
anns <- read_csv('/media/ag-cherrmann/cramirez/tcd8ExscSeq/data/annotations_cured_patients/annotation_after_treated_cells_hcv.csv')

colnames(anns) <- c('Cell', 'cell_type')
cell_types <- as.character(seurat$nina_annotations)
names(cell_types) <- colnames(seurat)
cell_types[anns$Cell] <- anns$cell_type
seurat$nina_annotations <- cell_types

## Adding information of the treatment
seurat$'condition' <- sapply(colnames(seurat), 
                    function(cell) ifelse(cell %in% anns$Cell, 'DAA_treated', 'Chronic'))



## Loading TOX signatures from Khan et al, 2019
tox.khan.ko.sign.url <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-019-1325-x/MediaObjects/41586_2019_1325_MOESM3_ESM.xlsx' 
tox.khan.ko.sign.df <- import(tox.khan.ko.sign.url, skip=2)
dim(tox.khan.ko.sign.df)
head(tox.khan.ko.sign.df)

n.top <- 10
tox.khan.ko.sign.df %>%
        mutate(rank=1:n())  %>%
        mutate(highlight=ifelse(rank<=n.top | 
                                        n() - n.top <= rank | 
                                        `Gene ID` == 'Tox',
                                TRUE, FALSE)) %>%
        mutate(gene_label=ifelse(highlight==TRUE,
                                 `Gene ID`, '')) %>%
        ggplot(aes(x=`Log2 Fold Change`, y=-log10(`Adjusted p-value`),
                   colour=highlight,
                   label=gene_label)) +
        geom_point() +
        geom_text_repel(max.overlaps = 100) +
        theme_classic() +
        labs(x='Log2 FC', y='-Log10(p-val adjusted)') +
        scale_color_manual(values = c('black', 'red')) +
        theme(legend.position = 'none')



## Loading stat5a signature
stat5a_sign <- read_lines('/media/ag-cherrmann/cramirez/tcd8ExscSeq/data/stat5_signature/stat5a_signature_kanai.txt', 
                          skip = 2)


human_mouse.mapping.url <- 'http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt'
mapping <- read.table(human_mouse.mapping.url,
                      sep = '\t')


## Mapping mouse to human genes
top.n <- 100
tox.sign.up <- head(tox.khan.ko.sign.df, n = top.n) %>% pull(`Gene ID`)
tox.sign.down <- tail(tox.khan.ko.sign.df, n = top.n) %>% pull(`Gene ID`)
tox.sign.list <- list(tox_sign_up=tox.sign.up, 
                      tox_sign_down=tox.sign.down)
sign.list <- lapply(tox.sign.list, function(tox_genes){
        genes <- filter(mapping, V3 %in% tox_genes) %>% pull(V1)
})
sign.list$'stat5a_sign' <- stat5a_sign



##---------------------------
## Scoring cells
seurat <- AddModuleScore(seurat, 
                         features = sign.list, 
                         name = names(sign.list), 
                         nbin = 10, 
                         assay = 'RNA')



#####################################################################
## Comparing Mem in chronic Vs cured cells
        
##--------------------------
## Visualizations
vln_plot <- seurat[[]] %>%
        filter(nina_annotations %in% c('mem')) %>%
        select(condition,
               nina_annotations,
               tox_sign_up1, 
               tox_sign_down2, 
               stat5a_sign3) %>%
        rename('TOX(-) Up Signature'=tox_sign_up1,
               'TOX(-) Down Signature'=tox_sign_down2,
               'STAT5a Signature'=stat5a_sign3) %>%
        melt() %>%
        ggplot(aes(x=condition, y=value)) +
                geom_violin(draw_quantiles = 0.5) +
                geom_jitter(alpha=0.3) +
                facet_wrap(~variable, scales = 'free') +
                theme_bw() +
                theme(legend.position = 'none',
                      panel.grid = element_blank()) +
                labs(x='', y='Signature Score')
        

scat_plot <- seurat[[]] %>%
        filter(nina_annotations %in% c('mem')) %>%
        select(condition,
               nina_annotations,
               tox_sign_up1, 
               tox_sign_down2, 
               stat5a_sign3) %>%
        ggplot(aes(x=stat5a_sign3, 
                   y=tox_sign_down2,
                   colour=condition)) +
                geom_point() +
                stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                       after_stat(rr.label), sep = "*\", \"*"))) +
                geom_smooth(method = 'lm') +
                theme_bw() +
                theme(panel.grid = element_blank()) +
                labs(x='STAT5a Signature', 
                     y='TOX(-) Down Signature',
                     colour='') 




pdf('/media/ag-cherrmann/cramirez/tcd8ExscSeq/figures/tox_vs_stat5a_mem.pdf',
    height = 5, width = 13)
vln_plot + scat_plot                
dev.off()


table(seurat$nina_annotations, seurat$condition)



#####################################################################
##  Comparing Memory Vs Intermediate

## Visualizations
vln_plot <- seurat[[]] %>%
        filter(nina_annotations %in% c('mem', 'trans')) %>%
        filter(condition == 'Chronic') %>%
        select(condition,
               nina_annotations,
               tox_sign_up1, 
               tox_sign_down2, 
               stat5a_sign3) %>%
        rename('TOX(-) Up Signature'=tox_sign_up1,
               'TOX(-) Down Signature'=tox_sign_down2,
               'STAT5a Signature'=stat5a_sign3) %>%
        melt() %>%
        ggplot(aes(x=nina_annotations, y=value)) +
        geom_violin(draw_quantiles = 0.5) +
        stat_compare_means(label.x = 1, label = 'p.format') +
        geom_jitter(alpha=0.3) +
        facet_wrap(~variable, scales = 'free') +
        theme_bw() +
        theme(legend.position = 'none',
              panel.grid = element_blank()) +
        labs(x='', y='Signature Score',
             title = 'Memory Vs Transitory',
             subtitle = 'Chronic')

scat_plot <- seurat[[]] %>%
        filter(nina_annotations %in% c('mem', 'trans')) %>%
        filter(condition == 'Chronic') %>%
        select(condition,
               nina_annotations,
               tox_sign_up1, 
               tox_sign_down2, 
               stat5a_sign3) %>%
        ggplot(aes(x=stat5a_sign3, 
                   y=tox_sign_down2,
                   colour=nina_annotations)) +
        geom_point() +
        stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                       after_stat(rr.label), sep = "*\", \"*"))) +
        geom_smooth(method = 'lm') +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        labs(x='STAT5a Signature', 
             y='TOX(-) Down Signature',
             colour='') 

vln_plot + scat_plot


pdf('/media/ag-cherrmann/cramirez/tcd8ExscSeq/figures/tox_vs_stat5a_memVsTrans_chronic.pdf',
    height = 5, width = 13)
vln_plot + scat_plot                
dev.off()






#####################################################################
##  Comparing Intermediate Vs Exhausted

## Visualizations
vln_plot <- seurat[[]] %>%
        filter(nina_annotations %in% c('trans', 'exh')) %>%
        filter(condition == 'Chronic') %>%
        select(condition,
               nina_annotations,
               tox_sign_up1, 
               tox_sign_down2, 
               stat5a_sign3) %>%
        rename('TOX(-) Up Signature'=tox_sign_up1,
               'TOX(-) Down Signature'=tox_sign_down2,
               'STAT5a Signature'=stat5a_sign3) %>%
        melt() %>%
        ggplot(aes(x=nina_annotations, y=value)) +
        geom_violin(draw_quantiles = 0.5) +
        stat_compare_means(label.x = 1, label = 'p.format') +
        geom_jitter(alpha=0.3) +
        facet_wrap(~variable, scales = 'free') +
        theme_bw() +
        theme(legend.position = 'none',
              panel.grid = element_blank()) +
        labs(x='', y='Signature Score',
             title = 'Intermediate Vs Exhausted',
             subtitle = 'Chronic')

scat_plot <- seurat[[]] %>%
        filter(nina_annotations %in% c('trans', 'exh')) %>%
        filter(condition == 'Chronic') %>%
        select(condition,
               nina_annotations,
               tox_sign_up1, 
               tox_sign_down2, 
               stat5a_sign3) %>%
        ggplot(aes(x=stat5a_sign3, 
                   y=tox_sign_down2,
                   colour=nina_annotations)) +
        geom_point() +
        stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                       after_stat(rr.label), sep = "*\", \"*"))) +
        geom_smooth(method = 'lm') +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        labs(x='STAT5a Signature', 
             y='TOX(-) Down Signature',
             colour='') 

vln_plot + scat_plot


pdf('/media/ag-cherrmann/cramirez/tcd8ExscSeq/figures/tox_vs_stat5a_TransVsExh_chronic.pdf',
    height = 5, width = 13)
vln_plot + scat_plot                
dev.off()



#####################################################################
##  Comparing Intermediate Vs Exhausted

## Visualizations
vln_plot <- seurat[[]] %>%
        filter(nina_annotations %in% c('mem', 'exh')) %>%
        filter(condition == 'Chronic') %>%
        select(condition,
               nina_annotations,
               tox_sign_up1, 
               tox_sign_down2, 
               stat5a_sign3) %>%
        rename('TOX(-) Up Signature'=tox_sign_up1,
               'TOX(-) Down Signature'=tox_sign_down2,
               'STAT5a Signature'=stat5a_sign3) %>%
        melt() %>%
        ggplot(aes(x=nina_annotations, y=value)) +
        geom_violin(draw_quantiles = 0.5) +
        stat_compare_means(label.x = 1, label = 'p.format') +
        geom_jitter(alpha=0.3) +
        facet_wrap(~variable, scales = 'free') +
        theme_bw() +
        theme(legend.position = 'none',
              panel.grid = element_blank()) +
        labs(x='', y='Signature Score',
             title = 'Memory Vs Exhausted',
             subtitle = 'Chronic')

scat_plot <- seurat[[]] %>%
        filter(nina_annotations %in% c('mem', 'exh')) %>%
        filter(condition == 'Chronic') %>%
        select(condition,
               nina_annotations,
               tox_sign_up1, 
               tox_sign_down2, 
               stat5a_sign3) %>%
        ggplot(aes(x=stat5a_sign3, 
                   y=tox_sign_down2,
                   colour=nina_annotations)) +
        geom_point() +
        stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                       after_stat(rr.label), sep = "*\", \"*"))) +
        geom_smooth(method = 'lm') +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        labs(x='STAT5a Signature', 
             y='TOX(-) Down Signature',
             colour='') 

vln_plot + scat_plot


pdf('/media/ag-cherrmann/cramirez/tcd8ExscSeq/figures/tox_vs_stat5a_MemVsExh_chronic.pdf',
    height = 5, width = 13)
vln_plot + scat_plot                
dev.off()









