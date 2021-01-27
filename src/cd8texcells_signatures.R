## Comparisson of the CD8 TEx cells signatures ins Maike Hofmann 
## cells data signatures and with previously reported signatures
library(readxl)
library(rio)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(viridis)

path2project <- '/media/ag-cherrmann/cramirez/tcd8ExscSeq/'
setwd(path2project)

miller_deg_path <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41590-019-0312-6/MediaObjects/41590_2019_312_MOESM4_ESM.xlsx'

## LCMV TEx Progenitor vs Terminally DEG
lcmv_deg <- import(miller_deg_path, which = 1)

## TIL TEx Progenitor vs Terminally DEG
til_deg <- import(miller_deg_path, which = 2)

df <- merge(lcmv_deg, til_deg, by = 'Gene', 
            suffixes = c('_lcmv', '_til'))
head(df)

pdf('figures/deg_lcmv_til.pdf')
gene_list <- 'Pdcd1'
th <- 5.0
df %>%
        filter(pvalue_lcmv < 0.05 & pvalue_til < 0.05) %>%
        mutate(highlight=ifelse( ( abs(log2FoldChange_lcmv) > th &
                                        abs(log2FoldChange_til) > th ) |
                                         Gene %in% gene_list,
                                TRUE, FALSE)) %>%
        mutate(gene_label=ifelse(highlight==TRUE, Gene, '')) %>%
        ggplot(aes(x=log2FoldChange_lcmv,
                   y=log2FoldChange_til, 
                   label=gene_label)) +
                geom_point(aes(colour=highlight)) +
                           theme_bw() +
                           theme(panel.grid = element_blank(),
                                 legend.position = 'none') +
                geom_text_repel()
dev.off()


#####################################################################
## Getting progenitor and terminal signatures in Hofmann data

seurat <- readRDS('data/integrated_miller_hofmann.rds')

## Subsetting to Exhausted cells
exhausted <- subset(seurat,
                    dataset == 'Hoffmann' & 
                            predicted %in% c('Progenitor Ex',
                                             'Terminally Ex'))
## DEA
Idents(exhausted) <- exhausted$predicted
deg <- FindMarkers(exhausted, 
                   ident.1 = 'Progenitor Ex', 
                   logfc.threshold = 0.1, 
                   min.pct = 0.01)

deg <- add_rownames(deg, var = 'Gene')

filtered_df <- filter(df, pvalue_lcmv < 0.05 & 
                              pvalue_til < 0.05)
shared_degs <- intersect(deg$Gene, toupper(filtered_df$Gene)) 
length(shared_degs)
shared_degs %>% sort

######################################################################


df <- mutate(df, Gene=toupper(Gene))

gene_list <- 'Pdcd1'
th <- 1.0
pdf('figures/intersection_lcmv_til_hofmann.pdf')
df.all %>%
        filter(pvalue_lcmv < 0.05 & pvalue_til < 0.05) %>%
        mutate(highlight=ifelse( ( abs(log2FoldChange_lcmv) > th &
                                           abs(log2FoldChange_til) > th ) &
                                         toupper(Gene) %in% shared_degs,
                                 'Hofmann', 'Non-shared')) %>%
        mutate(gene_label=ifelse(highlight=='Hofmann', Gene, '')) %>%
        ggplot(aes(x=log2FoldChange_lcmv,
                   y=log2FoldChange_til, 
                   label=gene_label)) +
        geom_point(aes(colour=pct.1,
                       size=abs(avg_logFC))) +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        geom_text_repel(colour='red') +
        scale_color_viridis() +
        xlab('Log FC (LCMV infected Pr vs Term DEG )') +
        ylab('Log FC (TIL infected Pr vs Term DEG )')
dev.off()

## Saving results
progenitor_signature <- filter(df.all, Gene %in% shared_degs)
write.table(progenitor_signature,
            file = 'analysis/progenitor_sign_lcmv_til_hofmann.tsv',
            sep = '\t',
            row.names = FALSE)
        
##########################################################################
## Yao C Signatures 
## https://www.nature.com/articles/s41590-019-0403-4

yao_progenitors_path <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41590-019-0403-4/MediaObjects/41590_2019_403_MOESM4_ESM.xlsx'
progenitors_deg <- import(yao_progenitors_path)

yao_memory_path <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41590-019-0403-4/MediaObjects/41590_2019_403_MOESM5_ESM.xlsx'
memory_deg <- import(yao_memory_path)

yao_clusters_degs_path <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41590-019-0403-4/MediaObjects/41590_2019_403_MOESM3_ESM.xlsx'
yao_clusters_degs <- import(yao_clusters_degs_path)

genes_intersection <- intersect(toupper(filter(yao_clusters_degs, 
                         `Cluster ID` ==3)$`Gene Symbol`), 
                                deg$Gene)

th <- 0.325
pdf('figures/yao_cluster3_deg.pdf')
yao_clusters_degs %>%
        filter(`Cluster ID`==3) %>%
        mutate(Gene=toupper(`Gene Symbol`)) %>%
        mutate(highlight=ifelse(Gene %in% genes_intersection,
                                'Hofmann', 'Non-shared'))  %>%
        mutate(gene_label=ifelse(highlight=='Hofmann' & 
                                    `Log Fold Change` > th &
                                         -log10(`Adjusted P Value`) > 20,
                        Gene, '')) %>%
        ggplot(aes(x=`Log Fold Change`,
                   y=-log10(`Adjusted P Value`),
                   colour=highlight, 
                   label=gene_label)) +
                geom_point() + 
                geom_text_repel() + 
                theme_bw() +
                theme(panel.grid = element_blank())
dev.off()
                
