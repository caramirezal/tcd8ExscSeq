## Matching gene expression and open chromatin regions
## in progenitor cells profiles from Miller and Hofmann data

## Dependencies
library(dplyr)
library(rio)
library(VennDetail)
library(ggrepel)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(org.Mm.eg.db)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

set.seed(333)

## Setting working directory
path2project <- '/Users/carlosramirez/sc/tcd8ExscSeq/'
setwd(path2project)

##############################################################
## Loading Footprinting
footprint <- read.table(
        'analysis/prog_vs_exh_lcmv_bindetect_results.txt', 
        header = TRUE
)
footprint <- filter(footprint,
                    progenitor_terminal_pvalue < 0.05)
head(footprint)
dim(footprint)

#############################################################
## Loading TF ranked by a random forest based in their 
## transcription factor activities and var importance to 
## predict Memory-like vs Terminally exhausted CD8 T cells
n_tfa <- 30
rforest <- read.table('analysis/tfa_rforest_ranked.tsv.gzip')
tfs.top <- rforest$gene[1:n_tfa]
tfs.top

############################################################
## Loading DEGs between (Progenitor) Memory-like vs 
## Exhausted cells

# degs_url <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41590-020-00817-w/MediaObjects/41590_2020_817_MOESM4_ESM.xlsx'
degs <- readRDS('analysis/mem_vs_exh_hofmann_degs.rds')
degs <- filter(degs, p_val_adj < 0.05)
dim(degs)
head(degs)

######################################################################
## Intersecting TFA and differential footprints

tfs.intersected <- intersect(
        toupper(footprint$name),
        toupper(tfs.top)
)


pdf('figures/footprinting_degs_intersection.pdf')
footprint %>%
        mutate(highlight=ifelse(toupper(name) %in%
                                       tfs.intersected,
                                TRUE, FALSE)
               ) %>%
        mutate(gene_label=ifelse(highlight==TRUE, 
                                 name, '')) %>%
        ggplot(aes(x=progenitor_terminal_change,
                   y=-log(progenitor_terminal_pvalue),
                   label=gene_label,
                   colour=highlight)) +
                geom_point() +
                geom_text_repel(colour='black', 
                                force = 0.1, 
                                max.overlaps = 500) +
                theme_bw() +
                theme(panel.grid = element_blank(),
                      legend.position = 'none') +
                scale_color_manual(values = c('steelblue',
                                              'red')) +
                xlab('Progenitor vs Terminal score change') +
                ylab('-Log(p-value)')
dev.off()

###################################################################
## Matching gene expression with open chromatin regions

url_peaks <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41590-019-0312-6/MediaObjects/41590_2019_312_MOESM5_ESM.xlsx'
peaks <- import(url_peaks, which = 2)

## Creating bed file
dir.create('tmp/')
peaks.bed <- select(peaks, 
                    Chr, Start, End,
                    `LCMV_R1_Stem-like`:`LCMV_R2_TerminallyExh`,
                    `log2FC: Group LCMV_Stem-like_Tet vs LCMV_TerminallyExh_Tet`,
                    `q-value: Group LCMV_Stem-like_Tet vs LCMV_TerminallyExh_Tet`)
colnames(peaks.bed) <- gsub(' ', '_', colnames(peaks.bed))
head(peaks.bed)
peaks.bed <- filter(
        peaks.bed, 
        `q-value:_Group_LCMV_Stem-like_Tet_vs_LCMV_TerminallyExh_Tet` < 0.05
)
dim(peaks.bed)
write.table(peaks.bed, 
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE,
            sep = '\t',
            gzfile('tmp/peaks.bed.gz'))


peak <- readPeakFile('tmp/peaks.bed.gz')
peak

peakAnno <- annotatePeak('tmp/peaks.bed.gz', 
                         tssRegion=c(-3000, 3000),
                         TxDb=txdb, 
                         annoDb="org.Mm.eg.db")
peak.ann <- as.data.frame(peakAnno) 
names(peak.ann)[6:11] <- names(peaks.bed)[4:9]
dim(peak.ann)        
head(peak.ann)
unlink('tmp/', recursive = TRUE)

## Saving annotated file
write.table(
        peak.ann, 
        gzfile('data/41590_2019_312_MOESM5_ESM_miller_ref_annotated_peaks.txt.gz'),
        sep = '\t',
        quote = FALSE, 
        row.names = FALSE
)

## Vulcano plot
length(intersect(toupper(peak.ann$SYMBOL), 
                 toupper(rownames(degs))))
peak.ann <- mutate(peak.ann, 
                   intersection=ifelse(toupper(SYMBOL) %in% 
                                        toupper(rownames(degs)),
                                'DP:DEG', 'non-shared')) %>%
            mutate(gene_label=ifelse(intersection=='DP:DEG' &
                                             abs(`log2FC:_Group_LCMV_Stem-like_Tet_vs_LCMV_TerminallyExh_Tet`) > 2.5,
                                     SYMBOL, ''))
head(peak.ann)
ggplot(peak.ann, aes(x=`log2FC:_Group_LCMV_Stem-like_Tet_vs_LCMV_TerminallyExh_Tet`,
                   y=-log10(`q-value:_Group_LCMV_Stem-like_Tet_vs_LCMV_TerminallyExh_Tet`),
                   colour=intersection,
                   label=gene_label)) +
                geom_point() +
                geom_point(aes(x=`log2FC:_Group_LCMV_Stem-like_Tet_vs_LCMV_TerminallyExh_Tet`,
                           y=-log10(`q-value:_Group_LCMV_Stem-like_Tet_vs_LCMV_TerminallyExh_Tet`),
                           colour=intersection), 
                           data=subset(peak.ann, 
                                       intersection=='DP:DEG')) +
                        scale_color_manual(values = c('red', 'steelblue')) +
                geom_text_repel(max.overlaps = 10000, 
                                force = 6, 
                                color='black') +
                theme_bw() +
                theme(panel.grid = element_blank()) +
                xlab('Log2 FC') +
                ylab('-Log10(p-value)') +
                ggtitle('LCMV infected Progenitor-like vs Terminally Exhausted CD8 T cells')
                        

