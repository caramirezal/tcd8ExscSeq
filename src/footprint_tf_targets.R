## Evaluation of targets of TFs in Progenitors-like TEx cells
## identified by differential footprint

## Dependencies
library(dplyr)
library(VennDetail)
library(ggrepel)

set.seed(333)

path2project <- '/Users/carlosramirez/sc/tcd8ExscSeq/'
setwd(path2project)

################################################################
## Reading peak scores
files <- list.files('analysis/tobias/progenitor_vs_terminal/', 
                    full.names = TRUE)
tfbs.list <- lapply(files,
                    function(x)
                            read.table(x, header = TRUE))
## Merging TFBS
tfbs.df  <- do.call(rbind, tfbs.list)
tfbs.df <- filter(tfbs.df, progenitor_bound == 1)
tfbs.df <- mutate(tfbs.df, tf=gsub('_.*', '', TFBS_name))
dim(tfbs.df)
head(tfbs.df)
rm(tfbs.list)

###############################################################
## Loading Footprint
footprint <- read.table(
        'analysis/prog_vs_exh_lcmv_bindetect_results.txt', 
        header = TRUE
)
footprint <- filter(footprint,
                    progenitor_terminal_pvalue < 0.05)

head(footprint)
dim(footprint)


#############################################################
## Read DEGs 
degs <- read.table('analysis/mem_vs_exh_hofmann_degs.tsv.gz')
degs <- filter(degs, p_val_adj < 0.05)
head(degs)
dim(degs)

###############################################################
## Reading SCENIC regulons
regulons.scenic <- read.csv(
        'analysis/scenic/scenic_hgv_6000/reg.csv',
        skip = 2, header = TRUE
) 
dim(regulons.scenic)
head(regulons.scenic, 3)

##############################################################
## i) From ATAC Footprint selecting top 50 FC TFs 
## & targets in DEGs
top50TFs <- footprint %>%
        filter(!grepl(':|\\(', name)) %>%
        arrange(desc(progenitor_terminal_change)) %>%
        head(50) %>%
        select(name) %>%
        unlist %>%
        unname() %>%
        toupper()
FP.selected <- filter(tfbs.df, toupper(tf) %in% top50TFs)
FP.selected <- filter(FP.selected, toupper(gene_name) %in% rownames(degs))
FP.selected <- select(FP.selected, tf, gene_name)
FP.selected <- FP.selected[!duplicated(FP.selected), ]
dim(FP.selected)
head(FP.selected)
## [1] 3165    3
## There are allmost 20,000 

## ii) From TFA SCENIC select TFs with targets in DEGs
# function to extract TFs in column 7 from SCENIC results
extract_targets <- function(x){
        targets <- regmatches(x, gregexpr("'[^']*'", x))[[1]] 
        gsub("'", '', targets)
}
## Extracting targets having as targets at least one of the 
## genes in the progenitor signature
df <- data.frame()
i <- 1
for (i in 1:nrow(regulons.scenic)){
        targets <- extract_targets(regulons.scenic$X.7[i])
        tf.df <- data.frame(TF=rep(regulons.scenic$TF[i], 
                                   length(targets)),
                            target=targets)
        df <- rbind(df, tf.df)
}
## Removing duplications
df <- df[!duplicated(df), ]
df.scenic <- filter(df, target %in% rownames(degs))
## dim(df.scenic)
## head(df.scenic)
## [1] 756   3

## iii) Intersect interactions
df.scenic <- mutate(df.scenic, 
                  interaction=paste(TF, target, sep = ':'))
FP.selected <- mutate(FP.selected, 
                      interaction=paste(tf, toupper(gene_name), 
                                        sep = ':'))

## Visualizing intersection
int2.vp <- venndetail(list(SCENIC=df.scenic$interaction,
                           TFBS=FP.selected$interaction))
plot(int2.vp)
dev.off()
plot(int2.vp, type = "vennpie")

inter.filtered <- intersect(df.scenic$interaction,
                            FP.selected$interaction)

df.filtered <- filter(df.scenic, interaction %in% inter.filtered)
dim(df.filtered)
df.filtered

write.table(df.filtered,
            file = 'analysis/filtered_interactions_scenic_dif_peaks.tsv',
            sep = '\t',
            quote = FALSE,
            row.names = FALSE)

######################################################################

## 
tfs.ranked <- footprint %>%
        filter(!grepl(':|\\(', name)) %>%
        arrange(desc(progenitor_terminal_change)) %>%
        mutate(rank=1:n()) %>%
        select(name, rank)
dim(tfs.ranked)
head(tfs.ranked)


tfs.ranked$'n_Hofmann_DEG' <- 0
for (i in 1:nrow(tfs.ranked)){
        tfbs.sel <- filter(tfbs.df, tf==tfs.ranked$name[i])
        if (nrow(tfbs.sel) > 0 ) {
                targets <- toupper(unique(tfbs.sel$gene_name))
                tfs.ranked$n_Hofmann_DEG[i] <- length(
                            unique(intersect(targets,
                                   rownames(degs)))
                            )
        }
}
tfs.ranked <- mutate(tfs.ranked, 
                     percentage_degs=100*n_Hofmann_DEG/nrow(degs))      
head(tfs.ranked)

tfs.ranked <- arrange(tfs.ranked, desc(percentage_degs)) %>%
        mutate(rank2=1:n()) %>%
        mutate(highlight=ifelse(rank<30 |
                                        name %in% 'PRDM1', 
                                TRUE, FALSE)) 
tfs.ranked <- mutate(tfs.ranked,
                     tf_name=ifelse(highlight==TRUE, name, ''))
tfs.highlight <- subset(tfs.ranked, highlight==TRUE)

pdf('figures/tf_percentage_of_target_degs.pdf')
ggplot(tfs.ranked, aes(x=rank2, 
                   y=percentage_degs,
                   label=tf_name)) +
                geom_point() +
                geom_point(data=tfs.highlight,
                           aes(x=rank2,
                               y=percentage_degs),
                               color='red') +
                geom_text_repel(max.overlaps = 100, 
                                force = 50) +
                xlab('Rank') +
                ylab('% DEGs in targets') +
                theme_bw() +
                theme(panel.grid = element_blank())

tfs.ranked %>%
        select(rank, percentage_degs) %>%
        ggplot(aes(x=rank, 
                   y=percentage_degs)) +
                geom_bar(position="stack", 
                         stat="identity",
                         fill='steelblue') +
                theme_bw() +
                theme(axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      panel.grid = element_blank()) +
                xlab('Rank (Footprint score)') +
                ylab('% DEGs in targets')

dev.off()
