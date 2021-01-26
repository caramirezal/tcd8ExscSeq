## Matching gene expression and open chromatin regions
## in progenitor cells profiles from Miller and Hofmann data

## Dependencies
library(dplyr)
library(rio)
library(VennDetail)
library(ggrepel)

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

############################################################
## Loading DEGs between (Progenitor) Memory-like vs 
## Exhausted cells

# degs_url <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41590-020-00817-w/MediaObjects/41590_2020_817_MOESM4_ESM.xlsx'
degs <- readRDS('analysis/mem_vs_exh_hofmann_degs.rds')
degs <- filter(degs, p_val_adj < 0.05)
dim(degs)
head(degs)

tfs.intersected <- intersect(
        toupper(footprint$name),
        toupper(rownames(degs))
)


ven <- venndetail(list(ATAC = toupper(footprint$name), 
                       RNA = toupper(rownames(degs))))

pdf('figures/venn_footprinting_degs_intersection.pdf')
plot(ven)
dev.off()

pdf('figures/footprinting_degs_intersection.pdf')
footprint %>%
        mutate(highlight=ifelse(toupper(name) %in%
                                        toupper(rownames(degs)),
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

files <- list.files('analysis/tobias/progenitor_vs_terminal/', 
                    full.names = TRUE)
tfbs.list <- lapply(files,
                    function(x)
                            read.table(x, 
                                       header = TRUE))

## Merging TFBS
tfbs.df  <- do.call(rbind, tfbs.list)
tfbs.df <- filter(tfbs.df, progenitor_bound == 1 )
tfbs.df <- mutate(tfbs.df, tfbs=gsub('_.*', '', TFBS_name))
tfbs.df <- mutate(tfbs.df, tfbs=toupper(tfbs))
dim(tfbs.df)
head(tfbs.df)
rm(tfbs.list)
length(unique(tfbs.df$tfbs))

intersect(toupper(tfbs.df$gene_name), rownames(degs))

any_promoter <- tfbs.df %>% 
        filter(name=='any_promoter') %>% 
        select(gene_name) %>%
        unlist() %>%
        unique() %>%
        toupper()
distal_enhancer <- tfbs.df %>% 
        filter(name=='distal_enhancer') %>% 
        select(gene_name) %>%
        unlist() %>%
        unique() %>%
        toupper()
any_internal <- tfbs.df %>% 
        filter(name=='any_internal') %>% 
        select(gene_name) %>%
        unlist() %>%
        unique() %>%
        toupper()

tfs.intersected <- intersect(
        toupper(tfbs.df$gene_name),
        toupper(rownames(degs))
) %>% sort
tfs.intersected

ven <- venndetail(list(distal_enhancer = distal_enhancer,
                       any_promoter = any_promoter, 
                       any_internal = any_internal,
                       DEGS = toupper(rownames(degs))))
plot(ven, type = "vennpie")
plot.new()
plot(ven)

#####################################################################

hofmann <- readRDS('data/maike2020/hofmann_hcv_seu.rds')
hofmann
