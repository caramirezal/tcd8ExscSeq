library(GEOquery)

## Downloading Miller Data
if ( ! dir.exists('data/miller2019/')){
    dir.create('data/miller2019/')
}
barcodes_url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE122712&format=file&file=GSE122712%5Fbarcodes%2Etsv%2Egz'
genes_url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE122712&format=file&file=GSE122712%5Fgenes%2Etsv%2Egz'
counts_url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE122712&format=file&file=GSE122712%5Fmatrix%2Emtx%2Egz' 
if (! file.exists('data/miller2019/barcodes.tsv')){
        download.file(barcodes_url, destfile = 'data/miller2019/GSE122712_barcodes.tsv.gz')
        gunzip('data/miller2019/GSE122712_barcodes.tsv.gz',
               destname = 'data/miller2019/barcodes.tsv', 
               remove = FALSE)
}
if (! file.exists('data/miller2019/genes.tsv')) {
        download.file(genes_url, destfile = 'data/miller2019/GSE122712_genes.tsv.gz')
        gunzip('data/miller2019/GSE122712_genes.tsv.gz', 
               destname = 'data/miller2019/genes.tsv', 
               remove = FALSE)
}
if (! file.exists('data/miller2019/matrix.mtx.gz')) {
        download.file(counts_url, destfile = 'data/miller2019/matrix.mtx.gz')
        gunzip('data/miller2019/matrix.mtx.gz', 
               destname = 'data/miller2019/matrix.mtx', 
               remove = FALSE)
}


