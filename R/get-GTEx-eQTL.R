
library(readr)

gtex.files <- list.files('data/base/GTEx/', full.names = T)
gtex.tissues <- sapply(strsplit(basename(gtex.files), '.', fixed=T), `[`, 1)
gtex.df <- Map(function(f, tissue) {

  x <- read_delim(f, '\t')
  x$eQTL.tissue <- tissue
  x$eQTL.source <- 'GTEx'
  x

  }, gtex.files, gtex.tissues)

names(gtex.df) <- NULL
gtex.df <- do.call(rbind, gtex.df)
saveRDS(gtex.df, 'data/processed/GTEx.Rds')
