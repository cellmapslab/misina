
library(readr)

gtex.files <- list.files('data/base/GTEx', full.names = T)
gtex.tissues <- sapply(strsplit(basename(gtex.files), '_', fixed=T), function(x){paste(x[-length(x)], collapse = ' ')})
gtex.df <- Map(function(f, tissue) {

  print(tissue)
  x <- read_delim(f, '\t')
  x$eQTL.tissue <- tissue
  x$eQTL.source <- 'GTEx'
  x <- x[, c('rs_id_dbSNP142_GRCh37p13','beta', 't_stat', 'p_value', 'gene_name', 'eQTL.source', 'eQTL.tissue')]
  x

  }, gtex.files, gtex.tissues)

names(gtex.df) <- NULL
gtex.df <- do.call(rbind, gtex.df)

names(gtex.df)[1] <- 'SNP'
my_db <- src_sqlite("data/processed/GTExv6.sqlite", create = T)
sq <- copy_to(my_db, gtex.df, name='GTExv6', temporary = FALSE, indexes = list('SNP'))

#saveRDS(gtex.df, 'data/processed/GTExv6.Rds')