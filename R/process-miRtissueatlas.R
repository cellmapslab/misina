
library(tidyr)
library(dplyr)
library(jsonlite)

mat <- read.delim('data/base/mirTissueAtlas/data_matrix_quantile.txt', 
                  stringsAsFactors = F, dec=',', check.names = F)
mat$mir <- rownames(mat)
longmat <- gather(mat, 'tissue', 'expression', -mir)
tissue.groups <- trimws(gsub('[\\.\\d_]', ' ', longmat$tissue, perl = T))
longmat$tissue <- tissue.groups


mat2 <- read.delim('data/base/mirTissueAtlas/final_GSE11879_series_matrix_quantile.txt', 
                   stringsAsFactors = F, check.names = F)
mat2$mir <- rownames(mat2)
longmat2 <- gather(mat2, 'tissue', 'expression', -mir)
facs2 <- tolower(gsub(' Total RNA. Replicate *\\d*', '', longmat2$tissue, perl = T))
facs2 <- gsub('human ', '', facs2, perl = T)
longmat2$tissue <- facs2

final.matrix <- bind_rows(longmat, longmat2)
final.matrix.grouped <- aggregate(final.matrix$expression, list(mir=final.matrix$mir, tissue=final.matrix$tissue), median)
colnames(final.matrix.grouped)[3] <- 'expression'
#x = spread(final.matrix.grouped, tissue, expression)

final.matrix.json <- by(final.matrix.grouped, final.matrix.grouped$mir, function(df){ 
  toJSON(df[, c('expression', 'tissue')])
  }, simplify = F)
final.matrix.json <- data.frame(do.call(rbind, final.matrix.json), stringsAsFactors = F)
final.matrix.json <- data.frame(mir=rownames(final.matrix.json), mirtissueatlas.median.expression=final.matrix.json[[1]], stringsAsFactors = F)
saveRDS(final.matrix.json, 'data/processed/mirtissueatlas_median_json_2016_may.Rds')
saveRDS(final.matrix, 'data/processed/mirtissueatlas_grouped_2016_may.Rds')

# human readable ----------------------------------------------------------

result.human <- by(final.matrix.grouped, final.matrix.grouped$mir, function(df){ 
  paste(df$tissue, round(df$expression, 1), sep = ':', collapse = ',')
  }, simplify = F)

df.human <- data.frame(do.call(rbind, result.human), stringsAsFactors = F)
result.human <- data.frame(mir=rownames(df.human), mirtissueatlas.median.expression=df.human[[1]], stringsAsFactors = F)
saveRDS(result.human, 'data/processed/mirtissueatlas_median_humanreadable_2016_april.Rds')

