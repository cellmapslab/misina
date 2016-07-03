
library(tidyr)
library(jsonlite)

mirmine <- read.csv('data/base/miRmine - Human miRNA Expression Database.csv',
                    header = T, stringsAsFactors = F, check.names = F)

mirmine <- gather(mirmine, 'tissue', 'expression', -`miRNA ID`)
colnames(mirmine)[1] <- 'mir'

tissues <- as.factor(trimws(unlist(lapply(regmatches(mirmine$tissue, regexec('\\((.*?)\\)', mirmine$tissue)), `[`, 2))))

mirmine.grouped <- aggregate(mirmine$expression, list(mir=mirmine$mir, tissue=tissues), median)
colnames(mirmine.grouped)[3] <- 'expression'

result <- by(mirmine.grouped, mirmine.grouped$mir, function(df){ 
  toJSON(df[, c('expression', 'tissue')])
  }, simplify = F)

df <- data.frame(do.call(rbind, result), stringsAsFactors = F)
result <- data.frame(mir=rownames(df), mirmine.median.expression=df[[1]], stringsAsFactors = F)

saveRDS(result, 'data/processed/mirmine_median_json_2016_april.Rds')

mirmine$group <- as.factor(tissues)
saveRDS(mirmine, 'data/processed/mirmine_grouped_2016_april.Rds')

result <- by(mirmine.grouped, mirmine.grouped$mir, function(df){ 
  toJSON(df[, c('expression', 'tissue')])
  }, simplify = F)

df <- data.frame(do.call(rbind, result), stringsAsFactors = F)
result <- data.frame(mir=rownames(df), mirmine.median.expression=df[[1]], stringsAsFactors = F)

saveRDS(result, 'data/processed/mirmine_median_json_2016_april.Rds')



# human readable grouping -------------------------------------------------

result.human <- by(mirmine.grouped, mirmine.grouped$mir, function(df){ 
  paste(df$tissue, round(df$expression,1), sep = ':', collapse = ',')
  }, simplify = F)

df.human <- data.frame(do.call(rbind, result.human), stringsAsFactors = F)
result.human <- data.frame(mir=rownames(df.human), mirmine.median.expression=df.human[[1]], stringsAsFactors = F)
saveRDS(result.human, 'data/processed/mirmine_median_humanreadable_2016_april.Rds')
