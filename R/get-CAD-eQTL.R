
library(readr)

eqtl <- read_delim('data/processed/eQTL_CAD_SNPs0.05_monomacro.csv', '\t')
names(eqtl)[1] <- 'SNPs'
extended.eqtl <- extend.with.LD(eqtl)
extended.eqtl <- extended.eqtl[,c(7:9, 1:6)]
extended.eqtl$eQTL.source <- 'DHM'
extended.eqtl$eQTL.tissue <- 'Heart'
saveRDS(extended.eqtl, 'data/processed/eQTL-with-LD.Rds')
