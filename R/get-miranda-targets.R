
library(readr)

mc <- read_delim('data/base/miranda-conserved-human_predictions_S_C_aug2010.txt', '\t')
mc2 <- read_delim('data/base/miranda-nonconserved-human_predictions_S_0_aug2010.txt', '\t')

nmc <- nrow(mc)
nmc2 <- nrow(mc2)

mc <- rbind(mc, mc2)
rm(mc2)

pos.df <- as.data.frame(t(sapply(strsplit(mc$genome_coordinates, ':', fixed=T),
                                 `[`, 2:4)),
                        stringsAsFactors=F)

start.end <- as.data.frame(do.call(rbind,
                                   strsplit(sapply(strsplit(pos.df[,2], ',', fixed=T),
                                                   `[`, 1), '-', fixed=T)),
                           stringsAsFactors=F)

mir.gr <- GRanges(seqnames=pos.df[,1],
                  ranges = IRanges(start=as.numeric(start.end[,1]), end=as.numeric(start.end[,2])),
                  strand = substr(pos.df[,3], 1, 1))

mc$conserved <- as.factor(rep(c('conserved', 'nonconserved'), times=c(nmc, nmc2)))
mcols(mir.gr) <- mc[,c('#mirbase_acc', 'mirna_name', 'gene_symbol', 'mirsvr_score', 'conserved', 'seed_cat')]
names(mcols(mir.gr)) <- c('mirbase_acc', 'mir', 'gene', 'score', 'miranda.conserved', 'seed.category')
seqlevelsStyle(mir.gr) <- 'UCSC'
mir.gr <- keepStandardChromosomes(mir.gr)
genome(mir.gr) <- 'hg19'

mer.map = c("0" = 'NA', "6"="6mer", "7"="7mer")
mcols(mir.gr)$seed.category <- mer.map[as.character(mir.gr$seed.category)]

saveRDS(mir.gr, 'data/processed/miranda.Rds')
