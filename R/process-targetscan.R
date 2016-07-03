
library(readr)
library(dplyr)
library(rtracklayer)

mir.target.file <- 'data/base/Predicted_Targets.hg19.bed'
mir.info.file <- 'data/base/Conserved_Family_Info.txt'
mir.family.file <- 'data/base/miR_Family_Info.txt'


# Read the bed file as a tsv ----------------------------------------------

mir.targets.gr <- read_delim(mir.target.file, delim='\t', col_names = F)
mir.targets.gr$geneName <- sapply(strsplit(mir.targets.gr$X4, ':'), `[[`, 1)
mir.targets.gr$mir <- sapply(strsplit(mir.targets.gr$X4, ':'), `[[`, 2)
mir.targets.gr$X4 <- NULL
mir.targets.gr$X7 <- NULL
mir.targets.gr$X8 <- NULL
mir.targets.gr$X9 <- NULL
mir.targets.gr$X10 <- NULL
mir.targets.gr$X11 <- NULL
mir.targets.gr$X12 <- NULL
colnames(mir.targets.gr) <- c('seqname', 'start', 'end', 'score', 'strand', 'gene', 'mir.family')

#get seed type info from TargetScan info file
info <- read_delim(mir.info.file, delim='\t', na = 'NULL')
info <- info[info$`Species ID` == 9606,]
info <- info[, c("Gene Symbol", "miR Family", 'Seed match')]
colnames(info) <- c('gene', 'mir.family', 'seed.match')

x <- full_join(mir.targets.gr, info)
x <- group_by(x, seqname, start, end, score, strand, gene, mir.family) %>% distinct()
last <- left_join(mir.targets.gr, x)

mir.family <- read_delim(mir.family.file, delim='\t', na = c('', '-'))
mir.family <- mir.family[mir.family$`Species ID` == 9606,]
mir.family <- mir.family[, c(1,4,7)]
colnames(mir.family) <- c('mir.family', 'mir', 'mirbase_acc')

last <- left_join(last, mir.family, by='mir.family')

gr <- GRanges(seqnames = last$seqname,
              ranges = IRanges(start = last$start, last$end),
              strand = last$strand)
colnames(last)[8] <- 'seed.category'
mcols(gr) <- last[, c('score', 'gene', 'mir', 'mir.family', 'seed.category', 'mirbase_acc')]

saveRDS(gr, 'data/processed/targetscan.Rds')
