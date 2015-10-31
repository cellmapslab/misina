
library(readr)
library(dplyr)

mir.target.file <- 'data/processed/TargetScan_hg19Cons_ALL_CHRS.bed'
mir.info.file <- 'data/base/Predicted_Targets_Info.txt'

mir.targets.gr <- import(mir.target.file, genome = 'hg19')
mir.targets.gr$geneName <- sapply(strsplit(mir.targets.gr$name, ':'), `[[`, 1)
mir.targets.gr$mir <- sapply(strsplit(mir.targets.gr$name, ':'), `[[`, 2)
mir.targets.gr$name <- NULL
mcols(mir.targets.gr) <- mcols(mir.targets.gr)[, c('score', 'geneName', 'mir')]
names(mcols(mir.targets.gr)) <- c('score', 'gene', 'mir')
mir.targets.gr <- keepStandardChromosomes(mir.targets.gr)
mir.targets.gr


#get seed type info from TargetScan info file
info <- read_delim(mir.info.file, delim='\t', na = 'NULL')
info <- info[info$`Species ID` == 9606,]
info <- info[, c("Gene Symbol", "miR Family", 'Seed match')]
colnames(info)[colnames(info) == 'Seed match'] <- "seed.match"

target.tbl <- tbl_df(as.data.frame(mcols(mir.targets.gr)))

x <- full_join(target.tbl, info, 
              by=c('gene'='Gene Symbol', 'mir'='miR Family'))
x <- group_by(x, gene, mir) %>% distinct()
last <- left_join(target.tbl, x, by=c('score', 'gene', 'mir'))
last <- rename(last, seed.category=seed.match)

mcols(mir.targets.gr) <- last
saveRDS(mir.targets.gr, 'data/processed/targetscan.Rds')
