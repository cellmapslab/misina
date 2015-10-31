source('R/parse-SNP.R')
source('R/snipe.R')
source('R/snap.R')

library(dplyr)

# input data files
dbsnp.file <- 'data/processed/dbSNP.GRCh37.p13.build142.sqlite'
dbsnp.file2 <- 'data/processed/dbsnp.leveldb'

starbase.gr <- readRDS('data/processed/starbase.Rds')
targetscan.gr <- readRDS('data/processed/targetscan.Rds')
miranda.gr <- readRDS('data/processed/miranda.Rds')
mir.targets.gr <- merge.granges.aggressively(meta.columns=list(mir.target.db=c('Starbase', 'TargetScan', 'miranda')),
                                             starbase.gr, targetscan.gr, miranda.gr)

############### Bipolar SNPs
bp.file <- 'data/base/BP-GWAS-risk-SNPs.csv'
bp.snps <- read.table(bp.file, stringsAsFactors = F, sep = ';',
                      header = T, quote = '"', strip.white = T)
bp.snps <- bp.snps[bp.snps$P.GC < 5e-8,]
colnames(bp.snps)[1] <- 'SNPs'

CAD.SNP.df <- extend.with.LD(bp.snps)
CAD.SNP.gr <- get.hg19.positions(CAD.SNP.df, dbSNP.file = dbsnp.file)

snp.mir.overlap.matrix <- as.matrix(findOverlaps(CAD.SNP.gr, mir.targets.gr))
if(length(snp.mir.overlap.matrix) == 0)
  stop('There is no overlap between risk snps and miR targets.')

result.table <- generate.final.table(mir.targets.gr,
                                     CAD.SNP.gr,
                                     snp.mir.overlap.matrix)
