
source('R/parse-SNP.R')
source('R/snipe.R')
source('R/snap.R')

library(dplyr)

# input data files
dbsnp.file <- 'data/processed/dbSNP.GRCh37.p13.build142.sqlite'

starbase.gr <- readRDS('data/processed/starbase.Rds')
targetscan.gr <- readRDS('data/processed/targetscan.Rds')
miranda.gr <- readRDS('data/processed/miranda.Rds')
mir.targets.gr <- merge.granges.aggressively(meta.columns=list(mir.target.db=c('Starbase', 'TargetScan', 'miranda')),
                                             starbase.gr, targetscan.gr, miranda.gr)


# input data files
scz.file <- 'data/base/SCZ-GWAS-risk-SNPs.csv'

############## Schizophrenia risk SNPs
scz.snps <- read.table(scz.file, stringsAsFactors = F,
                      header = T, quote = '"', sep=';', strip.white = T)

colnames(scz.snps)[colnames(scz.snps) == 'Index.SNP'] <- 'SNPs'

CAD.SNP.df <- extend.with.LD(scz.snps)
CAD.SNP.gr <- get.hg19.positions(CAD.SNP.df, dbSNP.file = dbsnp.file)

snp.mir.overlap.matrix <- as.matrix(findOverlaps(CAD.SNP.gr, mir.targets.gr))
if(length(snp.mir.overlap.matrix) == 0)
  stop('There is no overlap between risk snps and miR targets.')

result.table <- generate.final.table(mir.targets.gr,
                                     CAD.SNP.gr,
                                     snp.mir.overlap.matrix)

write.table(result.table, 'snp-mir-SCZ-starbase-newpipeline.tsv', row.names=F, sep='\t')
