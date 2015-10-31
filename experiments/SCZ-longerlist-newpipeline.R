
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
scz.file <- 'data/base/scz_longer_list_from_fabian.csv'

############## Schizophrenia risk SNPs
scz.snps <- read.table(scz.file, stringsAsFactors = F,
                      header = F, quote = '"', strip.white = T)

colnames(scz.snps) <- c('Status', 'SNPs', 'CHR', 'POS', 'Effect.Allele',
                        'Other.Allele', 'Directions', 'Effect', 'SE', 'CI', 'pval')

CAD.SNP.df <- extend.with.LD(scz.snps, self.snp.label = 'risk.snp')
CAD.SNP.gr <- get.hg19.positions(CAD.SNP.df, dbSNP.file = dbsnp.file)

snp.mir.overlap.matrix <- as.matrix(findOverlaps(CAD.SNP.gr, mir.targets.gr))
if(length(snp.mir.overlap.matrix) == 0)
  stop('There is no overlap between risk snps and miR targets.')

result.table <- generate.final.table(mir.targets.gr,
                                     CAD.SNP.gr,
                                     snp.mir.overlap.matrix)

#write.table(result.table, 'snp-mir-SCZ-longerlist-newpipeline.tsv', row.names=F, sep='\t')

gtex.eqtl <- readRDS('data/processed/GTEx.Rds')
gtex.eqtl <- gtex.eqtl[, c('SNP', 'T_Stat', 'P_Val', 'Gene_Name', 'eQTL.source', 'eQTL.tissue')]
names(gtex.eqtl) <- c('SNP', 'eQTL.tstat', 'eQTL.pvalue', 'eQTL.Gene', 'eQTL.Source', 'eQTL.Tissue')
gtex.eqtl[] <- lapply(gtex.eqtl, as.character)

final.eqtl <- bind_rows(gtex.eqtl)

ultimate <- merge(result.table, final.eqtl, all.x=T, by='SNP')
ucol <- ncol(ultimate)
#move eqtl columns towards the beginning
ultimate <- ultimate[, c(1:8, (ucol-6):ucol, 9:(ucol-7))]
ultimate <- aggregate(ultimate,
                      list(SNP=ultimate$SNP, MIRR=ultimate$mir),
                      function(x)paste(unique(x[x!='']), collapse=','))
ultimate <- ultimate[,-(1:2)]

ultimate <- ultimate[order(ultimate$SNP),]
write.table(ultimate, '/tmp/mir-snp-bp.tsv', row.names = F, quote = F, sep = '\t')
