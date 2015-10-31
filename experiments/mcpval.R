library(data.table)
library(ggplot2)
library(GGally)
library(Cairo)
library(RColorBrewer)
library(mgcv)
source('R/parse-SNP.R')

all.snps <- fread('data/base/dbSNP.GRCh37.p13.build142.bed')
all.snps.gr <- GRanges(seqnames = all.snps[,V1],
                       ranges = IRanges(start = all.snps[,V2]+1, width = 1),
                       SNPid = all.snps[,V4])



starbase.file <- 'data/base/starBase_Human_Pan-Cancer_MiRNA-TargetInteractions2015-01-13_22-20.csv'
dbsnp.file <- 'data/processed/dbSNP.GRCh37.p13.build142.sqlite'

res <- build.mir.target.gr(starbase.file)
starbase.df <- res$starbase
mir.targets.gr <- res$mir.targets.gr
rm(res)

##########################  CAD risk snps
CAD.SNP.df <- read.csv2('data/base/knownLociCAD_Jan2015_CW.csv') #new list from Jeanette
####CAD.SNP.df <- read.csv2('CAD_risk_SNPs_Nov2014.csv')

#WARNING!!!
#dbSNP and Starbase genome assemblies must be exactly same
CAD.SNP.gr <- build.CAD.SNP.gr(CAD.SNP.df$Marker_ID, dbsnp.file)
CAD.SNP.gr <- unique(CAD.SNP.gr)
CAD.SNP.gr <- subsetByOverlaps(CAD.SNP.gr, all.snps.gr)

mcols(all.snps.gr)$inMir <- F
all.vs.mir.overlap <- findOverlaps(mir.targets.gr, all.snps.gr)
mcols(all.snps.gr)$inMir[all.vs.mir.overlap@subjectHits] <- T

mcols(all.snps.gr)$isRisk.withLD <- F
all.vs.risk.overlap <- findOverlaps(CAD.SNP.gr, all.snps.gr)
mcols(all.snps.gr)$isRisk.withLD[all.vs.risk.overlap@subjectHits] <- T
tbl.withLD <- table(all.snps.gr$inMir, all.snps.gr$isRisk.withLD)
fisher.test(tbl.withLD)

mcols(all.snps.gr)$isRisk.woLD <- F
all.vs.risk.overlap2 <- findOverlaps(CAD.SNP.gr[CAD.SNP.gr$IsProxyOf == ''], all.snps.gr)
mcols(all.snps.gr)$isRisk.woLD[all.vs.risk.overlap2@subjectHits] <- T
tbl.woLD <- table(all.snps.gr$inMir, all.snps.gr$isRisk.woLD)
fisher.test(tbl.woLD)
