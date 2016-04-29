
source('R/parse-SNP.R')
source('R/snipe.R')
source('R/snap.R')

library(dplyr)
library(RcppLevelDB)

# input data files
dbsnp.file <- 'data/processed/dbsnp.leveldb'

starbase.gr <- readRDS('data/processed/starbase.Rds')
targetscan.gr <- readRDS('data/processed/targetscan.Rds')
miranda.gr <- readRDS('data/processed/miranda.Rds')
mir.targets.gr <- merge.granges.aggressively(meta.columns=list(mir.target.db=c('Starbase', 'TargetScan', 'miranda')),
                                             list(starbase.gr, targetscan.gr, miranda.gr))


##########################  CAD risk snps
CAD.SNP.df <- read.csv2('data/base/knownLociCAD_Jan2015_CW.csv', stringsAsFactors = F, strip.white = T) #new list from Jeanette
conf.snps <- data.frame(SNPs=CAD.SNP.df$Marker_ID,
                        Phenotype='Coronary artery disease',
                        Risk.SNP.Source='Confidential', stringsAsFactors = F)

####CAD.SNP.df <- read.csv2('CAD_risk_SNPs_Nov2014.csv')
stroke.snps <- read.csv('data/processed/stroke-snps.csv', stringsAsFactors = F)

nature.snps.raw <- read.csv('data/base/nature-cardiogram.csv', stringsAsFactors = F, strip.white = T)
nature.snps <- data.frame(SNPs=nature.snps.raw$rs_number,
                          Phenotype='Coronary artery disease',
                          Risk.SNP.Source='CARDIOGRAMplusC4D',
                          stringsAsFactors = F)

gwas.cat.snps.raw <- readRDS('data/processed/gwas-heart-diseases.Rds')
gwas.cat.snps <- data.frame(SNPs=gwas.cat.snps.raw$SNPs,
                            Phenotype=gwas.cat.snps.raw$Disease.Trait,
                            Risk.SNP.Source='GWAS catalog',
                            stringsAsFactors = F)
total.snps <- rbind(conf.snps, stroke.snps, nature.snps, gwas.cat.snps)
total.snps <- aggregate(total.snps,
                         list(SNPs=total.snps$SNPs),
                         function(x)paste0(unique(x), collapse=','))
total.snps <- total.snps[,-1]

#WARNING!!!
#dbSNP and Starbase genome assemblies must match
CAD.SNP.df <- extend.with.LD(total.snps, rsquare=0.78, self.snp.label = 'risk.snp')
CAD.SNP.gr <- get.hg19.positions2(CAD.SNP.df, dbSNP.file = dbsnp.file)

#get p-values from a separate file
# CAD.ALL.pvalue.df <- fread('data/processed/1000G.MA.cad.add.fixrand.subset.Feb2015.CW.csv', #same file with "sort -uV"
#                            sep=';', header=T)
# snp.p <- as.data.frame(CAD.ALL.pvalue.df[,c('legendrs', 'p_dgc', 'n_studies'), with=F])
# mcols(CAD.SNP.gr)$pval <- snp.p[ match(CAD.SNP.gr$SNP, snp.p$legendrs),2]
# mcols(CAD.SNP.gr)$n.studies <- snp.p[ match(CAD.SNP.gr$SNP, snp.p$legendrs),3]

#find overlap between snps and miR BS+25nt flanking region
snp.mir.overlap.hits <- findOverlaps(CAD.SNP.gr, unique(mir.targets.gr)+25) #take flanking region (25ntx2) into acc.
snp.mir.overlap.matrix <- as.matrix(snp.mir.overlap.hits)
if(length(snp.mir.overlap.matrix) == 0)
  stop('There is no overlap between risk snps and miR targets.')

result.table <- generate.final.table(unique(mir.targets.gr),
                                     CAD.SNP.gr,
                                     snp.mir.overlap.matrix, annotate=T)


# visualize the results ---------------------------------------------------

visualize(unique(mir.targets.gr), CAD.SNP.gr, snp.mir.overlap.matrix)


# eQTL --------------------------------------------------------------------

extended.eqtl <- readRDS('data/processed/eQTL-mono-macro-with-LD.Rds')
extended.eqtl <- extended.eqtl[,c('SNP', 'IsProxyOf', 't.stat', 'p.value', 'FDR', 'gene', 'eQTL.source', 'eQTL.tissue')]
names(extended.eqtl) <- c('SNP', 'eQTL.IsProxyOf', 'eQTL.tstat', 'eQTL.pvalue', 'eQTL.FDR', 'eQTL.Gene', 'eQTL.Source', 'eQTL.Tissue')

gtex.eqtl <- readRDS('data/processed/GTEx.Rds')
gtex.eqtl <- gtex.eqtl[, c('SNP', 'T_Stat', 'P_Val', 'Gene_Name', 'eQTL.source', 'eQTL.tissue')]
names(gtex.eqtl) <- c('SNP', 'eQTL.tstat', 'eQTL.pvalue', 'eQTL.Gene', 'eQTL.Source', 'eQTL.Tissue')
gtex.eqtl[] <- lapply(gtex.eqtl, as.character)

final.eqtl <- bind_rows(extended.eqtl, gtex.eqtl)
final.eqtl$eQTL.IsProxyOf[is.na(final.eqtl$eQTL.IsProxyOf)] <- 'direct'
final.eqtl$eQTL.IsProxyOf[final.eqtl$eQTL.IsProxyOf == ''] <- 'direct'

ultimate <- merge(result.table, final.eqtl, all.x=T, by='SNP')
ucol <- ncol(ultimate)
#move eqtl columns towards the beginning
ultimate <- select(ultimate, SNP:mir.target.db, starts_with('eQTL'), everything())
ultimate <- aggregate(ultimate,
                      list(SNP=ultimate$SNP, MIRR=ultimate$mir),
                      function(x)paste(unique(x[x!='']), collapse=','))
ultimate <- ultimate[,-(1:2)]

ultimate <- ultimate[order(ultimate$SNP),]
#add one more column denoting if the target gene == eGene
ultimate$eQTL.Gene.Same.as.Target.gene <- simplify2array(Map(function(gene, egene){
  any(toupper(gene) == strsplit(toupper(egene), ',')[[1]])},
  ultimate$gene, ultimate$eQTL.Gene))

ultimate <- select(ultimate, SNP:mir.target.db, starts_with('eQTL'), everything())

write.table(ultimate, '~/Risk-SNPs-within-miR-BS-corrected-eQTL.tsv', row.names=F, sep='\t')


# Ghanbari section --------------------------------------------------------

l <- c('rs2048327', 'rs3217992', 'rs2306374', 'rs602633')
real.l <- c('rs1810126', 'rs3217992', 'rs9818870', 'rs12740374')

ghanbari <- c('rs3217992',
              'rs1063192',
              'rs7739181',
              'rs12740374',
              'rs6722332',
              'rs629301',
              'rs7528419',
              'rs12190287',
              'rs1810126',
              'rs3088442',
              'rs9818870',
              'rs3732837',
              'rs2036927')

plot.results <- function(result.table) {
  qplot(total.snps$Risk.SNP.Source) + coord_flip()
  ggsave('output/figs/risk-SNP-dist.png')
  qplot(result.table$gene) + coord_flip()
  ggsave('output/figs/result-gene-dist.png')
  qplot(result.table$Risk.SNP.Source) + coord_flip()
  ggsave('output/figs/result-risk-source.png')
  qplot(result.table$mir.target.db) + coord_flip()
  ggsave('output/figs/result-target-db.png')
}


# CairoPDF('output/figs/lddistance_kde.pdf')
#
# qplot(as.integer(result.table$LDDistance),
#       geom=c('density', 'rug'),
#       sides='x',
#       adjust=1/40,
#       xlab='LD Distance (bp)',
#       ylab='Density') + theme_minimal()
# dev.off()
#
#
#
# CairoPDF('output/figs/risksnp_count.pdf')
#
# snp.plot <- aggregate(result.table$mirTarget, list(snpid=result.table$SNPid), function(x)length(unique(x)))
# ggplot(snp.plot, aes(x=reorder(snpid, x), x)) +
#   geom_point() +
#   labs(x='Risk SNPs', y='Number of miRNA targets') +
#   scale_y_discrete(breaks=0:10) +
#   coord_flip() +
#   theme_minimal()
#
# dev.off()
#
#
# CairoPDF('output/figs/gene_count.pdf')
#
# gene.plot <- aggregate(result.table$SNPid, list(gene=result.table$genesymbol), function(x)length(unique(x)))
#
# ggplot(gene.plot, aes(x=reorder(as.character(gene), x), y=x)) +
#   geom_point() +
#   labs(x='Genes', y='SNP occurences') +
#   #scale_y_discrete(breaks=0:max(gene.plot$x)) +
#   coord_flip() +
#   theme_minimal()
#
# dev.off()
#
# library(stargazer)
#
# stargazer(result.table[,c(1, 4:7)],
#           summary=F,
#           rownames=F,
#           align=T,
#           font.size='scriptsize',
#           out='output/figs/result_table.tex')
