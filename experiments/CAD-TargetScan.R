
source('R/parse-SNP.R')
source('R/snipe.R')
source('R/snap.R')

# input data files --------------------------------------------------------

#Lines starting with 'track name...' and 'browser pack...' deleted
mir.target.file <- 'data/processed/TargetScan_hg19Cons_ALL_CHRS.bed'
dbsnp.file <- 'data/processed/dbSNP.GRCh37.p13.build142.sqlite'

mir.targets.gr <- import(mir.target.file, genome = 'hg19')
mir.targets.gr$geneName <- sapply(strsplit(mir.targets.gr$name, ':'), `[[`, 1)
mir.targets.gr$mir <- sapply(strsplit(mir.targets.gr$name, ':'), `[[`, 2)
mir.targets.gr$name <- NULL

##########################  CAD risk snps
CAD.SNP.df <- read.csv2('data/base/knownLociCAD_Jan2015_CW.csv', stringsAsFactors = F) #new list from Jeanette
####CAD.SNP.df <- read.csv2('CAD_risk_SNPs_Nov2014.csv')
stroke.snps <- read.csv('data/processed/stroke-snps.csv', stringsAsFactors = F)$SNPs

#WARNING!!!
#dbSNP and Starbase genome assemblies must be exactly same
CAD.SNP.gr <- build.CAD.SNP.gr(snps = c(CAD.SNP.df$Marker_ID, stroke.snps),
                               chr = c(CAD.SNP.df$CHR, rep('1', 11)),
                               pos = c(CAD.SNP.df$POS, rep(1, 11)),
                               dbSNP.file = dbsnp.file)

#get p-values from a separate file
CAD.ALL.pvalue.df <- fread('data/processed/1000G.MA.cad.add.fixrand.subset.Feb2015.CW.csv', #same file with "sort -uV"
                           sep = ';', header = T)
snp.p <- as.data.frame(CAD.ALL.pvalue.df[,c('legendrs', 'p_dgc', 'n_studies'),
                                         with = F])

mcols(CAD.SNP.gr)$pval <- snp.p[ match(CAD.SNP.gr$SNPid, snp.p$legendrs),2]
mcols(CAD.SNP.gr)$n.studies <- snp.p[ match(CAD.SNP.gr$SNPid, snp.p$legendrs),3]


snp.mir.overlap.matrix <- as.matrix(findOverlaps(CAD.SNP.gr, mir.targets.gr))
if(length(snp.mir.overlap.matrix) == 0)
  stop('There is no overlap between risk snps and miR targets.')

result.table <- generate.final.table(mir.targets.gr,
                                     CAD.SNP.gr,
                                     snp.mir.overlap.matrix)

write.table(result.table, 'output/data/snp-mir-CAD-TargetScan-w-stroke.tsv',
            row.names = F, sep = '\t')

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
