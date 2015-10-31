source('R/parse-SNP.R')

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(ggplot2)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb, force=TRUE) <- paste0('chr', 1:22)

starbase.file <- 'data/base/starBase_Human_Pan-Cancer_MiRNA-TargetInteractions2015-01-13_22-20.csv'

starbase.gr <- build.mir.target.gr(starbase.file, strand.info=T)$mir.targets.gr
genome(starbase.gr) <- 'hg19'
names(starbase.gr) <- mcols(starbase.gr)$mir.name


mir <- locateVariants(starbase.gr, txdb, AllVariants())
annot.count <- as.data.frame(table(mir$LOCATION))
colnames(annot.count) <- c('Annotation', 'Freq')
annot.count$Perc = paste0(round(annot.count$Freq / sum(annot.count$Freq) * 100, 3), '%')


ggplot(annot.count, aes(x=Annotation, y=Freq)) +
  geom_bar(stat='identity') +
  geom_text(aes(label = Perc), vjust=-0.5) +
  scale_y_continuous(breaks=seq(0, 1.7e6, 1e5)) +
  labs(y='Frequency', title='miR binding site annotations')
ggsave('output/figs/mir-annotation.png', width=12, height=10)
