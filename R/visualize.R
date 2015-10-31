
library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)


visualize <- function(mir.targets.gr, snps, hit.matrix) {

  hit <- snps[hit.matrix[10,1]]
  hit$id <- hit$SNP
  range <- 100
  from <- start(hit) - range
  to <- end(hit) + range

  chr <- as.character(unique(seqnames(hit)))
  gen <- unique(genome(hit))

  overlapping.mirs <- subsetByOverlaps(mir.targets.gr, hit+range)
  overlapping.mirs$id <- overlapping.mirs$mir
  #overlapping.mirs$group <- as.character(overlapping.mirs$mir.target.db)
  op <- overlapping.mirs[strand(overlapping.mirs) == '+']
  strand(op) <- '-'
  on <- overlapping.mirs[strand(overlapping.mirs) == '-']
  strand(on) <- '+'
  overlapping.mirs <- c(op, on)

  atrack <- AnnotationTrack(overlapping.mirs,
                            #groupAnnotation='group',
                            #just.group='above',
                            #collapse=F,
                            #mergeGroups=F,
                            cex.feature = 0.7,
                            fontcolor.feature = "white",
                            name = "miR Targets")
  #feature(atrack) <- overlapping.mirs$group

  gtrack <- GenomeAxisTrack()
  itrack <- IdeogramTrack(genome = gen, chromosome = chr)

  #TODOL use BSGenomeView for specific range
  strack <- SequenceTrack(Hsapiens, chromosome = chr, genome=gen, add53=TRUE,
                          complement=T)

  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  grtrack <- GeneRegionTrack(txdb, genome = gen,
                             chromosome = chr, name = "UCSC known genes",
                             transcriptAnnotation = "symbol",
                             fontcolor.feature = "black",
                             fontcolor= "black",
                             showFeatureId=T,
                             showId=T,
                            rotation=45,
                            collapse=T,
                            geneSymbols=T)
                             #shape='arrow')

  snpLocations <- UcscTrack(genome=gen, chromosome=chr,
                            track="snp142", from=from, to=to, trackType="AnnotationTrack",
                            start="chromStart", end="chromEnd", id="name", feature="func",
                            strand="strand", shape="box", stacking="dense",
                            name="SNPs", showFeatureId=F, fill='black')

  #workaround for a bug: resize SNPs to length of 1
  snpLocations@range <- resize(snpLocations@range, 1, fix="end", ignore.strand=T)


  plotTracks(list(itrack,
                  gtrack,
                  atrack,
                  snpLocations,
                  strack,
                  AnnotationTrack(hit, name='SNP Hit')),
                  #grtrack),
             from = from,
             to = to,
             featureAnnotation = "id",
             fontcolor.feature = "black",
             #Starbase="darkred",
             #miranda="darkgreen",
             cex.feature = 0.7,
             showOverplotting=TRUE,
             #showId=T,
             collapse=F,
             mergeGroups=F
             #collapseTranscripts = TRUE,
             #shape = "arrow")
  )

}
