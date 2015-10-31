
build.starbase <- function(target.file, strand.info=T) {

  starbase <- fread(target.file)

  #TODO: use strsplit to speed up
  starbase.positions <- do.call(rbind,
                                lapply(regmatches(starbase$position,
                                                  regexec('(chr\\w+):(\\d+)-(\\d+)\\[(.)\\]',
                                                          starbase$position)),
                                       `[`, c(2L,3L,4L,5L)))

  #check if all rows matched
  stopifnot(nrow(starbase) == nrow(starbase.positions))

  #we can also simply use mir families
  #mir.names <- sub('(hsa)-(\\w+)-(\\d+).*', '\\1-\\2-\\3', starbase$name, perl=T)

  if (strand.info)
    st <- starbase.positions[,4]
  else
    st <- rep('*', nrow(starbase.positions))

  #build a GenomicRanges object
  mir.targets.gr <- GRanges(seqnames=starbase.positions[,1],
                            IRanges(start=as.integer(starbase.positions[,2]),
                                    end=as.integer(starbase.positions[,3])),
                            strand=st,
                            mir=as.character(starbase$name),
                            gene=as.character(starbase$geneName))
  genome(mir.targets.gr) <- 'hg19'
  #autoplot(mir.targets.gr) + coord_flip()

  #incorporate the remaining miR target metadata
  #mcols(mir.targets.gr) <- cbind(mcols(mir.targets.gr), DataFrame(starbase[,-c(1,2,3),with=F]))

  #filter out long miR binding sites
  mir.targets.gr <- mir.targets.gr[width(mir.targets.gr) <= 31]
  mir.targets.gr <- keepStandardChromosomes(mir.targets.gr)
  mir.targets.gr$seed.category <- NA_character_

  return(mir.targets.gr)
}

s <- build.starbase('data/base/starBase_Human_Pan-Cancer_MiRNA-TargetInteractions2015-01-13_22-20.csv')
saveRDS(s, 'data/processed/starbase.Rds')
