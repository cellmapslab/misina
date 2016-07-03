
library(ggplot2)
library(ggbio)
library(GenomicRanges)
library(data.table)
library(RCurl)
library(Cairo)
library(NCBI2R)
library(gwascat)
library(rtracklayer)
library(RSQLite)
library(dplyr)

within.seed <- function(snps, mir.bs) {
  findOverlaps(snps, flank(resize(mir.bs, 1, fix='end'), 6, start=T))
}

#we use fix='end' because these mir BSs are on the mRNA side, not in the miR side
snp.relative.pos <- function(snps, mir.bs) {
  d <- start(resize(mir.bs, 1, fix='end')) - start(snps)
  d <- d * ifelse(strand(mir.bs) == '+', 1, -1)
  d + ifelse(d >= 0, 1, 0) #1nt upstream of 5' end of miR is called -1, not zero
}

merge.granges.aggressively <- function(meta.columns, grs) {
  gr <- unname(grs)
  gr.lengths <- sapply(gr, length)

  #remove existing metadata columns temporarily and combine GRanges
  res <- do.call(c, lapply(gr, function(g){mcols(g) <- NULL; g}))

  #rowbind metadata columns using nice dplyr function
  mcols(res) <- bind_rows(lapply(gr, function(g){as.data.frame(mcols(g))}))

  #define additional metadata columns given in argument "meta.columns"
  meta.cols <- data.frame(sapply(names(meta.columns),
                            function(col){rep(meta.columns[[col]], times=gr.lengths)}),
                          stringsAsFactors = F)
  #add new metadata columns
  mcols(res) <- cbind(mcols(res), meta.cols)

  res
}

extend.with.LD <- function(snps.df, dataset=c('1kg', 'hapmap'),
                           aggregate.results=T, self.snp.label=NULL, ...) {

  dataset <- match.arg(dataset)

  CAD.SNP.ids <- do.call(rbind, lapply(regmatches(snps.df$SNPs,
                                                  regexec('.*(rs\\d+).*',
                                                          snps.df$SNPs)),
                                       `[`, c(2L)))

  CAD.SNP.ids <- unique(c(na.exclude(CAD.SNP.ids)))

  if (dataset == '1kg') {
    #use snipa to get LD information
    CAD.proxy.df <- snipa.get.ld.by.snp(CAD.SNP.ids, ...)
    colnames(CAD.proxy.df)[match(c('QRSID', 'RSID', 'DIST'), colnames(CAD.proxy.df))] <- c('SNP', 'Proxy', 'Distance')

  } else {
    #use hapmap through SNAP
    CAD.proxy.df <- SNAP.LD(CAD.SNP.ids, ...)
  }

  #filter out snps in LD with themselves
  #CAD.proxy.df <- subset(CAD.proxy.df, SNP != Proxy)

  #filter out unnecessary cols
  CAD.proxy.df <- CAD.proxy.df[, c('SNP', 'Proxy', 'Distance', 'R2')]

  #add SNPs that are excluded from the query result by SNiPA
  excluded.SNPs <- setdiff(CAD.SNP.ids, CAD.proxy.df$SNP)
  if (length(excluded.SNPs) > 0) {
    CAD.proxy.df <- rbind(CAD.proxy.df,
                          data.frame(SNP=excluded.SNPs,
                                     Proxy=excluded.SNPs,
                                     Distance=0,
                                     R2=1, stringsAsFactors = F))
  }

  CAD.proxy.df <- merge(snps.df, CAD.proxy.df, by.x='SNPs', by.y='SNP', all.x=T)

  na.rows <- apply(CAD.proxy.df[,c('Proxy', 'Distance', 'R2')], 1, function(r)sum(is.na(r))==3)
  CAD.proxy.df[na.rows,'Proxy']=CAD.proxy.df[na.rows,'SNPs']
  CAD.proxy.df[na.rows,'Distance']=0
  CAD.proxy.df[na.rows,'R2']=1

  if(!is.null(self.snp.label)) {
    #do not set original snp id to empty, otherwise it would vanish in
    #the aggregation step, and # of snps in IsProxyOf and those in t.stat column
    #would not match
    CAD.proxy.df[CAD.proxy.df$SNPs==CAD.proxy.df$Proxy,'SNPs']= self.snp.label
  }

  if (aggregate.results) {
    #group and collapse on proxies
    CAD.proxy.df <- aggregate(CAD.proxy.df,
                              list(Proxy=CAD.proxy.df$Proxy),
                              function(x)paste(unique(x[x!='']), collapse=','))
    CAD.proxy.df <- CAD.proxy.df[,-1]
  }

  colnames(CAD.proxy.df)[match(c('SNPs', 'Proxy'), colnames(CAD.proxy.df))] <- c('IsProxyOf', 'SNP')

  CAD.proxy.df
}

get.hg19.positions2 <- function(CAD.SNP.df, dbSNP.file) {
  message('Loading dbSNP...')

  ldb <- new(LevelDB, dbSNP.file)
  query.snps <- CAD.SNP.df$SNP
  dbsnp <- ldb$Get(query.snps)
  rm(ldb) #close connection properly
  gc()

  m <- sapply(dbsnp, anyNA)
  dbsnp <- as.data.frame(t(sapply(dbsnp[!m], rbind)), stringsAsFactors=F)
  dbsnp$SNP <- query.snps[!m]

  colnames(dbsnp) <- c('chr', 'start', 'strand', 'SNP')
  dbsnp$chr <- as.character(dbsnp$chr)
  dbsnp$strand <- as.character(dbsnp$strand)
  dbsnp$start <- as.integer(dbsnp$start)

  #no need to add +1 since we took 3rd column from the bed file

  #all.y makes sure that all SNPs in CAD.SNP.df will be in returned data.table
  #if some dbSNP ids are not found in dbSNP file, NAs will be added
  SNPs.final <- merge(dbsnp, CAD.SNP.df, by='SNP', all.y=T)
  na.index <- which(rowSums(is.na(SNPs.final)) > 0)
  if (length(na.index) > 0) {
    print('!!! Some SNPs are not available in dbSNP:')
    print(SNPs.final$SNP[na.index])
    SNPs.final <- SNPs.final[-na.index,]
  }

  CAD.snps.gr <- GRanges(seqnames=SNPs.final$chr,
                         ranges = IRanges(start=SNPs.final$start, width=1))
  mcols(CAD.snps.gr) <- SNPs.final[,-c(2,3,4)]
  genome(CAD.snps.gr) <- 'hg19'

    #also parse SNPs without any dbSNP ids
#   SNPs.noid.index <- substr(CAD.SNP.df, 1, 2) != 'rs'
#   if(any(SNPs.noid.index)) {
#     SNPs.noid <- data.frame(SNP=snps[SNPs.noid.index],
#                             chr=chr[SNPs.noid.index],
#                             start=as.integer(pos[SNPs.noid.index]))
#
#     SNPs.noid.gr <- GRanges(seqnames=SNPs.noid$chr,
#                             ranges = IRanges(start=SNPs.noid$start,width=1),
#                             SNPid=as.character(SNPs.noid$SNP),
#                             IsProxyOf=rep('', nrow(SNPs.noid)),
#                             LDDistance=rep(0, nrow(SNPs.noid)))
#     genome(SNPs.noid.gr) <- 'hg19'
#     CAD.snps.gr <- c(CAD.snps.gr, SNPs.noid.gr)
#   }

  return(CAD.snps.gr)


}

get.hg19.positions <- function(CAD.SNP.df, dbSNP.file) {
  message('Loading dbSNP...')

  conn <- dbConnect(RSQLite::SQLite(), dbSNP.file)
  query.snps <- paste0('"', paste0(CAD.SNP.df$SNP, collapse='","'), '"')

  #TODO: Use dbbind() when new RSQLite (> 1.0) is published
  query <- paste0('select V1,V2,V3,V4 from f where V4 in (', query.snps,')')
  #check if the query length is longer than sqlite limit
  if(nchar(query)>1e6) {
    #TODO: split snps into smaller pieces
    stop('dbSNP SQL query too long...')
  }

  dbsnp <- dbGetQuery(conn, query)
  dbDisconnect(conn)

  colnames(dbsnp) <- c('chr', 'start', 'end', 'SNP')

  #all.y makes sure that all SNPs in CAD.SNP.df will be in returned data.table
  #if some dbSNP ids are not found in dbSNP file, NAs will be added
  SNPs.final <- merge(dbsnp, data.table(CAD.SNP.df), by='SNP', all.y=T)
  na.index <- which(rowSums(is.na(SNPs.final)) > 0)
  if (length(na.index) > 0) {
    print('!!! Some SNPs are not available in dbSNP:')
    print(SNPs.final$SNP[na.index])
    SNPs.final <- SNPs.final[-na.index,]
  }

  #dbSNP file is off-by-one
  #see https://github.com/ropensci/rsnps/issues/16
  SNPs.final$start <- as.integer(SNPs.final$start) + 1

  CAD.snps.gr <- GRanges(seqnames=SNPs.final$chr,
                         ranges = IRanges(start=SNPs.final$start, width=1))
  mcols(CAD.snps.gr) <- SNPs.final[,-c(2,3,4)]
  genome(CAD.snps.gr) <- 'hg19'

  #also parse SNPs without any dbSNP ids
#   SNPs.noid.index <- substr(CAD.SNP.df, 1, 2) != 'rs'
#   if(any(SNPs.noid.index)) {
#     SNPs.noid <- data.frame(SNP=snps[SNPs.noid.index],
#                             chr=chr[SNPs.noid.index],
#                             start=as.integer(pos[SNPs.noid.index]))
#
#     SNPs.noid.gr <- GRanges(seqnames=SNPs.noid$chr,
#                             ranges = IRanges(start=SNPs.noid$start,width=1),
#                             SNPid=as.character(SNPs.noid$SNP),
#                             IsProxyOf=rep('', nrow(SNPs.noid)),
#                             LDDistance=rep(0, nrow(SNPs.noid)))
#     genome(SNPs.noid.gr) <- 'hg19'
#     CAD.snps.gr <- c(CAD.snps.gr, SNPs.noid.gr)
#   }

  return(CAD.snps.gr)


}

build.CAD.SNP.gr <- function(snps,
                             chr,
                             pos,
                             dbSNP.file, sep1=':', sep2=':',
                             dataset=c('1kg', 'hapmap')) {

  dataset <- match.arg(dataset)

  stopifnot(length(snps) == length(chr))
  stopifnot(length(snps) == length(pos))
  no.chr <- substr(chr, 1, 3) != 'chr'
  chr[no.chr] <- paste0('chr', chr[no.chr])

  CAD.SNP.ids <- do.call(rbind, lapply(regmatches(snps,
                                                  regexec('.*(rs\\d+).*',
                                                          snps)),
                                       `[`, c(2L)))

  CAD.SNP.ids <- unique(c(na.exclude(CAD.SNP.ids)))

  if (dataset == '1kg') {
    #use snipa to get LD information
    CAD.proxy.df <- snipa.get.ld.by.snp(CAD.SNP.ids)
    colnames(CAD.proxy.df)[match(c('QRSID', 'RSID', 'DIST'), colnames(CAD.proxy.df))] <- c('SNP', 'Proxy', 'Distance')

  } else {
    #use hapmap through SNAP
    CAD.proxy.df <- SNAP.LD(CAD.SNP.ids)
  }

  #filter out snps in LD with themselves
  CAD.proxy.df <- subset(CAD.proxy.df, SNP != Proxy)

  #group and collapse on proxies
  CAD.proxy.df <- aggregate(CAD.proxy.df[,c('SNP', 'Distance')],
                            list(Proxy=CAD.proxy.df$Proxy),
                            paste, collapse=',')

  CAD.SNP.df <- rbind(CAD.proxy.df, data.frame(list(Proxy=unique(CAD.SNP.ids), SNP='', Distance=0)))
  colnames(CAD.SNP.df) <- c('SNP', 'IsProxyOf', 'Distance')

  #use dbsnp sqlite instance to get SNP positions
  message('Loading dbSNP...')

  conn <- dbConnect(RSQLite::SQLite(), dbSNP.file)
  query.snps <- paste0('"', paste0(CAD.SNP.df$SNP, collapse='","'), '"')

  #TODO: Use dbbind() when new RSQLite (> 1.0) is published
  query <- paste0('select V1,V2,V3,V4 from f where V4 in (', query.snps,')')
  #check if the query length is longer than sqlite limit
  if(nchar(query)>1e6) {
    #TODO: split snps into smaller pieces
    stop('dbSNP SQL query too long...')
  }

  dbsnp <- dbGetQuery(conn, query)
  dbDisconnect(conn)

  colnames(dbsnp) <- c('chr', 'start', 'end', 'SNP')

  #all.y makes sure that all SNPs in CAD.SNP.df will be in returned data.table
  #if some dbSNP ids are not found in dbSNP file, NAs will be added
  SNPs.final <- merge(dbsnp, data.table(CAD.SNP.df), by='SNP', all.y=T)
  na.index <- which(rowSums(is.na(SNPs.final)) > 0)
  if (length(na.index) > 0) {
    print('!!! Some SNPs are not available in dbSNP:')
    print(SNPs.final$SNP[na.index])
    SNPs.final <- SNPs.final[-na.index,]
  }

  #dbSNP file is off-by-one
  #see https://github.com/ropensci/rsnps/issues/16
  SNPs.final$start <- as.integer(SNPs.final$start) + 1

  CAD.snps.gr <- GRanges(seqnames=SNPs.final$chr,
                         ranges = IRanges(start=SNPs.final$start, width=1),
                         SNP=as.character(SNPs.final$SNP),
                         IsProxyOf=SNPs.final$IsProxyOf,
                         LDDistance=SNPs.final$Distance)
  genome(CAD.snps.gr) <- 'hg19'

  #also parse SNPs without any dbSNP ids
  SNPs.noid.index <- substr(snps, 1, 2) != 'rs'
  if(any(SNPs.noid.index)) {
    SNPs.noid <- data.frame(SNP=snps[SNPs.noid.index],
                            chr=chr[SNPs.noid.index],
                            start=as.integer(pos[SNPs.noid.index]))

    SNPs.noid.gr <- GRanges(seqnames=SNPs.noid$chr,
                            ranges = IRanges(start=SNPs.noid$start,width=1),
                            SNP=as.character(SNPs.noid$SNP),
                            IsProxyOf=rep('', nrow(SNPs.noid)),
                            LDDistance=rep(0, nrow(SNPs.noid)))
    genome(SNPs.noid.gr) <- 'hg19'
    CAD.snps.gr <- c(CAD.snps.gr, SNPs.noid.gr)
  }

  return(CAD.snps.gr)

}

generate.final.table <- function(targets.gr, CAD.snps.gr, snp.mir.overlap.matrix,
                                 annotate=T, aggregate.results=T) {
  overlapping.targets.gr <- targets.gr[snp.mir.overlap.matrix[,2]]
  overlapping.snps <- CAD.snps.gr[snp.mir.overlap.matrix[,1]]

  final.table <- mcols(overlapping.targets.gr)
  final.table$mir.target.pos <- paste(seqnames(overlapping.targets.gr),
                                      start(overlapping.targets.gr),
                                      end(overlapping.targets.gr),
                                      strand(overlapping.targets.gr),
                                      sep = ':')

#   final.table$SNPid  <- mcols(overlapping.snps)$SNP
#   final.table$IsProxyOf <- mcols(overlapping.snps)$IsProxyOf
#   final.table$LDDistance <- mcols(overlapping.snps)$LDDistance
#   final.table$pval <- mcols(overlapping.snps)$pval
#   final.table$nstudies <- mcols(overlapping.snps)$n.studies

  final.table <- cbind(final.table, mcols(overlapping.snps))

  final.table$SNP.position  <- paste0(seqnames(overlapping.snps),
                                      ':',
                                      start(overlapping.snps))
  sbc <- ncol(mcols(overlapping.targets.gr))
  final.table <- final.table[,c((sbc+1):ncol(final.table), 1:sbc)]

  #add gwas catalog info

  data(gwrngs19) #from gwascat package
  gwas.df <- data.frame(mcols(gwrngs19)[c('SNPs', 'Disease.Trait', 'PUBMEDID', 'Reported.Gene.s.')])
  gwas.df.snp <- gwas.df[,-4]
  colnames(gwas.df.snp)[colnames(gwas.df.snp)=='SNPs'] <- 'SNP'

  final.table <- merge(final.table, gwas.df.snp,
                       all.x = T, sort = F, by='SNP')

  ncf <- ncol(final.table)
  colnames(final.table)[c(ncf-1, ncf)] <- c('DiseaseBySNP', 'DiseaseBySNP.PUBMEDID')

  #add disease info using Gene name
  #TODO: also consider Mapped_gene column of gwas.df
  #TODO: some rows have genes separated by commas
  #TODO: consider multiple matches:
  # "If there is more than one match, all possible matches contribute one row each."
  gwas.df.gene <- gwas.df[,-1]
  colnames(gwas.df.gene)[colnames(gwas.df.gene)=='Reported.Gene.s.'] <- 'gene'

  final.table <- merge(final.table, gwas.df.gene,
                       all.x = T, sort = F, by='gene')
  ncf <- ncol(final.table)
  colnames(final.table)[c(ncf-1, ncf)] <- c('DiseaseByGene', 'DiseaseByGene.PUBMEDID')
  #final.table <- final.table[,c(2,1, 3:ncf), with=F]
  final.table <- final.table[,c(2,1, 3:ncf)]

  if(annotate) {
    #add gene annotation info through NCBI2R package

    #annot <- data.table(AnnotateSNPList(unique(as.character(final.table$SNPid))))
    #annot <- annot[, c('marker', 'fxn_class', 'genesymbol', 'genename', 'genesummary'), with=F]
    annot.snps <- unique(as.character(final.table$SNP))
    annot.snps <- annot.snps[substr(annot.snps, 1, 2) == 'rs']
    annot <- AnnotateSNPList(annot.snps)
    annot <- annot[, c('marker', 'fxn_class', 'genesymbol', 'genename', 'genesummary')]
    #annot <- annot[, c(1, 6, 2, 14, 13)]
    colnames(annot)[colnames(annot)=='marker'] <- 'SNP'

    final.table <- merge(final.table, annot, all.x=T,
                         sort=F, by='SNP')

    #final.table <- final.table[, -which(names(final.table) == 'geneName'), with=F]
    final.table <- final.table[, -which(names(final.table) == 'genesymbol')]
  }

  #add one column for SNPs in miR seed regions
  snps.in.seed <- within.seed(overlapping.snps, overlapping.targets.gr)
  snp.seed.df <- data.frame(SNP=overlapping.snps[queryHits(snps.in.seed)]$SNP,
                            mir=overlapping.targets.gr[subjectHits(snps.in.seed)]$mir,
                            SNP.bwn.2.7=rep(T, length(snps.in.seed)),
                            stringsAsFactors = F)
  final.table <- merge(final.table, snp.seed.df, by.x=c('SNP', 'mir'),
                       by.y=c('SNP', 'mir'), all.x=T)
  #move new column next to mir column
  final.table <- final.table[, c(1:12, length(final.table), 13:(length(final.table)-1))]
  final.table$SNP.bwn.2.7[is.na(final.table$SNP.bwn.2.7)] <- F


  #add one more column for SNP pos
  snps.rel.pos <- snp.relative.pos(overlapping.snps, overlapping.targets.gr)
  snp.pos.df <- data.frame(SNP=overlapping.snps$SNP,
                            mir=overlapping.targets.gr$mir,
                            SNP.position.in.miR=snps.rel.pos,
                            stringsAsFactors = F)
  final.table <- merge(final.table, snp.pos.df, by.x=c('SNP', 'mir'),
                       by.y=c('SNP', 'mir'), all.x=T)
  #move new column next to mir column
  final.table <- final.table[, c(1:13, length(final.table), 14:(length(final.table)-1))]

  if (aggregate.results) {
    #merge duplicate columns through paste0(..., collapse=',')
    final.table <- aggregate(final.table,
                             list(Group=final.table$SNP, Group2=final.table$mir),
                             function(x)paste0(unique(x), collapse=','))
    final.table <- final.table[,-(1:2)] #exclude group col.
  }

  as.data.frame(final.table, stringsAsFactors=F)
}

get.mir.expression <- function(df, json=T) {
  if (json) {
    mirmine <- readRDS('data/processed/mirmine_median_json_2016_april.Rds')
  } else {
    mirmine <- readRDS('data/processed/mirmine_median_humanreadable_2016_april.Rds')
  }
  result <- merge(df, mirmine, all.x=T, sort=F)
  result$mir.split <- NULL

  if (json) {
    mirtissue <- readRDS('data/processed/mirtissueatlas_median_json_2016_may.Rds')
  } else {
    mirtissue <- readRDS('data/processed/mirtissueatlas_median_humanreadable_2016_april.Rds')
  }
  result <- merge(result, mirtissue, all.x=T, sort=F)
  result
}