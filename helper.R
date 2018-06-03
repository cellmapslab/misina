suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(grasp2db))

suppressPackageStartupMessages(source('R/parse-SNP.R'))
suppressPackageStartupMessages(source('R/snipe.R'))

# input data files
dbsnp.file <- 'data/processed/dbSNP.GRCh37.p13.build142.sqlite'

#GTEx db
gtexdb.file <- 'data/processed/GTExv6.sqlite'

extract.snp.df <- function(inputs) {

  if (!is.null(inputs$grasp.pheno) && inputs$grasp.pheno != '') {
    g <- GRASP2()
    dbDisconnect(g$con) #workaround for ugly rsqlite error
    g$con <- dbConnect(g$con)
    s <- g %>% tbl(., 'study')
    v <- g %>% tbl(., 'variant')

    ph <- inputs$grasp.pheno
    Encoding(ph) <- 'utf-8' #revert back to original encoding
    #required to use %in% statement in dplyr
    ph <- as.list(ph)

    ph.pmid <- s %>% dplyr::select(PMID, PaperPhenotypeDescription) %>%
      filter(PaperPhenotypeDescription %in% ph) %>% dplyr::select(PMID) %>%
      as.data.frame %>% `[[`(., 1)

    #required to use %in% statement in dplyr
    ph.pmid <- as.list(ph.pmid)

    ret <- v %>% dplyr::select(SNPidInPaper, Phenotype, dbSNPfxn, PMID, Pvalue, NegativeLog10PBin) %>%
      filter(PMID %in% ph.pmid) %>% filter(NegativeLog10PBin >= 8) %>% dplyr::select(-PMID) %>% as.data.frame
    colnames(ret)[colnames(ret) == 'SNPidInPaper'] <- 'SNPs'

    ret <- aggregate(ret,
                     list(SNPs=ret$SNPs),
                     function(x)paste0(unique(x), collapse=','))
    ret <- ret[,-1]

  } else {
    ret <- inputs$snp.df
  }

  return(ret)

}

run.pipeline <- function(inputs) {

  cat('Starting pipeline...')
  total.snps <- extract.snp.df(inputs)
  original.snp.count <- length(unique(total.snps$SNPs))
  additional.columns <- colnames(total.snps)
  additional.columns <- additional.columns[additional.columns != 'SNPs']

  ld.cutoff <- as.numeric(inputs$ld.cutoff)
  ld.population <- inputs$ld.population
  mir.target.requested <- inputs$mir.target.db

  cat('Loading mir target datasets...')
  mir.targets.gr <- readRDS('data/processed/mir-all-targets.Rds')
  mir.targets.gr <- mir.targets.gr[mir.targets.gr$mir.target.db %in% mir.target.requested]
  cat('Done')

  cat('Performing LD imputation...')
  SNP.df <- extend.with.LD(total.snps,
                           rsquare = ld.cutoff,
                           self.snp.label = 'risk.snp',
                           population = ld.population,
                           aggregate.results=F)
  cat('Done')
  total.snp.count <- length(unique(SNP.df$SNP))

  cat('Getting hg19 positions of all SNPs...')
  SNP.gr <- get.hg19.positions(SNP.df, dbSNP.file = dbsnp.file)
  cat('Done')

  cat('Generating report now...')
  snp.mir.overlap.hits <- findOverlaps(SNP.gr, unique(mir.targets.gr))
  snp.mir.overlap.matrix <- as.matrix(snp.mir.overlap.hits)
  if (length(snp.mir.overlap.matrix) == 0)
    stop('There is no overlap between risk snps and miR targets.')

  result.table <- generate.final.table(unique(mir.targets.gr),
                                       SNP.gr,
                                       snp.mir.overlap.matrix,
                                       annotate = F,
                                       aggregate.results = F)
  cat('Done')

  # eQTL enrichment analysis ------------------------------------------------

  cat('Performing eQTL enrichment...')
  gtex <- src_sqlite(gtexdb.file)
  gtex.eqtl <- tbl(gtex, 'GTExv6')

  # rename variables separately: https://github.com/tidyverse/dplyr/issues/2943
  gtex.eqtl <- dplyr::rename(gtex.eqtl, eQTL.beta=beta)
  gtex.eqtl <- dplyr::rename(gtex.eqtl, eQTL.tstat=t_stat)
  gtex.eqtl <- dplyr::rename(gtex.eqtl, eQTL.pvalue=p_value)
  gtex.eqtl <- dplyr::rename(gtex.eqtl, eQTL.Gene=gene_name)
  gtex.eqtl <- dplyr::rename(gtex.eqtl, eQTL.Source=eQTL.source)
  gtex.eqtl <- dplyr::rename(gtex.eqtl, eQTL.Tissue=eQTL.tissue)

  tmp.snps <- as.list(result.table$SNP)
  gtex.eqtl <- dplyr::collect(gtex.eqtl %>% dplyr::filter(SNP %in% tmp.snps))

  ultimate <- dplyr::left_join(result.table, gtex.eqtl, by='SNP', copy=T)
  #move eqtl columns towards the beginning
  ultimate <- dplyr::select(ultimate, SNP:mir.target.db, starts_with('eQTL'), everything())

  ultimate <- ultimate[order(ultimate$SNP),]
  #add one more column denoting if the target gene == eGene
  ultimate$eQTL.identical.target <- simplify2array(Map(function(gene, egene){
    any(toupper(gene) == strsplit(toupper(egene), ',')[[1]])},
    ultimate$gene, ultimate$eQTL.Gene))

  ultimate <- dplyr::select(ultimate, SNP:mir.target.db, starts_with('eQTL'), everything())
  cat('Done')

  cat('Finished pipeline...')

  #save names of additional columns
  attr(ultimate, 'additional.columns') <- additional.columns
  attr(ultimate, 'original.snp.count') <- original.snp.count
  attr(ultimate, 'total.snp.count') <- total.snp.count

  return(ultimate)
}


