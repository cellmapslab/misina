suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(RcppLevelDB))
suppressPackageStartupMessages(library(grasp2db))

suppressPackageStartupMessages(source('R/parse-SNP.R'))
suppressPackageStartupMessages(source('R/snipe.R'))

# input data files
dbsnp.file <- 'data/processed/dbsnp.leveldb'

starbase.gr <- readRDS('data/processed/starbase.Rds')
targetscan.gr <- readRDS('data/processed/targetscan.Rds')
miranda.gr <- readRDS('data/processed/miranda.Rds')
mir.target.avail <- list(Starbase = starbase.gr,
                         miranda = miranda.gr,
                         TargetScan = targetscan.gr)

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
    
    ph.pmid <- s %>% select(PMID, PaperPhenotypeDescription) %>%
      filter(PaperPhenotypeDescription %in% ph) %>% select(PMID) %>% 
      as.data.frame %>% `[[`(., 1)
    
    #required to use %in% statement in dplyr
    ph.pmid <- as.list(ph.pmid)
    
    ret <- v %>% select(SNPidInPaper, Phenotype, chr_hg19, pos_hg19, dbSNPfxn, PMID, Pvalue) %>% 
      filter(PMID %in% ph.pmid) %>% select(-PMID) %>% as.data.frame 
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
  ld.cutoff <- as.numeric(inputs$ld.cutoff)
  ld.population <- inputs$ld.population
  mir.target.requested <- inputs$mir.target.db
  
  mir.targets.gr <- merge.granges.aggressively(meta.columns=list(mir.target.db=mir.target.requested),
                                               mir.target.avail[mir.target.requested])
  
  SNP.df <- extend.with.LD(total.snps, 
                           rsquare = ld.cutoff, 
                           self.snp.label = 'risk.snp',
                           population = ld.population)
  
  SNP.gr <- get.hg19.positions2(SNP.df, dbSNP.file = dbsnp.file)
  
  snp.mir.overlap.hits <- findOverlaps(SNP.gr, unique(mir.targets.gr))
  snp.mir.overlap.matrix <- as.matrix(snp.mir.overlap.hits)
  if (length(snp.mir.overlap.matrix) == 0)
    stop('There is no overlap between risk snps and miR targets.')
  
  result.table <- generate.final.table(unique(mir.targets.gr),
                                       SNP.gr,
                                       snp.mir.overlap.matrix, 
                                       annotate=F,
                                       aggregate=F)
  
  cat('Finished pipeline...')
  return(result.table)
}
