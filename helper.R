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

extract.snp.df <- function(inputs) {
  
  if (!is.null(inputs$grasp.pheno) && inputs$grasp.pheno != '') {
    s <- GRASP2() %>% tbl(., 'study')
    v <- GRASP2() %>% tbl(., 'variant')
    
    ph <- inputs$grasp.pheno
    Encoding(ph) <- 'utf-8' #revert back to original encoding
    
    ph.pmid <- s %>% select(PMID, PaperPhenotypeDescription) %>%
      filter(PaperPhenotypeDescription %in% ph) %>% select(PMID) %>% 
      as.data.frame %>% `[[`(., 1)
    
    ret <- v %>% select(SNPidInPaper, Phenotype, chr_hg19, pos_hg19, dbSNPfxn, PMID, Pvalue) %>% 
      filter(PMID %in% ph.pmid) %>% select(-PMID) %>%as.data.frame %>% rename(SNPs=SNPidInPaper) %>% as.data.frame
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
  mir.target.db <- inputs$mir.target.db
  
  mir.targets.gr <- merge.granges.aggressively(meta.columns=list(mir.target.db=c('Starbase', 'TargetScan', 'miranda')),
                                               starbase.gr, targetscan.gr, miranda.gr)
  
  SNP.df <- extend.with.LD(total.snps, rsquare = ld.cutoff, self.snp.label = 'risk.snp')
  SNP.gr <- get.hg19.positions2(SNP.df, dbSNP.file = dbsnp.file)
  
  snp.mir.overlap.hits <- findOverlaps(SNP.gr, unique(mir.targets.gr))
  snp.mir.overlap.matrix <- as.matrix(snp.mir.overlap.hits)
  if (length(snp.mir.overlap.matrix) == 0)
    stop('There is no overlap between risk snps and miR targets.')
  
  result.table <- generate.final.table(unique(mir.targets.gr),
                                       CAD.SNP.gr,
                                       snp.mir.overlap.matrix, annotate=T)
  
  cat('Finished pipeline...')
  return(result.table)
}
