suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(RcppLevelDB))
suppressPackageStartupMessages(library(grasp2db))

suppressPackageStartupMessages(source('R/parse-SNP.R'))
suppressPackageStartupMessages(source('R/snipe.R'))

# input data files
#dbsnp.file <- 'data/processed/dbsnp.leveldb'
dbsnp.file <- 'data/processed/dbSNP.GRCh37.p13.build142.sqlite'
#dbsnp.file <- '/storage/cmbstore/projects/misina/dbSNP.GRCh37.p13.build142.sqlite'

#GTEx db
dbsnp.file <- 'data/processed/GTExv6.sqlite'
#dbsnp.file <- '/storage/cmbstore/projects/misina/GTExv6.sqlite'

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
    
    ph.pmid <- s %>% dplyr::select(PMID, PaperPhenotypeDescription) %>%
      filter(PaperPhenotypeDescription %in% ph) %>% dplyr::select(PMID) %>% 
      as.data.frame %>% `[[`(., 1)
    
    #required to use %in% statement in dplyr
    ph.pmid <- as.list(ph.pmid)
    
    ret <- v %>% dplyr::select(SNPidInPaper, Phenotype, dbSNPfxn, PMID, Pvalue) %>% 
      filter(PMID %in% ph.pmid) %>% dplyr::select(-PMID) %>% as.data.frame 
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
  additional.columns <- colnames(total.snps)
  additional.columns <- additional.columns[additional.columns != 'SNPs']
  
  ld.cutoff <- as.numeric(inputs$ld.cutoff)
  ld.population <- inputs$ld.population
  mir.target.requested <- inputs$mir.target.db
  
  cat('Merging mir target datasets...')
  mir.targets.gr <- merge.granges.aggressively(meta.columns=list(mir.target.db=mir.target.requested),
                                               mir.target.avail[mir.target.requested])
  cat('Done')
  
  cat('Performing LD imputation...')
  SNP.df <- extend.with.LD(total.snps, 
                           rsquare = ld.cutoff, 
                           self.snp.label = 'risk.snp',
                           population = ld.population,
                           aggregate.results=F)
  cat('Done')
  
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
                                       annotate=F,
                                       aggregate=F)
  cat('Done')
  
  # eQTL enrichment analysis ------------------------------------------------
  
  cat('Performing eQTL enrichment...')
  gtex.eqtl <- readRDS('data/processed/GTExv6.Rds')
  names(gtex.eqtl) <- c('SNP', 'eQTL.beta', 'eQTL.tstat', 'eQTL.pvalue', 'eQTL.Gene', 'eQTL.Source', 'eQTL.Tissue')
  gtex.eqtl[] <- lapply(gtex.eqtl, as.character)
  
  ultimate <- merge(result.table, gtex.eqtl, all.x=T, by='SNP')
  #move eqtl columns towards the beginning
  ultimate <- dplyr::select(ultimate, SNP:mir.target.db, starts_with('eQTL'), everything())
  
  ultimate <- ultimate[order(ultimate$SNP),]
  #add one more column denoting if the target gene == eGene
  ultimate$eQTL.Gene.Same.as.Target.gene <- simplify2array(Map(function(gene, egene){
    any(toupper(gene) == strsplit(toupper(egene), ',')[[1]])},
    ultimate$gene, ultimate$eQTL.Gene))
  
  ultimate <- dplyr::select(ultimate, SNP:mir.target.db, starts_with('eQTL'), everything())
  cat('Done')
  
  cat('Finished pipeline...')
  
  #save names of additional columns
  attr(ultimate, 'additional.columns') <- additional.columns
  
  return(ultimate)
}


