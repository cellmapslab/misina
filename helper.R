library(dplyr)
library(RcppLevelDB)

source('R/parse-SNP.R')
source('R/snipe.R')

# input data files
dbsnp.file <- 'data/processed/dbsnp.leveldb'

starbase.gr <- readRDS('data/processed/starbase.Rds')
targetscan.gr <- readRDS('data/processed/targetscan.Rds')
miranda.gr <- readRDS('data/processed/miranda.Rds')

run.pipeline <- function(inputs) {
  
  total.snps <- inputs$snp.list
  ld.cutoff <- as.numeric(input$ld.cutoff)
  ld.population <- input$ld.population
  mir.target.db <- input$mir.target.db
  
  mir.targets.gr <- merge.granges.aggressively(meta.columns=list(mir.target.db=c('Starbase', 'TargetScan', 'miranda')),
                                               starbase.gr, targetscan.gr, miranda.gr)
  
  SNP.df <- extend.with.LD(total.snps, rsquare = ld.cutoff, self.snp.label = 'risk.snp')
  
#   SNP.gr <- get.hg19.positions2(SNP.df, dbSNP.file = dbsnp.file)
#   
#   snp.mir.overlap.hits <- findOverlaps(SNP.gr, unique(mir.targets.gr))
#   snp.mir.overlap.matrix <- as.matrix(snp.mir.overlap.hits)
#   if(length(snp.mir.overlap.matrix) == 0)
#     stop('There is no overlap between risk snps and miR targets.')
#   
#   result.table <- generate.final.table(unique(mir.targets.gr),
#                                        CAD.SNP.gr,
#                                        snp.mir.overlap.matrix, annotate=T)
#   
#   return(result.table)
  
  return(SNP.df)
}
