

library(gwascat)

data("gwrngs19")

phenotypes <- c("Coronary artery disease",
                "Coronary artery disease or ischemic stroke",
                "Coronary artery disease or large artery stroke",
                "Coronary heart disease",
                "Stroke",
                "Myocardial infarction")

gwas.risk <- as(subsetByTraits(gwrngs19, tr=phenotypes), 'GRanges')
mcols(gwas.risk)$SNPs <- trimws(mcols(gwas.risk)$SNPs)
saveRDS(gwas.risk, 'data/processed/gwas-heart-diseases.Rds')
