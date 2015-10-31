
##########################  CAD risk snps
CAD.SNP.df <- read.csv2('data/base/knownLociCAD_Jan2015_CW.csv', stringsAsFactors = F, strip.white = T) #new list from Jeanette
conf.snps <- data.frame(SNPs=CAD.SNP.df$Marker_ID,
                        Phenotype='Coronary artery disease',
                        Risk.SNP.Source='Confidential', stringsAsFactors = F)

####CAD.SNP.df <- read.csv2('CAD_risk_SNPs_Nov2014.csv')
stroke.snps <- read.csv('data/processed/stroke-snps.csv', stringsAsFactors = F)

nature.snps.raw <- read.csv('data/base/nature-cardiogram.csv', stringsAsFactors = F, strip.white = T)
nature.snps <- data.frame(SNPs=nature.snps.raw$rs_number,
                          Phenotype='Coronary artery disease',
                          Risk.SNP.Source='CARDIOGRAMplusC4D',
                          stringsAsFactors = F)

gwas.cat.snps.raw <- readRDS('data/processed/gwas-heart-diseases.Rds')
gwas.cat.snps <- data.frame(SNPs=gwas.cat.snps.raw$SNPs,
                            Phenotype=gwas.cat.snps.raw$Disease.Trait,
                            Risk.SNP.Source='GWAS catalog',
                            stringsAsFactors = F)
total.snps <- rbind(conf.snps, stroke.snps, nature.snps, gwas.cat.snps)
total.snps <- aggregate(total.snps,
                         list(SNPs=total.snps$SNPs),
                         function(x)paste0(unique(x), collapse=','))
total.snps <- total.snps[,-1]

write.table(total.snps, 'data/processed/all-heart-risk-SNPs.tsv', row.names = F,
          quote = F, sep='\t')
