source('helper.R')

df <- data.frame(SNPs=c('rs10476052', 'rs1048920', 'rs1049633', 'rs113767110',
                        'rs11739062', 'rs11749762', 'rs1931895'),
                 stringsAsFactors=F)

results <- run.pipeline(list(snp.df=df,
                            ld.cutoff=0.8,
                            ld.population='eur',
                            mir.target.db=c('TargetScan','miranda','Starbase')))
