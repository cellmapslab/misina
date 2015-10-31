
library(data.table)
library(ggplot2)
library(GGally)
library(Cairo)
library(RColorBrewer)
library(mgcv)
source('R/parse-SNP.R')

all.snps <- fread('data/base/SCZ-GWAS-all-SNPs.tsv')
all.snps.gr <- GRanges(seqnames = all.snps$hg19chrc,
                       ranges = IRanges(start = all.snps$bp, width = 1),
                       SNPid = all.snps$snpid,
                       pval = -log10(all.snps$p))

risk.snps <-  read.table('output/data/snp-mir-SCZ-longerlist-newpipeline.tsv',
                         header=T, stringsAsFactors = F)$SNP

flanking.region <- 250000
ld.data <- LDSearch2(risk.snps,
                     RSquaredLimit = 0,
                     distanceLimit = flanking.region/1e3)
result <- list()

for(risk.snp in risk.snps) {

  risk.snp.gr <- unique(all.snps.gr[mcols(all.snps.gr)$SNPid==risk.snp])
  flanking.snps <- unique(subsetByOverlaps(all.snps.gr,risk.snp.gr + flanking.region))
  mcols(flanking.snps)$Distance <- start(risk.snp.gr) - start(flanking.snps)
  flanking.snps.df <- data.frame(mcols(flanking.snps))
  ld.data.risk.snp <- ld.data[ld.data$SNP == risk.snp,
                              c('Proxy', 'RSquared', 'DPrime', 'MAF')]
  colnames(ld.data.risk.snp)[1] <- 'SNPid'
  ld.data.risk.snp$RSquared <- as.numeric(ld.data.risk.snp$RSquared)
  ld.data.risk.snp$DPrime <- as.numeric(ld.data.risk.snp$DPrime)
  ld.data.risk.snp$MAF <- as.numeric(ld.data.risk.snp$MAF)

  res <- merge(flanking.snps.df, ld.data.risk.snp)
  res$Sentinel <- risk.snp
  result[[risk.snp]] <- res

  #generate locuszoom-like plot P-VALUES
  ggplot(res, aes(x=Distance,
                  y=pval,
                  color=RSquared)) +
    geom_point(size=4) +
    geom_point(data=subset(res, SNPid == risk.snp),
              aes(label=SNPid), color='red', size=6, shape=18) +
    #scale_color_gradientn(colours = brewer.pal(n=11, name='RdYlBu')) +
    labs(x='Distance', y='-log10(pvalue)', title=risk.snp) + theme_minimal()
  ggsave(paste0('output/figs/SCZ/mylocuszoom-pval/', risk.snp, '.pdf'),
         width = 14, height=12)


  #generate scatterplot matrix
  pdf(paste0('output/figs/SCZ/splom/splom-color-', risk.snp, '.pdf'), 16, 12)
  res$RSquaredDiscrete <- cut(res$RSquared, 5)
  ggplot <- function(...) ggplot2::ggplot(...) +
    scale_color_manual(values=colorRampPalette(c("darkblue", "red"))(5))
  unlockBinding("ggplot",parent.env(asNamespace("GGally")))
  assign("ggplot",ggplot,parent.env(asNamespace("GGally")))

  print(ggpairs(res,
                columns = 2:(ncol(res)-2),
                title=risk.snp,
                color='RSquaredDiscrete',
                group=1,
                params=c(shape=1),
                lower = list(continuous = "smooth")))
  dev.off()
  ggplot <- function(...) ggplot2::ggplot(...)
  unlockBinding("ggplot",parent.env(asNamespace("GGally")))
  assign("ggplot",ggplot,parent.env(asNamespace("GGally")))


}

final <- do.call(rbind, result)
final$Sentinel <- factor(final$Sentinel)


# plot splom for all SNPs
# workaround for color palette problem
ggplot <- function(...) ggplot2::ggplot(...) + scale_color_brewer(palette="Paired")
unlockBinding("ggplot",parent.env(asNamespace("GGally")))
assign("ggplot",ggplot,parent.env(asNamespace("GGally")))
CairoPNG('output/figs/SCZ/splom/allsnps.png', 1920, 1080)
ggpairs(final,
        columns = 2:(ncol(res)-2),
        title='All SNPs',
        #params=c(alpha=0.5),
        color='Sentinel')
dev.off()
# restore original ggplot
ggplot <- function(...) ggplot2::ggplot(...)
unlockBinding("ggplot",parent.env(asNamespace("GGally")))
assign("ggplot",ggplot,parent.env(asNamespace("GGally")))


ggplot(final, aes(x=Distance,
                  y=pval,
                  color=RSquared)) +
  geom_point(size=4, alpha=0.6) +
  geom_point(data=subset(final, SNPid %in% risk.snps),
             aes(label=SNPid), color='red', size=6, shape=18) +
  facet_wrap(~Sentinel, scales='free') +
  #scale_color_gradientn(colours = brewer.pal(n=11, name='RdYlBu')) +
  labs(x='Distance', y='-log10(pvalue)', title='SNPs in miR binding sites') +
  theme_minimal()

ggsave('output/figs/SCZ/mylocuszoom-pval/all-snps-pval.pdf',
       width = 18, height=12)

