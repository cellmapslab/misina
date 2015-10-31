
library(data.table)
library(ggplot2)
library(GGally)
library(Cairo)
library(RColorBrewer)
library(mgcv)
source('R/parse-SNP.R')

all.snps <- fread('data/processed/1000G.MA.cad.add.fixrand.subset.Feb2015.CW.csv')
all.snps.gr <- GRanges(seqnames = paste0('chr', all.snps$chr),
                       ranges = IRanges(start = all.snps$bp_hg19, width = 1),
                       SNPid = all.snps$legendrs,
                       pval = -log10(all.snps$p_dgc),
                       qval = -log10(all.snps$q_pvalue),
                       n.studies = all.snps$n_studies)

risk.snps <-  read.table('output/data/Risk-SNPs-within-miR-BS-corrected-eQTL.tsv',
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
  ggsave(paste0('output/figs/mylocuszoom-pval/', risk.snp, '.pdf'),
         width = 14, height=12)

  #generate locuszoom-like plot Q-VALUES
  ggplot(res, aes(x=Distance,
                  y=qval,
                  color=RSquared)) +
    geom_point(size=4) +
    geom_point(data=subset(res, SNPid == risk.snp),
               aes(label=SNPid), color='red', size=6, shape=18) +
    #scale_color_gradientn(colours = brewer.pal(n=11, name='RdYlBu')) +
    labs(x='Distance', y='-log10(qvalue)', title=risk.snp) + theme_minimal()
  ggsave(paste0('output/figs/mylocuszoom-qval/', risk.snp, '.pdf'),
         width = 14, height=12)



  #generate scatterplot matrix
  pdf(paste0('output/figs/splom/splom-color-', risk.snp, '.pdf'), 16, 12)
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
CairoPNG('output/figs/splom/allsnps.png', 1920, 1080)
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


# regional assoc. plot with SENTINEL facet
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

ggsave('output/figs/mylocuszoom-pval/all-snps-pval.pdf',
       width = 18, height=12)


# regional assoc. plot with p values + SENTINEL + n.studies facet
ggplot(final, aes(x=Distance,
                  y=pval,
                  color=RSquared)) +
  geom_point(size=2, alpha=0.6) +
  geom_point(data=subset(final, SNPid %in% risk.snps),
             aes(label=SNPid), color='red', size=3, shape=18) +
  facet_grid(Sentinel~n.studies, scales='free') +
  #scale_color_gradientn(colours = brewer.pal(n=11, name='RdYlBu')) +
  labs(x='Distance', y='-log10(pvalue)', title='SNPs in miR binding sites')

ggsave('output/figs/mylocuszoom-pval/all-snps-pval-with-nstudies.pdf',
       width = 18, height=12)

# regional assoc. plot with q values + SENTINEL facet

ggplot(final, aes(x=Distance,
                  y=qval,
                  color=RSquared)) +
  geom_point(size=4, alpha=0.6) +
  geom_point(data=subset(final, SNPid %in% risk.snps),
             aes(label=SNPid), color='red', size=6, shape=18) +
  facet_wrap(~Sentinel, scales='free') +
  #scale_color_gradientn(colours = brewer.pal(n=11, name='RdYlBu')) +
  labs(x='Distance', y='-log10(qvalue)', title='SNPs in miR binding sites') +
  theme_minimal()

ggsave('output/figs/mylocuszoom-qval/all-snps-qval.pdf',
       width = 18, height=12)
