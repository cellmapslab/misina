
mac.cis    <- read.delim('data/base/eqtl_CAD_mac_adjust_batch_covariates_age_gender_cis.res.txt',
                         stringsAsFactors = F) ; mac.cis$cell <- 'macrophage'; mac.cis$type <- 'cis'

mac.trans  <- read.delim('data/base/eqtl_CAD_mac_adjust_batch_covariates_age_gender_trans.res.txt',
                         stringsAsFactors = F) ; mac.trans$cell <- 'macrophage'; mac.trans$type <- 'trans'

mono.cis   <- read.delim('data/base/eqtl_CAD_mon_adjust_batch_covariates_age_gender_cis.res.txt',
                         stringsAsFactors = F) ; mono.cis$cell <- 'monocytes'; mono.cis$type <- 'cis'

mono.trans <- read.delim('data/base/eqtl_CAD_mon_adjust_batch_covariates_age_gender_trans.res.txt',
                         stringsAsFactors = F) ; mono.trans$cell <- 'monocytes'; mono.trans$type <- 'trans'

eqtl_cad <- rbind(mac.cis, mac.trans, mono.cis, mono.trans)

#library(BiocInstaller)
#biocLite('illuminaHumanv3.db')
library(illuminaHumanv3.db)

#filter out bad probes
probe.qual <- unlist(as.list(illuminaHumanv3PROBEQUALITY[unique(eqtl_cad$gene)]))
stopifnot(length(probe.qual) == length(unique(eqtl_cad$gene)))
bad.probes <- names(probe.qual[probe.qual %in% c('Bad', 'No match')])
eqtl_cad <- eqtl_cad[!(eqtl_cad$gene %in% bad.probes),]

#library(lumi)
#wow, nice bug: https://support.bioconductor.org/p/70706/
#do not forget to unique() the list before calling probeID2nuID function
# eqtl_ids  <- lumi::probeID2nuID(unique(eqtl_cad$gene),
#              lib.mapping = 'lumiHumanIDMapping',
#              chipVersion = 'HumanRef8_V3_0_R3_11282963_A')

#use re-annotated "romoat" gene symbols from illuminaHumanv3.db package
probe2gene <- as.data.frame(illuminaHumanv3SYMBOLREANNOTATED[unique(eqtl_cad$gene)])

#warn if there are unmapped probe IDs
if (nrow(probe2gene) != length(unique(eqtl_cad$gene))) {
  warning(sprintf('%d probe IDs are not mapped...', length(unique(eqtl_cad$gene))-nrow(probe2gene)))
}

eqtl_cad <- merge(eqtl_cad, probe2gene, by.x = 'gene', by.y = 'IlluminaID', all.y = T) #merge only those are mapped

library(dplyr)
library(magrittr)
eqtl_cad_tbl <- eqtl_cad %>% tbl_df %>% dplyr::select(-gene) %>% rename(gene=SymbolReannotated)

#select only significant associations
eqtl_cad_tbl %<>% filter(FDR < 0.05)

#filter out non-rs SNP ids
eqtl_cad_tbl %<>% filter(substr(SNP, 1, 2) == 'rs')

#filter out multiple p-values for SNP-celltype-gene triplets and
#pick the one with best pvalue
eqtl_cad_tbl %<>%
  group_by(SNP, cell, gene) %>%
  slice(which.min(FDR)) %>%
  ungroup

eqtl_cad_tbl %<>% rename(SNPs=SNP)

extended.eqtl <- extend.with.LD(as.data.frame(eqtl_cad_tbl), self.snp.label = 'direct', aggregate.results = F)
extended.eqtl.tbl <- extended.eqtl %>% tbl_df %>% dplyr::select(SNP:R2, IsProxyOf:gene)
extended.eqtl.tbl %<>% rename(eQTL.tissue=cell)
extended.eqtl.tbl$eQTL.source <- 'DHM'

saveRDS(as.data.frame(extended.eqtl.tbl), 'data/processed/eQTL-mono-macro-with-LD-nonaggregated.Rds')


# probeID2nuID(c('ILMN_2271149', 'ILMN_1709590', 'ILMN_1781388'),
#              lib.mapping = 'lumiHumanIDMapping',
#              chipVersion = 'HumanRef8_V3_0_R3_11282963_A')
#
# nuID2EntrezID(c('cFXlpCSLtEnp.XR_oA', '9ZuoXv7t7.n1_oGkS0', 'QSZuoXv7t7.n1_oGkQ'),
#               lib.mapping = 'lumiHumanIDMapping')