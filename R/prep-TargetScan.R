

# targetscan <- read.delim('data/base/Predicted_Targets_Info.txt',
#                         stringsAsFactors = F,
#                         na.strings = 'NULL')
#select only Hsapiens
#targetscan <- targetscan[targetscan$Species.ID==9606,]

targetscan <- read.delim('data/base/Conserved_Site_Context_Scores.txt',
                        stringsAsFactors = F,
                        na.strings = 'NULL')
#select only Hsapiens
targetscan <- targetscan[targetscan$Gene.Tax.ID==9606,]



#refGenes table from UCSC to map NM_ ids to genomic coordinates
ref <- read.delim('data/base/refGenes_UCSC.tsv', stringsAsFactors = F)
#ref <- ref[substr(ref$name, 1, 2) == 'NM',]
plus.ind <- ref$strand == '+'
ref.plus <- ref[plus.ind,]
neg.ind <- ref$strand == '-'
ref.neg <- ref[neg.ind,]

#Use txEnd
ref.df <- data.frame(id = c(ref.plus$name, ref.neg$name),
                     chr = c(ref.plus$chrom, ref.neg$chrom),
                     start = c(ref.plus$cdsEnd, ref.neg$txStart),
                     end = c(ref.plus$txEnd, ref.neg$cdsStart),
                     strand = c(ref.plus$strand, ref.neg$strand),
                     symbol = c(ref.plus$name2, ref.neg$name2), stringsAsFactors = F)

ref.df <- ref.df[ref.df$end - ref.df$start != 0,] #filter out ~800 0-UTR transcripts
ref.df$length <- ref.df$end - ref.df$start


targetscan.narrow <- targetscan[, c('Transcript.ID', 'miRNA', 'UTR_start', 'UTR.end')]

merged.mir <- merge(targetscan.narrow, ref.df, by.x='Transcript.ID', by.y='id')

w <- which((merged.mir$UTR_start > merged.mir$length) | (merged.mir$UTR.end > merged.mir$length))
sum((merged.mir$UTR_start > merged.mir$length) | (merged.mir$UTR.end > merged.mir$length))

# NM_001003716 is a nice example that is outside the UTR_start UTR_end range :'(