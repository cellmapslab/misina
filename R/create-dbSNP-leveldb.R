
library(RcppLevelDB)
ldb = new(LevelDB, "/tmp/dbsnp.leveldb")

#iterate over files downloaded from dbSNP
for (f in list.files('~/dbsnp', '*.bed', full.names = T)) {

  cat('Processing file', f)
  bed = read.table(f, skip=1, stringsAsFactors = F)
  x = mapply(function(rsid, chr, pos, str){
    ldb$Put(rsid, list(chr, pos, str))
    }, bed$V4, bed$V1, bed$V3, bed$V6)
  gc()
}
